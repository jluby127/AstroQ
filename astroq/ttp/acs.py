"""Ant Colony System heuristic for the Traveling Telescope Problem.

A fast primal heuristic for the TTP that operates on the same precomputed
``nodes`` / ``arcs`` data produced by :class:`astroq.ttp.model.TTPModel`. It is
an adaptation of the time-dependent OPTW Ant Colony System of Verbeeck,
Vansteenwegen & Aghezzaf (2017, Ann. Oper. Res. 254:481-505) to the TTP:

* Time-dependent "travel" (slew) is **piecewise-constant per window**: the slew
  from node ``i`` to ``j`` departing at minute ``t`` is the precomputed
  worst-case value ``arcs[(i, j, window_of(t))]`` -- simpler than the paper's
  piecewise-linear ``mu*w + nu``.
* The objective mirrors the MILP (Handley+ 2024 eq. 10 as implemented in
  ``TTPModel.build_model``): maximize ``sum(priority) - slew_penalty*t_slew
  - slew_penalty*idle_ratio*t_idle_between``.
* Time windows map as ``o_i -> t_early``, ``c_i -> t_late`` on the visit
  *completion* time ``ti`` (TTP convention: ``t_early + t_visit <= ti <=
  t_late``); the first visit is pinned to ``t_start`` like the MILP's
  ``first_exposure_at_start`` constraint.

For correctness at the small TTP scale (``N`` ~ tens of nodes) every candidate
move is checked by an O(N) forward simulation (:meth:`_evaluate`) rather than
the paper's incremental ``max_shift`` bookkeeping; this trades a little speed
for a formulation that is obviously faithful to the MILP feasibility region.

The solver is state-agnostic for construction: with ``n_states > 1`` the arc
cost used during search is the minimum over feasible wrap-state pairs, and a
per-node wrap state is assigned afterwards by a shortest-path pass over the
fixed visit order (:meth:`assign_states`).
"""

import logging
import time

import numpy as np

logs = logging.getLogger(__name__)


def _solve_with_seed(solver, seed):
    """Top-level worker for parallel multi-start (must be picklable).

    Reseeds the (already-pickled) solver and runs one full ACS search.
    """
    solver.rng = np.random.default_rng(int(seed))
    return solver.solve()


class ACSSolver:
    """Ant Colony System heuristic over a built TTP node/arc catalog.

    Args:
        nodes (pd.DataFrame): ``TTPModel.nodes`` (anchor, real..., anchor).
        arcs_lookup (dict): ``{(i, j, m): t_slew}`` for single state, or
            ``{(i, j, m, si, sj): t_slew}`` for ``n_states > 1``.
        w (np.ndarray): slot boundaries in minutes, length ``M + 1``.
        dur_min (float): night duration in minutes.
        multi_visit_groups (dict): ``{unique_id: [node_id, ...]}`` in visit_seq
            order for requests with more than one visit.
        N (int): number of nodes (including the two anchors).
        M (int): number of slew windows.
        slew_penalty (float): objective weight on total slew minutes.
        idle_ratio (float): idle weight as a fraction of ``slew_penalty``.

    Keyword Args:
        node_states (dict | None): ``{node_id: [feasible_state, ...]}`` when
            ``n_states > 1``.
        n_states (int): number of cable-wrap states.
        rng (np.random.Generator | None): random source.
        params (dict | None): ACS hyperparameters (see :attr:`_DEFAULT_PARAMS`).
    """

    _DEFAULT_PARAMS = {
        "alpha": 1.0,        # pheromone weight
        "beta": 2.0,         # greedy (heuristic) weight
        "rho": 0.1,          # local pheromone evaporation
        "rho_global": 0.1,   # global pheromone reinforcement
        "max_ants": 10,      # solutions constructed per iteration
        "nb_max": 25,        # neighbor-list cap per node
        "tau_init": 1.0,
        "time_limit_s": 3.0,
        "max_iter": 200,
        "reset_no_improve": 30,
        "swap_window": 4,    # only swap positions within this span
    }

    def __init__(
        self,
        nodes,
        arcs_lookup,
        w,
        dur_min,
        multi_visit_groups,
        N,
        M,
        slew_penalty,
        idle_ratio,
        *,
        node_states=None,
        n_states=1,
        rng=None,
        params=None,
    ):
        self.nodes = nodes
        self.arcs_lookup = arcs_lookup
        self.w = np.asarray(w, dtype=float)
        self.dur_min = float(dur_min)
        self.multi_visit_groups = multi_visit_groups or {}
        self.N = int(N)
        self.M = int(M)
        self.slew_penalty = float(slew_penalty)
        self.idle_ratio = float(idle_ratio)
        self.n_states = int(n_states)
        self.node_states = node_states or {}
        self.rng = rng if rng is not None else np.random.default_rng(0)
        self.params = {**self._DEFAULT_PARAMS, **(params or {})}

        # Real nodes are 1..N-2; 0 and N-1 are the start / end anchors.
        self.real_nodes = list(range(1, self.N - 1))
        self.t_early = nodes["t_early"].to_numpy(dtype=float)
        self.t_late = nodes["t_late"].to_numpy(dtype=float)
        self.t_visit = nodes["t_visit"].to_numpy(dtype=float)
        self.tau_intra = nodes["tau_intra"].to_numpy(dtype=float)
        self.priority = nodes["priority"].to_numpy(dtype=float)
        self.uid = nodes["unique_id"].to_numpy()

        # First exposure is pinned to the start (or earliest feasible) time.
        if self.real_nodes:
            self.t_start = max(0.0, float(self.t_early[self.real_nodes].min()))
        else:
            self.t_start = 0.0

        # node_id -> (group_list, position) for tau_intra ordering checks.
        self._group_of = {}
        for uid, members in self.multi_visit_groups.items():
            for pos, nid in enumerate(members):
                self._group_of[nid] = (members, pos)

        self._deadline = None
        self._build_arc_cost()
        self._build_neighbors()
        self._build_rev_neighbors()

    # ------------------------------------------------------------------ setup
    def _window_of(self, t_depart):
        """Slot index ``m`` containing departure minute ``t_depart``."""
        m = int(np.searchsorted(self.w, t_depart, side="right") - 1)
        return min(max(m, 0), self.M - 1)

    def _build_arc_cost(self):
        """Collapse the arc catalog to a single cost per ``(i, j, m)``.

        For ``n_states > 1`` the search cost is the minimum slew over feasible
        wrap-state pairs (states are assigned later on the fixed order).
        """
        cost = {}
        if self.n_states > 1:
            for (i, j, m, si, sj), val in self.arcs_lookup.items():
                key = (i, j, m)
                if key not in cost or val < cost[key]:
                    cost[key] = val
        else:
            cost = dict(self.arcs_lookup)
        self._cost = cost

    def _arc_cost(self, i, j, m):
        """Slew minutes for arc ``(i, j)`` departing in window ``m``.

        Anchor arcs (``i == 0`` or ``j == N-1``) are free; missing real arcs
        fall back to ``0.0`` (kept consistent with the MILP's ``.get(.., 0.0)``).
        """
        if i == 0 or j == self.N - 1:
            return 0.0
        return self._cost.get((i, j, m), 0.0)

    def _build_neighbors(self):
        """Reward-ranked, reachability-pruned neighbor lists per node.

        ``j`` is a neighbor of ``i`` when it is time-window-reachable from
        ``i`` at the earliest departure, ranked by ``priority_j / (t_visit_j +
        slew)`` (Verbeeck eq. 5). Capped at ``nb_max``.
        """
        nb_max = int(self.params["nb_max"])
        neighbors = {}
        # window at the earliest possible departure, used only for ranking.
        for i in [0, *self.real_nodes]:
            t_dep = self.t_start if i == 0 else max(
                self.t_early[i] + self.t_visit[i], self.t_start
            )
            m = self._window_of(t_dep)
            cands = []
            for j in self.real_nodes:
                if j == i:
                    continue
                slew = self._arc_cost(i, j, m)
                arrive = max(t_dep + slew, self.t_early[j])
                comp = arrive + self.t_visit[j]
                if comp > min(self.t_late[j], self.dur_min) + 1e-9:
                    continue
                ratio = self.priority[j] / max(self.t_visit[j] + slew, 1e-6)
                cands.append((ratio, j))
            cands.sort(reverse=True)
            neighbors[i] = [j for _, j in cands[:nb_max]]
        self.neighbors = neighbors

    def _build_rev_neighbors(self):
        """Inverse of :attr:`neighbors`: ``rev[j]`` = nodes that can precede j.

        Used to restrict insert / replace candidate positions to geometric
        neighbors so local search visits a small candidate set per move.
        """
        rev = {j: set() for j in range(self.N)}
        for i, lst in self.neighbors.items():
            for j in lst:
                rev[j].add(i)
        self.rev_neighbors = rev

    # ------------------------------------------------------------- evaluation
    def _resolve(self, order):
        """Forward-simulate a visit order into a *realizable* schedule.

        Walks ``order`` and, for ``n_states > 1``, picks each node's wrap state
        to maintain continuity with the previous kept node (the in- and out-arc
        of a visit share its state). A node is **dropped** when no feasible
        wrap state / arc / time-window placement exists, so the returned
        ``order`` is always physically realizable.

        Returns a dict ``{order, states, ti, objective, feasible}`` where
        ``order`` is the kept subsequence and ``ti`` its completion minutes.
        """
        kept, states, ti = [], {}, []
        total_slew = 0.0
        group_last = {}
        prev, prev_state, last_comp = 0, 0, None

        for v in order:
            dep = self.t_start if prev == 0 else last_comp
            m = self._window_of(dep)
            best = None  # (state, slew, comp)
            for sv in (self.node_states.get(v, [0]) if self.n_states > 1 else [0]):
                if prev == 0:
                    slew = 0.0  # start-anchor arcs are free in every state
                elif self.n_states > 1:
                    val = self.arcs_lookup.get((prev, v, m, prev_state, sv))
                    if val is None:
                        continue
                    slew = val
                else:
                    slew = self._arc_cost(prev, v, m)
                arrive = max(dep + slew, self.t_early[v])
                comp = arrive + self.t_visit[v]
                if comp > min(self.t_late[v], self.dur_min) + 1e-9:
                    continue
                if not self._group_ready(v, comp, group_last):
                    continue
                if best is None or slew < best[1]:
                    best = (sv, slew, comp)
            if best is None:
                continue  # drop v: no feasible state/time placement
            sv, slew, comp = best
            kept.append(v)
            states[v] = sv
            ti.append(comp)
            total_slew += slew
            self._group_ready(v, comp, group_last, commit=True)
            prev, prev_state, last_comp = v, sv, comp

        if not kept:
            return {"order": [], "states": {}, "ti": [], "objective": 0.0,
                    "feasible": True}
        sum_visit = float(self.t_visit[kept].sum())
        idle = ti[-1] - sum_visit - total_slew
        obj = (
            float(self.priority[kept].sum())
            - self.slew_penalty * total_slew
            - self.slew_penalty * self.idle_ratio * idle
        )
        return {"order": kept, "states": states, "ti": ti, "objective": obj,
                "feasible": True}

    def _evaluate(self, order):
        """Objective of a (resolved) visit order; thin wrapper over _resolve."""
        res = self._resolve(order)
        return res["feasible"], res["objective"], res

    # ----------------------------------------------------------- construction
    def _construct(self, tau):
        """Build one solution by pheromone/greedy roulette construction."""
        alpha, beta = self.params["alpha"], self.params["beta"]
        order = []
        used = np.zeros(self.N, dtype=bool)
        ti_last = self.t_start
        u = 0  # start anchor
        group_last = {}

        while True:
            cands, weights = [], []
            m = self._window_of(ti_last)
            for j in self.neighbors[u]:
                if used[j]:
                    continue
                slew = self._arc_cost(u, j, m)
                arrive = max(ti_last + slew, self.t_early[j])
                comp = arrive + self.t_visit[j]
                if comp > min(self.t_late[j], self.dur_min) + 1e-9:
                    continue
                if not self._group_ready(j, comp, group_last):
                    continue
                eta = self.priority[j] / max(slew + self.t_visit[j], 1e-6)
                # discourage long waits (Verbeeck eq. 11)
                wait = max(0.0, self.t_early[j] - (ti_last + slew))
                li = max(1e-6, 1.0 - wait / max(self.dur_min, 1e-6))
                weights.append(li * (tau[u, j] ** alpha) * (eta ** beta))
                cands.append((j, comp))
            if not cands:
                break
            probs = np.asarray(weights, dtype=float)
            probs /= probs.sum()
            pick = int(self.rng.choice(len(cands), p=probs))
            j, comp = cands[pick]
            order.append(j)
            used[j] = True
            self._group_ready(j, comp, group_last, commit=True)
            tau[u, j] *= (1.0 - self.params["rho"])  # local evaporation
            ti_last = comp
            u = j
        return order

    def _group_ready(self, nid, comp, group_last, *, commit=False):
        """Like :meth:`_check_group` but non-mutating unless ``commit``."""
        info = self._group_of.get(nid)
        if info is None:
            return True
        members, pos = info
        if pos > 0:
            prev = members[pos - 1]
            if prev not in group_last:
                return False
            if comp < group_last[prev] + self.tau_intra[nid] - 1e-9:
                return False
        if commit:
            group_last[nid] = comp
        return True

    # ------------------------------------------------------------ local search
    def _expired(self):
        """True once the solve deadline has passed (guards hot loops)."""
        return self._deadline is not None and time.time() > self._deadline

    def _profile(self, order):
        """Forward profile + Verbeeck-style ``max_shift`` slack for a tour.

        ``order`` must be a *resolved* (realizable) sequence. Returns per-position
        completion times, the wait absorbed before each node, and ``max_shift``
        (how far each node may be delayed before some downstream node or the end
        of night becomes infeasible). With ``n_states > 1`` slews use the
        min-over-states cost, so ``max_shift`` is an optimistic lower bound and
        the filter never rejects a feasible insertion (exact check is _resolve).
        """
        L = len(order)
        comp = [0.0] * L
        wait = [0.0] * L
        prev, prev_state, last_comp = 0, 0, self.t_start
        for k, v in enumerate(order):
            dep = self.t_start if prev == 0 else last_comp
            m = self._window_of(dep)
            best = None
            for sv in (self.node_states.get(v, [0]) if self.n_states > 1 else [0]):
                if prev == 0:
                    slew = 0.0
                elif self.n_states > 1:
                    val = self.arcs_lookup.get((prev, v, m, prev_state, sv))
                    if val is None:
                        continue
                    slew = val
                else:
                    slew = self._arc_cost(prev, v, m)
                if best is None or slew < best[1]:
                    best = (sv, slew)
            sv, slew = best if best is not None else (0, 0.0)
            arrive = max(dep + slew, self.t_early[v])
            comp[k] = arrive + self.t_visit[v]
            wait[k] = max(0.0, self.t_early[v] - (dep + slew))
            prev, prev_state, last_comp = v, sv, comp[k]

        ub = [min(self.t_late[v], self.dur_min) for v in order]
        max_shift = [0.0] * L
        if L:
            max_shift[L - 1] = ub[L - 1] - comp[L - 1]
            for k in range(L - 2, -1, -1):
                max_shift[k] = min(ub[k] - comp[k], wait[k + 1] + max_shift[k + 1])
        return {"comp": comp, "wait": wait, "max_shift": max_shift, "ub": ub}

    def _local_search(self, order):
        """Insert / swap / replace moves until no feasible improvement.

        Operates on the *resolved* order (so candidate positions index a
        realizable tour) and uses neighbor-restricted candidate generation plus
        the ``max_shift`` insertion filter to keep each pass cheap.
        """
        res = self._resolve(order)
        if not res["feasible"]:
            return list(order), -np.inf
        best_order, best_obj = res["order"], res["objective"]
        improved = True
        while improved and not self._expired():
            improved = False
            order2, obj2 = self._try_insert(best_order, best_obj)
            if obj2 > best_obj + 1e-9:
                best_order, best_obj, improved = order2, obj2, True
                continue
            order2, obj2 = self._try_swap(best_order, best_obj)
            if obj2 > best_obj + 1e-9:
                best_order, best_obj, improved = order2, obj2, True
                continue
            order2, obj2 = self._try_replace(best_order, best_obj)
            if obj2 > best_obj + 1e-9:
                best_order, best_obj, improved = order2, obj2, True
        return best_order, best_obj

    def _insert_positions(self, order, v):
        """Candidate insertion indices for ``v`` (neighbor-restricted).

        A position is allowed at the tour boundaries, or where the preceding
        node can reach ``v`` (``prev in rev_neighbors[v]``) or ``v`` can reach
        the following node (``next in neighbors[v]``).
        """
        L = len(order)
        nbset = set(self.neighbors.get(v, []))
        rev = self.rev_neighbors.get(v, set())
        positions = []
        for pos in range(L + 1):
            prev = order[pos - 1] if pos > 0 else 0
            nxt = order[pos] if pos < L else self.N - 1
            if prev == 0 or nxt == self.N - 1 or prev in rev or nxt in nbset:
                positions.append(pos)
        return positions

    def _insert_passes_filter(self, order, prof, v, pos):
        """O(1) ``max_shift`` feasibility pre-screen for inserting ``v`` at pos."""
        L = len(order)
        prev = order[pos - 1] if pos > 0 else 0
        dep = self.t_start if pos == 0 else prof["comp"][pos - 1]
        m = self._window_of(dep)
        slew_in = 0.0 if prev == 0 else self._arc_cost(prev, v, m)
        arrive = max(dep + slew_in, self.t_early[v])
        comp_v = arrive + self.t_visit[v]
        if comp_v > min(self.t_late[v], self.dur_min) + 1e-9:
            return False
        if pos == L:
            return comp_v <= self.dur_min + 1e-9
        nxt = order[pos]
        slew_out = self._arc_cost(v, nxt, self._window_of(comp_v))
        old_slew = 0.0 if prev == 0 else self._arc_cost(prev, nxt, m)
        shift = (comp_v + slew_out) - (dep + old_slew)
        slack = prof["wait"][pos] + prof["max_shift"][pos]
        return shift <= slack + 1e-9

    def _try_insert(self, order, base_obj):
        """Best-improvement insertion, neighbor-restricted + max_shift filtered."""
        in_sol = set(order)
        prof = self._profile(order)
        best_order, best_obj = order, base_obj
        for v in self.real_nodes:
            if v in in_sol:
                continue
            if self._expired():
                break
            for pos in self._insert_positions(order, v):
                if not self._insert_passes_filter(order, prof, v, pos):
                    continue
                cand = order[:pos] + [v] + order[pos:]
                ok, obj, _ = self._evaluate(cand)
                if ok and obj > best_obj + 1e-9:
                    best_order, best_obj = cand, obj
        return best_order, best_obj

    def _try_swap(self, order, base_obj):
        """First-improvement swap of nearby scheduled positions (window)."""
        n = len(order)
        win = int(self.params["swap_window"])
        for a in range(n):
            if self._expired():
                break
            for b in range(a + 1, min(a + 1 + win, n)):
                cand = list(order)
                cand[a], cand[b] = cand[b], cand[a]
                ok, obj, _ = self._evaluate(cand)
                if ok and obj > base_obj + 1e-9:
                    return cand, obj
        return order, base_obj

    def _try_replace(self, order, base_obj):
        """Best-improvement exchange of a node for an unused neighbor of prev."""
        in_sol = set(order)
        best_order, best_obj = order, base_obj
        for pos in range(len(order)):
            if self._expired():
                break
            prev = order[pos - 1] if pos > 0 else 0
            cand_v = self.neighbors.get(prev, []) if prev != 0 else self.real_nodes
            for v in cand_v:
                if v in in_sol:
                    continue
                cand = list(order)
                cand[pos] = v
                ok, obj, _ = self._evaluate(cand)
                if ok and obj > best_obj + 1e-9:
                    best_order, best_obj = cand, obj
        return best_order, best_obj

    # -------------------------------------------------------------- main loop
    def solve_multistart(self, n_starts=1, *, base_seed=0, parallel=True):
        """Run ``n_starts`` independent ACS searches; return the best tour.

        ACS is stochastic, so several short independent runs from different
        seeds typically beat a single long run and parallelize trivially. Each
        start uses the full ``time_limit_s`` budget; with ``parallel=True`` they
        run in separate processes so wall time stays ~``time_limit_s``.
        """
        if n_starts <= 1:
            return self.solve()
        seeds = [int(base_seed) + 7919 * i for i in range(n_starts)]
        results = None
        if parallel:
            try:
                from concurrent.futures import ProcessPoolExecutor

                with ProcessPoolExecutor(max_workers=n_starts) as ex:
                    results = list(ex.map(_solve_with_seed, [self] * n_starts, seeds))
            except Exception as exc:  # fall back to sequential on any pool error
                logs.warning("ACS parallel multistart failed (%s); sequential", exc)
                results = None
        if results is None:
            results = []
            for s in seeds:
                self.rng = np.random.default_rng(s)
                results.append(self.solve())

        def _key(r):
            return r["objective"] if r.get("order") else -np.inf

        best = max(results, key=_key)
        logs.info(
            "ACS multistart: %d starts, best obj=%.2f (all: %s)",
            n_starts,
            _key(best),
            ", ".join(f"{_key(r):.1f}" for r in results),
        )
        return best

    def solve(self):
        """Run the ACS and return the best tour found.

        Returns a dict with ``order`` (list of real node ids), ``ti``
        (completion minutes), ``objective``, and ``feasible``.
        """
        tau = np.full((self.N, self.N), self.params["tau_init"], dtype=float)
        best_order, best_obj = [], -np.inf
        no_improve = 0
        t0 = time.time()
        self._deadline = t0 + float(self.params["time_limit_s"])

        if not self.real_nodes:
            return {"order": [], "ti": [], "objective": 0.0, "feasible": True}

        for it in range(int(self.params["max_iter"])):
            iter_best_order, iter_best_obj = None, -np.inf
            for _ in range(int(self.params["max_ants"])):
                order = self._construct(tau)
                order, obj = self._local_search(order)
                if obj > iter_best_obj:
                    iter_best_order, iter_best_obj = order, obj
                if self._expired():
                    break
            if iter_best_obj > best_obj + 1e-9:
                best_order, best_obj = iter_best_order, iter_best_obj
                no_improve = 0
            else:
                no_improve += 1
            # global pheromone reinforcement on the iteration best.
            if iter_best_order:
                rg = self.params["rho_global"]
                prev = 0
                for nid in iter_best_order:
                    tau[prev, nid] = (1 - rg) * tau[prev, nid] + rg * max(
                        iter_best_obj, 0.0
                    )
                    prev = nid
            if no_improve >= int(self.params["reset_no_improve"]):
                tau[:] = self.params["tau_init"]
                no_improve = 0
            if time.time() - t0 > self.params["time_limit_s"]:
                break

        res = self._resolve(best_order)
        logs.info(
            "ACS: %d/%d visits, objective=%.2f, %.2fs",
            len(res["order"]),
            len(self.real_nodes),
            res["objective"],
            time.time() - t0,
        )
        res["feasible"] = bool(res["order"]) or not self.real_nodes
        return res
