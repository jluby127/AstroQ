"""Traveling Telescope Problem (TTP) solver

The model implements the MILP of Handley et al. 2024 (arXiv:2310.18497).
"""

# Standard library imports
import logging
import time

# Third-party imports
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.table import QTable
import gurobipy as gp
from gurobipy import GRB

logs = logging.getLogger(__name__)

class TTPModel:
    """MILP solver for the Traveling Telescope Problem (Handley+ 2024).

    Args:
        requests (astropy.table.QTable): one row per request, with native
            astropy-typed columns:

            unique_id        str                            primary key
            coord            SkyCoord column                ICRS
            first_available  Time column                    earliest start of accessibility window
            last_available   Time column                    latest end of accessibility window
            t_visit          Quantity column (time)         per-visit duration
            n_intra_max      int column                     max visits per night
            tau_intra        Quantity column (time)         min spacing between visits within a night
            weight           float column                   objective weight; higher = more important

        night_start (astropy.time.Time): start of the observing interval.
        night_end (astropy.time.Time): end of the observing interval.

    Keyword Args:
        slew_fn (callable): ``slew_fn(coord_a, coord_b, window_start,
            window_end) -> np.ndarray`` of worst-case slew minutes per
            (pair, window). ``coord_a`` and ``coord_b`` are pair-aligned
            1-D ``SkyCoord`` arrays of length ``P``; ``window_start`` and
            ``window_end`` are 1-D ``Time`` arrays of length ``M`` giving
            the bounds of each slew slot. The return has shape ``(P, M)``
            and its ``[k, m]`` entry is the slew time from ``coord_a[k]``
            to ``coord_b[k]`` over ``[window_start[m], window_end[m]]``.
            The implementation owns the per-window sampling cadence and
            reduction (max, mean, percentile, ...).
        n_slots (int): number of TTP slew slots ``M`` (Handley+ 2024 §2.2).
            ``n_slots=1`` is the recommended default.

    The input QTable ``self.requests`` is the single source of truth for
    the per-request inputs (native astropy / numpy types). Downstream
    relational state lives in three DataFrames built by the solver:

    * ``self.nodes`` -- one row per MILP node (:meth:`build_nodes`).
      Scalar-only dtypes (no object columns) so it is round-trip-safe
      through ``to_hdf``/``to_csv``.
    * ``self.arcs`` -- precomputed arc catalog ``(i, j, m)`` (:meth:`build_arcs`).
    * ``self.schedule`` -- post-solve output parallel to ``nodes``
      (:meth:`build_schedule`); ``None`` if Gurobi finds no incumbent
      within the time limit.

    Notes:
        Internal naming is aligned with Handley+ 2024 (``N``, ``M``, ``Yi``,
        ``Xijm``, ``arcs``) and AstroQ vocabulary (``t_visit``, ``tau_intra``).

    Usage:
        tm = TTPModel(...)
        tm.build_nodes() # builds nodes from requests
        tm.build_arcs()  # builds arcs from nodes
        tm.build_model() # builds the MILP

        # set the Gurobi parameters
        tm.model.params.TimeLimit = 300  # set the time limit
        tm.model.params.MIPGap = 0.05  # set the MIP gap
        tm.model.params.OutputFlag = 0 # set the output flag
        tm.model.update() # update the model

        tm.run_model()      # run the model
        tm.build_schedule() # builds the schedule
        tm.to_string()      # summarizes the solve for logging.
    """

    #: Slew tie-breaker weight in the objective. Same constant is used to
    #: recover total slew time inside :meth:`build_schedule`.

    # Constant that balances slew time vs. number of targets.
    # Interpretation: if dropping the highest-weight target saves this many
    # minutes of slew time, drop it
    _SLEW_MINUTES_FOR_TOP_TARGET = 30
    _SLEW_IDLE_PENALTY_RATIO = 0.5  # idle-between weight, as fraction of slew penalty

    #: Per-column spec for the ``requests`` QTable. Value is either an astropy
    #: type the column must be an instance of, or a time-equivalent unit the
    #: column's ``.unit`` must be convertible to. ``None`` skips type checking.
    _COLUMN_SPECS = {
        "unique_id":       None,      # str, primary key
        "coord":           SkyCoord,  # ICRS
        "first_available": Time,
        "last_available":  Time,
        "t_visit":         u.s,       # Quantity with time units
        "n_intra_max":     None,      # int
        "tau_intra":       u.s,       # Quantity with time units
        "weight":          None,      # float
    }

    def __init__(
        self,
        requests,
        night_start,
        night_end,
        *,
        slew_fn,
        n_slots=1,
        n_states=1,
    ):
        if not isinstance(requests, QTable):
            raise TypeError(
                "`requests` must be an astropy.table.QTable; got "
                f"{type(requests).__name__}"
            )
        missing = [c for c in self._COLUMN_SPECS if c not in requests.colnames]
        if missing:
            raise ValueError(f"`requests` missing required columns: {missing}")
        for col, spec in self._COLUMN_SPECS.items():
            if spec is None or len(requests) == 0:
                continue
            column = requests[col]
            if isinstance(spec, u.UnitBase):
                unit = column.unit
                if unit is None or not unit.is_equivalent(spec):
                    raise TypeError(
                        f"`{col}` column must be a Quantity with units "
                        f"equivalent to {spec}; got unit={unit}"
                    )
            elif not isinstance(column, spec):
                raise TypeError(
                    f"`{col}` column must be an astropy {spec.__name__}; "
                    f"got {type(column).__name__}"
                )

        self.requests = requests
        self.night_start = night_start
        self.night_end = night_end

        if len(self.requests) > 0:
            dfa = self.requests["first_available"] - night_start
            if dfa.min().to_value(u.s) > 0:
                logs.warning(
                    "min(first_available) is {:.1f} after night_start".format(
                        dfa.min().to(u.min)
                    )
                )

        self.slew_fn = slew_fn
        self.n_slots = n_slots
        #: Number of cable-wrap states. ``1`` selects the legacy single-cut
        #: path; ``>1`` enables the state-aware MILP (``slew_fn`` must then
        #: return a ``(P, M, S, S)`` tensor with ``NaN`` for impossible arcs).
        self.n_states = int(n_states)
        self.S = self.n_states
        self.dur_min = float(
            np.round((night_end - night_start).to_value(u.min), 0)
        )
        self.M = self.n_slots
        self.schedule = None
        self.stats = {}

    def build_nodes(self):
        """Expand the request frame into the ``self.nodes`` DataFrame.

        Sets attributes:
            N, M, nodes, multi_visit_groups.

        ``self.nodes`` carries only scalar dtypes. All times are in minutes
        from ``night_start``; sky coordinates appear as ``ra`` / ``dec``
        degrees. The :class:`~astropy.coordinates.SkyCoord` array used for
        slew evaluation is read from ``self.requests["coord"]`` and indexed
        by ``request_idx``.
        """
        # Local scratch DataFrame: scalar-only projection of self.requests
        # for the relational joins below. Built fresh here so the QTable
        # remains the single source of truth on the model.
        r = self.requests
        rdf = self.requests["unique_id","n_intra_max","weight"].to_pandas()
        rdf["t_early"] = (r["first_available"] - self.night_start).to_value(u.min).astype(float)
        rdf["t_late"] =  (r["last_available"] - self.night_start).to_value(u.min).astype(float)
        rdf["tau_intra"] = r["tau_intra"].to_value(u.min).astype(float)
        rdf["t_visit"] = r["t_visit"].to_value(u.min).astype(float)
        rdf["ra"] = r["coord"].ra.deg.astype(float)
        rdf["dec"] = r["coord"].dec.deg.astype(float)
        rdf["request_idx"] = np.arange(len(r), dtype=np.int64)
        rdf["is_anchor"] = False

        # Attach visit_seq via a simple cross-join + filter.
        max_intra = int(rdf.n_intra_max.max()) if len(rdf) else 0
        visit_seq_table = pd.DataFrame(
            {
                "visit_seq": np.arange(max_intra, dtype=np.int64),
            }
        )
        visits = (
            rdf.merge(visit_seq_table, how="cross")
            .query("visit_seq < n_intra_max")
            .sort_values(["request_idx", "visit_seq"], kind="stable")
            .reset_index(drop=True)
        )

        anchor_template = {
            "unique_id": "",
            "request_idx": -1,
            "visit_seq": 0,
            "is_anchor": True,
            "t_early": 0.0,
            "t_late": self.dur_min,
            "t_visit": 0.0,
            "tau_intra": 0.0,
            "weight": 0.0,
            "n_intra_max": 0,
            "ra": np.nan,
            "dec": np.nan,
        }
        anchor_df = pd.DataFrame([anchor_template, anchor_template])

        self.nodes = pd.concat(
            [anchor_df.iloc[[0]], visits, anchor_df.iloc[[1]]],
            ignore_index=True,
        )
        self.N = len(self.nodes)

        # Multi-visit groups: keyed by unique_id, value is the list of node
        # ids (in visit_seq order) for requests with n_intra_max > 1.
        gb = self.nodes[~self.nodes["is_anchor"]].groupby("unique_id", sort=False)
        sizes = gb.size()
        self.multi_visit_groups = {
            uid: gb.get_group(uid).sort_values("visit_seq").index.tolist()
            for uid in sizes.index[sizes > 1]
        }

    def build_arcs(self):
        """Build the precomputed, state-indexed arc catalog ``self.arcs``.

        Splits ``[night_start, night_end]`` into ``M`` equal-length windows
        and asks ``self.slew_fn`` for the worst-case slew (minutes) on each
        ordered pair of real nodes within each window. The per-window
        sampling cadence and reduction policy live inside ``slew_fn``; see
        the class docstring for the callable contract.

        The catalog is *always* keyed ``(i, j, m, si, sj)`` (departure window
        ``m``, wrap states ``si``/``sj``). A legacy single-state ``slew_fn`` may
        return a 2-D ``(P, M)`` array; it is promoted to ``(P, M, 1, 1)`` so the
        single-state model is just the ``S = 1`` degenerate case (one state,
        ``si = sj = 0``). Arcs with a non-finite slew (an unreachable winding)
        are dropped, and ``self.node_states`` records the wrap states each real
        node can actually be observed in (``{n: [0]}`` when ``S == 1``).
        """
        if not hasattr(self, "nodes"):
            raise RuntimeError("call build_nodes() before build_arcs()")

        # Equal-length window bounds across the night.
        fractions = np.linspace(0.0, 1.0, self.M + 1)
        slot_bounds = self.night_start + (self.night_end - self.night_start) * fractions
        window_start = slot_bounds[:-1]
        window_end = slot_bounds[1:]

        # Index the per-request SkyCoord array by the real-node request_idx.
        non_anchor = self.nodes[~self.nodes.is_anchor]
        node_coords = self.requests["coord"][non_anchor["request_idx"].to_numpy()]
        n_real = len(node_coords)

        # Pair list: ordered (i, j) with i != j over real-node positions
        # (0-based). `id` in self.nodes for real nodes is `position + 1`.
        ii, jj = np.meshgrid(np.arange(n_real), np.arange(n_real), indexing="ij")
        mask = ii != jj
        i_pos = ii[mask]  # starting node of arc
        j_pos = jj[mask]  # ending node of arc
        n_pair = len(i_pos)

        # slew_fn returns (n_pair, M, S, S) with NaN where a winding is
        # unreachable; a legacy single-state slew_fn returns (n_pair, M), which
        # we promote to (n_pair, M, 1, 1). Flatten, then drop the NaN
        # (impossible) arcs so they can never be selected.
        tau = np.asarray(
            self.slew_fn(
                node_coords[i_pos], node_coords[j_pos], window_start, window_end
            )
        )
        if tau.ndim == 2:
            tau = tau[:, :, None, None]
        if tau.shape != (n_pair, self.M, self.S, self.S):
            raise ValueError(
                f"slew_fn returned shape {tau.shape}; expected "
                f"{(n_pair, self.M, self.S, self.S)} for n_states={self.S}"
            )
        p_g, m_g, si_g, sj_g = np.mgrid[
            0:n_pair, 0 : self.M, 0 : self.S, 0 : self.S
        ]
        i_id = i_pos[p_g.ravel()] + 1
        j_id = j_pos[p_g.ravel()] + 1
        m_lev = m_g.ravel()
        si = si_g.ravel()
        sj = sj_g.ravel()
        t_slew = tau.ravel()
        keep = np.isfinite(t_slew)
        self.arcs = pd.DataFrame(
            {
                "t_slew": t_slew[keep],
            },
            index=pd.MultiIndex.from_arrays(
                [i_id[keep], j_id[keep], m_lev[keep], si[keep], sj[keep]],
                names=("i", "j", "m", "si", "sj"),
            ),
        )
        # Feasible states per real node: states in which the node appears as
        # either endpoint of a surviving arc (``[0]`` for every node if S == 1).
        node_states = {n: set() for n in range(1, self.N - 1)}
        for n, s in zip(i_id[keep], si[keep]):
            node_states[n].add(int(s))
        for n, s in zip(j_id[keep], sj[keep]):
            node_states[n].add(int(s))
        self.node_states = {n: sorted(s) for n, s in node_states.items()}

        # Slot bounds (minutes from night_start) used by build_model.
        self.w = self.dur_min * fractions

    # ------------------------------------------------------------- MILP build
    def build_model(self):
        """Construct the (state-indexed) TTP MILP.

        This is the Handley+ 2024 formulation (eqs. 2-9, B3, 10) generalized
        with a cable-wrap state index. The single-state model is the degenerate
        ``S == 1`` case (every node has one state ``0``), so there is a single
        code path for both.

        Paper map (Handley+ 2024 -> this code):
          * ``X_{ijm}`` (arc i->j departing in window m) -> ``Xijm[(i, j, m, si,
            sj)]`` with the extra wrap states ``si`` (at i) and ``sj`` (at j).
          * ``Y_i`` (node i visited) -> ``Yi[(i, s)]`` (node i visited in wrap
            state s). ``one_state`` caps ``sum_s Yi[(i, s)] <= 1``.
          * eqs. 2/3 (one arc out of start / into end), eq. 4 (visit once),
            eq. 5 (flow conservation), eq. 6 (``ti`` <-> ``tijm``), eq. 7
            (exposure/slew timing), eq. 8 (slot bounds), eq. 9 (rise/set),
            eq. B3 (intra-night separation), eq. 10 (objective) all carry their
            paper labels in the comments below.

        Divergence from Handley+ 2024:
          * Wrap states: arcs/visits are state-indexed and flow conservation is
            enforced per ``(node, state)`` so the telescope holds its winding
            through a visit (state continuity). Impossible windings were dropped
            (NaN) in :meth:`build_arcs`, so only feasible arcs exist. Anchors
            (start ``0`` / end ``N-1``) are stateless and their arcs cost 0.
          * ``Yi`` is left continuous in ``[0, 1]``: ``visit_once`` pins it to a
            sum of binary in-arcs (and ``one_state`` caps it), so it is 0/1 at
            every integer-feasible point. This removes redundant binaries and
            attacks the wrap-state branching symmetry.
          * ``first_exposure`` pins the first visit to the night start (or the
            earliest feasible time), preventing leading idle in underfilled
            schedules (not in Handley+ 2024).

        Variables are built over the sparse set of feasible arc keys via
        adjacency tables (``out_by_node`` / ``in_by_ns`` / ...) rather than dense
        ``range(N) x range(N) x range(M)`` loops, which keeps the model small
        when many windings are infeasible.
        """
        if not hasattr(self, "arcs"):
            raise RuntimeError("call build_arcs() before build_model()")
        self.model = gp.Model("TTP")
        N, M = self.N, self.M
        nodes = self.nodes
        real_nodes = range(1, N - 1)

        # Per-node feasible states (computed in build_arcs from surviving arcs).
        node_states = self.node_states

        # Slew cost for the internal (real-real) arcs that survived NaN drop.
        arc_cost = self.arcs["t_slew"].to_dict()  # keys (i, j, m, si, sj)

        # ---- Enumerate every arc key used by the model + adjacency tables.
        arc_keys = list(arc_cost.keys())          # internal feasible arcs
        cost = dict(arc_cost)                      # anchor arcs default to 0
        # start-anchor -> real node-state (one per window); end real -> end anchor.
        for j in real_nodes:
            for sj in node_states[j]:
                for m in range(M):
                    key = (0, j, m, 0, sj)
                    arc_keys.append(key)
                    cost[key] = 0.0
        for i in real_nodes:
            for si in node_states[i]:
                for m in range(M):
                    key = (i, N - 1, m, si, 0)
                    arc_keys.append(key)
                    cost[key] = 0.0

        out_by_node = {i: [] for i in range(N)}    # keyed by source node id
        in_by_node = {j: [] for j in range(N)}     # keyed by dest node id
        out_by_ns = {}                             # keyed by (node, state)
        in_by_ns = {}
        for key in arc_keys:
            i, j, m, si, sj = key
            out_by_node[i].append(key)
            in_by_node[j].append(key)
            out_by_ns.setdefault((i, si), []).append(key)
            in_by_ns.setdefault((j, sj), []).append(key)

        # Real node-states that actually exist (have a Y variable).
        ns_keys = [(i, s) for i in real_nodes for s in node_states[i]]

        # ---- Variables.
        # Yi is implied-integer: visit_once pins Yi[(j,s)] = sum of binary in-arcs
        # (<=1 via one_state), so it is 0/1 at every integer-feasible point. Leaving
        # it continuous keeps integrality on the arcs (Xijm) and removes ~|states|
        # redundant binaries from the branch set -- this attacks the wrap-state
        # branching symmetry directly. (build_schedule reads only Xijm/ti.)
        self.Yi = self.model.addVars(
            ns_keys, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name="Yi"
        )
        self.Xijm = self.model.addVars(arc_keys, vtype=GRB.BINARY, name="Xijm")
        self.tijm = self.model.addVars(arc_keys, vtype=GRB.CONTINUOUS, name="tijm")
        self.ti = self.model.addVars(range(N), vtype=GRB.CONTINUOUS, lb=0, name="ti")

        self.model.addConstr(self.ti[0] == 0.0, "anchor_start_time")

        # First exposure pinned to start (or earliest feasible) -- same policy
        # as the single-state model, now per start-anchor arc.
        t_earliest = nodes[~nodes.is_anchor].t_early.min()
        t_start = max(0.0, t_earliest)
        for key in out_by_node[0]:
            _, j, _, _, _ = key
            t_visit_j = float(nodes.at[j, "t_visit"])
            self.model.addGenConstrIndicator(
                self.Xijm[key], 1, self.ti[j], GRB.EQUAL,
                t_start + t_visit_j, name=f"first_exposure_{key}",
            )

        # eq. 2 - exactly one arc out of the start anchor.
        self.model.addConstr(
            gp.quicksum(self.Xijm[k] for k in out_by_node[0]) == 1, "start_anchor"
        )
        # eq. 3 - exactly one arc into the end anchor.
        self.model.addConstr(
            gp.quicksum(self.Xijm[k] for k in in_by_node[N - 1]) == 1, "end_anchor"
        )

        # eq. 4 - visit indicator per real node-state; at most one state used.
        for j in real_nodes:
            for sj in node_states[j]:
                self.model.addConstr(
                    gp.quicksum(self.Xijm[k] for k in in_by_ns.get((j, sj), []))
                    == self.Yi[(j, sj)],
                    f"visit_once_{j}_{sj}",
                )
            self.model.addConstr(
                gp.quicksum(self.Yi[(j, s)] for s in node_states[j]) <= 1,
                f"one_state_{j}",
            )

        # eq. 5 - flow conservation per (real node, state): in == out. Keying by
        # the same state on both sides enforces wrap-state continuity.
        for k in real_nodes:
            for sk in node_states[k]:
                self.model.addConstr(
                    gp.quicksum(self.Xijm[a] for a in in_by_ns.get((k, sk), []))
                    - gp.quicksum(self.Xijm[a] for a in out_by_ns.get((k, sk), []))
                    == 0,
                    f"flow_{k}_{sk}",
                )

        # eq. 6 - link ti to tijm (per non-end node; sum over its out-arcs).
        for i in range(N - 1):
            self.model.addConstr(
                self.ti[i] == gp.quicksum(self.tijm[k] for k in out_by_node[i]),
                f"tijm_def_{i}",
            )

        # eq. 7 - exposure/slew time linking (per non-start node).
        for j in range(1, N):
            t_visit_j = float(nodes.at[j, "t_visit"])
            self.model.addConstr(
                self.ti[j]
                >= gp.quicksum(
                    self.tijm[k] + (cost[k] + t_visit_j) * self.Xijm[k]
                    for k in in_by_node[j]
                ),
                f"exp_constr_{j}",
            )

        # eq. 8 - slot bounds on tijm (per arc).
        for k in arc_keys:
            m = k[2]
            self.model.addConstr(
                self.tijm[k] >= self.w[m] * self.Xijm[k], f"t_min_{k}"
            )
            self.model.addConstr(
                self.tijm[k] <= self.w[m + 1] * self.Xijm[k], f"t_max_{k}"
            )

        # eq. 9 - node accessibility (per real node; visited = sum over states).
        for i in real_nodes:
            row = nodes.loc[i]
            visited = gp.quicksum(self.Yi[(i, s)] for s in node_states[i])
            self.model.addConstr(
                self.ti[i] >= (row.t_early + row.t_visit) * visited,
                f"rise_constr_{i}",
            )
            self.model.addConstr(
                self.ti[i] <= row.t_late * visited, f"set_constr_{i}"
            )

        # eq. B3 - intra-night separation (multi-visit only).
        for indices in self.multi_visit_groups.values():
            for k in range(1, len(indices)):
                cur, prev = indices[k], indices[k - 1]
                self.model.addConstr(
                    gp.quicksum(self.tijm[a] for a in out_by_node[cur])
                    >= gp.quicksum(self.tijm[a] for a in out_by_node[prev])
                    + gp.quicksum(self.Yi[(cur, s)] for s in node_states[cur])
                    * nodes.at[cur, "tau_intra"],
                    f"intra_sep_constr_{prev}_{cur}",
                )

        # Total slew over internal (real-real) arcs only.
        self.t_slew = self.model.addVar(lb=0.0, name="t_slew")
        self.model.addConstr(
            self.t_slew
            == gp.quicksum(arc_cost[k] * self.Xijm[k] for k in arc_cost),
            "t_slew_def",
        )

        self.t_visit = self.model.addVar(lb=0.0, name="t_visit")
        self.model.addConstr(
            self.t_visit
            == gp.quicksum(
                nodes.at[j, "t_visit"] * self.Yi[(j, s)]
                for j in real_nodes
                for s in node_states[j]
            ),
            "t_visit_def",
        )

        self.t_idle_between = self.model.addVar(lb=0.0, name="t_idle_between")
        self.model.addConstr(
            self.t_idle_between == self.ti[N - 1] - self.t_visit - self.t_slew,
            name="t_idle_between_def",
        )

        # eq. 10 - objective.
        W_max = float(self.nodes.loc[1 : N - 1, "weight"].max())
        slew_penalty = W_max / self._SLEW_MINUTES_FOR_TOP_TARGET
        self.model.setObjective(
            gp.quicksum(
                nodes.at[j, "weight"] * self.Yi[(j, s)]
                for j in real_nodes
                for s in node_states[j]
            )
            - slew_penalty * self.t_slew
            - slew_penalty * self._SLEW_IDLE_PENALTY_RATIO * self.t_idle_between,
            GRB.MAXIMIZE,
        )
        self.model.update()

    def run_model(self):
        """Solve the MILP and build ``schedule`` / ``stats``.

        On success, ``schedule`` is a DataFrame and ``stats`` is populated.
        If Gurobi has no incumbent (``SolCount == 0``), logs a warning and
        leaves ``schedule`` as ``None``.
        """
        if not hasattr(self, "model"):
            raise RuntimeError("call build_model() before run_model()")
        logs.info(f"Solving TTP for {self.N - 2} visits")
        t0 = time.time()
        self.model.optimize()
        if self.model.Status == GRB.INFEASIBLE:
            logs.critical("TTP infeasible; computing IIS.")
            self.model.computeIIS()
            for c in self.model.getConstrs():
                if c.IISConstr:
                    logs.critical(c.ConstrName)

        logs.info(f"TTP solve finished in {time.time() - t0:.3f}s")
        if self.model.SolCount == 0:
            logs.warning(
                "No incumbent TTP solution within time limit. "
                "Try raising ``self.model.params.TimeLimit`` before solving."
            )

    # ---------------------------------------------------------- post-process

    def build_schedule(self):
        """Walk the chosen path and populate ``schedule``.

        ``schedule`` is parallel to ``nodes`` (same index = node id ``i``),
        with ``scheduled``, ``t_slew``, and solve-time columns added.
        """
        if not hasattr(self, "model"):
            raise RuntimeError("call build_model() before build_schedule()")

        arcs_selected = []
        for key, var in self.Xijm.items():
            i, j, m, si, sj = key
            if var.X > 0.5 and j != 0 and i != self.N - 1:
                arcs_selected.append(
                    {"i": i, "j": j, "m": m, "si": si, "sj": sj, "ti": self.ti[i].X}
                )
        arcs_selected = pd.DataFrame(arcs_selected)
        self._finalize_schedule(arcs_selected)

    def _finalize_schedule(self, arcs_selected):
        """Assemble ``self.schedule`` / ``self.stats`` from selected arcs.

        ``arcs_selected`` has one row per scheduled real node (its outgoing arc)
        with columns ``i, j, m, si, sj, ti``. ``wrap_state`` is always present
        (``0`` for every node when ``S == 1``).
        """
        # Guarantee the merge keys exist even when nothing was scheduled.
        cols = ["i", "j", "m", "si", "sj", "ti"]
        if arcs_selected is None or len(arcs_selected) == 0:
            arcs_selected = pd.DataFrame(columns=cols)

        # merge selected arcs with nodes, unvisited nodes will have NaN for ti
        schedule = pd.merge(
            self.nodes.query("~is_anchor"),
            arcs_selected,
            left_index=True,
            right_on=["i"],
            how="left",
        )
        schedule["t_start"] = schedule["ti"] - schedule["t_visit"]
        schedule["t_end"] = schedule["ti"]
        schedule["scheduled"] = ~schedule["ti"].isna()
        # wrap_state = the winding node ``i`` is observed in (0 when S == 1).
        schedule["wrap_state"] = schedule["si"]

        merge_on = ["i", "j", "m", "si", "sj"]
        schedule = pd.merge(
            schedule,
            self.arcs["t_slew"],
            left_on=merge_on,
            right_index=True,
            how="left",
        ).sort_values(by="t_start", na_position="last")
        schedule["order"] = range(len(schedule))

        # self.nodes has only scalar columns (ra/dec as floats, no SkyCoord);
        # schedule inherits that and is round-trip-safe through to_csv/to_hdf.
        self.schedule = schedule
        scheduled = self.schedule[self.schedule["scheduled"]]
        stats = {
            "dur_min": self.dur_min,
            "n_requested": self.N - 2,
            "n_scheduled": len(scheduled),
            "t_first_start": scheduled["t_start"].min(),
            "t_last_end": scheduled["t_end"].max(),
            "t_visit_sum": scheduled["t_visit"].sum(),
            "t_slew_sum": scheduled["t_slew"].sum(),
            "t_idle_sum": self.dur_min
            - scheduled["t_visit"].sum()
            - scheduled["t_slew"].sum(),
        }
        stats["t_idle_after_last"] = self.dur_min - stats["t_last_end"]
        stats["t_idle_before_last"] = stats["t_idle_sum"] - stats["t_idle_after_last"]
        self.stats = stats

    # ------------------------------------------------------- ACS warm start
    def _window_of(self, t_depart):
        """Slot index ``m`` (in ``0..M-1``) containing departure minute ``t``."""
        m = int(np.searchsorted(self.w, t_depart, side="right") - 1)
        return min(max(m, 0), self.M - 1)

    def _tour_arcs_dataframe(self, result):
        """Build an ``arcs_selected`` frame from an ACS/heuristic tour dict."""
        order = result.get("order", [])
        if not order:
            return pd.DataFrame(columns=["i", "j", "m", "si", "sj", "ti"])
        states = result.get("states") or {n: 0 for n in order}
        ti = list(result["ti"])
        seq = [0, *order, self.N - 1]
        ti_full = [0.0, *ti, float(ti[-1])]
        rows = []
        for k in range(len(seq) - 2):
            i = seq[k + 1]
            j = seq[k + 2]
            m = self._window_of(ti_full[k])
            si = states.get(i, 0)
            sj = 0 if j == self.N - 1 else states.get(j, 0)
            rows.append(
                {"i": i, "j": j, "m": m, "si": si, "sj": sj, "ti": ti_full[k + 1]}
            )
        return pd.DataFrame(rows)

    def schedule_from_tour(self, result):
        """Populate ``schedule`` / ``stats`` from a tour dict without Gurobi.

        Requires ``build_nodes`` and ``build_arcs`` (not ``build_model``).
        """
        if not hasattr(self, "nodes"):
            raise RuntimeError("call build_nodes() before schedule_from_tour()")
        self._finalize_schedule(self._tour_arcs_dataframe(result))

    def seed_from_tour(self, result):
        """Set a Gurobi MIPStart from an ACS warm-start ``result``.

        Call after :meth:`build_model`. ``result`` is the dict returned by
        :func:`astroq.ttp.acs.acs_warm_start` (``order`` / ``states`` / ``ti`` /
        ``feasible``). Sets ``.Start`` on the tour's arc / visit variables so
        Gurobi begins from a strong incumbent; unset variables are completed by
        Gurobi. Safe no-op if the tour is empty or infeasible.
        """
        if not hasattr(self, "model"):
            raise RuntimeError("call build_model() before seed_from_tour()")
        order = result.get("order", [])
        if not order or not result.get("feasible", False):
            return
        states = result.get("states") or {n: 0 for n in order}
        ti = list(result["ti"])
        seq = [0, *order, self.N - 1]
        ti_full = [0.0, *ti, float(ti[-1])]

        # Clear any stale starts, then set the tour's arcs / visits to 1.
        for var in self.Xijm.values():
            var.Start = 0.0
        for var in self.Yi.values():
            var.Start = 0.0

        for k in range(len(seq) - 1):
            i, j = seq[k], seq[k + 1]
            m = self._window_of(ti_full[k])
            si = 0 if i == 0 else states.get(i, 0)
            sj = 0 if j == self.N - 1 else states.get(j, 0)
            key = (i, j, m, si, sj)
            if key in self.Xijm:
                self.Xijm[key].Start = 1.0
        for nid in order:
            ykey = (nid, states.get(nid, 0))
            if ykey in self.Yi:
                self.Yi[ykey].Start = 1.0
        self.model.update()

    def to_string(self, *, header="Stats for TTP Solution"):
        """Return a human-readable summary of the solve from ``self.stats``."""
        if not self.stats:
            raise RuntimeError("call build_schedule() before to_string()")
        s = self.stats
        rows = [
            ("Observations Requested:", s["n_requested"], "d"),
            ("Observations Scheduled:", s["n_scheduled"], "d"),
            ("Observing Duration (min):", s["dur_min"], ".1f"),
            ("First Exposure Start (min):", s["t_first_start"], ".1f"),
            ("Last Exposure End (min):", s["t_last_end"], ".1f"),
            ("Visit Time (min):", s["t_visit_sum"], ".1f"),
            ("Slew Time (min):", s["t_slew_sum"], ".1f"),
            ("Idle Time (min):", s["t_idle_sum"], ".1f"),
            ("Idle After Last (min):", s["t_idle_after_last"], ".1f"),
            ("Idle Before Last (min):", s["t_idle_before_last"], ".1f"),
        ]
        label_w = max(len(label) for label, _, _ in rows)
        value_w = 7
        divider = "-" * (2 + label_w + 1 + value_w)
        fmt = f"  {{:<{label_w}}} {{:>{value_w}{{spec}}}}"

        lines = [header, divider]
        lines.extend(fmt.format(label, value, spec=spec) for label, value, spec in rows)
        lines.append(divider)
        return "\n".join(lines) + "\n"
