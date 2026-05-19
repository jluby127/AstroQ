"""Traveling Telescope Problem (TTP) solver

The model implements the MILP of Handley et al. 2024 (arXiv:2310.18497).
"""

# Standard library imports
import logging
import time
from collections import defaultdict
from itertools import permutations, product

# Third-party imports
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time, TimeDelta
import gurobipy as gp
from gurobipy import GRB

logs = logging.getLogger(__name__)

REQUIRED_COLUMNS = [
    "unique_id", "ra", "dec",
    "exptime",
    "n_exp",
    "n_intra_max",
    "tau_intra",
    "priority",
    "first_available",
    "last_available",
]

class TTPModel:
    """MILP solver for the Traveling Telescope Problem (Handley+ 2024).

    Args:
        requests_frame (pd.DataFrame): one row per request, required columns

            unique_id        str            primary key
            ra, dec          float, deg
            exptime          float, seconds (single shot exposure time)
            n_exp            int            shots per visit (multi-shot exposures)
            n_intra_max      int            visits per night
            tau_intra        float, hours   minimum spacing between visits within a night
            priority         int            objective weight; higher = more important
            first_available  str            ISO-8601 (caller computes; e.g. via Access)
            last_available   str            ISO-8601

        night_start (astropy.time.Time): start of the observing interval.
        night_end (astropy.time.Time): end of the observing interval.
        outdir (str): directory in which to write ``TTPstatistics.txt``.

    Keyword Args:
        observer (astroplan.Observer): site fixture for alt/az lookups. Used by
            :meth:`_compute_tau_slew` to grid the slew tensor and by
            :meth:`digest_gurobi` to evaluate alt/az at solver-chosen times.
        slew_rate (float): mean slew rate, degrees/second.
        wrap_limit (float | None): azimuth wrap limit, degrees. ``None`` means
            no wrap (any az difference > 180 wraps the short way).
        readout_time (float): detector readout between shots of a visit, sec.
        n_slots (int): number of TTP slew slots ``M`` (Handley+ 2024 §2.2).
            ``n_slots=1`` is the recommended default.
        inaccessible_zones (list[tuple]): retained on the model only as a
            convenience for the plotter (``astroq.ttp.plot``); the solver
            does not use them.
        runtime (float): Gurobi time limit, seconds.
        optgap (float): Gurobi MIP gap.
        slew_sample_cadence_min (int): max spacing in minutes at which to
            sample the slew tensor within a slot.
        output_flag (bool): pass-through to ``gurobipy.Model.Params.OutputFlag``.

    The class is structured around three DataFrames:

    * ``self.requests_frame`` -- one row per request
    * ``self.nodes`` -- one row per MILP node Built by :meth:`_build_nodes`
    * ``self.tau_slew`` -- per-arc slew tensor, indexed by ``(i, j, m)``.
       Built by :meth:`_compute_tau_slew`.

    Internal naming is aligned with the Handley 2024 paper (``N``, ``M``, ``Yi``,
    ``Xijm``, ``tau_slew``) and with AstroQ vocabulary elsewhere (``t_visit``,
    ``tau_intra``).
    """

    #: Slew tie-breaker weight in the objective. Same constant is used to
    #: recover total slew time inside :meth:`optimization_status`.
    _SLEW_PENALTY = 1 / 100

    def __init__(
        self,
        requests_frame,
        night_start,
        night_end,
        outdir,
        *,
        observer,
        slew_rate,
        wrap_limit=None,
        readout_time=0.0,
        n_slots=1,
        inaccessible_zones=None,
        runtime=300,
        optgap=0.01,
        slew_sample_cadence_min=30,
        output_flag=True,
    ):
        self.requests_frame = requests_frame.reset_index(drop=True).copy()
        self.night_start = night_start
        self.night_end = night_end
        self.outdir = outdir
        self.observer = observer
        self.slew_rate = slew_rate
        self.wrap_limit = wrap_limit
        self.readout_time = readout_time
        self.n_slots = n_slots
        self.inaccessible_zones = list(inaccessible_zones or [])
        self.runtime = runtime
        self.optgap = optgap
        self.slew_sample_cadence_min = slew_sample_cadence_min
        self.output_flag = output_flag

        # Validate required columns.
        missing = [c for c in REQUIRED_COLUMNS if c not in self.requests_frame.columns]
        if missing:
            raise ValueError(
                f"TTPModel.requests_frame missing required columns: {missing}. "
                f"Required: {REQUIRED_COLUMNS}"
            )

        # Attach per-row SkyCoord (``coord`` column). Stored as a list of
        # scalar SkyCoord instances so Observer.altaz(..., grid_times_targets)
        # accepts it directly. Derived from ra/dec; not serialization-safe.
        self.requests_frame["coord"] = list(
            SkyCoord(
                self.requests_frame.ra.values * u.deg,
                self.requests_frame.dec.values * u.deg,
                frame="icrs",
            )
        )

    # ---------------------------------------------------------- helpers

    def _visit_duration(self, exptime_s, n_shots):
        """Total visit duration in *minutes* including readout between shots.

        Same canonical formula as ``Queue.visit_duration``; duplicated here so
        the model is queue-free.
        """
        return (exptime_s * n_shots + self.readout_time * (n_shots - 1)) / 60.0

    def _wrap_az(self, angle_deg):
        """Vectorized wrap-frame shift. ``wrap_limit=None`` means no shift."""
        if not self.wrap_limit:
            return np.asarray(angle_deg)
        a = np.asarray(angle_deg) + (360 - self.wrap_limit)
        return np.where(a > 360, a - 360, a)

    def _short_az_sep(self, az_sep):
        """If telescope has no wrap, az slews never exceed 180 deg."""
        az_sep = np.asarray(az_sep)
        if self.wrap_limit:
            return az_sep
        return np.where(az_sep > 180, 360 - az_sep, az_sep)

    def _slew_minutes(self, alt_a, alt_b, az_a, az_b):
        """Worst-case slew time (minutes) for an alt/az pair, broadcast-friendly."""
        az_sep = self._short_az_sep(
            np.abs(self._wrap_az(az_a) - self._wrap_az(az_b))
        )
        alt_sep = np.abs(np.asarray(alt_a) - np.asarray(alt_b))
        return np.maximum(az_sep, alt_sep) / (60.0 * self.slew_rate)

    def _minutes_from_start(self, t):
        """``astropy.time.Time`` -> minutes since ``self.night_start``."""
        return (t.jd - self.night_start.jd) * 24 * 60

    def _time_from_minutes(self, mins):
        """Minutes since ``self.night_start`` -> ``Time`` (broadcast-friendly)."""
        return self.night_start + TimeDelta(np.asarray(mins) * 60, format="sec")

    # ------------------------------------------------------------- node setup

    def _build_nodes(self):
        """Expand the request frame into the ``self.nodes`` DataFrame.

        Sets attributes:
            N, M, dur, nodes, multi_visit_groups.

        All times in ``self.nodes`` are in minutes from ``night_start``.
        """
        self.dur = float(np.round(self._minutes_from_start(self.night_end), 0))

        anchor_row = {
            "unique_id": "",
            "request_idx": -1,
            "visit_seq": 0,
            "is_anchor": True,
            "t_early": 0.0,
            "t_late": self.dur,
            "t_visit": 0.0,
            "tau_intra": 0.0,
            "priority": 0,
            "coord": None,
        }

        rows = [dict(anchor_row)]  # start anchor at index 0
        for row_idx, row in self.requests_frame.iterrows():
            visits = int(row.n_intra_max)
            for v in range(visits):
                rows.append({
                    "unique_id": row.unique_id,
                    "request_idx": int(row_idx),
                    "visit_seq": v,
                    "is_anchor": False,
                    "t_early": self._minutes_from_start(Time(row.first_available)),
                    "t_late": self._minutes_from_start(Time(row.last_available)),
                    "t_visit": self._visit_duration(float(row.exptime), int(row.n_exp)),
                    "tau_intra": float(row.tau_intra) * 60.0,  # hours -> minutes
                    "priority": int(row.priority),
                    "coord": row.coord,
                })
        rows.append(dict(anchor_row))  # end anchor at index N-1

        self.nodes = pd.DataFrame(rows)
        self.N = len(self.nodes)
        self.M = self.n_slots

        # Multi-visit groups: keyed by unique_id, value is the list of node
        # ids (in visit_seq order) for requests with n_intra_max > 1.
        self.multi_visit_groups = {}
        real = self.nodes[~self.nodes["is_anchor"]]
        for uid, group in real.groupby("unique_id", sort=False):
            if len(group) > 1:
                self.multi_visit_groups[uid] = (
                    group.sort_values("visit_seq").index.tolist()
                )

    # -------------------------------------------------------------- slew grid

    def _compute_tau_slew(self):
        """Build the per-slot worst-case slew tensor.

        Within slot ``m``, sample alt/az for every real (non-anchor) node at
        least every ``slew_sample_cadence_min`` minutes (and at least 3 times
        per slot); the worst-case slew between any two real nodes is the
        maximum of ``max|delta_alt|`` and ``max|delta_az|`` (wrap-aware),
        divided by ``slew_rate``.

        Sets attributes ``w`` and ``tau_slew``.

        ``tau_slew`` is a DataFrame with MultiIndex ``(i, j, m)`` and columns
        ``i_id``, ``j_id`` (the ``unique_id`` of each endpoint, for human
        introspection) and ``tau`` (minutes). Inner arcs only (``i, j`` in
        ``range(1, N - 1)``); anchor arcs are MILP bookkeeping and default to
        ``0.0`` via ``.get(..., 0.0)`` inside :meth:`build_model`.
        """
        M = self.M
        samples_per_slot = int(max(self.dur / (M * self.slew_sample_cadence_min), 3))

        slot_bounds = Time(
            np.linspace(self.night_start.jd, self.night_end.jd, M + 1, endpoint=True),
            format="jd",
        )
        # Sample grid: shape (M, samples_per_slot) -> flatten for one altaz call.
        sample_jd = np.array([
            np.linspace(slot_bounds[m].jd, slot_bounds[m + 1].jd, samples_per_slot)
            for m in range(M)
        ])
        times = Time(sample_jd.ravel(), format="jd")

        real_ids = range(1, self.N - 1)
        n_real = self.N - 2
        real_coords = [self.nodes.at[i, "coord"] for i in real_ids]

        altaz = self.observer.altaz(
            times, real_coords, grid_times_targets=True
        )
        # Reshape to (n_real, M, samples_per_slot) so the slot dimension is explicit.
        alts = altaz.alt.deg.reshape(n_real, M, samples_per_slot)
        azs = altaz.az.deg.reshape(n_real, M, samples_per_slot)

        records = []
        for m in range(M):
            for i, j in permutations(real_ids, 2):
                # Real node ids are 1..N-2; subtract 1 to index alts/azs.
                tau = float(np.round(
                    self._slew_minutes(
                        alts[i - 1, m], alts[j - 1, m],
                        azs[i - 1, m], azs[j - 1, m],
                    ).max(),
                    3,
                ))
                records.append({
                    "i": i,
                    "j": j,
                    "m": m,
                    "i_id": self.nodes.at[i, "unique_id"],
                    "j_id": self.nodes.at[j, "unique_id"],
                    "tau": tau,
                })

        self.tau_slew = (
            pd.DataFrame.from_records(
                records,
                columns=["i", "j", "m", "i_id", "j_id", "tau"],
            ).set_index(["i", "j", "m"])
        )

        # Slot bounds expressed as minutes from start.
        self.w = (slot_bounds.jd - slot_bounds[0].jd) * 24 * 60

    # ------------------------------------------------------------- MILP build

    def build_model(self):
        """Construct the TTP MILP (Handley+ 2024, eqs. 2-9, B3, 10).

        Constraints are added in the order they appear in Handley+ 2024 §2.4.
        Each constraint family iterates one of three slices of ``range(N)``
        that line up with the paper:

        * ``range(1, N - 1)`` -- real visit nodes (paper's ``1..N``).
        * ``range(N - 1)`` -- arc sources, i.e. start anchor + real
          (paper's ``0..N``).
        * ``range(1, N)`` -- arc targets, i.e. real + end anchor
          (paper's ``1..N+1``).

        Per-node scalars (``t_early``, ``t_late``, ``t_visit``, ``tau_intra``,
        ``priority``) are pulled directly from ``self.nodes`` at constraint
        construction. The only hot cube loop is eq. 8 which reads only
        ``self.w``, so no pre-extraction is needed.
        """
        self.model = gp.Model("TTP")
        self.model.Params.OutputFlag = self.output_flag

        N, M = self.N, self.M
        nodes = self.nodes

        # O(1) slew lookup for hot loops (anchor arcs absent => .get(..., 0.0)).
        tau_slew_lookup = self.tau_slew["tau"].to_dict()

        self.Yi = self.model.addVars(
            range(N), vtype=GRB.BINARY, name="Yi",
        )
        self.Xijm = self.model.addVars(
            range(N), range(N), range(M), vtype=GRB.BINARY, name="Xijm",
        )
        self.tijm = self.model.addVars(
            range(N), range(N), range(M), vtype=GRB.CONTINUOUS, name="tijm",
        )
        self.ti = self.model.addVars(
            range(N), vtype=GRB.CONTINUOUS, lb=0, name="ti",
        )

        # eq. 2 - exactly one arc out of the start anchor.
        self.model.addConstr(
            gp.quicksum(
                self.Xijm[0, j, m]
                for j in range(1, N)
                for m in range(M)
            ) == 1,
            "start_anchor",
        )

        # eq. 3 - exactly one arc into the end anchor.
        self.model.addConstr(
            gp.quicksum(
                self.Xijm[i, N - 1, m]
                for i in range(N - 1)
                for m in range(M)
            ) == 1,
            "end_anchor",
        )

        # eq. 4 - visit indicator (one constraint per non-start node).
        for j in range(1, N):
            self.model.addConstr(
                gp.quicksum(
                    self.Xijm[i, j, m]
                    for i in range(N - 1)
                    for m in range(M)
                ) == self.Yi[j],
                f"visit_once_{j}",
            )

        # eq. 5 - flow conservation (one constraint per real node).
        for k in range(1, N - 1):
            self.model.addConstr(
                gp.quicksum(
                    self.Xijm[i, k, m]
                    for i in range(N - 1)
                    for m in range(M)
                )
                - gp.quicksum(
                    self.Xijm[k, j, m]
                    for j in range(1, N)
                    for m in range(M)
                ) == 0,
                f"flow_constr_{k}",
            )

        # eq. 6 - link ti to tijm (one constraint per non-end node).
        for i in range(N - 1):
            self.model.addConstr(
                self.ti[i] == gp.quicksum(
                    self.tijm[i, j, m]
                    for j in range(1, N)
                    for m in range(M)
                ),
                f"tijm_def_{i}",
            )

        # eq. 7 - exposure/slew time linking (one constraint per non-start
        # node). Anchor arcs (i=0) are MILP bookkeeping with tau_slew = 0; the
        # .get(..., 0.0) fallback keeps tau_slew restricted to real arcs.
        for j in range(1, N):
            t_visit_j = nodes.at[j, "t_visit"]
            self.model.addConstr(
                self.ti[j] >= gp.quicksum(
                    self.tijm[i, j, m]
                    + (tau_slew_lookup.get((i, j, m), 0.0) + t_visit_j) * self.Xijm[i, j, m]
                    for i in range(N - 1)
                    for m in range(M)
                ),
                f"exp_constr_{j}",
            )

        # eq. 8 - slot bounds on tijm (one pair per (i, j, m)).
        for i, j, m in product(range(N), range(N), range(M)):
            self.model.addConstr(
                self.tijm[i, j, m] >= self.w[m] * self.Xijm[i, j, m],
                f"t_min_{i}_{j}_{m}",
            )
            self.model.addConstr(
                self.tijm[i, j, m] <= self.w[m + 1] * self.Xijm[i, j, m],
                f"t_max_{i}_{j}_{m}",
            )

        # eq. 9 - node accessibility (one pair per node).
        for i in range(N):
            row = nodes.loc[i]
            self.model.addConstr(
                self.ti[i] >= (row.t_early + row.t_visit) * self.Yi[i],
                f"rise_constr_{i}",
            )
            self.model.addConstr(
                self.ti[i] <= row.t_late * self.Yi[i],
                f"set_constr_{i}",
            )

        # eq. B3 - intra-night separation (multi-visit only).
        for indices in self.multi_visit_groups.values():
            for k in range(1, len(indices)):
                self.model.addConstr(
                    gp.quicksum(
                        self.tijm[indices[k], j, m]
                        for j in range(1, N)
                        for m in range(M)
                    )
                    >= gp.quicksum(
                        self.tijm[indices[k - 1], j, m]
                        for j in range(1, N)
                        for m in range(M)
                    )
                    + self.Yi[indices[k]] * nodes.at[indices[k], "tau_intra"],
                    f"intra_sep_constr_{indices[k - 1]}_{indices[k]}",
                )

        # eq. 10 - objective: priority-weighted visit count minus slew tie-breaker.
        self.model.setObjective(
            gp.quicksum(
                nodes.at[j, "priority"] * self.Yi[j]
                for j in range(1, N - 1)
            )
            - self._SLEW_PENALTY * gp.quicksum(
                tau_slew_lookup.get((i, j, m), 0.0) * self.Xijm[i, j, m]
                for i in range(1, N - 1)
                for j in range(1, N - 1)
                for m in range(M)
            ),
            GRB.MAXIMIZE,
        )

        self.model.params.TimeLimit = self.runtime
        self.model.params.MIPGap = self.optgap
        self.model.update()

    def optimize_model(self):
        """Optimize ``self.model`` and report infeasibility via IIS.

        Mirrors :meth:`astroq.splan.SemesterPlanner.optimize_model`.
        """
        logs.debug("Begin TTP solve.")
        t0 = time.time()
        self.model.optimize()
        if self.model.Status == GRB.INFEASIBLE:
            logs.critical("TTP infeasible; computing IIS.")
            self.model.computeIIS()
            for c in self.model.getConstrs():
                if c.IISConstr:
                    logs.critical(c.ConstrName)
        logs.info(f"TTP solve finished in {time.time() - t0:.3f}s")

    def run_model(self):
        """Build and solve the MILP, then post-process the solution.

        Pipeline:

            _build_nodes -> _compute_tau_slew -> build_model
              -> optimize_model -> digest_gurobi -> optimization_status

        Callers may invoke the individual stages directly to inspect
        intermediate state (e.g. ``self.nodes``, ``self.tau_slew``) without
        paying the solve cost.
        """
        self._build_nodes()
        logs.info(f"Solving TTP for {self.N - 2} visits with Gurobi")
        self._compute_tau_slew()
        self.build_model()
        self.optimize_model()
        if self.model.SolCount > 0:
            self.digest_gurobi()
            self.optimization_status()
        else:
            logs.warning(
                "No incumbent TTP solution within time limit. "
                "Try increasing the ``runtime`` parameter."
            )

    # ---------------------------------------------------------- post-process

    def _extract_tour(self):
        """Walk the chosen tour 0 -> ... -> N-1, populating ``self.tour``.

        Reads solution values directly off the Gurobi ``Var.X`` attribute of
        ``self.Yi`` / ``self.Xijm`` / ``self.ti`` / ``self.tijm``. Stores:

        * ``self.tour`` (DataFrame): one row per real node in traversal order.
          Columns: ``node_id``, ``unique_id``, ``t_start``, ``t_end``,
          ``arc_slot``, ``arc_from``.
        * ``self._Y`` (np.ndarray): per-node visit indicators (length ``N``).
        * ``self._chosen_out`` (dict[int, list[tuple]]):
          ``{i: [(j, m, tijm_val), ...]}`` for every chosen arc i -> j (used
          by :meth:`optimization_status`).
        """
        Y_arr = np.array([self.Yi[i].X for i in range(self.N)])

        # Filter out arcs that exist in the variable space but are
        # unconstrained in the MILP (and thus may be set to 1 by Gurobi
        # without any feasibility penalty):
        #   * j == 0          arcs into the start anchor
        #   * i == self.N - 1 arcs out of the end anchor
        chosen_out = defaultdict(list)
        for (i, j, m), var in self.Xijm.items():
            if var.X > 0.5 and j != 0 and i != self.N - 1:
                chosen_out[i].append((j, m, self.tijm[i, j, m].X))

        tour_rows = []
        node = 0
        order = 0
        # Flow conservation guarantees exactly one outgoing arc per node up to
        # the end anchor; the loop terminates when the end anchor is reached.
        while True:
            succs = chosen_out[node]
            assert len(succs) == 1, (
                f"TTP tour broken: node {node} has {len(succs)} chosen "
                f"outgoing arcs (expected 1). chosen_out={dict(chosen_out)}"
            )
            j, m, _tval = succs[0]
            if j == self.N - 1:
                break
            t_end = self.ti[j].X
            tour_rows.append({
                "order": order,
                "node_id": j,
                "unique_id": self.nodes.at[j, "unique_id"],
                "t_start": t_end - self.nodes.at[j, "t_visit"],
                "t_end": t_end,
                "arc_slot": m,
                "arc_from": node,
            })
            node = j
            order += 1

        if tour_rows:
            self.tour = pd.DataFrame(tour_rows).set_index("order")
        else:
            self.tour = pd.DataFrame(
                columns=[
                    "node_id", "unique_id", "t_start", "t_end",
                    "arc_slot", "arc_from",
                ]
            )
            self.tour.index.name = "order"

        self._Y = Y_arr
        self._chosen_out = chosen_out

    def digest_gurobi(self):
        """Interpret the Gurobi solution; populate schedule / plotly / paths.

        Walks the chosen tour (via :meth:`_extract_tour`) rather than sorting
        visits by float ``ti``; the tour is the canonical TTP output and is
        guaranteed to be in execution order even under numerically degenerate
        solves.

        Sets attributes ``extras``, ``num_scheduled``, ``schedule``,
        ``plotly``, ``times``, ``az_path``, ``alt_path``, ``tour``.
        """
        self._extract_tour()
        self.num_scheduled = len(self.tour)

        # Sanity check: tour order should agree with argsort(ti). Disagreement
        # is loud, not silent -- usually a numerically degenerate solve.
        t_end = self.tour["t_end"].to_numpy()
        if len(t_end) and not np.array_equal(
            np.argsort(t_end, kind="stable"), np.arange(len(t_end))
        ):
            logs.warning(
                "TTP tour order disagrees with argsort(ti); using tour order. "
                "Likely a numerically degenerate Gurobi solution."
            )

        self._tour_to_legacy_dicts()

    def _tour_to_legacy_dicts(self):
        """Map ``self.tour`` -> legacy dict outputs.

        The wire format of ``schedule`` / ``plotly`` / ``extras`` (dict keys,
        value dtypes, ``Time`` in JD where applicable) is the public contract
        consumed by :mod:`astroq.nplan`, :mod:`astroq.ttp.plot`, and
        :mod:`astroq.plot`; do not change it without coordinating with those
        modules.
        """
        tour = self.tour
        tour_nodes = tour["node_id"].tolist()

        # ----- extras: real nodes not in the tour ------------------------
        scheduled = set(tour_nodes)
        unscheduled = [i for i in range(1, self.N - 1) if i not in scheduled]
        if unscheduled:
            self.extras = {
                "Starname": self.nodes.loc[unscheduled, "unique_id"].tolist(),
                "First Available": [
                    self._time_from_minutes(self.nodes.at[i, "t_early"]).isot[11:16]
                    for i in unscheduled
                ],
                "Last Available": [
                    self._time_from_minutes(self.nodes.at[i, "t_late"]).isot[11:16]
                    for i in unscheduled
                ],
            }
        else:
            self.extras = {"Starname": [], "First Available": [], "Last Available": []}

        # ----- per-visit arrays -----------------------------------------
        if tour_nodes:
            request_indices = self.nodes.loc[tour_nodes, "request_idx"].to_numpy()
            rows = self.requests_frame.iloc[request_indices]
            visit_dur = self.nodes.loc[tour_nodes, "t_visit"].to_numpy()
            t_start = tour["t_start"].to_numpy()
            t_end = tour["t_end"].to_numpy()

            t_start_time = self._time_from_minutes(t_start)
            t_end_time = self._time_from_minutes(t_end)

            # Two batched altaz calls cover the start and end of every visit.
            coords = [self.nodes.at[i, "coord"] for i in tour_nodes]
            aa_start = self.observer.altaz(t_start_time, coords)
            aa_end = self.observer.altaz(t_end_time, coords)
            az_start = np.atleast_1d(aa_start.az.deg)
            alt_start = np.atleast_1d(aa_start.alt.deg)
            az_end = np.atleast_1d(aa_end.az.deg)
            alt_end = np.atleast_1d(aa_end.alt.deg)
            az_path = np.empty(2 * len(tour_nodes))
            alt_path = np.empty(2 * len(tour_nodes))
            az_path[0::2], az_path[1::2] = az_start, az_end
            alt_path[0::2], alt_path[1::2] = alt_start, alt_end
            times_list = [
                t for pair in zip(list(t_start_time), list(t_end_time)) for t in pair
            ]

            self.schedule = {
                "Order": list(range(len(tour_nodes))),
                "Starname": tour["unique_id"].tolist(),
                "Time": [self.night_start.jd + t / (24 * 60) for t in t_start],
            }
            self.plotly = {
                "Starname": tour["unique_id"].tolist(),
                "First Available": self.nodes.loc[tour_nodes, "t_early"].to_numpy(),
                "Last Available": self.nodes.loc[tour_nodes, "t_late"].to_numpy(),
                "Start Exposure": t_start,
                "Minutes the from Start of the Night": (t_start + t_end) / 2,
                "Stop Exposure": t_end,
                "N_shots": rows["n_exp"].astype(int).to_numpy(),
                "Exposure Time (min)": rows["exptime"].astype(float).to_numpy(),
                "Total Exp Time (min)": visit_dur,
                "Priority": rows["priority"].astype(int).tolist(),
            }
            self.times = times_list
            self.az_path = az_path
            self.alt_path = alt_path
        else:
            self.schedule = {"Order": [], "Starname": [], "Time": []}
            self.plotly = {
                "Starname": [],
                "First Available": np.array([]),
                "Last Available": np.array([]),
                "Start Exposure": np.array([]),
                "Minutes the from Start of the Night": np.array([]),
                "Stop Exposure": np.array([]),
                "N_shots": np.array([], dtype=int),
                "Exposure Time (min)": np.array([]),
                "Total Exp Time (min)": np.array([]),
                "Priority": [],
            }
            self.times = []
            self.az_path = np.array([])
            self.alt_path = np.array([])

    def optimization_status(self):
        """Write ``TTPstatistics.txt`` and log solver wall-clock stats.

        ``time_slewing`` is computed directly from the chosen-arc tau_slew
        values rather than from the fractional part of the objective bound;
        this is correct under non-uniform priorities.
        """
        # Sum tau_slew over inner arcs (exclude anchor arcs) by joining the
        # chosen-arc keys against self.tau_slew.
        chosen_arcs = [
            (i, j, m)
            for i, succs in self._chosen_out.items()
            for (j, m, _tval) in succs
            if 0 < i < self.N - 1 and 0 < j < self.N - 1
        ]
        if chosen_arcs:
            chosen_idx = pd.MultiIndex.from_tuples(
                chosen_arcs, names=["i", "j", "m"]
            )
            time_slewing = float(
                self.tau_slew.loc[chosen_idx, "tau"].sum()
            )
        else:
            time_slewing = 0.0

        t_visit = self.nodes["t_visit"].to_numpy()
        time_exposing = float(np.dot(self._Y, t_visit))
        time_idle = self.dur - time_exposing - time_slewing

        self.time_idle = time_idle
        self.time_slewing = time_slewing
        self.time_exposing = time_exposing
        self.solve_time = self.model.Runtime

        lines = [
            "Stats for TTP Solution",
            "------------------------------------",
            f"    Model ran for {self.solve_time:.2f} seconds",
            f"     Observations Requested: {self.N - 2}",
            f"     Observations Scheduled: {self.num_scheduled}",
            "------------------------------------",
            f"   Observing Duration (min): {self.dur:.2f}",
            f"  Time Spent Exposing (min): {self.time_exposing:.2f}",
            f"      Time Spent Idle (min): {self.time_idle:.2f}",
            f"   Time Spent Slewing (min): {self.time_slewing:.2f}",
            "------------------------------------",
        ]
        block = "\n".join(lines) + "\n"

        with open(self.outdir + "/TTPstatistics.txt", "w") as fh:
            fh.write(block)
        logs.info("\n" + block)
