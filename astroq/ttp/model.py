"""Traveling Telescope Problem (TTP) solver

The model implements the MILP of Handley et al. 2024 (arXiv:2310.18497).

Required columns of ``requests_frame`` (see :data:`REQUIRED_COLUMNS`)::

    unique_id        str            primary key
    ra, dec          float, deg
    exptime          float, seconds (single shot exposure time)
    n_exp            int            shots per visit (multi-shot exposures)
    n_intra_max      int            visits per night
    tau_intra        float, hours   minimum spacing between visits within a night
    priority         int            objective weight; higher = more important
    first_available  str            ISO-8601 (caller computes; e.g. via Access)
    last_available   str            ISO-8601

Internal naming is aligned with the Handley 2024 paper (``N``, ``M``, ``Yi``,
``Xijm``, ``tau_slew``) and with AstroQ vocabulary elsewhere (``t_visit``,
``tau_intra``).
"""

# Standard library imports
import logging
import time
from collections import defaultdict
from itertools import permutations

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

    The class follows the same layout as ``astroq.splan.SemesterPlanner``:
    ``self.model`` is the Gurobi model, Gurobi variables are named class
    attributes (``self.Yi`` / ``self.Xijm`` / ``self.ti`` / ``self.tijm``),
    and three entry points -- ``build_model``, ``optimize_model``,
    ``run_model`` -- mirror the splan equivalents.

    Args:
        requests_frame (pd.DataFrame): one row per request, AstroQ-vocab columns
            (see :data:`REQUIRED_COLUMNS`). A derived ``coord`` column
            (``SkyCoord``) is attached on init and is NOT serialization-safe.
        night_start (astropy.time.Time): start of the observing interval.
        night_end (astropy.time.Time): end of the observing interval.
        outdir (str): directory in which to write ``TTPstatistics.txt``.

    Keyword Args:
        observer (astroplan.Observer): site fixture for alt/az lookups. Used by
            :meth:`compute_tau_slew` to grid the slew tensor and by
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

    Attributes:
        requests_frame, night_start, night_end, outdir, observer, slew_rate,
        wrap_limit, readout_time, n_slots, inaccessible_zones, runtime, optgap,
        slew_sample_cadence_min, output_flag: as above.
        N (int): total node count (= visit count + 2 anchor nodes).
        M (int): equal to ``n_slots``.
        dur (float): observing-interval duration, minutes.
        node_to_request (dict[int, int]): node index -> row index in
            ``requests_frame``.
        multi_visit_ind (dict[int, list[int]]): row index -> list of node
            indices for multi-visit requests.
        t_visit (np.ndarray): per-node visit duration (minutes, including
            readout).
        tau_intra (np.ndarray): per-node intra-night cadence (minutes).
        t_early, t_late (np.ndarray): per-node first/last allowed completion
            times (minutes from ``night_start``).
        priorities (np.ndarray): per-node objective weights.
        tau_slew (dict): slew tensor ``{(i,j,m): minutes}``.
        w (np.ndarray): slot bound times (minutes from ``night_start``).
        model (gurobipy.Model): underlying Gurobi model.
        Yi, Xijm, ti, tijm (gurobipy tupledict): MILP decision variables.
        schedule, plotly, extras, times, az_path, alt_path,
        num_scheduled, time_idle, time_slewing, time_exposing, solve_time:
        populated by :meth:`digest_gurobi` / :meth:`optimization_status`.
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

        self.create_nodes()
        self.compute_tau_slew()
        self.run_model()

    # ------------------------------------------------------------------ setup

    def _visit_duration(self, exptime_s, n_shots):
        """Total visit duration in *minutes* including readout between shots.

        Same canonical formula as ``Queue.visit_duration``; duplicated here so
        the model is queue-free.
        """
        return (exptime_s * n_shots + self.readout_time * (n_shots - 1)) / 60.0

    # --------------------------------------------------------- geometry helpers

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

    # ------------------------------------------------------------ time helpers

    def _minutes_from_start(self, t):
        """``astropy.time.Time`` -> minutes since ``self.night_start``."""
        return (t.jd - self.night_start.jd) * 24 * 60

    def _time_from_minutes(self, mins):
        """Minutes since ``self.night_start`` -> ``Time`` (broadcast-friendly)."""
        return self.night_start + TimeDelta(np.asarray(mins) * 60, format="sec")

    # ------------------------------------------------------------- node setup

    def create_nodes(self):
        """Expand the request frame into TTP nodes (one per visit + 2 anchors).

        For each request with ``n_intra_max`` visits per night we create that
        many consecutive nodes; the multi-visit ordering constraint is layered
        on later in :meth:`build_model`.

        Sets attributes:
            N, dur, t_early, t_late, t_visit, tau_intra, priorities,
            node_to_request, multi_visit_ind, _node_coords.

        All times are in minutes from ``night_start``.
        """
        self.dur = np.round(self._minutes_from_start(self.night_end), 0)

        node_to_request = defaultdict(int)
        multi_visit_ind = defaultdict(list)

        i = 1
        for row_idx, row in self.requests_frame.iterrows():
            visits = int(row.n_intra_max)
            for _ in range(visits):
                if visits > 1:
                    multi_visit_ind[row_idx].append(i)
                node_to_request[i] = row_idx
                i += 1

        N = int(self.requests_frame.n_intra_max.sum()) + 2
        self.N = N

        t_early, t_late = [], []
        t_visit, tau_intra, priorities = [], [], []

        for i in range(N):
            if i == 0 or i == N - 1:
                t_early.append(0)
                t_late.append(self.dur)
                t_visit.append(0)
                tau_intra.append(0)
                priorities.append(0)
            else:
                row = self.requests_frame.iloc[node_to_request[i]]
                priorities.append(int(row.priority))
                t_visit.append(self._visit_duration(float(row.exptime), int(row.n_exp)))
                tau_intra.append(float(row.tau_intra) * 60.0)        # hours -> minutes
                t_early.append(self._minutes_from_start(Time(row.first_available)))
                t_late.append(self._minutes_from_start(Time(row.last_available)))

        self.t_early = np.array(t_early)
        self.t_late = np.array(t_late)
        self.priorities = np.array(priorities)
        self.t_visit = np.array(t_visit)
        self.tau_intra = np.array(tau_intra)
        self.node_to_request = node_to_request
        self.multi_visit_ind = multi_visit_ind

        # Cache per-node SkyCoord (real nodes only, length N-2; index by node-1).
        self._node_coords = [
            self.requests_frame.iloc[node_to_request[i]]["coord"]
            for i in range(1, N - 1)
        ]

    # -------------------------------------------------------------- slew grid

    def compute_tau_slew(self):
        """Build the per-slot worst-case slew tensor ``tau_slew[i,j,m]``.

        Within slot ``m``, sample alt/az for every target at least every
        ``slew_sample_cadence_min`` minutes (and at least 3 times per slot);
        the worst-case slew between any two targets is the maximum of
        ``max|delta_alt|`` and ``max|delta_az|`` (wrap-aware), divided by
        ``slew_rate``.

        Sets attributes ``M``, ``w``, ``tau_slew``.
        """
        self.M = self.n_slots
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

        altaz = self.observer.altaz(
            times, self._node_coords, grid_times_targets=True
        )
        # Reshape to (N-2, M, samples_per_slot) so the slot dimension is explicit.
        alts = altaz.alt.deg.reshape(self.N - 2, M, samples_per_slot)
        azs = altaz.az.deg.reshape(self.N - 2, M, samples_per_slot)

        tau_slew = defaultdict(float)  # minutes
        for m in range(M):
            for i, j in permutations(range(1, self.N - 1), 2):
                tau_slew[(i, j, m)] = np.round(
                    self._slew_minutes(
                        alts[i - 1, m], alts[j - 1, m],
                        azs[i - 1, m], azs[j - 1, m],
                    ).max(),
                    3,
                )

        # Slot bounds expressed as minutes from start.
        self.w = (slot_bounds.jd - slot_bounds[0].jd) * 24 * 60
        self.tau_slew = tau_slew
        import pdb; pdb.set_trace()

    # ------------------------------------------------------------- MILP build

    def build_model(self):
        """Construct the TTP MILP (Handley+ 2024, eqs. 2-9, B3, 10).

        Constraints are added in the order they appear in Handley+ 2024 §2.4.
        Each constraint family uses an explicit Python ``for`` loop so the
        index structure is visible at the indentation level. ``gp.quicksum``
        generators are used only for the *summation* dimensions.
        """
        self.model = gp.Model("TTP")
        self.model.Params.OutputFlag = self.output_flag

        N, M = self.N, self.M
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

        # eqs. 2-3 - anchor flow: exactly one arc out of node 0, exactly one
        # arc into node N-1.
        self.model.addConstr(
            gp.quicksum(
                self.Xijm[0, j, m]
                for j in range(N)
                for m in range(M)
            ) == 1,
            "start_anchor",
        )
        self.model.addConstr(
            gp.quicksum(
                self.Xijm[i, N - 1, m]
                for i in range(N)
                for m in range(M)
            ) == 1,
            "end_anchor",
        )

        # eq. 4 - visit indicator (one constraint per j)
        for j in range(N)[1:]:
            self.model.addConstr(
                gp.quicksum(
                    self.Xijm[i, j, m]
                    for i in range(N)[:-1]
                    for m in range(M)
                ) == self.Yi[j],
                f"visit_once_{j}",
            )

        # eq. 5 - flow conservation (one constraint per internal node k)
        for k in range(N)[1:-1]:
            self.model.addConstr(
                gp.quicksum(
                    self.Xijm[i, k, m]
                    for i in range(N)[:-1]
                    for m in range(M)
                )
                - gp.quicksum(
                    self.Xijm[k, j, m]
                    for j in range(N)[1:]
                    for m in range(M)
                ) == 0,
                f"flow_constr_{k}",
            )

        # eq. 6 - link ti to tijm (one constraint per i)
        for i in range(N)[:-1]:
            self.model.addConstr(
                self.ti[i] == gp.quicksum(
                    self.tijm[i, j, m]
                    for j in range(N)[1:]
                    for m in range(M)
                ),
                f"tijm_def_{i}",
            )

        # eq. 7 - exposure/slew time linking (one constraint per j)
        for j in range(N)[1:]:
            self.model.addConstr(
                self.ti[j] >= gp.quicksum(
                    self.tijm[i, j, m]
                    + (self.tau_slew[(i, j, m)] + self.t_visit[j]) * self.Xijm[i, j, m]
                    for i in range(N)[:-1]
                    for m in range(M)
                ),
                f"exp_constr_{j}",
            )

        # eq. 8 - slot bounds on tijm (one pair per (i, j, m))
        for i in range(N):
            for j in range(N):
                for m in range(M):
                    self.model.addConstr(
                        self.tijm[i, j, m] >= self.w[m] * self.Xijm[i, j, m],
                        f"t_min_{i}_{j}_{m}",
                    )
                    self.model.addConstr(
                        self.tijm[i, j, m] <= self.w[m + 1] * self.Xijm[i, j, m],
                        f"t_max_{i}_{j}_{m}",
                    )

        # eq. 9 - node accessibility (one pair per i)
        for i in range(N):
            self.model.addConstr(
                self.ti[i] >= (self.t_early[i] + self.t_visit[i]) * self.Yi[i],
                f"rise_constr_{i}",
            )
            self.model.addConstr(
                self.ti[i] <= self.t_late[i] * self.Yi[i],
                f"set_constr_{i}",
            )

        # eq. B3 - intra-night separation (multi-visit only)
        for indices in self.multi_visit_ind.values():
            for k in range(1, len(indices)):
                self.model.addConstr(
                    gp.quicksum(
                        self.tijm[indices[k], j, m]
                        for j in range(N)[1:]
                        for m in range(M)
                    )
                    >= gp.quicksum(
                        self.tijm[indices[k - 1], j, m]
                        for j in range(N)[1:]
                        for m in range(M)
                    )
                    + self.Yi[indices[k]] * self.tau_intra[indices[k]],
                    f"intra_sep_constr_{indices[k - 1]}_{indices[k]}",
                )

        # eq. 10 - objective: priority-weighted visit count minus slew tie-breaker
        self.model.setObjective(
            gp.quicksum(
                self.priorities[j] * self.Yi[j]
                for j in range(N)[1:-1]
            )
            - self._SLEW_PENALTY * gp.quicksum(
                self.tau_slew[(i, j, m)] * self.Xijm[i, j, m]
                for i in range(N)[1:-1]
                for j in range(N)[1:-1]
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
        """Build, optimize, and post-process the TTP MILP.

        Mirrors :meth:`astroq.splan.SemesterPlanner.run_model`.
        """
        logs.info(f"Solving TTP for {self.N - 2} exposures with Gurobi")
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
        """Walk the chosen tour 0 -> ... -> N-1 and read every Gurobi var once.

        Returns a dict with keys:
            tour_nodes (list[int]): real (non-anchor) node indices in
                traversal order.
            t_end (np.ndarray): completion time ``ti`` per tour node
                (minutes from ``night_start``).
            arc_slots (np.ndarray): slot ``m`` used for each arc *into*
                ``tour_nodes[k]``. ``arc_slots[0]`` is the 0 -> first-visit
                arc.
            Y (np.ndarray): per-node visit indicators (length ``N``).
            X_chosen (dict[int, list[tuple]]): ``{i: [(j, m, tijm_val), ...]}``
                for every chosen arc i -> j (used by ``optimization_status``).
        """
        Y_dict = self.model.getAttr("X", self.Yi)
        X_dict = self.model.getAttr("X", self.Xijm)
        ti_dict = self.model.getAttr("X", self.ti)
        tijm_dict = self.model.getAttr("X", self.tijm)

        Y_arr = np.array([Y_dict[i] for i in range(self.N)])

        # Filter out arcs that exist in the variable space but are
        # unconstrained in the MILP (and thus may be set to 1 by Gurobi
        # without any feasibility penalty):
        #   * j == 0          arcs into the start anchor
        #   * i == self.N - 1 arcs out of the end anchor
        # The original digest_gurobi avoided this bug by only iterating
        # real-internal i and j; we make the filter explicit here so
        # the tour walk is robust.
        chosen_out = defaultdict(list)
        for (i, j, m), v in X_dict.items():
            if v > 0.5 and j != 0 and i != self.N - 1:
                chosen_out[i].append((j, m, tijm_dict[(i, j, m)]))

        tour_nodes, arc_slots, arc_times = [], [], []
        node = 0
        # Flow conservation guarantees exactly one outgoing arc per node up to
        # the end anchor; the loop terminates when the end anchor is reached.
        while True:
            succs = chosen_out[node]
            assert len(succs) == 1, (
                f"TTP tour broken: node {node} has {len(succs)} chosen "
                f"outgoing arcs (expected 1). chosen_out={dict(chosen_out)}"
            )
            j, m, tval = succs[0]
            arc_slots.append(m)
            arc_times.append(tval)
            if j == self.N - 1:
                break
            tour_nodes.append(j)
            node = j

        t_end = np.array([ti_dict[i] for i in tour_nodes])
        return {
            "tour_nodes": tour_nodes,
            "t_end": t_end,
            "arc_slots": np.array(arc_slots),
            "arc_times": np.array(arc_times),
            "Y": Y_arr,
            "X_chosen": chosen_out,
        }

    def digest_gurobi(self):
        """Interpret the Gurobi solution; populate schedule / plotly / paths.

        Walks the chosen tour ``0 -> j_1 -> ... -> N-1`` rather than sorting
        visits by float ``ti``; the tour is the canonical TTP output and is
        guaranteed to be in execution order even under numerically degenerate
        solves.

        Sets attributes ``extras``, ``num_scheduled``, ``schedule``,
        ``plotly``, ``times``, ``az_path``, ``alt_path``, ``_tour``.
        """
        tour = self._extract_tour()
        self._tour = tour

        tour_nodes = tour["tour_nodes"]
        t_end = tour["t_end"]
        self.num_scheduled = len(tour_nodes)

        # Sanity check: tour order should agree with argsort(ti). Disagreement
        # is loud, not silent -- usually a numerically degenerate solve.
        if not np.array_equal(
            np.argsort(t_end, kind="stable"), np.arange(len(t_end))
        ):
            logs.warning(
                "TTP tour order disagrees with argsort(ti); using tour order. "
                "Likely a numerically degenerate Gurobi solution."
            )

        # ----- extras: real nodes not in the tour ------------------------
        real_nodes = range(1, self.N - 1)
        scheduled = set(tour_nodes)
        unscheduled = [i for i in real_nodes if i not in scheduled]
        if unscheduled:
            extras_rows = self.requests_frame.iloc[
                [self.node_to_request[i] for i in unscheduled]
            ]
            self.extras = {
                "Starname": list(extras_rows.unique_id),
                "First Available": [
                    self._time_from_minutes(self.t_early[i]).isot[11:16]
                    for i in unscheduled
                ],
                "Last Available": [
                    self._time_from_minutes(self.t_late[i]).isot[11:16]
                    for i in unscheduled
                ],
            }
        else:
            self.extras = {"Starname": [], "First Available": [], "Last Available": []}

        # ----- per-visit arrays -----------------------------------------
        rows = self.requests_frame.iloc[
            [self.node_to_request[i] for i in tour_nodes]
        ]
        visit_dur = self.t_visit[np.array(tour_nodes)] if tour_nodes else np.array([])
        t_start = t_end - visit_dur

        t_start_time = self._time_from_minutes(t_start) if len(t_start) else []
        t_end_time = self._time_from_minutes(t_end) if len(t_end) else []

        # Two batched altaz calls replace 2N scalar ones from the old loop.
        coords = [self._node_coords[i - 1] for i in tour_nodes]
        if coords:
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
        else:
            az_path = np.array([])
            alt_path = np.array([])
            times_list = []

        # ----- public outputs --------------------------------------------
        # schedule['Time']: start-of-exposure time in JD (legacy contract;
        # consumed by nplan.drop_gap_rows and the plotter).
        self.schedule = {
            "Order": list(range(len(tour_nodes))),
            "Starname": list(rows.unique_id),
            "Time": [self.night_start.jd + t / (24 * 60) for t in t_start],
        }

        self.plotly = {
            "Starname": list(rows.unique_id),
            "First Available": self.t_early[np.array(tour_nodes)] if tour_nodes else np.array([]),
            "Last Available": self.t_late[np.array(tour_nodes)] if tour_nodes else np.array([]),
            "Start Exposure": t_start,
            "Minutes the from Start of the Night": (t_start + t_end) / 2 if len(t_start) else np.array([]),
            "Stop Exposure": t_end,
            "N_shots": rows["n_exp"].astype(int).to_numpy() if tour_nodes else np.array([], dtype=int),
            "Exposure Time (min)": rows["exptime"].astype(float).to_numpy() if tour_nodes else np.array([]),
            "Total Exp Time (min)": visit_dur,
            "Priority": rows["priority"].astype(int).tolist(),
        }

        self.times = times_list
        self.az_path = az_path
        self.alt_path = alt_path

    def optimization_status(self):
        """Write ``TTPstatistics.txt`` and log solver wall-clock stats.

        ``time_slewing`` is computed directly from the chosen-arc tau_slew
        values rather than from the fractional part of the objective bound;
        this is correct under non-uniform priorities.
        """
        tour = self._tour

        # Sum tau_slew over inner arcs (exclude anchor arcs).
        time_slewing = 0.0
        for i, succs in tour["X_chosen"].items():
            if not (0 < i < self.N - 1):
                continue
            for (j, m, _tval) in succs:
                if 0 < j < self.N - 1:
                    time_slewing += self.tau_slew[(i, j, m)]

        time_exposing = float(np.dot(tour["Y"], self.t_visit))
        time_idle = self.dur - time_exposing - time_slewing

        self.time_idle = time_idle
        self.time_slewing = time_slewing
        self.time_exposing = time_exposing
        self.solve_time = self.model.Runtime

        # Build the stats block once, then write + log.
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
