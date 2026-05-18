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

    Attributes:
        requests_frame, night_start, night_end, outdir, observer, slew_rate,
        wrap_limit, readout_time, n_slots, inaccessible_zones, runtime, optgap,
        slew_sample_cadence_min: as above.
        N (int): total node count (= visit count + 2 anchor nodes).
        M (int): equal to ``n_slots``.
        dur (float): observing-interval duration, minutes.
        node_to_request (dict[int, int]): node index → row index in
            ``requests_frame``.
        multi_visit_ind (dict[int, list[int]]): row index → list of node
            indices for multi-visit requests.
        t_visit (np.ndarray): per-node visit duration (minutes, including
            readout).
        tau_intra (np.ndarray): per-node intra-night cadence (minutes).
        t_early, t_late (np.ndarray): per-node first/last allowed completion
            times (minutes from ``night_start``).
        priorities (np.ndarray): per-node objective weights.
        tau_slew (dict): slew tensor ``{(i,j,m): minutes}``.
        w (np.ndarray): slot bound times (minutes from ``night_start``).
        gurobi_model (gurobipy.Model): underlying Gurobi model.
        schedule, plotly, extras, times, az_path, alt_path,
        num_scheduled, estimated_slews, real_slews, time_idle, time_slewing,
        time_exposing, solve_time: populated by :meth:`digest_gurobi` /
        :meth:`optimization_status`.
    """

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

        self._validate_columns()
        self._attach_coords()
        self.create_nodes()
        self.compute_tau_slew()
        self.solve()

    # ------------------------------------------------------------------ setup

    def _validate_columns(self):
        missing = [c for c in REQUIRED_COLUMNS if c not in self.requests_frame.columns]
        if missing:
            raise ValueError(
                f"TTPModel.requests_frame missing required columns: {missing}. "
                f"Required: {REQUIRED_COLUMNS}"
            )

    def _attach_coords(self):
        """Materialize a per-row ``SkyCoord`` (``coord`` column).

        Stored as a list of scalar ``SkyCoord`` instances so
        ``Observer.altaz(time, list_of_coords, grid_times_targets=True)``
        accepts it directly. Derived from ``ra``/``dec``; should not be
        persisted to disk (rebuild on load).
        """
        df = self.requests_frame
        self.requests_frame["coord"] = list(
            SkyCoord(
                df.ra.values * u.deg,
                df.dec.values * u.deg,
                frame="icrs",
            )
        )

    def _visit_duration(self, exptime_s, n_shots):
        """Total visit duration in *minutes* including readout between shots.

        Same canonical formula as ``Queue.visit_duration``; duplicated here so
        the model is queue-free.
        """
        return (exptime_s * n_shots + self.readout_time * (n_shots - 1)) / 60.0

    # ------------------------------------------------------------- node setup

    def create_nodes(self):
        """Expand the request frame into TTP nodes (one per visit + 2 anchors).

        For each request with ``n_intra_max`` visits per night we create that
        many consecutive nodes; the multi-visit ordering constraint is layered
        on later in :meth:`to_gurobi_model`.

        Sets attributes:
            N, dur, t_early, t_late, t_visit, tau_intra, priorities,
            node_to_request, multi_visit_ind.

        All times are in minutes from ``night_start``.
        """
        self.dur = np.round((self.night_end - self.night_start).jd * 24 * 60, 0)

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
                first_jd = Time(row.first_available).jd
                last_jd = Time(row.last_available).jd
                t_early.append((first_jd - self.night_start.jd) * 24 * 60)
                t_late.append((last_jd - self.night_start.jd) * 24 * 60)

        self.t_early = np.array(t_early)
        self.t_late = np.array(t_late)
        self.priorities = np.array(priorities)
        self.t_visit = np.array(t_visit)
        self.tau_intra = np.array(tau_intra)
        self.node_to_request = node_to_request
        self.multi_visit_ind = multi_visit_ind

    # -------------------------------------------------------------- slew grid

    def compute_tau_slew(self):
        """Build the per-slot worst-case slew tensor ``tau_slew[i,j,m]``.

        Within slot ``m``, sample alt/az for every target at least every
        ``slew_sample_cadence_min`` minutes (and at least 3 times per slot);
        the worst-case slew between any two targets is the maximum of
        ``max|Δalt|`` and ``max|Δaz|`` (wrap-aware), divided by ``slew_rate``.

        Sets attributes ``M``, ``w``, ``tau_slew``.
        """
        M = self.n_slots
        self.M = M
        night_dur = (self.night_end.jd - self.night_start.jd) * 24 * 60  # minutes

        samples_per_slot = int(max(night_dur / (M * self.slew_sample_cadence_min), 3))

        slot_bounds = Time(
            np.linspace(self.night_start.jd, self.night_end.jd, M + 1, endpoint=True),
            format="jd",
        )

        t = []
        for m in range(M):
            for item in np.linspace(slot_bounds[m].jd, slot_bounds[m + 1].jd, samples_per_slot):
                t.append(item)
        times = Time(t, format="jd")

        N = self.N
        node_to_request = self.node_to_request

        # Reuse the cached per-row SkyCoord column. Each entry is a scalar
        # SkyCoord; astroplan's Observer.altaz accepts a list of them.
        target_coords = [
            self.requests_frame.iloc[node_to_request[i]]["coord"]
            for i in range(1, N - 1)
        ]
        coordinate_matrix = self.observer.altaz(
            times, target_coords, grid_times_targets=True
        )

        def to_wrap_frame(angle):
            if self.wrap_limit:
                angle = angle + (360 - self.wrap_limit)
                angle = np.array([x - 360 if x > 360 else x for x in angle])
            return angle

        def max_ang_sep(targ1, targ2, slot):
            slot_ind_start = slot * samples_per_slot
            slot_ind_end = (slot + 1) * samples_per_slot - 1

            # -1 below: coordinate_matrix excludes the anchor node 0.
            coords1 = coordinate_matrix[targ1 - 1, slot_ind_start:slot_ind_end + 1]
            coords2 = coordinate_matrix[targ2 - 1, slot_ind_start:slot_ind_end + 1]

            alt1 = coords1.alt.deg
            alt2 = coords2.alt.deg
            az1 = to_wrap_frame(coords1.az.deg)
            az2 = to_wrap_frame(coords2.az.deg)

            alt_sep = np.abs(alt1 - alt2)
            az_sep = np.abs(az1 - az2)

            # No-wrap telescopes: no slew greater than 180 deg can exist.
            if not self.wrap_limit:
                az_sep = [360 - x if x > 180 else x for x in az_sep]

            return max(max(alt_sep), max(az_sep))

        tau_slew = defaultdict(float)  # minutes
        for m in range(M):
            for targ1, targ2 in permutations(range(1, N - 1), 2):
                tau_slew[(targ1, targ2, m)] = np.round(
                    max_ang_sep(targ1, targ2, m) / (60 * self.slew_rate), 3
                )

        # Slot bounds expressed as minutes from start.
        self.w = (slot_bounds.jd - slot_bounds[0].jd) * 24 * 60
        self.tau_slew = tau_slew

    # ------------------------------------------------------------- MILP build

    def to_gurobi_model(self, output_flag=True):
        """Translate the TTP into a Gurobi MILP.

        Variable names match Handley 2024:
        ``Yi[i]`` (binary visit indicator), ``Xijm[i,j,m]`` (binary arc
        traversal), ``tijm[i,j,m]`` (continuous), ``ti[i]`` (continuous).
        """
        Mod = gp.Model("TTP")
        Mod.Params.OutputFlag = output_flag

        N = self.N
        M = self.M
        w = self.w
        tau_slew = self.tau_slew
        t_visit = self.t_visit
        t_early = self.t_early
        t_late = self.t_late
        tau_intra = self.tau_intra

        Yi = Mod.addVars(range(N), vtype=GRB.BINARY, name="Yi")
        Xijm = Mod.addVars(range(N), range(N), range(M), vtype=GRB.BINARY, name="Xijm")
        tijm = Mod.addVars(range(N), range(N), range(M), vtype=GRB.CONTINUOUS, name="tijm")
        ti = Mod.addVars(range(N), vtype=GRB.CONTINUOUS, lb=0, name="ti")

        Mod.addConstrs(
            (ti[i] == gp.quicksum(tijm[i, j, m] for j in range(N)[1:] for m in range(M))
             for i in range(N)[:-1]), "tijm_def",
        )
        Mod.addConstr(
            gp.quicksum(Xijm[0, j, m] for j in range(N) for m in range(M)) == 1,
            "start_origin",
        )
        Mod.addConstr(
            gp.quicksum(Xijm[i, N - 1, m] for i in range(N) for m in range(M)) == 1,
            "end_origin",
        )
        Mod.addConstrs(
            (gp.quicksum(Xijm[i, j, m] for i in range(N)[:-1] for m in range(M)) == Yi[j]
             for j in range(N)[1:]), "visit_once",
        )
        Mod.addConstrs(
            ((gp.quicksum(Xijm[i, k, m] for i in range(N)[:-1] for m in range(M))
              - gp.quicksum(Xijm[k, j, m] for j in range(N)[1:] for m in range(M)) == 0)
             for k in range(N)[:-1][1:]), "flow_constr",
        )
        Mod.addConstrs(
            (ti[j] >= gp.quicksum(
                tijm[i, j, m] + (tau_slew[(i, j, m)] + t_visit[j]) * Xijm[i, j, m]
                for i in range(N)[:-1] for m in range(M))
             for j in range(N)[1:]), "exp_constr",
        )
        Mod.addConstrs(
            ((tijm[i, j, m] >= w[m] * Xijm[i, j, m])
             for j in range(N) for m in range(M) for i in range(N)), "t_min",
        )
        Mod.addConstrs(
            (tijm[i, j, m] <= w[m + 1] * Xijm[i, j, m]
             for j in range(N) for m in range(M) for i in range(N)), "t_max",
        )
        Mod.addConstrs(
            (ti[i] >= (t_early[i] + t_visit[i]) * Yi[i] for i in range(N)), "rise_constr",
        )
        Mod.addConstrs((ti[i] <= t_late[i] * Yi[i] for i in range(N)), "set_constr")

        # Multi-visit intra-night cadence (Handley 2024 Appendix B.2 eq. B3).
        for targ in self.multi_visit_ind.keys():
            indices = self.multi_visit_ind[targ]
            Mod.addConstrs(
                ((gp.quicksum(tijm[indices[i], j, m] for j in range(N)[1:] for m in range(M))
                  >= (gp.quicksum(tijm[indices[i - 1], j, m] for j in range(N)[1:] for m in range(M))
                      + Yi[indices[i]] * tau_intra[indices[i]]))
                 for i in range(len(indices))[1:]), "intra_sep_constr",
            )

        priority_param = self.priorities
        slew_param = 1 / 100

        # Primary reward: visit count weighted by priority. Tie-break by slew time.
        Mod.setObjective(
            gp.quicksum(priority_param[j] * Yi[j] for j in range(N)[1:-1])
            - slew_param * gp.quicksum(
                tau_slew[(i, j, m)] * Xijm[i, j, m]
                for i in range(N)[1:-1] for j in range(N)[1:-1] for m in range(M)
            ),
            GRB.MAXIMIZE,
        )
        print("Building TTP")
        Mod.params.TimeLimit = self.runtime
        Mod.params.MIPGap = self.optgap
        Mod.update()

        self.gurobi_model = Mod

    def solve(self):
        """Build and optimize; if a feasible solution exists, post-process."""
        print(f"Solving TTP for {self.N - 2} exposures with Gurobi")
        self.to_gurobi_model()
        Mod = self.gurobi_model
        Mod.optimize()

        if Mod.SolCount > 0:
            self.digest_gurobi()
            self.optimization_status()
        else:
            print(
                "No incumbent solution in time allotted, aborting solve. "
                "Try increasing time_limit parameter."
            )

    # ---------------------------------------------------------- post-process

    def digest_gurobi(self):
        """Interpret the Gurobi solution; populate schedule / plotly / paths.

        Sets attributes ``extras``, ``num_scheduled``, ``estimated_slews``,
        ``real_slews``, ``schedule``, ``plotly``, ``times``, ``az_path``,
        ``alt_path``.
        """
        N = self.N
        M = self.M
        tau_slew = self.tau_slew
        Mod = self.gurobi_model

        num_scheduled = 0
        scheduled_targets = []
        extras, extra_rises, extra_sets = [], [], []
        for i in range(self.N)[1:-1]:
            Yvar = Mod.getVarByName(f"Yi[{i}]")
            if np.round(Yvar.X, 0) == 1:
                num_scheduled += 1
                v = Mod.getVarByName(f"ti[{i}]")
                scheduled_targets.append((i, v.X))
            else:
                row = self.requests_frame.iloc[self.node_to_request[i]]
                extras.append(row.unique_id)
                # Convert back to ISO times for the extras table.
                first_minutes = self.t_early[i]
                last_minutes = self.t_late[i]
                t1 = self.night_start + TimeDelta(first_minutes * 60, format="sec")
                t2 = self.night_start + TimeDelta(last_minutes * 60, format="sec")
                extra_rises.append(str(t1.isot)[11:16])
                extra_sets.append(str(t2.isot)[11:16])

        self.extras = {
            "Starname": extras,
            "First Available": extra_rises,
            "Last Available": extra_sets,
        }
        self.num_scheduled = num_scheduled

        def to_wrap_frame_scalar(angle):
            if self.wrap_limit:
                angle = angle + (360 - self.wrap_limit)
                if angle > 360:
                    angle -= 360
            return angle

        est_slews = []
        for i in range(N)[1:-1]:
            for j in range(N)[1:-1]:
                for m in range(M):
                    var = Mod.getVarByName(f"Xijm[{i},{j},{m}]").X
                    if np.round(var, 0) == 1:
                        est_slews.append(tau_slew[i, j, m])

        real_slews = []
        for i in range(N)[1:-1]:
            for j in range(N)[1:-1]:
                for m in range(M):
                    t = Mod.getVarByName(f"tijm[{i},{j},{m}]").X
                    if np.round(t, 1) != 0:
                        minutes = np.round(t, 1)
                        time_of_slew = self.night_start + TimeDelta(minutes * 60, format="sec")
                        coord_i = self.requests_frame.iloc[self.node_to_request[i]]["coord"]
                        coord_j = self.requests_frame.iloc[self.node_to_request[j]]["coord"]
                        altaz1 = self.observer.altaz(time_of_slew, coord_i)
                        altaz2 = self.observer.altaz(time_of_slew, coord_j)

                        alt1 = altaz1.alt.deg
                        alt2 = altaz2.alt.deg
                        az1 = to_wrap_frame_scalar(altaz1.az.deg)
                        az2 = to_wrap_frame_scalar(altaz2.az.deg)

                        alt_sep = np.abs(alt1 - alt2)
                        az_sep = np.abs(az1 - az2)

                        if not self.wrap_limit and az_sep >= 180:
                            az_sep = 360 - az_sep

                        separation = max(alt_sep, az_sep)
                        slew = separation / (60 * self.slew_rate)
                        real_slews.append(np.round(slew, 3))

        self.estimated_slews = est_slews
        self.real_slews = real_slews

        # Sort by true solver completion time (minutes from night start).
        # Never truncate to int before sorting: many visits can share the same
        # integer minute while their float ``ti`` ordering differs; truncating
        # ties them and yields visit sequences whose *start* times are not
        # monotonic, which makes ``plot_path_2D_interactive`` draw lines that
        # jump backward on the time axis.
        order = np.argsort([p[1] for p in scheduled_targets])
        scheduled_targets = [scheduled_targets[i] for i in order]

        starnames, orders = [], []
        t_starts, t_ends = [], []
        n_shots, priorities_out, exptimes = [], [], []
        ordered_target_nodes = []
        all_times, az_path, alt_path = [], [], []
        for order_idx, pair in enumerate(scheduled_targets):
            node_ind = pair[0]
            ordered_target_nodes.append(node_ind)
            row = self.requests_frame.iloc[self.node_to_request[node_ind]]
            starnames.append(row.unique_id)
            priorities_out.append(int(row.priority))
            n_shots.append(int(row.n_exp))
            exptimes.append(float(row.exptime))

            visit_dur_min = self.t_visit[node_ind]
            t1 = self.night_start + TimeDelta((pair[1] - visit_dur_min) * 60, format="sec")
            t2 = self.night_start + TimeDelta(pair[1] * 60, format="sec")
            t_starts.append(pair[1] - visit_dur_min)
            t_ends.append(pair[1])

            all_times.append(t1)
            all_times.append(t2)
            coord = row["coord"]
            coords_start = self.observer.altaz(t1, coord)
            coords_end = self.observer.altaz(t2, coord)
            az_path.append(coords_start.az.deg)
            alt_path.append(coords_start.alt.deg)
            az_path.append(coords_end.az.deg)
            alt_path.append(coords_end.alt.deg)
            orders.append(order_idx)

        t_starts = np.array(t_starts)
        t_ends = np.array(t_ends)
        rise_times = self.t_early[ordered_target_nodes]
        set_times = self.t_late[ordered_target_nodes]

        # ``schedule['Time']`` is the start-of-exposure time in JD (legacy
        # consumer contract: nplan.drop_gap_rows and plotter rely on it).
        self.schedule = {
            "Order": orders,
            "Starname": starnames,
            "Time": [self.night_start.jd + t / (24 * 60) for t in t_starts],
        }

        self.plotly = {
            "Starname": starnames,
            "First Available": rise_times,
            "Last Available": set_times,
            "Start Exposure": t_starts,
            "Minutes the from Start of the Night": (t_starts + t_ends) / 2,
            "Stop Exposure": t_ends,
            "N_shots": np.array(n_shots),
            "Exposure Time (min)": np.array(exptimes),
            "Total Exp Time (min)": self.t_visit[ordered_target_nodes],
            "Priority": priorities_out,
        }

        self.times = all_times
        self.az_path = az_path
        self.alt_path = alt_path

    def optimization_status(self):
        """Write ``TTPstatistics.txt`` and print solver wall-clock stats.

        NOTE: The objective-bound parse used to define ``obs_bound`` and
        ``slew_bound`` is brittle when priorities are non-uniform; those
        attributes were never read downstream and have been removed.
        """
        slew_param = 1 / 100
        Mod = self.gurobi_model

        # The MILP objective is (priority-weighted visit count) - slew_param * total_slew.
        # When priorities are uniform we can back out the slew time from the fractional
        # part of the objective; for non-uniform priorities we accept the looser estimate.
        obj = Mod.ObjVal
        lower = 1 + (obj // 1)
        slewtime = (lower - obj) * 1 / slew_param

        time_exposing = 0
        for i in range(self.N)[1:-1]:
            Y = Mod.getVarByName(f"Yi[{i}]").x
            time_exposing += Y * self.t_visit[i]
        time_idle = self.dur - time_exposing - slewtime

        self.time_idle = time_idle
        self.time_slewing = slewtime
        self.time_exposing = time_exposing
        self.solve_time = Mod.Runtime

        with open(self.outdir + "/TTPstatistics.txt", "w") as fh:
            fh.write("Stats for TTP Solution\n")
            fh.write("------------------------------------\n")
            fh.write(f"    Model ran for {self.solve_time:.2f} seconds\n")
            fh.write(f"     Observations Requested: {self.N - 2}\n")
            fh.write(f"     Observations Scheduled: {self.num_scheduled}\n")
            fh.write("------------------------------------\n")
            fh.write(f"   Observing Duration (min): {self.dur:.2f}\n")
            fh.write(f"  Time Spent Exposing (min): {self.time_exposing:.2f}\n")
            fh.write(f"      Time Spent Idle (min): {self.time_idle:.2f}\n")
            fh.write(f"   Time Spent Slewing (min): {self.time_slewing:.2f}\n")
            fh.write("------------------------------------\n")

        print("\n------------------------------------")
        print(f"    Model ran for {self.solve_time:.2f} seconds")
        print("------------------------------------")
        print(f"     Observations Requested: {self.N - 2}")
        print(f"     Observations Scheduled: {self.num_scheduled}")
        print("------------------------------------")
        print(f"   Observing Duration (min): {self.dur:.2f}")
        print(f"  Time Spent Exposing (min): {self.time_exposing:.2f}")
        print(f"      Time Spent Idle (min): {self.time_idle:.2f}")
        print(f"   Time Spent Slewing (min): {self.time_slewing:.2f}")
        print("------------------------------------")
