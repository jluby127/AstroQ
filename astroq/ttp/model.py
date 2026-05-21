"""Traveling Telescope Problem (TTP) solver

The model implements the MILP of Handley et al. 2024 (arXiv:2310.18497).
"""

# Standard library imports
import logging
import time
from itertools import product

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

    Keyword Args:
        observer (astroplan.Observer): site fixture for alt/az lookups in
            :meth:`build_arcs`.
        slew_rate (float): mean slew rate, degrees/second.
        wrap_limit (float | None): azimuth wrap limit, degrees. ``None`` means
            no wrap (any az difference > 180 wraps the short way).
        readout_time (float): detector readout between shots of a visit, sec.
        n_slots (int): number of TTP slew slots ``M`` (Handley+ 2024 §2.2).
            ``n_slots=1`` is the recommended default.
        slew_sample_cadence_min (int): max spacing in minutes at which to
            sample arcs within a slot.

    The class is structured around four DataFrames:

    * ``self.requests_frame`` -- one row per request
    * ``self.nodes`` -- one row per MILP node (:meth:`build_nodes`)
    * ``self.arcs`` -- precomputed arc catalog ``(i, j, m)`` (:meth:`build_arcs`)
    * ``self.schedule`` -- post-solve output parallel to ``nodes`` (:meth:`build_schedule`);
      ``None`` if Gurobi finds no incumbent within the time limit

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
    _SLEW_PENALTY = 1 / 100

    def __init__(
        self,
        requests_frame,
        night_start,
        night_end,
        *,
        observer,
        slew_rate,
        wrap_limit=None,
        readout_time=0.0,
        n_slots=1,
        slew_sample_cadence_min=30,
    ):

        # Validate required columns.
        for col in REQUIRED_COLUMNS:
            if col not in requests_frame.columns:
                raise ValueError(f"requests_frame missing required column: {col}")

        # Attach per-row SkyCoord 
        requests_frame["coord"] = SkyCoord(
            requests_frame.ra * u.deg, 
            requests_frame.dec.values * u.deg, 
            frame="icrs",
        )

        self.requests_frame = requests_frame.reset_index(drop=True).copy()
        self.night_start = night_start
        self.night_end = night_end
        self.observer = observer
        self.slew_rate = slew_rate
        self.wrap_limit = wrap_limit
        self.readout_time = readout_time
        self.n_slots = n_slots
        self.slew_sample_cadence_min = slew_sample_cadence_min
        self.M = self.n_slots
        self.schedule = None
        self.selected_arcs = None
        self.stats = {}

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
        """Slew time in minutes for aligned alt/az inputs (broadcasted).

        Returns an ndarray with the broadcast of the four inputs. Any
        reduction (max over time samples, etc.) is the caller's
        responsibility.
        """
        az_sep = self._short_az_sep(
            np.abs(self._wrap_az(az_a) - self._wrap_az(az_b))
        )
        alt_sep = np.abs(
            np.asarray(alt_a, dtype=np.float64)
            - np.asarray(alt_b, dtype=np.float64)
        )
        return np.maximum(az_sep, alt_sep) / (60.0 * float(self.slew_rate))

    def _minutes_from_start(self, t):
        """``astropy.time.Time`` -> minutes since ``self.night_start``."""
        return (t.jd - self.night_start.jd) * 24 * 60

    def _time_from_minutes(self, mins):
        """Minutes since ``self.night_start`` -> ``Time`` (broadcast-friendly)."""
        return self.night_start + TimeDelta(np.asarray(mins) * 60, format="sec")

    # ------------------------------------------------------------- node setup

    def build_nodes(self):
        """Expand the request frame into the ``self.nodes`` DataFrame.

        Sets attributes:
            N, M, dur, nodes, multi_visit_groups.

        All times in ``self.nodes`` are in minutes from ``night_start``.
        """
        self.dur = float(np.round(self._minutes_from_start(self.night_end), 0))

        # Start from a copy of the requests frame and only add/modify the
        reqs = self.requests_frame.copy()
        reqs["request_idx"] = np.arange(len(reqs), dtype=np.int64)
        reqs["t_early"] = self._minutes_from_start(Time(reqs.first_available.tolist()))
        reqs["t_late"] = self._minutes_from_start(Time(reqs.last_available.tolist()))
        reqs["t_visit"] = self._visit_duration(reqs.exptime, reqs.n_exp)
        reqs["tau_intra"] = reqs.tau_intra.astype(float) * 60.0  # hours -> minutes
        reqs["is_anchor"] = False

        # Attach visit_seq via a simple cross-join + filter.
        max_intra = int(reqs.n_intra_max.max())
        visit_seq_table = pd.DataFrame({
            "visit_seq": np.arange(max_intra, dtype=np.int64),
        })
        visits = (
            reqs.merge(visit_seq_table, how="cross")
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
            "t_late": self.dur,
            "t_visit": 0.0,
            "tau_intra": 0.0,
            "priority": 0,
            "coord": None,
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

    def build_arcs(self, *, samples_per_slot=3):
        """Build precomputed arc catalog ``self.arcs``.

        Within slot ``m``, sample alt/az for every real (non-anchor) node at
        least every ``slew_sample_cadence_min`` minutes (and at least
        ``samples_per_slot`` times per slot; default 3); the worst-case slew
        between any two real nodes is the maximum of ``max|delta_alt|`` and
        ``max|delta_az|`` (wrap-aware), divided by ``slew_rate``.

        Args:
            samples_per_slot (int): floor on the number of temporal samples per
                slew slot; the used count is ``int(max(...))`` of this and the
                value implied by ``dur``, ``self.M``, and
                ``slew_sample_cadence_min``.

        """
        n_samples = int(max(
            self.dur / (self.M * self.slew_sample_cadence_min),
            samples_per_slot,
        ))

        slot_bounds = Time(
            np.linspace(self.night_start.jd, self.night_end.jd, self.M + 1, endpoint=True),
            format="jd",
        )
        # Sample grid: shape (self.M, n_samples) -> flatten for one altaz call.
        sample_jd = np.array([
            np.linspace(slot_bounds[m].jd, slot_bounds[m + 1].jd, n_samples)
            for m in range(self.M)
        ])
        times = Time(sample_jd.ravel(), format="jd")
        altaz = self.observer.altaz(
            times, 
            self.nodes[~self.nodes.is_anchor].coord.to_list(),
            grid_times_targets=True
        )

        # Visit-sample table.
        id, m, sample = np.mgrid[
            1:self.N - 1, 0:self.M, 0:n_samples
        ]
        visit_samples = pd.DataFrame({
            "id": id.reshape(-1),
            "m":       m.reshape(-1),
            "sample":  sample.reshape(-1),
            "alt":     altaz.alt.deg.reshape(-1),
            "az":      altaz.az.deg.reshape(-1),
        })

        # Self-merge on (m, sample) to pair every visit with every other visit
        # drop self-pairs. 
        node_samples = (
            visit_samples
            .merge(visit_samples, on=["m", "sample"], suffixes=("_i", "_j"))
            .query("id_i != id_j")
            .rename(columns={"id_i": "i", "id_j": "j"})
            [["i", "j", "m", "sample", "alt_i", "az_i", "alt_j", "az_j"]]
        )

        node_samples["tau_sample"] = self._slew_minutes(
            node_samples["alt_i"].to_numpy(),
            node_samples["alt_j"].to_numpy(),
            node_samples["az_i"].to_numpy(),
            node_samples["az_j"].to_numpy(),
        )

        agg = node_samples.groupby(["i", "j", "m"], sort=False)["tau_sample"].max()

        uid = self.nodes["unique_id"].to_numpy()
        i_lev = agg.index.get_level_values("i").to_numpy(dtype=np.int64)
        j_lev = agg.index.get_level_values("j").to_numpy(dtype=np.int64)
        self.arcs = pd.DataFrame(
            {
                "i_id": uid[i_lev],
                "j_id": uid[j_lev],
                "t_slew": agg.to_numpy(dtype=np.float64),
            },
            index=agg.index,
        )

        # Slot bounds expressed as minutes from start.
        self.w = (slot_bounds.jd - slot_bounds[0].jd) * 24 * 60

    # ------------------------------------------------------------- MILP build
    def build_model(self):
        """Construct the TTP MILP (Handley+ 2024, eqs. 2-9, B3, 10).

        Constraints are added in the order they appear in Handley+ 2024 §2.4.
        """
        self.model = gp.Model("TTP")

        N, M = self.N, self.M
        nodes = self.nodes

        # O(1) slew lookup for hot loops (anchor arcs absent => .get(..., 0.0)).
        arcs_lookup = self.arcs["t_slew"].to_dict()

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
        # node). Anchor arcs (i=0) are MILP bookkeeping with t_slew = 0; the
        # .get(..., 0.0) fallback keeps arcs restricted to real arcs.
        for j in range(1, N):
            t_visit_j = nodes.at[j, "t_visit"]
            self.model.addConstr(
                self.ti[j] >= gp.quicksum(
                    self.tijm[i, j, m]
                    + (arcs_lookup.get((i, j, m), 0.0) + t_visit_j) * self.Xijm[i, j, m]
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
                arcs_lookup.get((i, j, m), 0.0) * self.Xijm[i, j, m]
                for i in range(1, N - 1)
                for j in range(1, N - 1)
                for m in range(M)
            ),
            GRB.MAXIMIZE,
        )

        self.model.update()

    def run_model(self):
        """Solve the MILP and build ``schedule`` / ``stats``.

        On success, ``schedule`` is a DataFrame and ``stats`` is populated.
        If Gurobi has no incumbent (``SolCount == 0``), logs a warning and
        leaves ``schedule`` as ``None``.
        """
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
        """Walk the chosen path and populate ``schedule`` + ``selected_arcs``.

        ``schedule`` is parallel to ``nodes`` (same index = node id ``i``),
        with ``scheduled``, ``t_slew``, and solve-time columns added.
        """
        # extract selected arcs from gurobi solution
        arcs_selected = []
        for (i, j, m), var in self.Xijm.items():
            if var.X > 0.5 and j != 0 and i != self.N - 1:
                arcs_selected.append({"i": i, "j": j, "m": m, "ti":self.ti[i].X})
        arcs_selected = pd.DataFrame(arcs_selected)

        # merge selected arcs with nodes, unvisited nodes will have NaN for ti
        schedule = pd.merge(
            self.nodes.query('~is_anchor'), 
            arcs_selected, 
            left_index=True, 
            right_on=["i"],
            how="left"
        )
        schedule['t_start'] = schedule['ti'] - schedule['t_visit']
        schedule['t_end'] = schedule['ti']
        schedule['scheduled'] = ~schedule['ti'].isna()

        schedule = pd.merge(
            schedule,
            self.arcs['t_slew'],
            left_on=["i","j","m"],
            right_index=True,
            how="left"
        ).sort_values(by='t_start',na_position='last')
        schedule['order'] = range(len(schedule))

        # Drop the SkyCoord cache so the schedule round-trips cleanly through
        # to_csv / to_hdf without object-dtype hazards. Plot adapters rebuild
        # coords from ra/dec on demand; the solver no longer needs `coord`.
        self.schedule = schedule.drop(columns=['coord'], errors='ignore')
        scheduled = self.schedule[self.schedule['scheduled']]
        self.stats = {
            "dur": self.dur,
            "n_requested": self.N - 2,
            "n_scheduled": len(scheduled),
            "t_visit_sum": scheduled['t_visit'].sum(),
            "t_slew_sum": scheduled['t_slew'].sum(),
            "t_idle_sum": self.dur - scheduled['t_visit'].sum() - scheduled['t_slew'].sum(),
        }

    def to_string(self, *, header="Stats for TTP Solution"):
        """Return a human-readable summary of the solve from ``self.stats``."""
        s = self.stats
        lines = [
            header,
            "------------------------------------",
            f"     Observations Requested: {s['n_requested']}",
            f"     Observations Scheduled: {s['n_scheduled']}",
            "------------------------------------",
            f"   Observing Duration (min): {s['dur']:.2f}",
            f"  Time Spent Exposing (min): {s['t_visit_sum']:.2f}",
            f"      Time Spent Idle (min): {s['t_idle_sum']:.2f}",
            f"   Time Spent Slewing (min): {s['t_slew_sum']:.2f}",
            "------------------------------------",
        ]
        return "\n".join(lines) + "\n"
