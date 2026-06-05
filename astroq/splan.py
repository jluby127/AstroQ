"""
Module that defines the SemesterPlanner class. This class is responsible for defining, building, and solving the
Gurobi model for semester-level observation planning. It is nearly completely agnostic to all astronomy knowledge.
"""

# Standard library imports
import json
import logging
import os
import time
import warnings
from configparser import ConfigParser
from datetime import datetime
from pathlib import Path

# Third-party imports
import gurobipy as gp
import h5py
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
from gurobipy import GRB
from jinja2 import Template

# Local imports
import astroq.access as ac
import astroq.history as hs
import astroq.queue

warnings.filterwarnings("ignore")

logs = logging.getLogger(__name__)

# Persistence schema for semester_planner.h5. Bump when the on-disk layout
# changes; from_hdf5 hard-fails on mismatch (mirrors nplan.py).
SEMESTER_PLANNER_H5_SCHEMA = 1


def _resolve_path(workdir, raw):
    """Resolve a possibly-relative config path against ``workdir``."""
    return raw if os.path.isabs(raw) else os.path.join(workdir, raw)


def _serialize_past_history(past_history):
    """Flatten ``{uid: StarHistory}`` to a JSON-safe dict."""
    return {
        uid: {
            "name": h.name,
            "date_last_observed": h.date_last_observed,
            "total_n_exposures": h.total_n_exposures,
            "total_n_visits": h.total_n_visits,
            "total_n_unique_nights": h.total_n_unique_nights,
            "total_open_shutter_time": h.total_open_shutter_time,
            "n_obs_on_nights": h.n_obs_on_nights,
            "n_visits_on_nights": h.n_visits_on_nights,
        }
        for uid, h in past_history.items()
    }


def _deserialize_past_history(data):
    """Inverse of :func:`_serialize_past_history`."""
    return {uid: hs.StarHistory(**fields) for uid, fields in data.items()}


class SemesterPlanner(object):
    """Semester-level scheduler: pick which targets get observed when.

    Formulates the cadenced-scheduling problem of Lubin et al. 2025
    (arXiv:2506.08195) as a Gurobi MILP over a (request, day, slot) grid
    and produces a sparse schedule for the entire semester. The night-level
    slew ordering is handled separately by :class:`astroq.nplan.NightPlanner`.

    Inputs (resolved relative to ``[global] workdir`` in the config):
        - ``request.csv``  — observing requests, one row per target.
        - ``allocation.csv`` — telescope time blocks for the semester.
        - ``past.csv`` — prior observations (caps future ``n_inter_max``).
        - ``custom.csv`` — PI-supplied per-target observability windows.
        - ``programs.csv`` — awarded nights per program (drives throttling).

    Key outputs (written to ``<workdir>/outputs``):
        - ``semester_plan.csv`` — sparse schedule with columns
          ``unique_id, d, s, starname``
        - ``request_selected.csv`` — tonight's targets, the handoff to
          :class:`astroq.nplan.NightPlanner`.
        - ``runReport.txt`` — human-readable per-program statistics.
        - ``semester_planner.h5`` — round-trippable snapshot consumed by
          :class:`astroq.nplan.NightPlanner` and the plotting layer.

    Lifecycle:

        >>> sp = SemesterPlanner("config.ini", run_band3=False)
        >>> sp.run_model()       # builds constraints, solves, writes outputs

    Persistence:
        :meth:`to_hdf5` stores the config text plus a handful of DataFrames
        and the precomputed ``access_record``; :meth:`from_hdf5` rehydrates
        a planner suitable for downstream consumers (no Gurobi state). The
        on-disk schema version is :data:`SEMESTER_PLANNER_H5_SCHEMA`.

    Args:
        cf (str): path to the ``config.ini`` file.
        run_band3 (bool): if True, schedule the band-3 filler program instead
            of the primary request list. Reads ``[data] filler_file`` and
            applies a 20-minute nautical-twilight buffer to ``allocation.csv``.
    """

    def __init__(self, cf, run_band3):
        """
        Initialize the SemesterPlanner object.

        Args:
            cf (str): the path to the config.ini file
            run_band3 (bool): whether to run the band 3 weather loss model
                (this will be unnecessary in the 2026A semester)
        """
        logs.debug("Building the SemesterPlanner.")
        self.start_the_clock = time.time()

        # Read config as text so we can persist it verbatim and recreate the
        # parser on from_hdf5. All scalar config values are exposed as
        # @property below; none are duplicated as instance attributes.
        self._config_ini_text = Path(cf).read_text()
        self.config = ConfigParser()
        self.config.read_string(self._config_ini_text)
        self.run_band3 = run_band3

        # Queue owns telescope/instrument knowledge; cached once.
        self.queue = astroq.queue.from_config(self.config)
        self.queue_name = self.config.get("global", "queue")

        os.makedirs(self.output_directory, exist_ok=True)

        # Band-3 fallback: extend allocation by 20 min into nautical twilight.
        if self.run_band3:
            self.add_twilights()

        # Active vs. all requests + data hygiene.
        self.requests_frame_all, self.requests_frame = self._load_requests_frame()

        # Per-request derived data.
        self.strategy = self._build_strategy()
        self.past_history = hs.process_star_history(
            _resolve_path(
                self.semester_directory, self.config.get("data", "past_file")
            )
        )
        self.slots_needed_for_exposure_dict = self._build_slots_required_dictionary()

        # Observability cube (single source of truth for which slots are valid).
        self.access_obj = ac.Access.from_planner(self)
        self.access_record = self.access_obj.build_access()
        self.observability = self.access_obj.observability(
            self.access_record.is_observable
        )

        # Pre-computed aggregations consumed by constraint methods. These exist
        # so the constraint expressions stay compact and aligned with Lubin
        # et al. notation -- keep them.
        self._build_constraint_lookups()

        # Gurobi model + variables.
        self.model = gp.Model("Semester_Scheduler")
        self._build_gurobi_vars()
        self._build_max_obs_dicts()

        logs.debug("Initializing complete.")

    # ------------------------------------------------------------------
    # Config-derived properties (externally consumed -- keep clean access).
    # ------------------------------------------------------------------

    @property
    def slot_size(self):
        return self.config.getint("semester", "slot_size")

    @property
    def current_day(self):
        return self.config.get("global", "current_day")

    @property
    def semester_directory(self):
        return self.config.get("global", "workdir")

    @property
    def output_directory(self):
        return os.path.join(self.semester_directory, "outputs")

    @property
    def semester_start_date(self):
        return self.config.get("global", "semester_start_day")

    @property
    def semester_length(self):
        start = datetime.strptime(self.semester_start_date, "%Y-%m-%d")
        end = datetime.strptime(
            self.config.get("global", "semester_end_day"), "%Y-%m-%d"
        )
        return int((end - start).days + 1)

    @property
    def semester_letter(self):
        return self.config.get("global", "semester")[-1]

    @property
    def hours_per_night(self):
        return self.config.getfloat("semester", "hours_per_night")

    @property
    def throttle_grace(self):
        return self.config.getfloat("semester", "throttle_grace")

    @property
    def allocation_file(self):
        return _resolve_path(
            self.semester_directory, self.config.get("data", "allocation_file")
        )

    @property
    def custom_file(self):
        return _resolve_path(
            self.semester_directory, self.config.get("data", "custom_file")
        )

    @property
    def programs_file(self):
        return _resolve_path(
            self.semester_directory, self.config.get("data", "programs_file")
        )

    @property
    def weather_loss_file(self):
        v = self.config.get("semester", "weather_loss_file", fallback=None)
        return v or None

    @property
    def run_weather_loss(self):
        return self.config.getboolean("semester", "run_weather_loss")

    # ------------------------------------------------------------------
    # Derived properties (computed from the above, not raw config).
    # ------------------------------------------------------------------

    @property
    def all_dates_array(self):
        return self.access_obj.all_dates_array

    @property
    def all_dates_dict(self):
        return self.access_obj.all_dates_dict

    @property
    def n_slots_in_night(self):
        return int(24 * 60 / self.slot_size)

    @property
    def n_nights_in_semester(self):
        return len(self.all_dates_dict) - self.all_dates_dict[self.current_day]

    @property
    def n_slots_in_semester(self):
        return self.n_slots_in_night * self.n_nights_in_semester

    @property
    def today_starting_night(self):
        return self.all_dates_dict[self.current_day]

    @property
    def today_starting_slot(self):
        return self.all_dates_dict[self.current_day] * self.n_slots_in_night

    # ------------------------------------------------------------------
    # Construction helpers (called from __init__ and from_hdf5).
    # ------------------------------------------------------------------

    def _load_requests_frame(self):
        """Read request.csv (or filler when run_band3), clean, validate."""
        key = "filler_file" if self.run_band3 else "request_file"
        request_file = _resolve_path(
            self.semester_directory, self.config.get("data", key)
        )
        if not os.path.exists(request_file):
            raise FileNotFoundError(f"Requests file not found: {request_file}")

        rfa = pd.read_csv(request_file)
        if "comments" not in rfa.columns:
            rfa["comments"] = ""
        return self._split_and_clean_requests(rfa, source=request_file)

    @staticmethod
    def _split_and_clean_requests(rfa, source="<hdf5>"):
        """Split into active/all frames, apply legacy data hygiene, validate.

        Idempotent: safe to call on already-cleaned frames (e.g. from HDF5).
        """
        mask = rfa["inactive"] == False  # noqa: E712 -- explicit bool match
        logs.warning(
            f"There are {len(rfa[~mask])} inactive of {len(rfa)} requests."
        )
        rf = rfa[mask].reset_index(drop=True).copy()

        # Tolerate "None" strings left over from the early-2025B webform.
        for col, default in (("n_intra_max", 1), ("n_intra_min", 1), ("tau_intra", 0)):
            rf[col] = rf[col].replace("None", np.nan).fillna(default)
        for band in (1, 2, 3):
            col = f"weather_band_{band}"
            if col in rf.columns:
                rf[col] = rf[col].replace("None", np.nan).fillna(False)
        rf["unique_id"] = rf["unique_id"].astype(str)
        rf["starname"] = rf["starname"].astype(str)

        # Fail loudly on duplicate unique_id -- otherwise gurobipy raises a
        # cryptic "Duplicate keys in Model.addVars()" later.
        dup_mask = rf["unique_id"].duplicated(keep=False)
        if dup_mask.any():
            dup_ids = sorted(rf.loc[dup_mask, "unique_id"].unique())
            raise ValueError(
                f"Duplicate unique_id among active requests in {source!r}: "
                f"{dup_ids}. Remove or merge duplicate rows so each active "
                f"request has one row."
            )
        return rfa, rf

    def _build_strategy(self):
        """Per-request strategy with seconds/hours converted to slot units."""
        rf = self.requests_frame
        strategy = rf[
            [
                "starname",
                "unique_id",
                "n_intra_min",
                "n_intra_max",
                "n_inter_max",
                "tau_inter",
            ]
        ].copy()
        strategy["t_visit"] = (
            (rf["exptime"] / 60 / self.slot_size).clip(lower=1).round().astype(int)
        )
        strategy["tau_intra"] = (
            (rf["tau_intra"] * 60 / self.slot_size).round().astype(int)
        )
        return strategy

    def _build_slots_required_dictionary(self):
        """Slots per visit per request, via :meth:`Queue.visit_seconds`."""
        slot_seconds = self.slot_size * 60.0
        out = {}
        for _, row in self.requests_frame.iterrows():
            total_s = self.queue.visit_seconds(
                float(row["exptime"]),
                int(row["n_exp"]),
                int(row["n_intra_max"]),
            )
            out[row["unique_id"]] = max(1, int(np.round(total_s / slot_seconds)))
        return out

    def _build_constraint_lookups(self):
        """Build aggregation tables consumed by the constraint methods."""
        self.observability_tuples = list(
            self.observability.itertuples(index=False, name=None)
        )
        self.joiner = pd.merge(self.strategy, self.observability, on=["unique_id"])

        self.observability_nights = (
            self.joiner.loc[self.joiner["n_intra_max"] > 1, ["unique_id", "d"]]
            .drop_duplicates()
            .copy()
        )
        self.multi_visit_requests = list(
            self.observability_nights["unique_id"].unique()
        )
        self.all_requests = list(self.requests_frame["unique_id"])
        self.schedulable_requests = list(self.joiner["unique_id"].unique())
        self.single_visit_requests = [
            uid for uid in self.schedulable_requests
            if uid not in self.multi_visit_requests
        ]

        warncount = sum(
            uid not in self.schedulable_requests for uid in self.all_requests
        )
        logs.warning(
            f"There are {warncount} targets out of {len(self.all_requests)} "
            f"that have no valid day/slot pairs and therefore are effectively "
            f"removed from the model."
        )

        self.all_valid_ds_for_request = (
            self.joiner.groupby(["unique_id"])[["d", "s"]].agg(list)
        )
        valid_s_for_rd = pd.merge(
            self.joiner.drop_duplicates(["unique_id", "d"]),
            self.joiner[["unique_id", "d", "s"]],
            suffixes=["", "3"],
            on=["unique_id"],
        ).query("d == d3")
        self.slots_on_day_for_r = (
            valid_s_for_rd.groupby(["unique_id", "d"])[["s3"]].agg(list)
        )

    def _build_gurobi_vars(self):
        """Build ``Yrds``, ``Wrd``, ``theta`` on ``self.model``."""
        self.Yrds = self.model.addVars(
            self.observability_tuples, vtype=GRB.BINARY, name="Requests_Slots"
        )
        if len(self.observability_nights) != 0:
            self.Wrd = self.model.addVars(
                list(self.observability_nights.itertuples(index=False, name=None)),
                vtype=GRB.BINARY,
                name="OnSky",
            )
        self.theta = self.model.addVars(self.all_requests, name="Shortfall")

    def _build_max_obs_dicts(self):
        """Past nights observed + per-request desired/absolute max-obs bounds."""
        max_bonus = self.config.getfloat("semester", "maximum_bonus_size")
        rf = self.requests_frame.set_index("unique_id")
        past_nights_observed_dict = {}
        desired_max_obs_allowed_dict = {}
        absolute_max_obs_allowed_dict = {}
        # absolute_max_obs is conditionally assigned inside the loop; this
        # initializer matches the pre-refactor behavior where a value would
        # have carried over from an earlier iteration.
        absolute_max_obs = 0
        for uid in self.all_requests:
            n_inter_max = int(rf.loc[uid, "n_inter_max"])
            past_obs = (
                self.past_history[uid].total_n_unique_nights
                if uid in self.past_history else 0
            )
            if past_obs > n_inter_max:
                desired_max_obs = past_obs
            else:
                desired_max_obs = n_inter_max - past_obs
                absolute_max_obs = desired_max_obs + int(n_inter_max * max_bonus)
                if past_obs > absolute_max_obs:
                    absolute_max_obs = past_obs
            past_nights_observed_dict[uid] = past_obs
            desired_max_obs_allowed_dict[uid] = desired_max_obs
            absolute_max_obs_allowed_dict[uid] = absolute_max_obs
        self.past_nights_observed_dict = past_nights_observed_dict
        self.desired_max_obs_allowed_dict = desired_max_obs_allowed_dict
        self.absolute_max_obs_allowed_dict = absolute_max_obs_allowed_dict

    # ==================================================================
    # Constraints (grouped by purpose).
    # ==================================================================

    # ---- variable/shortfall definition ----

    def constraint_build_theta_multivisit(self):
        """
        See Equation 3 in Lubin et al. 2025.

        Definition of the "shortfall" matrix, Theta. The shortfall is defined
        for each target, giving the difference between the requested number of
        nights and the sum of past and future scheduled observations.
        """
        logs.info("Constraint: Build theta variable")
        for starid in self.schedulable_requests:
            idx = self.requests_frame.index[
                self.requests_frame["unique_id"] == starid
            ][0]
            self.model.addConstr(
                self.theta[starid] >= 0, f"greater_than_zero_shortfall_{starid}"
            )
            available = list(
                zip(
                    self.all_valid_ds_for_request.loc[starid].d,
                    self.all_valid_ds_for_request.loc[starid].s,
                )
            )
            rhs = (
                self.requests_frame["n_inter_max"][idx]
                - self.past_nights_observed_dict[starid]
                - (
                    gp.quicksum(self.Yrds[starid, d, s] for d, s in available)
                ) / self.requests_frame["n_intra_max"][idx]
            )
            self.model.addConstr(
                self.theta[starid] >= rhs,
                f"greater_than_nobs_shortfall_{starid}",
            )

    # ---- per-slot reservation ----

    def constraint_reserve_multislot_exposures(self):
        """
        See Constraint 1 in Lubin et al. 2025.

        Reserve multiple time slots for exposures that require more than one
        time slot to complete, ensuring no other observations are scheduled
        during these slots.
        """
        logs.info("Constraint: Reserve slots for multi-slot exposures.")
        strategy = self.strategy
        max_t_visit = strategy.t_visit.max()
        R_ds = (
            self.observability.groupby(["d", "s"])["unique_id"].apply(set).to_dict()
        )
        R_geq_t_visit = {
            t: set(strategy.loc[strategy.t_visit >= t, "unique_id"])
            for t in range(1, max_t_visit + 1)
        }

        for d, s in self.observability.drop_duplicates(["d", "s"])[
            ["d", "s"]
        ].itertuples(index=False, name=None):
            rhs = []
            for delta in range(1, max_t_visit):
                s_shift = s - delta
                if (d, s_shift) in R_ds:
                    rhs.extend(
                        self.Yrds[r, d, s_shift]
                        for r in R_ds[d, s_shift] & R_geq_t_visit[delta + 1]
                    )
            lhs = 1 - gp.quicksum(self.Yrds[r, d, s] for r in R_ds[d, s])
            self.model.addConstr(
                lhs >= gp.quicksum(rhs), f"reserve_multislot_{d}d_{s}s"
            )

    # ---- cadence ----

    def constraint_enforce_internight_cadence(self):
        """
        See Constraint 3 in Lubin et al. 2025.

        Ensure that the minimum number of days pass between consecutive
        observations of a given target.
        """
        logs.info("Constraint: Enforce inter-night cadence.")
        # Defensive: coerce tau_inter to numeric without mutating self.joiner.
        tau_inter_numeric = pd.to_numeric(self.joiner["tau_inter"], errors="coerce")
        joiner_local = self.joiner.assign(tau_inter=tau_inter_numeric)

        intercadence = pd.merge(
            joiner_local.drop_duplicates(["unique_id", "d"]),
            joiner_local[["unique_id", "d", "s"]],
            suffixes=["", "3"],
            on=["unique_id"],
        ).query("d + 0 < d3 < d + tau_inter")
        intercadence_tracker = intercadence.groupby(["unique_id", "d"])[
            ["d3", "s3"]
        ].agg(list)

        # Inter-night cadence of 1 day has no forbidden future slots; skip
        # those rows and drop the duplicates-per-day rows up front.
        valid = self.joiner[self.joiner["tau_inter"] > 1].drop_duplicates(
            subset=["unique_id", "d"]
        )
        for _, row in valid.iterrows():
            constrained_slots_tonight = np.array(
                self.slots_on_day_for_r.loc[(row.unique_id, row.d)][0]
            )
            if (row.unique_id, row.d) not in intercadence_tracker.index:
                continue
            future = intercadence_tracker.loc[(row.unique_id, row.d)]
            ds_pairs = zip(
                np.array(future.d3).flatten(),
                np.array(future.s3).flatten(),
            )
            lhs = (
                gp.quicksum(
                    self.Yrds[row.unique_id, row.d, s2]
                    for s2 in constrained_slots_tonight
                )
                / row.n_intra_max
            )
            rhs = 1 - gp.quicksum(
                self.Yrds[row.unique_id, d3, s3] for d3, s3 in ds_pairs
            )
            self.model.addConstr(
                lhs <= rhs,
                f"enforce_internight_cadence_{row.unique_id}_{row.d}d_{row.s}s",
            )

    def constraint_build_enforce_intranight_cadence(self):
        """
        Constraint 4 in Lubin et al. 2025.

        Ensure that the minimum number of hours pass between consecutive
        observations of a given target on the same night.
        """
        logs.info("Constraint: Enforce intra-night cadence.")
        # Intra-night cadence of 0 has nothing to forbid; restrict up front.
        valid = self.joiner[self.joiner["n_intra_max"] > 1]
        intracadence_frame = pd.merge(
            valid.drop_duplicates(["unique_id", "d", "s"]),
            valid[["unique_id", "d", "s"]],
            suffixes=["", "3"],
            on=["unique_id", "d"],
        ).query("s + 0 < s3 < s + tau_intra")
        intracadence_frame = intracadence_frame.groupby(
            ["unique_id", "d", "s"]
        )[["s3"]].agg(list)

        for _, row in valid.iterrows():
            key = (row.unique_id, row.d, row.s)
            if key not in intracadence_frame.index:
                continue
            slots_to_constrain = list(intracadence_frame.loc[key][0])
            lhs = self.Yrds[row.unique_id, row.d, row.s]
            rhs = self.Wrd[row.unique_id, row.d] - gp.quicksum(
                self.Yrds[row.unique_id, row.d, s3] for s3 in slots_to_constrain
            )
            self.model.addConstr(
                lhs <= rhs,
                f"enforce_intranight_cadence_{row.unique_id}_{row.d}d_{row.s}s",
            )

    # ---- visit count bounds ----

    def constraint_set_max_desired_unique_nights_Wrd(self):
        """
        See Constraint 2 in Lubin et al. 2025.

        Limit the number of observations scheduled for a given target to the
        maximum value provided by the PI. This constraint may later be relaxed
        if Round 2 of scheduling is invoked.
        """
        logs.info("Constraint: Set desired maximum observations.")
        for name in self.multi_visit_requests:
            all_d = list(set(self.all_valid_ds_for_request.loc[name].d))
            self.model.addConstr(
                gp.quicksum(self.Wrd[name, d] for d in all_d)
                <= self.desired_max_obs_allowed_dict[name],
                f"max_desired_unique_nights_for_request_{name}",
            )
        for name in self.single_visit_requests:
            available = list(
                zip(
                    self.all_valid_ds_for_request.loc[name].d,
                    self.all_valid_ds_for_request.loc[name].s,
                )
            )
            self.model.addConstr(
                gp.quicksum(self.Yrds[name, d, s] for d, s in available)
                <= self.desired_max_obs_allowed_dict[name],
                f"max_desired_unique_nights_for_request_{name}",
            )

    def remove_constraint_set_max_desired_unique_nights_Wrd(self):
        """
        Bonus round: not in Lubin et al. 2025.

        Remove the maximum number of observations set by
        :meth:`constraint_set_max_desired_unique_nights_Wrd`.
        """
        logs.info("Constraint: Removing previous maximum observations constraint.")
        for name in self.multi_visit_requests:
            rm_const = self.model.getConstrByName(
                f"max_desired_unique_nights_for_request_{name}"
            )
            self.model.remove(rm_const)

    def constraint_set_max_absolute_unique_nights_Wrd(self):
        """
        Bonus round: not in Lubin et al. 2025.

        Set the maximum number of observations for a target to 150% of the
        original requested number.
        """
        logs.info("Constraint: Set absolute maximum observations.")
        for name in self.multi_visit_requests:
            all_d = list(set(self.all_valid_ds_for_request.loc[name].d))
            self.model.addConstr(
                gp.quicksum(self.Wrd[name, d] for d in all_d)
                <= self.absolute_max_obs_allowed_dict[name],
                f"max_absolute_unique_nights_for_request_{name}",
            )

    def constraint_set_min_max_visits_per_night(self):
        """
        See Constraint 5 in Lubin et al. 2025.

        Require that the number of scheduled visits to a target in a given
        night falls between the minimum and maximum values supplied by the PI.
        """
        logs.info("Constraint: Bound minimum and maximum visits per night.")
        per_day = self.joiner.drop_duplicates(subset=["unique_id", "d"])
        grouped_s = (
            self.joiner.groupby(["unique_id", "d"])["s"].unique().reset_index()
        )
        grouped_s.set_index(["unique_id", "d"], inplace=True)
        for _, row in per_day.iterrows():
            slots_tonight = list(grouped_s.loc[(row.unique_id, row.d)]["s"])
            name_tag = f"{row.unique_id}_{row.d}d_{row.s}s"
            visits_tonight = gp.quicksum(
                self.Yrds[row.unique_id, row.d, s3] for s3 in slots_tonight
            )
            if row.unique_id in self.multi_visit_requests:
                self.model.addConstr(
                    visits_tonight <= row.n_intra_max * self.Wrd[row.unique_id, row.d],
                    f"enforce_max_visits1_{name_tag}",
                )
                self.model.addConstr(
                    visits_tonight >= row.n_intra_min * self.Wrd[row.unique_id, row.d],
                    f"enforce_min_visits_{name_tag}",
                )
            else:
                self.model.addConstr(
                    visits_tonight <= row.n_intra_max,
                    f"enforce_min_visits_{name_tag}",
                )

    # ---- throttling & bonus round ----

    def constraint_throttle(self):
        """
        Not described in Lubin et al. 2025.

        Ensure that no program is scheduled for more time than they bring to
        the queue (within a grace amount).
        """
        logs.info("Constraint: Throttling over-requested programs.")
        program_frame = pd.read_csv(self.programs_file)

        # Map uid -> slots already burned in past observations.
        past_slots = {
            uid: h.total_n_exposures
            * self.slots_needed_for_exposure_dict.get(uid, 1)
            for uid, h in self.past_history.items()
        }

        merged_df = self.requests_frame[["unique_id", "program_code"]].copy()
        merged_df["past_slots_used"] = (
            merged_df["unique_id"].astype(str).map(past_slots).fillna(0)
        )
        merged_df = merged_df.merge(
            program_frame[["program", "nights"]],
            left_on="program_code",
            right_on="program",
            how="inner",
        )

        past_used_slots_by_program = (
            merged_df.groupby("program")["past_slots_used"].sum().to_dict()
        )
        program_requests_map = (
            merged_df.groupby("program")["unique_id"].apply(set).to_dict()
        )

        throttle_grace = self.throttle_grace
        for program in program_frame["program"]:
            row = program_frame.loc[program_frame["program"] == program].iloc[0]
            awarded_slots = (
                row["nights"] * self.hours_per_night * 60
            ) / self.slot_size
            awarded_slots_grace = int(awarded_slots * throttle_grace)

            uids_for_program = program_requests_map.get(program, set())
            schedulable_slots = gp.quicksum(
                self.Yrds[r, d, s] * self.slots_needed_for_exposure_dict[r]
                for r, d, s in self.observability_tuples
                if r in uids_for_program
            )

            past_used = past_used_slots_by_program.get(program, 0)
            if awarded_slots_grace < past_used:
                logs.warning(
                    f"Program {program} has already been over-observed. "
                    f"Setting award equal to past used."
                )
                logs.warning(
                    f"Therefore, Program {program}, will not be scheduled "
                    f"for any additional observations."
                )
                awarded_slots_grace = past_used

            self.model.addConstr(
                awarded_slots_grace - past_used >= schedulable_slots,
                f"throttle_program_{program}",
            )

    def constraint_fix_previous_objective(self, epsilon=0.03):
        """
        Bonus round: not in Lubin et al. 2025.

        Ensure that the Round-2 objective is within ``epsilon`` of Round-1.
        """
        logs.info("Constraint: Fixing the previous solution's objective value.")
        self.model.addConstr(
            gp.quicksum(self.theta[name] for name in self.requests_frame["unique_id"])
            <= self.model.objval + epsilon,
            "fix_previous_objective",
        )

    # ==================================================================
    # Objectives.
    # ==================================================================

    def set_objective_minimize_theta_time_normalized(self):
        """See Equation 1 in Lubin et al. 2025."""
        self.model.setObjective(
            gp.quicksum(
                self.theta[name] * self.slots_needed_for_exposure_dict[name]
                for name in self.schedulable_requests
            ),
            GRB.MINIMIZE,
        )

    def set_objective_maximize_slots_used(self):
        """Bonus round: maximize filled slots."""
        logs.info("Objective: Maximize the number of slots used.")
        self.model.setObjective(
            gp.quicksum(
                self.slots_needed_for_exposure_dict[uid] * self.Yrds[uid, d, s]
                for uid, d, s in self.observability_tuples
            ),
            GRB.MAXIMIZE,
        )

    # ==================================================================
    # Model orchestration.
    # ==================================================================

    def build_model_round1(self):
        """Round 1 constraints + objective per Lubin et al. 2025."""
        t1 = time.time()
        self.constraint_reserve_multislot_exposures()
        self.constraint_enforce_internight_cadence()
        self.constraint_set_max_desired_unique_nights_Wrd()
        self.constraint_build_enforce_intranight_cadence()
        self.constraint_set_min_max_visits_per_night()
        self.constraint_build_theta_multivisit()
        self.constraint_throttle()
        self.set_objective_minimize_theta_time_normalized()
        logs.info(f"Time to build constraints: {np.round(time.time() - t1, 3):.3f}")

    def build_model_round2(self):
        """Round 2 constraints + objective (bonus round)."""
        t1 = time.time()
        self.remove_constraint_set_max_desired_unique_nights_Wrd()
        self.constraint_set_max_absolute_unique_nights_Wrd()
        self.constraint_fix_previous_objective()
        self.set_objective_maximize_slots_used()
        logs.info(f"Time to build constraints: {np.round(time.time() - t1, 3):.3f}")

    def optimize_model(self):
        """Solve the Gurobi model (with IIS diagnostics on infeasibility)."""
        logs.debug("Begin model solve.")
        t1 = time.time()
        self.model.params.TimeLimit = self.config.getint("semester", "max_solve_time")
        self.model.Params.OutputFlag = self.config.getboolean(
            "semester", "show_gurobi_output"
        )
        # Allow stop at solver gap to prevent spending time on marginal gains.
        self.model.params.MIPGap = self.config.getfloat("semester", "max_solve_gap")
        self.model.params.Presolve = 2
        self.model.update()
        self.model.optimize()

        if self.model.Status == GRB.INFEASIBLE:
            logs.critical(
                "Model remains infeasible. Searching for invalid constraints."
            )
            self.model.computeIIS()
            logs.critical("Printing bad constraints:")
            for c in self.model.getConstrs():
                if c.IISConstr:
                    logs.critical("%s" % c.ConstrName)
            for c in self.model.getGenConstrs():
                if c.IISGenConstr:
                    logs.critical("%s" % c.GenConstrName)
        else:
            logs.debug("Model Successfully Solved.")
        logs.info("Time to finish solver: {:.3f}".format(time.time() - t1))

    def run_model(self):
        """Construct and solve the Gurobi model (with optional bonus round)."""
        self.round_info = "Round1"
        self.build_model_round1()
        self.optimize_model()
        self.serialize_results_csv()
        if self.config.getboolean("semester", "run_bonus_round"):
            self.round_info = "Round2"
            self.build_model_round2()
            self.optimize_model()
            self.serialize_results_csv()
        logs.info("Scheduling complete, clear skies!")

    # ==================================================================
    # Output.
    # ==================================================================

    def build_schedule(self):
        """
        Build the sparse schedule DataFrame from ``self.Yrds`` and write
        ``semester_plan.csv``.

        Sets ``self.schedule`` to a DataFrame with columns
        ``unique_id, d, s, starname`` -- one row per scheduled exposure start.
        """
        df = pd.DataFrame(self.Yrds.keys(), columns=["unique_id", "d", "s"])
        df["value"] = [self.Yrds[k].x for k in self.Yrds.keys()]
        sparse = df.query("value > 0").drop(columns=["value"]).copy()
        id_to_name = dict(
            zip(self.requests_frame["unique_id"], self.requests_frame["starname"])
        )
        sparse["starname"] = (
            sparse["unique_id"].map(id_to_name).fillna("NO MATCHING NAME")
        )
        sparse.to_csv(
            os.path.join(self.output_directory, "semester_plan.csv"),
            index=False,
            na_rep="",
        )
        self.schedule = sparse

    def _program_statistics(self):
        """Per-program awarded/requested/past/scheduled slot accounting."""
        slots_per_hour = 60 / self.slot_size
        slots_per_night = self.hours_per_night * slots_per_hour

        progs = pd.read_csv(self.programs_file).rename(
            columns={"program": "program_code", "nights": "awarded_nights"}
        )
        progs["awarded_hours"] = progs["awarded_nights"] * self.hours_per_night
        progs["awarded_slots"] = progs["awarded_nights"] * slots_per_night

        rf = self.requests_frame.copy()
        rf["t_slots_per_visit"] = (
            rf["unique_id"].map(self.slots_needed_for_exposure_dict).fillna(0)
        )
        rf["requested_slots"] = (
            rf["t_slots_per_visit"] * rf["n_intra_max"] * rf["n_inter_max"]
        )
        past_n_exp = pd.Series(
            {uid: h.total_n_exposures for uid, h in self.past_history.items()},
            dtype=float,
            name="past_n_exp",
        )
        rf = rf.merge(past_n_exp, left_on="unique_id", right_index=True, how="left")
        rf["past_n_exp"] = rf["past_n_exp"].fillna(0)
        rf["past_slots"] = rf["past_n_exp"] * rf["t_slots_per_visit"]

        per_program_req = (
            rf.groupby("program_code")[["requested_slots", "past_slots"]].sum()
        )

        sched = (
            self.schedule.copy()
            if getattr(self, "schedule", None) is not None
            else pd.DataFrame(columns=["unique_id"])
        )
        sched["sched_slots"] = (
            sched["unique_id"].map(self.slots_needed_for_exposure_dict).fillna(1)
        )
        sched = sched.merge(
            rf[["unique_id", "program_code"]], on="unique_id", how="inner"
        )
        per_program_sched = (
            sched.groupby("program_code")["sched_slots"].sum().rename("scheduled_slots")
        )

        stats = progs.merge(per_program_req, on="program_code", how="left").merge(
            per_program_sched, on="program_code", how="left"
        ).fillna(0)
        for col in ("requested_slots", "past_slots", "scheduled_slots"):
            stats[col.replace("slots", "hours")] = stats[col] / slots_per_hour
            stats[col.replace("slots", "nights")] = (
                stats[col] / slots_per_hour / self.hours_per_night
            )

        total = stats["scheduled_slots"] + stats["past_slots"]
        stats["fullnessA"] = np.where(
            stats["requested_slots"] > 0, 100 * total / stats["requested_slots"], 0.0
        )
        stats["fullnessB"] = np.where(
            stats["awarded_slots"] > 0, 100 * total / stats["awarded_slots"], 0.0
        )
        stats["fullnessC"] = np.where(
            stats["awarded_slots"] > 0,
            100 * stats["requested_slots"] / stats["awarded_slots"],
            0.0,
        )
        return stats.set_index("program_code")

    # Compiled once at module load so to_string stays a one-liner per program.
    _PROGRAM_REPORT_TEMPLATE = Template(
        """{{ program }}
 -- Awarded {{ "%.2f"|format(s['awarded_nights']) }} nights = {{ "%.2f"|format(s['awarded_hours']) }} hours = {{ "%.1f"|format(s['awarded_slots']) }} slots.
 -- Requested {{ "%.2f"|format(s['requested_nights']) }} nights = {{ "%.2f"|format(s['requested_hours']) }} hours = {{ "%.1f"|format(s['requested_slots']) }} slots
 ------ Fullness of requested to awarded: {{ "%.2f"|format(s['fullnessC']) }}%
 -- Past {{ "%.2f"|format(s['past_nights']) }} nights = {{ "%.2f"|format(s['past_hours']) }} hours = {{ "%.1f"|format(s['past_slots']) }} slots.
 -- Scheduled {{ "%.2f"|format(s['scheduled_nights']) }} nights = {{ "%.2f"|format(s['scheduled_hours']) }} hours = {{ "%.1f"|format(s['scheduled_slots']) }} slots
 ------ Fullness of past/scheduled to requested: {{ "%.2f"|format(s['fullnessA']) }}%
 ------ Fullness of past/scheduled to awarded: {{ "%.2f"|format(s['fullnessB']) }}%


"""
    )

    def to_string(self, *, header="Stats for Round1"):
        """Return a human-readable run report for the current schedule.

        Mirrors :meth:`astroq.ttp.model.TTPModel.to_string`. Requires that
        :meth:`build_schedule` has been called so ``self.schedule`` is set.
        """
        if getattr(self, "schedule", None) is None:
            raise RuntimeError("call build_schedule() before to_string()")
        sched = self.schedule

        is_alloc_2d = self.access_record["is_allocated"][0]
        total_slots_in_semester = is_alloc_2d.shape[0] * is_alloc_2d.shape[1]
        allocated_slots = int(np.sum(is_alloc_2d))

        scheduled_starting_slots = len(sched)
        reserved_per_row = (
            sched["unique_id"].map(self.slots_needed_for_exposure_dict).fillna(1) - 1
        )
        reserved_slots = int(reserved_per_row.clip(lower=0).sum())
        total_scheduled_slots = scheduled_starting_slots + reserved_slots
        empty_slots = allocated_slots - total_scheduled_slots

        rf = self.requests_frame
        slots_needed = (
            rf["unique_id"].map(self.slots_needed_for_exposure_dict).fillna(0)
        )
        total_slots_requested = int(
            (slots_needed * rf["n_intra_max"] * rf["n_inter_max"]).sum()
        )

        pct_available = (
            float(np.round(total_scheduled_slots * 100 / allocated_slots, 3))
            if allocated_slots > 0 else 0
        )
        pct_requested = (
            float(np.round(total_scheduled_slots * 100 / total_slots_requested, 3))
            if total_slots_requested > 0 else 0
        )

        lines = [
            header,
            "------------------------------------------------------",
            f"N slots in semester:{total_slots_in_semester}",
            f"N available slots:{allocated_slots}",
            f"N starting slots scheduled: {scheduled_starting_slots}",
            f"N reserved slots: {reserved_slots}",
            f"N total slots scheduled: {total_scheduled_slots}",
            f"N slots left empty: {empty_slots}",
            f"N slots requested (total): {total_slots_requested}",
            f"Utilization (% of available slots): {pct_available}%",
            f"Utilization (% of requested slots): {pct_requested}%",
            "",
            "Program Statistics:",
            "------------------------------------------------------",
        ]
        out = "\n".join(lines) + "\n"

        program_stats = self._program_statistics()
        if len(program_stats) == 0:
            return out + "  No program information available\n"
        for program in sorted(program_stats.index):
            out += self._PROGRAM_REPORT_TEMPLATE.render(
                program=program, s=program_stats.loc[program].to_dict()
            )
        return out

    def serialize_results_csv(self):
        """Write all per-run text artifacts and the HDF5 snapshot."""
        logs.debug("Building human readable schedule.")
        self.build_schedule()
        report = self.to_string()
        with open(os.path.join(self.output_directory, "runReport.txt"), "w") as fh:
            fh.write(report)

        today_idx = self.all_dates_dict[self.current_day]
        selected = list(
            {k[0] for k, v in self.Yrds.items() if v.x > 0 and k[1] == today_idx}
        )
        selected_df = self.requests_frame[
            self.requests_frame["unique_id"].isin(selected)
        ].copy()
        selected_df.to_csv(
            os.path.join(self.output_directory, "request_selected.csv"), index=False
        )
        self.to_hdf5()

    # ==================================================================
    # Persistence (mirrors astroq.nplan layout).
    # ==================================================================

    def to_hdf5(self, hdf5_path=None):
        """Persist ``config_ini_text`` + a few DataFrames + ``access_record``.

        Args:
            hdf5_path (str, optional): defaults to
                ``<output_directory>/semester_planner.h5``.
        """
        if hdf5_path is None:
            hdf5_path = os.path.join(self.output_directory, "semester_planner.h5")
        if os.path.exists(hdf5_path):
            os.remove(hdf5_path)

        self.requests_frame_all.to_hdf(
            hdf5_path, key="requests_frame_all", mode="a", format="table"
        )
        sched = getattr(self, "schedule", None)
        if sched is not None:
            fmt = "fixed" if sched.empty else "table"
            sched.to_hdf(hdf5_path, key="schedule", mode="a", format=fmt)

        with h5py.File(hdf5_path, "a") as f:
            f.attrs["schema_version"] = SEMESTER_PLANNER_H5_SCHEMA
            f.attrs["config_ini_text"] = self._config_ini_text
            f.attrs["run_band3"] = bool(self.run_band3)
            f.attrs["past_history_json"] = json.dumps(
                _serialize_past_history(self.past_history)
            )
            # Save each access_record field as a compressed dataset.
            for field_name in self.access_record.dtype.names:
                f.create_dataset(
                    f"access_record/{field_name}",
                    data=self.access_record[field_name],
                    compression="gzip",
                )

        logs.info(f"SemesterPlanner saved to HDF5: {hdf5_path}")
        return hdf5_path

    @classmethod
    def from_hdf5(cls, hdf5_path):
        """Rehydrate a SemesterPlanner from a snapshot written by :meth:`to_hdf5`.

        Skips Gurobi-only state (model, Yrds, Wrd, theta) and the constraint
        lookup tables -- those are only meaningful when solving. Downstream
        consumers (plot.py, nplan.py) read requests_frame*, schedule,
        access_record, past_history, and queue, all of which are restored.
        """
        with h5py.File(hdf5_path, "r") as f:
            schema = int(f.attrs.get("schema_version", 0))
            if schema != SEMESTER_PLANNER_H5_SCHEMA:
                raise ValueError(
                    f"semester_planner.h5 schema_version={schema} is unsupported "
                    f"(expected {SEMESTER_PLANNER_H5_SCHEMA}). Re-run plan-semester."
                )
            config_ini_text = f.attrs["config_ini_text"]
            if isinstance(config_ini_text, bytes):
                config_ini_text = config_ini_text.decode("utf-8")
            run_band3 = bool(f.attrs["run_band3"])
            past_history_data = json.loads(f.attrs["past_history_json"])

            ar_keys = list(f["access_record"].keys())
            field_data = {k: f[f"access_record/{k}"][:] for k in ar_keys}

        requests_frame_all = pd.read_hdf(hdf5_path, key="requests_frame_all")
        try:
            schedule = pd.read_hdf(hdf5_path, key="schedule")
        except KeyError:
            schedule = None

        instance = cls.__new__(cls)
        instance.start_the_clock = time.time()
        instance._config_ini_text = config_ini_text
        instance.config = ConfigParser()
        instance.config.read_string(config_ini_text)
        instance.run_band3 = run_band3
        instance.queue = astroq.queue.from_config(instance.config)
        instance.queue_name = instance.config.get("global", "queue")

        if "comments" not in requests_frame_all.columns:
            requests_frame_all = requests_frame_all.copy()
            requests_frame_all["comments"] = ""
        instance.requests_frame_all, instance.requests_frame = (
            cls._split_and_clean_requests(requests_frame_all)
        )

        instance.strategy = instance._build_strategy()
        instance.past_history = _deserialize_past_history(past_history_data)
        instance.slots_needed_for_exposure_dict = (
            instance._build_slots_required_dictionary()
        )
        instance.access_obj = ac.Access.from_planner(instance)
        instance.access_record = _recarray_from_fields(field_data, ar_keys)
        instance.schedule = schedule

        logs.info(f"SemesterPlanner loaded from HDF5: {hdf5_path}")
        return instance

    # ==================================================================
    # Misc.
    # ==================================================================

    def add_twilights(self):
        """Add 20-minute buffer to allocation times matching 12 deg twilight."""
        keck = self.queue.observatory
        allocation_df = pd.read_csv(self.allocation_file)

        for idx, row in allocation_df.iterrows():
            date_str = str(row["start"])[:10]
            day = Time(date_str, format="iso", scale="utc")

            evening_12 = keck.twilight_evening_nautical(day, which="next")
            morning_12 = keck.twilight_morning_nautical(day, which="next")

            # 10-minute window for the twilight-match heuristic.
            tol = TimeDelta(10, format="jd") / 1440
            start_time = Time(row["start"])
            if abs(start_time - evening_12) <= tol:
                adjusted_start = start_time - TimeDelta(20, format="jd") / 1440
                allocation_df.loc[idx, "start"] = adjusted_start.strftime(
                    "%Y-%m-%dT%H:%M"
                )
                logs.info(f"Adjusted start time for {date_str}: subtracted 20 min")

            stop_time = Time(row["stop"])
            if abs(stop_time - morning_12) <= tol:
                adjusted_stop = stop_time + TimeDelta(20, format="jd") / 1440
                allocation_df.loc[idx, "stop"] = adjusted_stop.strftime(
                    "%Y-%m-%dT%H:%M"
                )
                logs.info(f"Adjusted stop time for {date_str}: added 20 min")

        allocation_df.to_csv(self.allocation_file, index=False)
        logs.info("Allocation file updated with twilight adjustments")


def _recarray_from_fields(field_data, field_names):
    """Rebuild a recarray from per-field arrays loaded out of HDF5."""
    n_records = field_data[field_names[0]].shape[0]
    dtype_list = []
    for k in field_names:
        d = field_data[k]
        if d.ndim == 1:
            dtype_list.append((k, d.dtype))
        else:
            dtype_list.append((k, d.dtype, d.shape[1:]))
    arr = np.zeros(n_records, dtype=dtype_list)
    for k in field_names:
        arr[k] = field_data[k]
    return arr.view(np.recarray)
