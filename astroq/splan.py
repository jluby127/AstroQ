"""
Module that defines the SemesterPlanner class. This class is responsible for defining,
building, and solving the Gurobi model for semester-level observation planning. It is
nearly completely agnostic to all astronomy knowledge.
"""

import logging
import os
import time
from configparser import ConfigParser
from functools import cached_property
from pathlib import Path
from typing import Any

import gurobipy as gp
import h5py
import numpy as np
import pandas as pd
from astropy.time import Time
from gurobipy import GRB
import astroq.access as ac
import astroq.queue

logs = logging.getLogger(__name__)

# Schema for h5 serialization bump when the on-disk layout changes
SEMESTER_PLANNER_H5_SCHEMA = 4

# ---------------------------------------------------------------------------
# Input contracts. One spec per CSV: required column -> dtype that
# ``load_frame`` coerces to. ``Time`` marks ISO-datetime columns parsed to
# ``astropy.time.Time``. Repair (default filling, legacy "None" strings,
# junk filtering) is the prep stage's job (see astroq.queue.prep_common);
# ``load_frame`` validates and coerces only, and raises on anything else.
# Extra columns (``comments``, weather bands, ...) pass through untouched.
# ---------------------------------------------------------------------------

REQUEST_SCHEMA = {
    "unique_id": str,
    "target": str,
    "program_code": str,
    "ra": float,            # deg
    "dec": float,           # deg
    "exptime": float,       # seconds
    "n_exp": int,
    "n_inter_max": int,
    "tau_inter": int,       # days
    "n_intra_min": int,
    "n_intra_max": int,
    "tau_intra": float,     # hours
    "inactive": bool,
    "splan_weight": float,
}

PAST_SCHEMA = {
    "unique_id": str,
    "target": str,
    "timestamp": str,        # UT ISO
    "exposure_time": float,  # seconds
}

ALLOCATION_SCHEMA = {"start": Time, "stop": Time}

CUSTOM_SCHEMA = {"unique_id": str, "target": str, "start": Time, "stop": Time}

PROGRAMS_SCHEMA = {"program": str, "hours": float, "nights": float}

# Column-name lists kept for consumers that only need presence checks.
REQUEST_COLS = list(REQUEST_SCHEMA)
PAST_COLS = list(PAST_SCHEMA)


def load_frame(path, schema, name, *, key=None, empty_ok=False):
    """Read ``path`` and validate/coerce it against ``schema``.

    Args:
        path (str): CSV location.
        schema (dict): column -> target dtype (``str``/``int``/``float``/
            ``bool``) or ``astropy.time.Time`` for ISO-datetime columns.
        name (str): label used in error messages (e.g. ``"request.csv"``).
        key (str, optional): column whose values must be unique.
        empty_ok (bool): missing/zero-byte/header-only files return an empty
            frame with the schema's columns instead of raising.

    Returns:
        pandas.DataFrame with every schema column coerced; extra columns
        pass through untouched.

    Raises:
        FileNotFoundError: missing file when ``empty_ok=False``.
        ValueError: missing columns, null values in schema columns,
            non-integral values in int columns, unparseable values, or
            duplicate ``key`` values. No repair is attempted.
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        if empty_ok:
            return pd.DataFrame({c: pd.Series(dtype=object) for c in schema})
        raise FileNotFoundError(f"{name} not found: {path}")
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        if empty_ok:
            return pd.DataFrame({c: pd.Series(dtype=object) for c in schema})
        raise

    missing = [c for c in schema if c not in df.columns]
    if missing:
        raise ValueError(f"{name} missing required column(s): {missing}")

    nulls = [c for c in schema if df[c].isna().any()]
    if nulls:
        raise ValueError(
            f"{name} has null values in required column(s): {nulls}. "
            f"Inputs must arrive clean from the prep stage."
        )

    for col, dtype in schema.items():
        if dtype is bool:
            if df[col].dtype != bool:
                raise ValueError(
                    f"{name} column {col!r} must be boolean (True/False); "
                    f"got dtype {df[col].dtype}."
                )
        elif dtype is int:
            vals = pd.to_numeric(df[col])
            if (vals % 1 != 0).any():
                raise ValueError(f"{name} column {col!r} has non-integer values.")
            df[col] = vals.astype(int)
        elif dtype is float:
            df[col] = pd.to_numeric(df[col]).astype(float)
        elif dtype is Time:
            df[col] = df[col].apply(Time)
        else:
            df[col] = df[col].astype(str)

    if key is not None:
        dup = df[key].duplicated(keep=False)
        if dup.any():
            dup_vals = sorted(df.loc[dup, key].unique())
            raise ValueError(f"{name} has duplicate {key!r} values: {dup_vals}")

    return df.reset_index(drop=True)

_ROUND4_THROTTLE_GRACE = 2.0

# Default cap on Round-1-optimum multiplier for the upcoming-night round's
# global shortfall constraint. Referenced from both build_model_upcoming_night_round
# and log_report; keep as one constant so the two never drift apart.
_DEFAULT_GLOBAL_SHORTFALL_SLACK = 1.1

_ROUND_SPECS = {
    "Round1": (1, "Minimize time-weighted shortfall (Lubin et al.)"),
    "Round2": (2, "Maximize inter-program fill factors"),
    "Round3": (3, "Intra-program priorities (hold Round-2 fill)"),
    "Round4": (4, "Minimize empty slots (re-throttle)"),
    "UpcomingNight": (5, "Fill current night (cap global shortfall)"),
}

_PROGRAM_STATS_KEY = """\
** Key
aw     - awarded time (hr)
req    - requested time (hr)
past   - past executed time (hr)
proj   - projected time (past + scheduled future, hr)
miff%  - minimum fill factor (% of award); constrained lower bound active this round
maff%  - maximum fill factor (% of award); throttle ceiling active this round
past%  - past executed time (% of award)
proj%  - projected fill (% of award); should satisfy miff% <= proj% <= maff%
----------------------------------------------------------------"""


class SemesterPlanner:
    """Semester-level scheduler: pick which targets get observed when.

    Formulates the cadenced-scheduling problem of Lubin et al. 2025
    (arXiv:2506.08195) as a Gurobi MILP over a (request, day, slot) grid
    and produces a sparse schedule for the entire semester. The night-level
    slew ordering is handled separately by :class:`astroq.nplan.NightPlanner`.

    Inputs (resolved relative to ``[global] workdir`` in the config):
        - ``request.csv``  -- observing requests, one row per target.
        - ``allocation.csv`` -- telescope time blocks for the semester.
        - ``past.csv`` -- prior observations (caps future ``n_inter_max``).
        - ``custom.csv`` -- PI-supplied per-target observability windows.
        - ``programs.csv`` -- awarded nights per program (drives throttling).

    Key outputs (written to ``<workdir>/outputs``):
        - ``semester_plan.csv`` -- sparse schedule with columns
          ``unique_id, d, s, target``
        - ``request_selected.csv`` -- tonight's targets, the handoff to
          :class:`astroq.nplan.NightPlanner`.
        - ``semester_planner.h5`` -- round-trippable snapshot consumed by
          :class:`astroq.nplan.NightPlanner` and the plotting layer.

    The per-round run report is emitted via :meth:`log_report` (logged at
    INFO) rather than persisted to disk.

    Lifecycle:

        >>> sp = SemesterPlanner("config.ini")
        >>> sp.run_model()       # builds constraints, solves, writes outputs

    Persistence:
        :meth:`to_hdf5` stores the config text plus a handful of DataFrames
        and the precomputed ``access_record``; :meth:`from_hdf5` rehydrates
        a planner suitable for downstream consumers (no Gurobi state). The
        on-disk schema version is :data:`SEMESTER_PLANNER_H5_SCHEMA`.

    Args:
        cf (str): path to the ``config.ini`` file.
    """

    def __init__(self, cf, *, boost=None):
        """See class docstring."""
        logs.debug("Building the SemesterPlanner.")
        self.boost = boost

        # Read config as text so we can persist it verbatim and recreate the
        # parser on from_hdf5.
        self._config_ini_text = Path(cf).read_text()
        self.config = ConfigParser()
        self.config.read_string(self._config_ini_text)
        self.queue = astroq.queue.from_config(self.config)
        self.schedule = None
        self._round1_obj_val = None
        self._round1_weighted_theta = None
        self._round2_slots_by_program = None
        self._hold_fill_alpha = 0.0
        self._open_throttle_grace = None

        workdir = self.config.get("global", "workdir")
        self.output_directory = os.path.join(workdir, "outputs")
        self.allocation_file = self._resolve_path("allocation_file")
        self.custom_file = self._resolve_path("custom_file")
        self.programs_file = self._resolve_path("programs_file")
        os.makedirs(self.output_directory, exist_ok=True)

        self.requests_frame_all, self.requests_frame = self._load_requests_frame()
        self.past_df = self._load_past()

        # Per-request derived columns that depend on past_df live on
        # requests_frame (single source of truth, no parallel dict
        # attributes). Constraint methods derive `dict(zip(...))` adapters
        # locally where Gurobi's quicksum needs O(1) keyed lookup.
        self._attach_past_columns()

        # Observability cube (single source of truth for which slots are valid).
        self.access_obj = ac.Access.from_planner(self)
        self.access_record = self.access_obj.build_access()
        self.observability = self.access_obj.observability(
            self.access_record.is_observable
        )

        # Pre-computed aggregations consumed by the constraint methods. The
        # ones stored on self (joiner, observability_tuples,
        # all_valid_ds_for_request) are read from multiple constraints; the
        # rest live as locals at their call sites.
        self._build_constraint_lookups()
        self._log_boost_current_day_slots()

        self.build_gurobi_model()

        logs.debug("Initializing complete.")

    def _resolve_path(self, key):
        """Resolve a ``[data]`` config key against ``[global] workdir``."""
        raw = self.config.get("data", key)
        workdir = self.config.get("global", "workdir")
        return raw if os.path.isabs(raw) else os.path.join(workdir, raw)

    def _load_past(self):
        """Read ``past.csv`` validated against :data:`PAST_SCHEMA`.

        Empty/missing files yield an empty frame with the canonical schema.
        ``past.csv`` is expected to arrive clean from the prep stage: no
        junk-flagged visits and no ``junk`` column (any junk filtering belongs
        in the past-history producer, not here).
        """
        path = self._resolve_path("past_file")
        return load_frame(path, PAST_SCHEMA, "past.csv", empty_ok=True)

    # ------------------------------------------------------------------
    # Properties (date-derived; path attrs are set in __init__).
    # ------------------------------------------------------------------

    @cached_property
    def semester_length(self):
        """Inclusive semester span in nights (computed once)."""
        start = Time(
            self.config.get("global", "semester_start_day"),
            format="iso",
            scale="utc",
        )
        end = Time(
            self.config.get("global", "semester_end_day"),
            format="iso",
            scale="utc",
        )
        return int(round(end.jd - start.jd)) + 1

    @property
    def all_dates_array(self):
        return self.access_obj.all_dates_array

    @property
    def all_dates_dict(self):
        return self.access_obj.all_dates_dict

    @property
    def today_starting_night(self):
        return self.all_dates_dict[self.config.get("global", "current_day")]

    # ------------------------------------------------------------------
    # Construction helpers.
    # ------------------------------------------------------------------

    def _load_requests_frame(self):
        """Read + validate request.csv. Returns ``(all_frame, active_frame)``.

        Expects a clean request.csv from the prep stage: strategy defaults are
        already filled (see :func:`astroq.queue.prep_common.standardize_request_strategy`),
        so splan only validates/coerces against :data:`REQUEST_SCHEMA`, derives
        slot columns, and fails on duplicate active unique_id (which would
        otherwise surface as a cryptic Gurobi error).

        Slot columns appended by :meth:`_attach_slot_columns`:

        - ``t_visit_slots`` -- full per-visit duration in slots, from
          :meth:`astroq.queue.base.Queue.visit_seconds` (includes inter-shot
          readouts and slew overhead), rounded and clipped to >= 1. This is the
          slot reservation charged by every Gurobi consumer.
        - ``tau_intra_slots`` -- minimum intra-night spacing between visits, in
          slots.

        Original units of ``exptime`` (seconds) and ``tau_intra`` (hours) are
        left untouched.
        """
        request_file = self._resolve_path("request_file")
        rfa = load_frame(
            request_file, REQUEST_SCHEMA, f"request.csv ({request_file!r})"
        )
        logs.warning(
            f"There are {int(rfa['inactive'].sum())} inactive of {len(rfa)} requests."
        )

        self._attach_slot_columns(rfa)

        rf = rfa[~rfa["inactive"]].reset_index(drop=True).copy()

        dup_mask = rf["unique_id"].duplicated(keep=False)
        if dup_mask.any():
            dup_ids = sorted(rf.loc[dup_mask, "unique_id"].unique())
            raise ValueError(
                f"Duplicate unique_id among active requests in {request_file!r}: "
                f"{dup_ids}. Remove or merge duplicate rows so each active "
                f"request has one row."
            )

        return rfa, rf

    def _attach_slot_columns(self, rf):
        """Append ``t_visit_slots`` and ``tau_intra_slots`` columns to ``rf``.

        - ``t_visit_slots`` -- full per-visit duration in slots, from
          :meth:`Queue.visit_seconds` (includes inter-shot readouts and
          slew overhead). Rounded; clipped to >= 1.
        - ``tau_intra_slots`` -- minimum intra-night spacing between
          visits, in slots.

        Mutates and returns ``rf`` (idempotent).
        """
        slot_size = self.config.getfloat("semester", "slot_size")
        visit_s = self.queue.visit_seconds(rf["exptime"], rf["n_exp"])
        rf["t_visit_slots"] = (
            (visit_s / (slot_size * 60.0)).round().clip(lower=1).astype(int)
        )
        rf["tau_intra_slots"] = (
            (rf["tau_intra"] * 60 / slot_size).round().astype(int)
        )
        return rf

    def _build_constraint_lookups(self):
        """Build aggregation tables consumed by the constraint methods.

        Only the three multi-consumer tables (observability_tuples, joiner,
        all_valid_ds_for_request) are stored on self. Single-consumer
        derivations live at their call sites.
        """
        self.observability_tuples = list(
            self.observability.itertuples(index=False, name=None)
        )
        strategy_cols = [
            "unique_id",
            "target",
            "n_intra_min",
            "n_intra_max",
            "n_inter_max",
            "tau_inter",
            "t_visit_slots",
            "tau_intra_slots",
        ]
        self.joiner = pd.merge(
            self.requests_frame[strategy_cols], self.observability, on=["unique_id"]
        )

        schedulable_requests = set(self.joiner["unique_id"].unique())
        all_requests = list(self.requests_frame["unique_id"])
        missing = sum(uid not in schedulable_requests for uid in all_requests)
        logs.warning(
            f"There are {missing} targets out of {len(all_requests)} "
            f"that have no valid day/slot pairs and therefore are effectively "
            f"removed from the model."
        )

        self.all_valid_ds_for_request = (
            self.joiner.groupby(["unique_id"])[["d", "s"]].agg(list)
        )
        self.build_yrds_tuples()

    def build_yrds_tuples(self):
        """``yrds_tuples``: full observability."""
        self.yrds_tuples = list(self.observability_tuples)
        self._index_yrds_tuples()

    def _index_yrds_tuples(self):
        """Build ``_yrds_keys`` and ``_yrds_by_ds`` from ``yrds_tuples``."""
        self._yrds_keys = set(self.yrds_tuples)
        by_ds = {}
        for uid, d, s in self.yrds_tuples:
            by_ds.setdefault((int(d), int(s)), set()).add(uid)
        self._yrds_by_ds = by_ds

    def _log_boost_current_day_slots(self):
        """Report observable slot counts on current_day for each boosted target."""
        if self.boost is None:
            return
        current_day = self.config.get("global", "current_day")
        d_today = self.today_starting_night
        boost_by_uid = self._boost_by_uid
        uid_to_target = dict(
            zip(
                self.requests_frame_all["unique_id"],
                self.requests_frame_all["target"],
            )
        )
        factor = next(iter(boost_by_uid.values()))
        logs.info(
            "Boost on current_day=%s (d=%d), factor=%g:",
            current_day,
            d_today,
            factor,
        )
        joiner_uids = self.joiner["unique_id"]
        joiner_d = self.joiner["d"]
        for uid in boost_by_uid:
            n_slots = int(((joiner_uids == uid) & (joiner_d == d_today)).sum())
            target = uid_to_target.get(uid, "(unknown unique_id)")
            logs.info(
                "  %s (%s): %d observable slot(s) on current_day",
                uid,
                target,
                n_slots,
            )

    def build_gurobi_model(self):
        """Instantiate the Gurobi model and add ``Yrds``, ``Wrd``, ``theta``."""
        self.model = gp.Model("Semester_Scheduler")
        observability_nights = (
            self.joiner.loc[self.joiner["n_intra_max"] > 1, ["unique_id", "d"]]
            .drop_duplicates()
        )
        self.Yrds = self.model.addVars(
            self.yrds_tuples, vtype=GRB.BINARY, name="Requests_Slots"
        )
        if not observability_nights.empty:
            self.Wrd = self.model.addVars(
                list(observability_nights.itertuples(index=False, name=None)),
                vtype=GRB.BINARY,
                name="OnSky",
            )
        self.theta = self.model.addVars(
            list(self.requests_frame["unique_id"]), name="Shortfall"
        )

    def _attach_past_columns(self):
        """Attach past-history aggregates and max-obs caps to ``requests_frame``.

        Aggregates are indexed by ``unique_id`` over UT calendar nights
        (``timestamp[:10]``); missing uids default to 0 (or ``""``).
        ``desired_max_obs`` is the Round-1 night cap; ``absolute_max_obs``
        relaxes it by ``maximum_bonus_size`` for the bonus round. Both
        collapse to ``past_nights_observed`` when a target is over-observed
        so the model stays feasible.
        """
        rf = self.requests_frame
        uids = rf["unique_id"]

        if self.past_df.empty:
            agg = pd.DataFrame(
                {"nights": 0, "n_exp": 0, "last": ""}, index=uids,
            )
        else:
            night = self.past_df["timestamp"].str[:10]
            g = self.past_df.assign(_night=night).groupby("unique_id")
            agg = pd.DataFrame({
                "nights": g["_night"].nunique(),
                "n_exp": g.size(),
                "last": g["_night"].max(),
            }).reindex(uids).fillna({"nights": 0, "n_exp": 0, "last": ""})

        rf["past_nights_observed"] = agg["nights"].astype(int).to_numpy()
        rf["past_n_exposures"] = agg["n_exp"].astype(int).to_numpy()
        rf["past_date_last_observed"] = agg["last"].astype(str).to_numpy()

        bonus = self.config.getfloat("semester", "maximum_bonus_size")
        n_max = rf["n_inter_max"].to_numpy()
        past = rf["past_nights_observed"].to_numpy()
        over = past > n_max
        desired = np.where(over, past, n_max - past)
        absolute = np.where(
            over, past,
            np.maximum(desired + (n_max * bonus).astype(int), past),
        )
        rf["desired_max_obs"] = desired.astype(int)
        rf["absolute_max_obs"] = absolute.astype(int)

    # ==================================================================
    # Constraints 
    # ==================================================================

    def constraint_build_theta_multivisit(self):
        """Build the shortfall matrix, Theta.

        Notes:  
            Equation 3 in Lubin et al. 2025.
        """
        logs.info("Constraint: Build theta variable")
        rf_indexed = self.requests_frame.set_index("unique_id")
        for uid in self.joiner["unique_id"].unique():
            self.model.addConstr(
                self.theta[uid] >= 0, f"greater_than_zero_shortfall_{uid}"
            )
            ds_pairs = list(
                zip(
                    self.all_valid_ds_for_request.loc[uid].d,
                    self.all_valid_ds_for_request.loc[uid].s,
                )
            )
            row = rf_indexed.loc[uid]
            rhs = (
                row["n_inter_max"]
                - row["past_nights_observed"]
                - gp.quicksum(self.Yrds[uid, d, s] for d, s in ds_pairs)
                / row["n_intra_max"]
            )
            self.model.addConstr(
                self.theta[uid] >= rhs,
                f"greater_than_nobs_shortfall_{uid}",
            )

    def constraint_reserve_multislot_exposures(self):
        """
        See Constraint 1 in Lubin et al. 2025.

        Reserve multiple time slots for exposures that require more than one time slot
        to complete, ensuring no other observations are scheduled during these slots.
        """
        logs.info("Constraint: Reserve slots for multi-slot exposures.")
        rf = self.requests_frame
        max_t_visit = int(rf["t_visit_slots"].max())
        R_geq_t_visit = {
            t: set(rf.loc[rf["t_visit_slots"] >= t, "unique_id"])
            for t in range(1, max_t_visit + 1)
        }

        for d, s in self.observability.drop_duplicates(["d", "s"])[
            ["d", "s"]
        ].itertuples(index=False, name=None):
            uids_at = self._yrds_by_ds.get((int(d), int(s)), set())
            if not uids_at:
                continue
            rhs = []
            for delta in range(1, max_t_visit):
                s_shift = s - delta
                uids_shift = self._yrds_by_ds.get((int(d), int(s_shift)), set())
                for uid in uids_shift & R_geq_t_visit[delta + 1]:
                    rhs.append(self.Yrds[uid, d, s_shift])
            lhs = 1 - gp.quicksum(self.Yrds[uid, d, s] for uid in uids_at)
            self.model.addConstr(
                lhs >= gp.quicksum(rhs), f"reserve_multislot_{d}d_{s}s"
            )

    def constraint_enforce_internight_cadence(self):
        """
        See Constraint 3 in Lubin et al. 2025.

        Ensure that the minimum number of days pass between consecutive observations of
        a given target.
        """
        logs.info("Constraint: Enforce inter-night cadence.")
        joiner = self.joiner
        intercadence = pd.merge(
            joiner.drop_duplicates(["unique_id", "d"]),
            joiner[["unique_id", "d", "s"]],
            suffixes=["", "3"],
            on=["unique_id"],
        ).query("d + 0 < d3 < d + tau_inter")
        intercadence_tracker = intercadence.groupby(["unique_id", "d"])[
            ["d3", "s3"]
        ].agg(list)
        slots_on_day_for_r = (
            self.observability.groupby(["unique_id", "d"])["s"]
            .apply(list)
            .to_frame("s3")
        )

        # Inter-night cadence of 1 day has no forbidden future slots; skip
        # those rows and drop the duplicates-per-day rows up front.
        valid = joiner[joiner["tau_inter"] > 1].drop_duplicates(
            subset=["unique_id", "d"]
        )
        for _, row in valid.iterrows():
            constrained_slots_tonight = [
                int(s2)
                for s2 in slots_on_day_for_r.loc[(row.unique_id, row.d)][0]
                if (row.unique_id, int(row.d), int(s2)) in self._yrds_keys
            ]
            if not constrained_slots_tonight:
                continue
            if (row.unique_id, row.d) not in intercadence_tracker.index:
                continue
            future = intercadence_tracker.loc[(row.unique_id, row.d)]
            ds_pairs = [
                (int(d3), int(s3))
                for d3, s3 in zip(
                    np.array(future.d3).flatten(),
                    np.array(future.s3).flatten(),
                )
                if (row.unique_id, int(d3), int(s3)) in self._yrds_keys
            ]
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

        Ensure that the minimum number of hours pass between consecutive observations of
        a given target on the same night.
        """
        logs.info("Constraint: Enforce intra-night cadence.")
        valid = self.joiner[self.joiner["n_intra_max"] > 1]
        intracadence_frame = pd.merge(
            valid.drop_duplicates(["unique_id", "d", "s"]),
            valid[["unique_id", "d", "s"]],
            suffixes=["", "3"],
            on=["unique_id", "d"],
        ).query("s + 0 < s3 < s + tau_intra_slots")
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

    def constraint_set_max_desired_unique_nights_Wrd(self):
        """
        See Constraint 2 in Lubin et al. 2025.

        Limit the number of observations scheduled for a given target to the
        maximum value provided by the PI. This constraint may later be relaxed
        if Round 2 of scheduling is invoked.
        """
        logs.info("Constraint: Set desired maximum observations.")
        multi_visit_uids = self.multi_visit_uids
        schedulable_uids = set(self.joiner["unique_id"].unique())
        single_visit_uids = [
            uid for uid in schedulable_uids if uid not in multi_visit_uids
        ]
        desired_max_obs = self.requests_frame.set_index("unique_id")["desired_max_obs"]
        for uid in multi_visit_uids:
            all_d = list(set(self.all_valid_ds_for_request.loc[uid].d))
            self.model.addConstr(
                gp.quicksum(self.Wrd[uid, d] for d in all_d)
                <= desired_max_obs.loc[uid],
                f"max_desired_unique_nights_for_request_{uid}",
            )
        for uid in single_visit_uids:
            available = list(
                zip(
                    self.all_valid_ds_for_request.loc[uid].d,
                    self.all_valid_ds_for_request.loc[uid].s,
                )
            )
            self.model.addConstr(
                gp.quicksum(self.Yrds[uid, d, s] for d, s in available)
                <= desired_max_obs.loc[uid],
                f"max_desired_unique_nights_for_request_{uid}",
            )

    def remove_constraint_set_max_desired_unique_nights_Wrd(self):
        """
        Bonus round: not in Lubin et al. 2025.

        Remove the maximum number of observations set by
        :meth:`constraint_set_max_desired_unique_nights_Wrd`.
        """
        logs.info("Constraint: Removing previous maximum observations constraint.")
        for uid in self.multi_visit_uids:
            rm_const = self.model.getConstrByName(
                f"max_desired_unique_nights_for_request_{uid}"
            )
            self.model.remove(rm_const)

    def constraint_set_max_absolute_unique_nights_Wrd(self):
        """
        Bonus round: not in Lubin et al. 2025.

        Set the maximum number of observations for a target to 150% of the
        original requested number.
        """
        logs.info("Constraint: Set absolute maximum observations.")
        absolute_max_obs = self.requests_frame.set_index("unique_id")["absolute_max_obs"]
        for uid in self.multi_visit_uids:
            all_d = list(set(self.all_valid_ds_for_request.loc[uid].d))
            self.model.addConstr(
                gp.quicksum(self.Wrd[uid, d] for d in all_d)
                <= absolute_max_obs.loc[uid],
                f"max_absolute_unique_nights_for_request_{uid}",
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
        multi_visit_uids = self.multi_visit_uids
        for _, row in per_day.iterrows():
            slots_tonight = [
                int(s3)
                for s3 in grouped_s.loc[(row.unique_id, row.d)]["s"]
                if (row.unique_id, int(row.d), int(s3)) in self._yrds_keys
            ]
            if not slots_tonight:
                continue
            name_tag = f"{row.unique_id}_{row.d}d_{row.s}s"
            visits_tonight = gp.quicksum(
                self.Yrds[row.unique_id, row.d, s3] for s3 in slots_tonight
            )
            if row.unique_id in multi_visit_uids:
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
                    f"enforce_max_visits_{name_tag}",
                )

    @cached_property
    def multi_visit_uids(self):
        """uids that may receive >1 visit per night (Wrd is defined for these).

        Built once from ``joiner`` (itself constructed once in ``__init__``).
        """
        return set(
            self.joiner.loc[self.joiner["n_intra_max"] > 1, "unique_id"].unique()
        )

    @cached_property
    def _t_visit_slots_by_uid(self):
        """dict[unique_id -> t_visit_slots]. Built once; every constraint and
        objective method that needs per-request slot durations reads this
        instead of re-zipping ``requests_frame`` locally."""
        return dict(
            zip(self.requests_frame["unique_id"], self.requests_frame["t_visit_slots"])
        )

    @cached_property
    def _boost_by_uid(self):
        """dict[unique_id -> boost factor], or ``None`` if no boost was passed."""
        if self.boost is None:
            return None
        return dict(
            zip(self.boost["unique_id"].astype(str), self.boost["boost"].astype(float))
        )

    @cached_property
    def _programs_frame(self):
        """``programs.csv`` validated against :data:`PROGRAMS_SCHEMA`, indexed
        by program code. Read once and left unmutated (per-formula columns are
        derived on demand, e.g. via :attr:`_program_awarded_slots`, rather
        than written back onto this cached frame)."""
        df = load_frame(
            self.programs_file, PROGRAMS_SCHEMA, "programs.csv", key="program"
        )
        return df.set_index("program")

    @cached_property
    def _program_awarded_slots(self):
        """Series[program -> awarded slots this semester].

        Shared ``nights * hours_per_night * 60 / slot_size`` formula,
        previously re-derived independently in constraint_throttle,
        build_model_round2_priority_NEW, and to_string.
        """
        slot_size = self.config.getfloat("semester", "slot_size")
        hours_per_night = self.config.getfloat("semester", "hours_per_night")
        return self._programs_frame["nights"] * hours_per_night * 60 / slot_size

    @cached_property
    def _active_uids_by_program(self):
        """dict[program_code -> {unique_id}] over active requests, built via
        one groupby pass. Several priority-round methods previously rebuilt
        this via a per-program boolean-mask comprehension instead."""
        return (
            self.requests_frame.groupby("program_code")["unique_id"]
            .apply(set)
            .to_dict()
        )

    @cached_property
    def _yrds_frame(self):
        """(unique_id, d, s, program_code, t_visit_slots, var) for every
        schedulable triple, merged once. Several methods previously
        recomputed a program's scheduled-slot total by rescanning
        ``self.yrds_tuples`` once per program (O(programs x tuples)); this
        lets that collapse to a single groupby pass (O(tuples))."""
        df = pd.DataFrame(self.yrds_tuples, columns=["unique_id", "d", "s"])
        df = df.merge(
            self.requests_frame[["unique_id", "program_code"]],
            on="unique_id",
            how="left",
        )
        if df["program_code"].isna().any():
            raise ValueError(
                "_yrds_frame: some yrds_tuples reference a unique_id missing "
                "from requests_frame; every schedulable request must be active."
            )
        df["t_visit_slots"] = df["unique_id"].map(self._t_visit_slots_by_uid)
        df["var"] = [
            self.Yrds[r, d, s] for r, d, s in zip(df["unique_id"], df["d"], df["s"])
        ]
        return df

    @cached_property
    def _program_slot_expr(self):
        """dict[program_code -> gp.LinExpr] of scheduled slot-time
        (symbolic). Round-invariant -- Yrds vars and program membership
        don't change across rounds -- so this is safe to cache."""
        return {
            p: gp.quicksum(v * n for v, n in zip(g["var"], g["t_visit_slots"]))
            for p, g in self._yrds_frame.groupby("program_code")
        }

    def _program_slot_value(self):
        """dict[program_code -> float] of scheduled slot-time at the
        *current* Gurobi solution. Not cached -- ``.X`` changes every round."""
        return {
            p: sum(v.X * n for v, n in zip(g["var"], g["t_visit_slots"]))
            for p, g in self._yrds_frame.groupby("program_code")
        }

    # ---- throttling & bonus round ----

    def _past_slots_by_program(self):
        """Past slots consumed per program, summed over ALL request rows.

        Counts both active and inactive requests: inactive targets can never
        be scheduled, but their past observations still consume the program's
        throttle budget, so a PI cannot reclaim time by flipping a target
        inactive. Returns ``dict[program_code -> int past_slots]``.
        """
        rfa = self.requests_frame_all
        t_visit = dict(zip(rfa["unique_id"], rfa["t_visit_slots"]))
        if self.past_df.empty:
            past_n = pd.Series(dtype="int64")
        else:
            past_n = self.past_df.groupby("unique_id").size()

        out = {}
        for uid, prog in zip(rfa["unique_id"], rfa["program_code"]):
            slots = int(past_n.get(uid, 0)) * t_visit.get(uid, 0)
            out[prog] = out.get(prog, 0) + slots
        return out

    def constraint_throttle(self, throttle_grace=1.0):
        """
        Not described in Lubin et al. 2025.

        Ensure that no program is scheduled for more time than they bring to
        the queue (within a grace amount). Past usage is counted over ALL
        request rows (active and inactive) via
        :meth:`_past_slots_by_program`, while only active targets contribute
        schedulable slots (inactive targets have no ``Yrds`` variables).
        """
        logs.info("Constraint: Throttling over-requested programs.")
        awarded_slots_grace_by_program = (
            self._program_awarded_slots * throttle_grace
        ).astype(int)

        # Past budget: ALL rows (active + inactive).
        past_used_slots_by_program = self._past_slots_by_program()

        # Schedulable budget: only ACTIVE targets get Yrds variables, so
        # _program_slot_expr (built from requests_frame/yrds_tuples, both
        # active-only) stays restricted to active uids.
        program_slot_expr = self._program_slot_expr

        clamped = []
        for program, awarded_slots_grace in awarded_slots_grace_by_program.items():
            awarded_slots_grace = int(awarded_slots_grace)
            schedulable_slots = program_slot_expr.get(program, 0)
            past_used = past_used_slots_by_program.get(program, 0)
            if awarded_slots_grace < past_used:
                clamped.append(program)
                awarded_slots_grace = past_used

            self.model.addConstr(
                awarded_slots_grace - past_used >= schedulable_slots,
                f"throttle_program_{program}",
            )

        if clamped:
            logs.warning(
                "Throttle: %d program(s) at/over grace from past alone "
                "(no new scheduling allowed): %s",
                len(clamped),
                ", ".join(sorted(str(p) for p in clamped)),
            )

    def constraint_fix_previous_objective(self, epsilon=0.03):
        """
        Bonus round: not in Lubin et al. 2025.

        Ensure that the Round-2 objective is within ``epsilon`` of Round-1.
        """
        logs.info("Constraint: Fixing the previous solution's objective value.")
        self.model.addConstr(
            gp.quicksum(self.theta[uid] for uid in self.requests_frame["unique_id"])
            <= self.model.objval + epsilon,
            "fix_previous_objective",
        )

    def constraint_fix_global_shortfall(self, slack_factor):
        """Upcoming-night round: cap weighted shortfall at Round-1 optimum * slack."""
        if self._round1_weighted_theta is None:
            raise RuntimeError(
                "Round 1 must be solved before fixing global shortfall."
            )
        cap = self._round1_weighted_theta * slack_factor
        logs.info(
            "Constraint: weighted shortfall <= Round-1 optimum * %g "
            "(cap=%.3f from Round-1 weighted shortfall=%.3f)",
            slack_factor,
            cap,
            self._round1_weighted_theta,
        )
        self.model.addConstr(
            self._weighted_theta_expr() <= cap,
            "fix_global_shortfall_upcoming_night",
        )

    # ==================================================================
    # Objectives.
    # ==================================================================

    def _weighted_theta_expr(self):
        """Time-weighted global shortfall (Round 1 objective without boost)."""
        schedulable_uids = list(self.joiner["unique_id"].unique())
        t_visit_slots = self._t_visit_slots_by_uid
        return gp.quicksum(
            self.theta[uid] * t_visit_slots[uid] for uid in schedulable_uids
        )

    def _eval_weighted_theta(self):
        """Evaluate weighted shortfall at the current Gurobi solution."""
        schedulable_uids = list(self.joiner["unique_id"].unique())
        t_visit_slots = self._t_visit_slots_by_uid
        return sum(
            self.theta[uid].X * t_visit_slots[uid] for uid in schedulable_uids
        )

    def set_objective_minimize_theta_time_normalized(self):
        """See Equation 1 in Lubin et al. 2025."""
        theta_obj = self._weighted_theta_expr()
        if self.boost is not None:
            boost_by_uid = self._boost_by_uid
            d_today = self.today_starting_night
            boost_terms = [
                boost_by_uid[uid] * self.Yrds[uid, d, s]
                for uid, d, s in self.observability_tuples
                if d == d_today and uid in boost_by_uid
            ]
            if boost_terms:
                logs.info(
                    "Objective: boost term for %d unique_id(s) on current_day=%s.",
                    len(boost_by_uid),
                    self.config.get("global", "current_day"),
                )
                theta_obj -= gp.quicksum(boost_terms)
        self.model.setObjective(theta_obj, GRB.MINIMIZE)

    def set_objective_maximize_slots_used(self):
        """Bonus round: maximize filled slots."""
        logs.info("Objective: Maximize the number of slots used.")
        t_visit_slots = self._t_visit_slots_by_uid
        self.model.setObjective(
            gp.quicksum(
                t_visit_slots[uid] * self.Yrds[uid, d, s]
                for uid, d, s in self.yrds_tuples
            ),
            GRB.MAXIMIZE,
        )

    def set_objective_maximize_slots_used_tonight(self):
        """Upcoming-night round: maximize filled slots on ``current_day``."""
        d_today = self.today_starting_night
        current_day = self.config.get("global", "current_day")
        logs.info(
            "Objective: Maximize slot usage on upcoming night current_day=%s (d=%d).",
            current_day,
            d_today,
        )
        t_visit_slots = self._t_visit_slots_by_uid
        self.model.setObjective(
            gp.quicksum(
                t_visit_slots[uid] * self.Yrds[uid, d, s]
                for uid, d, s in self.yrds_tuples
                if d == d_today
            ),
            GRB.MAXIMIZE,
        )

    def set_objective_minimize_empty_slots(self):
        """Bonus round: minimize empty slots."""
        logs.info("Objective: Minimize the number of empty slots.")
        t_visit_slots = self._t_visit_slots_by_uid
        total_slots = self.semester_length * self.access_obj.nslots
        self.model.setObjective(
            (total_slots - gp.quicksum(
                t_visit_slots[uid] * self.Yrds[uid, d, s]
                for uid, d, s in self.yrds_tuples
            )),
            GRB.MINIMIZE,
        )

    def build_model_round2_priority_NEW(self):
        """
        New round 2 objective.

        maximize the completion rate of all programs
        """
        # self.constraint_fix_previous_objective()

        awarded_slots_by_program = self._program_awarded_slots.to_dict()
        program_slot_expr = self._program_slot_expr
        program_slot_value = self._program_slot_value()

        self.program_awarded_slots = {}
        self.program_fill_factor = {}
        for p in self._active_uids_by_program:
            awarded_slots = awarded_slots_by_program.get(p)
            if awarded_slots is None or awarded_slots <= 0:
                logs.warning(
                    "Program %s missing or has non-positive awarded slots; skipping fill factor.",
                    p,
                )
                continue
            slots_used = program_slot_expr.get(p, 0)
            self.program_awarded_slots[p] = float(awarded_slots)
            self.program_fill_factor[p] = slots_used / awarded_slots

            slots_used_val = program_slot_value.get(p, 0.0)
            fill_factor_val = slots_used_val / awarded_slots

            logs.info(
                "Round2 priority: program %s awarded_slots=%.0f slots_used=%.0f fill_factor=%.3f",
                p,
                awarded_slots,
                slots_used_val,
                fill_factor_val,
            )

            self.model.addConstr(
                slots_used_val/awarded_slots <= self.program_fill_factor[p],
                "maintain_fill_factor_for_program_" + p,
            )

        self.model.setObjective(
            gp.quicksum(self.program_fill_factor[p] for p in self.program_fill_factor.keys()),
            GRB.MAXIMIZE,
        )

    def remove_constraint_throttle(self):
        """Remove throttle constraints from a prior round."""
        logs.info("Constraint: Removing previous throttle constraints.")
        for program in self._programs_frame.index:
            rm_const = self.model.getConstrByName(f"throttle_program_{program}")
            if rm_const is not None:
                self.model.remove(rm_const)

    def capture_theta_prior(self):
        """Fix Gurobi ``theta_prior`` vars to each schedulable request's solved shortfall."""
        schedulable = set(self.joiner["unique_id"].unique())
        self.theta_prior = {}
        for uid in schedulable:
            prior_val = float(self.theta[uid].X)
            self.theta_prior[uid] = float(prior_val)

    def remove_constraint_fix_previous_completion_rates(self):
        """Drop per-target ``theta <= theta_prior`` constraints from an earlier round."""
        for constr in list(self.model.getConstrs()):
            if constr.ConstrName.startswith("theta_le_prior_"):
                self.model.remove(constr)

    def constraint_fix_previous_completion_rates(self, epsilon=0.0):
        """Keep each target's shortfall no worse than the prior round (lower ``theta`` is better)."""
        self.capture_theta_prior()
        logs.info(
            "Constraint: Per-target shortfall must stay at or below prior round "
            "(epsilon=%g).",
            epsilon,
        )
        for uid, prior_var in self.theta_prior.items():
            self.model.addConstr(
                self.theta[uid] <= prior_var,
                f"theta_le_prior_{uid}",
            )


    def build_model_round4_priority_NEW(self):
        """
        New round 4 objective.

        Minimize the number of empty slots.
        """
        self.constraint_fix_previous_completion_rates()
        # grace = self.config.getfloat("semester", "throttle_grace")
        self.remove_constraint_throttle()
        self.constraint_throttle(throttle_grace=_ROUND4_THROTTLE_GRACE)
        self._open_throttle_grace = _ROUND4_THROTTLE_GRACE
        self.set_objective_minimize_empty_slots()

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
        self.constraint_throttle(throttle_grace=1.0)
        self.set_objective_minimize_theta_time_normalized()
        logs.debug(f"Time to build constraints: {np.round(time.time() - t1, 3):.3f}")

    def build_model_round2(self):
        """Round 2 constraints + objective (bonus round)."""
        t1 = time.time()
        self.remove_constraint_set_max_desired_unique_nights_Wrd()
        self.constraint_set_max_absolute_unique_nights_Wrd()
        self.constraint_fix_previous_objective()
        self.set_objective_maximize_slots_used()
        logs.debug(f"Time to build constraints: {np.round(time.time() - t1, 3):.3f}")

    def build_model_upcoming_night_round(self):
        """Upcoming-night round: cap global shortfall, fill ``current_day``."""
        t1 = time.time()
        if self._open_throttle_grace is not None:
            self.remove_constraint_throttle()
            self.constraint_throttle(throttle_grace=self._open_throttle_grace)
        slack = self.config.getfloat(
            "semester", "global_shortfall_slack", fallback=_DEFAULT_GLOBAL_SHORTFALL_SLACK
        )
        self.constraint_fix_global_shortfall(slack_factor=slack)
        self.set_objective_maximize_slots_used_tonight()
        logs.debug(f"Time to build constraints: {np.round(time.time() - t1, 3):.3f}")

    def _planned_round_count(self):
        """Total rounds that ``run_model`` will execute."""
        n = 1
        if self.config.getboolean("semester", "run_bonus_round"):
            n += 3
        if self.config.getboolean("semester", "run_upcoming_night_round", fallback=False):
            n += 1
        return n

    def _begin_round(self, round_label, *, round_num, round_total):
        """Log a visible banner at the start of each scheduling round."""
        _, desc = _ROUND_SPECS[round_label]
        logs.info(
            "===== Semester Round %d/%d: %s =====", round_num, round_total, desc
        )

    def _log_solver_config_once(self):
        """Emit semester solver settings once per ``run_model`` invocation."""
        if not self.config.has_option("semester", "method"):
            raise ValueError(
                "[semester] method is required; expected milp or norel+milp"
            )
        method = self.config.get("semester", "method").strip().lower()
        if method not in ("milp", "norel+milp"):
            raise ValueError(
                f"[semester] method={method!r} invalid; expected milp or norel+milp"
            )

        max_solve_time = self.config.getfloat("semester", "max_solve_time")
        warmstart_time = self.config.getfloat("semester", "warmstart_time", fallback=0.0)
        if method == "norel+milp":
            if warmstart_time >= max_solve_time:
                raise ValueError(
                    f"[semester] max_solve_time={max_solve_time} must exceed "
                    f"warmstart_time={warmstart_time} for method=norel+milp"
                )
            norel_budget = warmstart_time
            milp_time = max_solve_time - warmstart_time
        else:
            if warmstart_time > 0:
                logs.warning(
                    "[semester] warmstart_time=%g ignored for method=milp",
                    warmstart_time,
                )
            norel_budget = 0.0
            milp_time = max_solve_time

        logs.info(
            "Semester: method=%s warmstart=%gs milp=%gs rounds=%d",
            method,
            norel_budget,
            milp_time,
            self._planned_round_count(),
        )
        return method, norel_budget, milp_time

    def _throttle_grace_for_round(self, round_label):
        """Throttle grace factor in effect when reporting program statistics."""
        if round_label == "Round4":
            return _ROUND4_THROTTLE_GRACE
        if round_label == "UpcomingNight" and self._open_throttle_grace is not None:
            return self._open_throttle_grace
        return self.config.getfloat("semester", "throttle_grace", fallback=1.0)

    def _capture_round2_slots_by_program(self):
        """Record per-program scheduled slot-time at end of Round 2."""
        self._round2_slots_by_program = self._program_slot_value()

    def optimize_model(self, round_label="Round1", *, build_secs=0.0):
        """Solve the Gurobi model (with IIS diagnostics on infeasibility)."""
        logs.debug("Begin model solve for %s.", round_label)
        # method/norel_budget/milp_time are parsed and validated once per
        # run_model() call by _log_solver_config_once, not re-derived here.
        norel_budget = self._solver_norel_budget
        milp_time = self._solver_milp_time

        is_first_round = round_label == "Round1"
        show_gurobi = self.config.getboolean("semester", "show_gurobi_output")
        if not is_first_round and self.config.has_option(
            "semester", "show_gurobi_output_later_rounds"
        ):
            show_gurobi = self.config.getboolean(
                "semester", "show_gurobi_output_later_rounds"
            )
        if is_first_round:
            presolve = 2
            mip_focus = 1
        else:
            presolve = self.config.getint(
                "semester", "presolve_later_rounds", fallback=1
            )
            mip_focus = 2

        self.model.params.TimeLimit = milp_time
        self.model.Params.OutputFlag = int(show_gurobi)
        self.model.params.MIPGap = self.config.getfloat("semester", "max_solve_gap")
        self.model.params.NoRelHeurTime = norel_budget if is_first_round else 0.0
        self.model.params.Presolve = presolve
        self.model.params.MIPFocus = mip_focus
        self.model.params.DegenMoves = -1
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
                    logs.critical("%s", c.ConstrName)
            for c in self.model.getGenConstrs():
                if c.IISGenConstr:
                    logs.critical("%s", c.GenConstrName)
        else:
            logs.debug("Model Successfully Solved.")

        runtime = float(self.model.Runtime)
        nodes = float(self.model.NodeCount)
        gap_pct = 100.0 * float(self.model.MIPGap)
        if self.model.SolCount > 0:
            obj_str = f"{float(self.model.ObjVal):.4g}"
        else:
            obj_str = "n/a"

        logs.info(
            "%s solve: status=%d runtime=%.1fs nodes=%.0f gap=%.2f%% "
            "obj=%s build=%.1fs presolve=%d",
            round_label,
            int(self.model.Status),
            runtime,
            nodes,
            gap_pct,
            obj_str,
            build_secs,
            int(presolve),
        )

    def run_model(self):
        """Construct and solve the Gurobi model (with optional bonus round)."""
        self._round1_obj_val = None
        self._round1_weighted_theta = None
        self._round2_slots_by_program = None
        self._hold_fill_alpha = 0.0
        self._open_throttle_grace = None
        self._solver_method, self._solver_norel_budget, self._solver_milp_time = (
            self._log_solver_config_once()
        )

        round_steps = [("Round1", self.build_model_round1)]
        if self.config.getboolean("semester", "run_bonus_round"):
            round_steps.extend(
                [
                    ("Round2", self.build_model_round2_priority_NEW),
                    ("Round3", self.build_model_round3_priority),
                    ("Round4", self.build_model_round4_priority_NEW),
                ]
            )
        if self.config.getboolean("semester", "run_upcoming_night_round", fallback=False):
            round_steps.append(("UpcomingNight", self.build_model_upcoming_night_round))

        round_total = len(round_steps)
        for round_num, (round_label, build_fn) in enumerate(round_steps, start=1):
            self._begin_round(round_label, round_num=round_num, round_total=round_total)
            t_build = time.time()
            build_fn()
            self.optimize_model(round_label, build_secs=time.time() - t_build)
            if round_label == "Round1":
                self._round1_obj_val = self.model.objVal
                self._round1_weighted_theta = self._eval_weighted_theta()
            self._finalize_round(round_label)

        logs.info("Scheduling complete, clear skies!")

    def _finalize_round(self, round_label):
        """Build schedule, log report, persist per-night handoff + snapshot."""
        self.build_schedule()
        if round_label == "Round2":
            self._capture_round2_slots_by_program()
        self.log_report(round_label)
        self.write_request_selected()
        self.to_hdf5()

    # ==================================================================
    # Output.
    # ==================================================================

    def build_schedule(self):
        """
        Build the sparse schedule DataFrame from ``self.Yrds`` and write
        ``semester_plan.csv``.

        Sets ``self.schedule`` to a DataFrame with columns
        ``unique_id, d, s, target`` -- one row per scheduled exposure start.
        """
        df = pd.DataFrame(self.Yrds.keys(), columns=["unique_id", "d", "s"])
        df["value"] = [self.Yrds[k].x for k in self.Yrds.keys()]
        sparse = df.query("value > 0").drop(columns=["value"]).copy()
        sparse = sparse.merge(
            self.requests_frame[["unique_id", "target"]],
            on="unique_id",
            how="left",
        )
        sparse["target"] = sparse["target"].fillna("NO MATCHING NAME")
        sparse.to_csv(
            os.path.join(self.output_directory, "semester_plan.csv"),
            index=False,
            na_rep="",
        )
        self.schedule = sparse

    def to_string(self, round_label="Round1", *, header="Semester Planner Statistics"):
        """Run report: top-level summary Series + per-program hours DataFrame.

        Requires that :meth:`build_schedule` has been called so
        ``self.schedule`` is set.
        """
        if self.schedule is None:
            raise RuntimeError("call build_schedule() before to_string()")

        slot_size = self.config.getfloat("semester", "slot_size")
        hours_per_night = self.config.getfloat("semester", "hours_per_night")
        slots_per_hour = 60 / slot_size

        def slot_demand_slots(frame):
            return int(
                (
                    frame["t_visit_slots"]
                    * frame["n_intra_max"]
                    * frame["n_inter_max"]
                ).sum()
            )

        # ---- top-level summary as a Series ----
        today_idx = self.today_starting_night
        is_alloc_2d = self.access_record["is_allocated"][0]
        allocated = int(is_alloc_2d.sum())
        allocated_future = int(is_alloc_2d[today_idx:].sum())
        allocated_today = int(is_alloc_2d[today_idx].sum())

        sched = self.schedule
        sched_future = sched[sched["d"] >= today_idx]
        sched_today = sched[sched["d"] == today_idx]
        t_visit_slots = self.requests_frame.set_index("unique_id")["t_visit_slots"]
        slots_per_visit_future = sched_future["unique_id"].map(t_visit_slots).fillna(1)
        slots_per_visit_today = sched_today["unique_id"].map(t_visit_slots).fillna(1)
        future_reserved = int(slots_per_visit_future.sum())
        today_reserved = int(slots_per_visit_today.sum())

        summary = pd.Series(
            {
                "Total requests": len(self.requests_frame_all),
                "Total requests (active)": len(self.requests_frame),
                "Total allocated slots": allocated,
                "Total slots requested": slot_demand_slots(self.requests_frame_all),
                "Total slots requested (active)": slot_demand_slots(
                    self.requests_frame
                ),
                "Future allocated slots": allocated_future,
                "Future reserved slots": future_reserved,
                "Future fill factor": (
                    100 * future_reserved / allocated_future
                    if allocated_future
                    else 0.0
                ),
                "Current day allocated slots": allocated_today,
                "Current day reserved slots": today_reserved,
                "Current day fill factor": (
                    100 * today_reserved / allocated_today
                    if allocated_today
                    else 0.0
                ),
            }
        )

        # ---- per-program table (hours) ----
        awarded = self._programs_frame["nights"] * hours_per_night
        awarded_slots = self._program_awarded_slots

        rf = self.requests_frame.copy()
        rf["requested_h"] = (
            rf["t_visit_slots"] * rf["n_intra_max"] * rf["n_inter_max"]
        ) / slots_per_hour
        requested_by_prog = rf.groupby("program_code")["requested_h"].sum()
        past_slots_by_prog = self._past_slots_by_program()
        past_by_prog = (
            pd.Series(past_slots_by_prog, dtype="float64") / slots_per_hour
        )

        sched_with_prog = sched.merge(
            self.requests_frame[["unique_id", "program_code", "t_visit_slots"]],
            on="unique_id",
            how="left",
        )
        sched_with_prog["scheduled_h"] = (
            sched_with_prog["t_visit_slots"].fillna(1) / slots_per_hour
        )
        scheduled_h = sched_with_prog.groupby("program_code")["scheduled_h"].sum()

        table = (
            pd.DataFrame({"aw": awarded})
            .join(requested_by_prog.rename("req"), how="left")
            .join(past_by_prog.rename("past"), how="left")
            .join(scheduled_h.rename("sched"), how="left")
            .fillna(0.0)
        )
        table["proj"] = table["past"] + table["sched"]

        throttle_grace = self._throttle_grace_for_round(round_label)
        hold_active = round_label in ("Round3", "Round4", "UpcomingNight")

        miff_pct = []
        maff_pct = []
        for prog, row in table.iterrows():
            aw_h = float(row["aw"])
            aw_slots = float(awarded_slots.get(prog, 0.0))
            past_slots = float(past_slots_by_prog.get(prog, 0))

            if aw_slots > 0:
                grace_slots = int(aw_slots * throttle_grace)
                if grace_slots < past_slots:
                    grace_slots = int(past_slots)
                maff_pct.append(100.0 * grace_slots / aw_slots)
            else:
                maff_pct.append(0.0)

            if hold_active and self._round2_slots_by_program is not None:
                r2_slots = float(self._round2_slots_by_program.get(prog, 0.0))
                min_slots = past_slots + r2_slots * (1.0 - self._hold_fill_alpha)
                if aw_slots > 0:
                    miff_pct.append(100.0 * min_slots / aw_slots)
                else:
                    miff_pct.append(0.0)
            else:
                miff_pct.append(0.0)

        table["miff%"] = miff_pct
        table["maff%"] = maff_pct
        aw_col = table["aw"]
        has_aw = aw_col > 0
        table["past%"] = np.where(has_aw, 100.0 * table["past"] / aw_col, 0.0)
        table["proj%"] = np.where(has_aw, 100.0 * table["proj"] / aw_col, 0.0)
        table = table.sort_index()

        for prog, row in table.iterrows():
            miff = float(row["miff%"])
            maff = float(row["maff%"])
            proj = float(row["proj%"])
            if proj < miff - 0.05 or proj > maff + 0.05:
                logs.warning(
                    "Program %s: proj%%=%.1f outside [miff%%=%.1f, maff%%=%.1f] "
                    "in %s report",
                    prog,
                    proj,
                    miff,
                    maff,
                    round_label,
                )

        display = table[["aw", "req", "past", "proj", "miff%", "maff%", "past%", "proj%"]].copy()
        program_table = display.copy()
        for col in ("aw", "req", "past", "proj"):
            program_table[col] = program_table[col].map(lambda x: f"{x:.1f}")
        for col in ("miff%", "maff%", "past%", "proj%"):
            program_table[col] = program_table[col].map(lambda x: f"{x:.1f}%")

        divider = "-" * 54
        stats_divider = "-" * 19 + " Program Statistics " + "-" * 19
        parts = [
            header,
            divider,
            summary.to_string(float_format=lambda x: f"{int(round(x))}"),
            "",
            stats_divider,
            program_table.to_string(),
            "",
            _PROGRAM_STATS_KEY,
            "",
        ]
        return "\n".join(parts) + "\n"

    def log_report(self, round_label):
        """Emit the run-report text to stdout (no log prefix on table lines)."""
        if round_label == "UpcomingNight" and self._round1_weighted_theta is not None:
            slack = self.config.getfloat(
                "semester", "global_shortfall_slack",
                fallback=_DEFAULT_GLOBAL_SHORTFALL_SLACK,
            )
            cap = self._round1_weighted_theta * slack
            current_theta = self._eval_weighted_theta()
            logs.info(
                "UpcomingNight: Round-1 weighted shortfall=%.3f cap=%.3f "
                "post-round weighted shortfall=%.3f",
                self._round1_weighted_theta,
                cap,
                current_theta,
            )
        report = self.to_string(round_label)
        if report:
            logs.info("Run report (%s):", round_label)
            print(report.rstrip(), flush=True)

    def write_request_selected(self):
        """Write ``request_selected.csv`` -- the handoff to ``NightPlanner``."""
        today_idx = self.all_dates_dict[self.config.get("global", "current_day")]
        selected = {
            k[0] for k, v in self.Yrds.items() if v.x > 0 and k[1] == today_idx
        }
        selected_df = self.requests_frame[
            self.requests_frame["unique_id"].isin(selected)
        ].copy()
        selected_df["nplan_weight"] = 1.0
        selected_df.to_csv(
            os.path.join(self.output_directory, "request_selected.csv"),
            index=False,
        )

    # ==================================================================
    # Serialization
    # ==================================================================

    def to_hdf5(self, hdf5_path=None):
        """Persist ``config_ini_text`` + a few DataFrames + ``access_record``.

        Write is atomic: the snapshot lands at a ``.tmp`` sibling and is
        renamed into place once all writes succeed. A crash midway never
        clobbers the previous snapshot.

        Args:
            hdf5_path (str, optional): defaults to
                ``<output_directory>/semester_planner.h5``.
        """
        if hdf5_path is None:
            hdf5_path = os.path.join(self.output_directory, "semester_planner.h5")
        tmp_path = hdf5_path + ".tmp"
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

        self.requests_frame_all.to_hdf(
            tmp_path, key="requests_frame_all", mode="a", format="table"
        )
        past_fmt = "fixed" if self.past_df.empty else "table"
        self.past_df.to_hdf(tmp_path, key="past_df", mode="a", format=past_fmt)
        if self.schedule is not None:
            fmt = "fixed" if self.schedule.empty else "table"
            self.schedule.to_hdf(tmp_path, key="schedule", mode="a", format=fmt)

        with h5py.File(tmp_path, "a") as f:
            f.attrs["schema_version"] = SEMESTER_PLANNER_H5_SCHEMA
            f.attrs["config_ini_text"] = self._config_ini_text
            f.create_dataset(
                "access_record", data=self.access_record, compression="gzip"
            )

        os.replace(tmp_path, hdf5_path)
        logs.info(f"SemesterPlanner saved to HDF5: {hdf5_path}")
        return hdf5_path

    @classmethod
    def from_hdf5(cls, hdf5_path):
        """Rehydrate a SemesterPlanner from a snapshot written by :meth:`to_hdf5`.

        Skips Gurobi-only state (model, Yrds, Wrd, theta) and the constraint
        lookup tables -- those are only meaningful when solving. Downstream
        consumers (plot.py, nplan.py) read requests_frame*, schedule,
        access_record, past_df, and queue, all of which are restored.
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
            access_record = f["access_record"][:].view(np.recarray)

        requests_frame_all = pd.read_hdf(hdf5_path, key="requests_frame_all")
        try:
            past_df = pd.read_hdf(hdf5_path, key="past_df")
        except KeyError:
            past_df = pd.DataFrame(columns=PAST_COLS)
        try:
            schedule = pd.read_hdf(hdf5_path, key="schedule")
        except KeyError:
            schedule = None

        instance = cls.__new__(cls)
        instance._config_ini_text = config_ini_text
        instance.config = ConfigParser()
        instance.config.read_string(config_ini_text)
        instance.queue = astroq.queue.from_config(instance.config)

        workdir = instance.config.get("global", "workdir")
        instance.output_directory = os.path.join(workdir, "outputs")
        instance.allocation_file = instance._resolve_path("allocation_file")
        instance.custom_file = instance._resolve_path("custom_file")
        instance.programs_file = instance._resolve_path("programs_file")

        # Re-derive slot columns on the full frame (active + inactive) so the
        # throttle can count past usage on every row. Columns are pure
        # functions of the persisted data, so we don't ship them on disk. The
        # persisted frame was already validated when written; no repair here.
        instance._attach_slot_columns(requests_frame_all)
        instance.requests_frame_all = requests_frame_all
        instance.requests_frame = (
            requests_frame_all[~requests_frame_all["inactive"]]
            .reset_index(drop=True)
            .copy()
        )

        instance.past_df = past_df
        instance._attach_past_columns()
        instance.access_obj = ac.Access.from_planner(instance)
        instance.access_record = access_record
        instance.schedule = schedule

        logs.info(f"SemesterPlanner loaded from HDF5: {hdf5_path}")
        return instance

    def constraint_hold_program_fill_factors(self, alpha=0.0):
        """
        Bonus round constraint: not featured in Lubin et al. 2025.

        After Round 2, record each program's total scheduled slot-time
        (``sum(Yrds * t_visit_slots)``). In Round 3, constrain the program's
        scheduled slot-time to stay within alpha (relative tolerance) of that
        Round 2 value.
        """
        self._hold_fill_alpha = float(alpha)
        logs.info("Constraint: Holding program fill factors.")

        program_slot_value = self._program_slot_value()
        program_slot_expr = self._program_slot_expr

        for p in self._active_uids_by_program:
            r2_slots_awarded = program_slot_value.get(p, 0.0)
            logs.info(
                f"Holding program {p} to at least {alpha * 100:.1f}% less than "
                f"Round 2 scheduled slots: {r2_slots_awarded:.0f}"
            )
            future_slots = program_slot_expr.get(p, 0)
            self.model.addConstr(
                r2_slots_awarded * (1-alpha)  <= future_slots,
                'hold_program_fill_factors_lower_' + p,
            )
            # Removed because we want to fill as much as possible
            # self.model.addConstr(
            #     r2_slots_awarded * (1+alpha)  >= future_slots,
            #     'hold_program_fill_factors_upper_' + p,
            # )

    def build_model_round2_priority(self):
        """
        Implement the constraints and objective function for Round 2. Not described in Lubin et al. 2025.

        Returns:
            None
        """
        t1 = time.time()
        self.constraint_fix_previous_objective()
        self.set_objective_priorities_INTER()
        logs.debug(f"Time to build constraints: {np.round(time.time()-t1,3):.3f}")

    def build_model_round3_priority(self):
        """
        Implement the constraints and objective function for Round 2. Not described in Lubin et al. 2025.

        Returns:
            None
        """
        t1 = time.time()
        self.constraint_hold_program_fill_factors()
        self.set_objective_priorities_INTRA()
        logs.debug(f"Time to build constraints: {np.round(time.time()-t1,3):.3f}")

    def set_objective_priorities_INTER(self):
        """
        Set inter-program priority objective:


        where N_p = sum over r in R_p of 2^{w_r} * t_{visit,r} * n_{intra,max,r} * n_{inter,max,r}
        (normalization per program), w_r = weight of request r, Y_{r,d,s} = binary schedule variable.
        """
        logs.info("Objective: Inter-program priorities.")

        program_request_ids = self._active_uids_by_program

        program_frame = self._programs_frame
        if "priority" not in program_frame.columns:
            logs.warning(
                f"{self.programs_file} has no 'priority' column; using priority 1.0 for all programs."
            )
            N_p = {p: 1.0 for p in program_request_ids}
        else:
            priority_by_program = program_frame["priority"].astype(float).to_dict()
            N_p = {p: priority_by_program[p] for p in program_request_ids}

        self.model.setObjective(
            gp.quicksum(
                N_p[p] * gp.quicksum(
                    self.Yrds[r, d, s]
                    for r, d, s in self.observability_tuples
                    if r in program_request_ids[p]
                )
                for p in program_request_ids
            ),
            GRB.MAXIMIZE
        )

    def set_objective_priorities_INTRA(self):
        """
        Set intra-program priority objective:


        """
        logs.info("Objective: Intra-program priorities.")

        weight_by_id = self.requests_frame.set_index('unique_id')['splan_weight']
        t_visit_slots = self._t_visit_slots_by_uid
        program_request_ids = self._active_uids_by_program

        self.model.setObjective(
            gp.quicksum(
                gp.quicksum(
                    (1.0 / weight_by_id.loc[r])
                    * self.Yrds[r, d, s]
                    * t_visit_slots[r]
                    for r, d, s in self.yrds_tuples
                    if r in program_request_ids[p]
                )
                for p in program_request_ids
            ),
            GRB.MAXIMIZE
        )