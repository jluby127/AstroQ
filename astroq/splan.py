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
import astroq.io
import astroq.queue
from astroq.io import PAST_COLS

logs = logging.getLogger(__name__)

# Schema for h5 serialization bump when the on-disk layout changes
SEMESTER_PLANNER_H5_SCHEMA = 6

# Request columns denormalized onto the observability grid to form
# ``request_slots`` -- the single relational table build_model is defined over.
STRATEGY_COLS = [
    "unique_id",
    "target",
    "program_code",
    "n_intra_min",
    "n_intra_max",
    "n_inter_max",
    "tau_inter",
    "t_visit_slots",
    "tau_intra_slots",
]


_ROUND4_THROTTLE_GRACE = 2.0

# Default cap on Round-1-optimum multiplier for the upcoming-night round's
# global shortfall constraint. Referenced from both build_model_upcoming_night
# and log_report; keep as one constant so the two never drift apart.
_DEFAULT_GLOBAL_SHORTFALL_SLACK = 1.1

# label -> (banner description, per-round build method). Each build method
# layers round-specific constraint deltas + an objective on the structural
# model from build_model().
_ROUND_SPECS = {
    "Round1": (
        "Minimize time-weighted shortfall (Lubin et al.)",
        "build_model_round1",
    ),
    "Round2": (
        "Maximize inter-program fill factors",
        "build_model_round2",
    ),
    "Round3": (
        "Intra-program priorities (hold Round-2 fill)",
        "build_model_round3",
    ),
    "Round4": (
        "Minimize empty slots (re-throttle)",
        "build_model_round4",
    ),
    "UpcomingNight": (
        "Fill current night (cap global shortfall)",
        "build_model_upcoming_night",
    ),
}

# The two supported scheduling modes, selected by ``[semester] mode``.
MODE_SEQUENCES = {
    "round1": ("Round1",),
    "full": ("Round1", "Round2", "Round3", "Round4", "UpcomingNight"),
}


def _initial_round_state():
    """Cross-round values captured as rounds solve. ``throttle_grace``
    mirrors the grace factor currently enforced on the model."""
    return {
        "round1_weighted_theta": None,
        "round2_slots_by_program": None,
        "hold_fill_alpha": 0.0,
        "throttle_grace": 1.0,
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
        self.round_state = _initial_round_state()

        self.requests = self._load_frame("request")
        self._attach_slot_columns(self.requests)

        self.past = self._load_frame("past")
        self.allocation = self._load_frame("allocation")
        self.custom = self._load_frame("custom")
        self.programs = self._load_frame("programs")

        # Per-request derived columns that depend on past live on
        # requests (single source of truth, no parallel dict
        # attributes). Per-program budget/past aggregates live on
        # ``self.programs`` via :meth:`_attach_program_columns`.
        self._attach_past_columns()
        self._attach_program_columns()

        # Observability cube (single source of truth for which slots are valid).
        self.access_obj = ac.Access.from_planner(self)
        self.access_record = self.access_obj.build_access()
        self.observability = self.access_obj.observability(
            self.access_record.is_observable
        )

        self._log_boost_current_day_slots()

        self.build_model()

        logs.debug("Initializing complete.")

    def _load_frame(self, kind):
        """Load a validated CSV frame. ``kind`` maps to ``{kind}_file`` in config."""
        key = f"{kind}_file"
        raw = self.config.get("semester", key)
        path = raw if os.path.isabs(raw) else os.path.join(self.workdir, raw)
        return astroq.io.read_csv(path, kind)

    def _ensure_output_dir(self):
        os.makedirs(self.output_directory, exist_ok=True)

    # ------------------------------------------------------------------
    # Properties (date-derived; paths derived from config).
    # ------------------------------------------------------------------

    @property
    def workdir(self):
        return self.config.get("global", "workdir")

    @property
    def output_directory(self):
        return os.path.join(self.workdir, "outputs")

    @property
    def requests_active(self):
        return self.requests[~self.requests["inactive"]]

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

    def _log_boost_current_day_slots(self):
        """Report observable slot counts on current_day for each boosted target."""
        if self.boost is None:
            return
        current_day = self.config.get("global", "current_day")
        d_today = self.today_starting_night
        boost_by_uid = self._boost_by_uid
        uid_to_target = dict(
            zip(
                self.requests["unique_id"],
                self.requests["target"],
            )
        )
        factor = next(iter(boost_by_uid.values()))
        logs.info(
            "Boost on current_day=%s (d=%d), factor=%g:",
            current_day,
            d_today,
            factor,
        )
        obs_uids = self.observability["unique_id"]
        obs_d = self.observability["d"]
        for uid in boost_by_uid:
            n_slots = int(((obs_uids == uid) & (obs_d == d_today)).sum())
            target = uid_to_target.get(uid, "(unknown unique_id)")
            logs.info(
                "  %s (%s): %d observable slot(s) on current_day",
                uid,
                target,
                n_slots,
            )

    def build_model(self):
        """Gurobi variables plus every structural constraint, inline.

        First builds the relational tables the model is defined over:

        - ``request_slots`` -- one row per observable ``(unique_id, d, s)``: the
          observability grid with the request's :data:`STRATEGY_COLS` joined on,
          plus an ``rds`` column holding that row's ``(unique_id, d, s)`` key
          into ``self.Yrds``. The single table the structural constraints,
          throttle, and priority rounds are defined over. Gurobi variables are
          looked up as ``self.Yrds[k]`` at constraint-build time, never stored
          in the frame.
        - ``yrds_tuples`` -- the same ``(unique_id, d, s)`` set as plain tuples
          (the ``Yrds`` variable keys).
        - ``schedulable_uids`` -- unique_ids with at least one observable slot,
          in first-appearance order.

        Everything here is required by every scheduling mode; rounds layer
        objectives and round-specific deltas on top (``build_model_round*``).
        Paper map (Lubin et al. 2025 -> code): ``Y_{r,d,s}`` -> ``Yrds``,
        ``W_{r,d}`` -> ``Wrd``, shortfall -> ``theta``; constraint numbers
        below refer to that paper. Set-building is relational (merges /
        groupbys over ``request_slots``); Python loops only emit ``addConstr``.
        """
        t0 = time.time()
        self.model = gp.Model("Semester_Scheduler")

        # ---- relational tables the model is defined over ----
        self.request_slots = self.observability.merge(
            self.requests_active[STRATEGY_COLS], on="unique_id"
        )
        # rds: each row's (unique_id, d, s) key into self.Yrds. Pure data, so the
        # frame never holds Gurobi objects; constraints look up self.Yrds[k].
        self.request_slots["rds"] = list(
            zip(
                self.request_slots["unique_id"],
                self.request_slots["d"],
                self.request_slots["s"],
            )
        )
        self.yrds_tuples = list(
            self.observability.itertuples(index=False, name=None)
        )
        self.schedulable_uids = list(self.observability["unique_id"].unique())

        # diagnostics: requests with no observable slot are absent from the model
        schedulable = set(self.schedulable_uids)
        all_requests = list(self.requests_active["unique_id"])
        missing = sum(uid not in schedulable for uid in all_requests)
        logs.warning(
            f"There are {missing} targets out of {len(all_requests)} that have "
            f"no valid day/slot pairs and therefore are effectively removed "
            f"from the model."
        )

        rs = self.request_slots

        # ---- variables ----
        self.Yrds = self.model.addVars(
            self.yrds_tuples, vtype=GRB.BINARY, name="Requests_Slots"
        )
        wrd_keys = rs.query("n_intra_max > 1")[["unique_id", "d"]].drop_duplicates()
        if not wrd_keys.empty:
            self.Wrd = self.model.addVars(
                list(wrd_keys.itertuples(index=False, name=None)),
                vtype=GRB.BINARY,
                name="OnSky",
            )
        self.theta = self.model.addVars(
            list(self.requests_active["unique_id"]), name="Shortfall"
        )

        # ---- Eq. 3: shortfall definition. theta_r >= remaining nights owed
        # minus scheduled nights (visits / n_intra_max); lb=0 from addVars. ----
        logs.info("Constraint: Build theta variable")
        rf_indexed = self.requests_active.set_index("unique_id")
        for uid, grp_keys in rs.groupby("unique_id", sort=False)["rds"]:
            row = rf_indexed.loc[uid]
            self.model.addConstr(
                self.theta[uid]
                >= row["n_inter_max"]
                - row["past_nights_observed"]
                - gp.quicksum(self.Yrds[k] for k in grp_keys) / row["n_intra_max"],
                f"greater_than_nobs_shortfall_{uid}",
            )

        # ---- Reserve slots for multi-slot exposures. Constraint 1 in Lubin et al.
        # (2026)
        #
        # We forbid any request from being scheduled at (d,s) if there is a *previous*
        # request that started within t_visit_slots of (d,s). This prevents two two
        # requests from overlpaping. 
        #
        # Note on implementation: the same behavior can be achieved by requiring that no
        # slots be schedule t_visit_slots after a multi-slot exposure. However, since
        # most exposures are multi-slot, nearly every (r,d,s) results in a seperate
        # constraint. In earlier testing, this resulted a long presolve. This
        # implementation introduces a constraint for every unique (d,s).
        #
        # The rs_multislot has one row per (r,d,s,s_future) where s_future is a a slot
        # within t_visit_slots of r,d,s. 
        logs.info("Constraint: Reserve slots for multi-slot exposures.")
        max_t_visit = int(rs["t_visit_slots"].max())
        deltas = pd.DataFrame({"delta": np.arange(1, max_t_visit)})  # [] if all t == 1
        rs_multislot = (
            rs.query("t_visit_slots > 1")[["d", "s", "t_visit_slots", "rds"]]
            .merge(deltas, how="cross") 
            .query("delta < t_visit_slots")
            .assign(s_future=lambda f: f["s"] + f["delta"])
            .set_index(["d", "s_future"])
            .sort_index()
        )

        for ds_on, rds_on in rs.groupby(["d", "s"], sort=False)["rds"]:
            # rds_off: r,d,s of visits that cover ds_on. Empty when no exposure reaches
            # this slot -- the constraint then collapses to one start per (d, s), still
            # required so two requests can't share a slot.
            rds_off = rs_multislot["rds"].loc[[ds_on]] if ds_on in rs_multislot.index else []
            self.model.addConstr(
                gp.quicksum(self.Yrds[rds] for rds in rds_on)
                + gp.quicksum(self.Yrds[rds] for rds in rds_off)
                <= 1,
                "reserve_multislot_{0}d_{1}s".format(*ds_on),
            )

        # ---- Enforce desired maximum unique nights. Constraint 2 in Lubin et al.
        # (2026).
        #
        # Cap scheduled observations at desired_max_obs per target. Single/multi-visit
        # requests are treated differently. For single visit requests. We require the
        # sum of the future scheduled visits not exceed the max value using the Yrds
        # variable (Lubin et al. 2026 Eq. 7). For multi-visit requests, Yrds is replaced
        # by Wrd since there may be multiple Yrds on a single night (Lubin et al. 2026
        # Eq. 8). We handle these cases seperately since there are usually few multi-
        # visit requrests, and thus we can keep the Wrd varaible small.
        logs.info("Constraint: Enforce maximum n_inter_max")
        desired_max_obs = rf_indexed["desired_max_obs"]
        r_maxnights_single = (
            rs.query("n_intra_max == 1")
            .groupby("unique_id", sort=False)["rds"]
            .agg(list)
        )
        for r, rds_keys in r_maxnights_single.items():
            self.model.addConstr(
                gp.quicksum(self.Yrds[rds] for rds in rds_keys) <= desired_max_obs.loc[r],
                "max_desired_unique_nights_for_request_{0}".format(r),
            )
        r_maxnights_multi = (
            rs.query("n_intra_max > 1")
            .drop_duplicates(["unique_id", "d"])
            .groupby("unique_id", sort=False)["d"]
            .agg(list)
        )
        for r, days in r_maxnights_multi.items():
            self.model.addConstr(
                gp.quicksum(self.Wrd[r, d] for d in days) <= desired_max_obs.loc[r],
                "max_desired_unique_nights_for_request_{0}".format(r),
            )

        # ---- Enforce inter-night cadence Constraint 3 in Lubin et al. (2026).
        #
        # If a request is observed on night d, it must be turned off on every future
        # night d_future with d < d_future < d + tau_inter. rd_internight maps
        # (unique_id, d) -> the rds that are forbidden when the request is observed on
        # night d.
        logs.info("Constraint: Enforce inter-night cadence.")
        rd_internight = (
            rs.query("tau_inter > 1")[["unique_id", "d", "tau_inter"]]
            .drop_duplicates(["unique_id", "d"])
            .merge(rs[["unique_id", "d", "rds"]],suffixes=["", "_future"],on="unique_id")
            .query("d < d_future < d + tau_inter")
            .groupby(["unique_id", "d"], sort=False)["rds"]
            .agg(list)
        )
        for rd_on, grp in rs.groupby(["unique_id", "d"], sort=False):
            if rd_on not in rd_internight.index:
                continue
            rds_on = grp["rds"] # (r,d,s) of request r on night d
            rds_off = rd_internight.loc[rd_on] # forbidden future (r,d,s)
            n_intra_max = grp["n_intra_max"].iloc[0]

            # term1: is request r observed on night d? Divide by n_intra_max so a
            # fully-fired multi-shot night reads as exactly 1 (a single visit as
            # 1/n_intra_max) while still permitting all n_intra_max shots. term2:
            # is r observed anywhere in the forbidden window? In any integer
            # solution term1 >= 1/n_intra_max forces term2 = 0.
            term1 = gp.quicksum(self.Yrds[rds] for rds in rds_on) / n_intra_max
            term2 = gp.quicksum(self.Yrds[rds] for rds in rds_off)
            self.model.addConstr(
                term1 + term2 <= 1,
                "enforce_internight_cadence_{0}_{1}".format(*rd_on),
            )

        # ---- Enforce intra-night cadence. Constraint 4 Lubin et al. (2026). For every
        # (r,d,s) of a multi-visit request, find all future (r,d,s_future) where s <
        # s_future < s + tau_intra_slots. If (r,d,s) is on then all (r,d,s_future) must
        # be off. rs_intranight contains one row per (r,d,s) with a list all forbidden
        # (r,d,s_future) 
        logs.info("Constraint: Enforce intra-night cadence.")
        rs_intranight = (
            pd.merge(
                rs.query("n_intra_max > 1")[["unique_id", "d", "s", "tau_intra_slots"]],
                rs.query("n_intra_max > 1")[["unique_id", "d", "s", "rds"]],
                suffixes=["", "_future"],   # left s -> s, right s -> s_future
                on=["unique_id", "d"],
            )
            .query("s < s_future < s + tau_intra_slots")
            .groupby(["unique_id", "d", "s"], sort=False)["rds"]
            .agg(list)
        )
        for rds_on, rds_off in rs_intranight.items():
            self.model.addConstr(
                self.Yrds[rds_on] + gp.quicksum(self.Yrds[rds] for rds in rds_off)
                <= self.Wrd[rds_on[:2]],
                "enforce_intranight_cadence_{0}_{1}d_{2}s".format(*rds_on),
            )

        # ---- Enforce min/max visits per night. Constraint 5 in Lubin et al. (2026). On
        # days when a multi-visit request is scheduled Wrds = 1, the sum of the visits
        # must be between n_intra_min and n_intra_max. There is an omission in Lubin et
        # al (2026) where there is no explicit constraint applying this rule to
        # single-visit requests. We implement that here. 
        logs.info("Constraint: Enforce min/max visits per night.")
        rd_multi = (
            rs.query("n_intra_max > 1")
            .groupby(["unique_id", "d"], sort=False)
            .agg(
                rds=("rds", list),
                n_intra_min=("n_intra_min", "first"),
                n_intra_max=("n_intra_max", "first"),
            )
        )
        for rd_on, row in rd_multi.iterrows():
            r, d = rd_on
            n_visits_intra = gp.quicksum(self.Yrds[rds] for rds in row.rds)
            self.model.addConstr(
                n_visits_intra <= row.n_intra_max * self.Wrd[r, d],
                "enforce_max_visits1_{0}_{1}d".format(*rd_on),
            )
            self.model.addConstr(
                n_visits_intra >= row.n_intra_min * self.Wrd[r, d],
                "enforce_min_visits_{0}_{1}d".format(*rd_on),
            )

        rd_single = (
            rs.query("n_intra_max == 1")
            .groupby(["unique_id", "d"], sort=False)["rds"]
            .agg(list)
        )
        for rd_on, rds_on in rd_single.items():
            self.model.addConstr(
                gp.quicksum(self.Yrds[rds] for rds in rds_on) <= 1,
                "enforce_max_visits_{0}_{1}d".format(*rd_on),
            )

        # ---- Throttle: structural, but re-parameterized by later rounds
        # (Round 4 re-adds it with a wider grace), so it stays a method. ----
        self.constraint_throttle(throttle_grace=1.0)

        logs.info(f"Structural model built in {time.time() - t0:.3f}s")

    def _attach_past_columns(self):
        """Attach past-history aggregates and the max-obs cap to ``requests_frame``.

        Aggregates are indexed by ``unique_id`` over UT calendar nights
        (``timestamp[:10]``); missing uids default to 0 (or ``""``).
        ``desired_max_obs`` is the per-target night cap (Constraint 2); it
        collapses to ``past_nights_observed`` when a target is over-observed
        so the model stays feasible.
        """
        rf = self.requests_active
        uids = rf["unique_id"]

        if self.past.empty:
            agg = pd.DataFrame(
                {"nights": 0, "n_exp": 0, "last": ""}, index=uids,
            )
        else:
            night = self.past["timestamp"].str[:10]
            g = self.past.assign(_night=night).groupby("unique_id")
            agg = pd.DataFrame({
                "nights": g["_night"].nunique(),
                "n_exp": g.size(),
                "last": g["_night"].max(),
            }).reindex(uids).fillna({"nights": 0, "n_exp": 0, "last": ""})

        rf["past_nights_observed"] = agg["nights"].astype(int).to_numpy()
        rf["past_n_exposures"] = agg["n_exp"].astype(int).to_numpy()
        rf["past_date_last_observed"] = agg["last"].astype(str).to_numpy()

        n_max = rf["n_inter_max"].to_numpy()
        past = rf["past_nights_observed"].to_numpy()
        over = past > n_max
        rf["desired_max_obs"] = np.where(over, past, n_max - past).astype(int)

    def _attach_program_columns(self):
        """Attach semester budget and past-slot aggregates to ``self.programs``.

        ``past_slots`` counts ALL request rows (active + inactive), matching
        throttle semantics.
        """
        slot_size = self.config.getfloat("semester", "slot_size")
        hours_per_night = self.config.getfloat("semester", "hours_per_night")

        self.programs["awarded_hours"] = self.programs["nights"] * hours_per_night
        self.programs["awarded_slots"] = (
            self.programs["awarded_hours"] * 60 / slot_size
        )

        rfa = self.requests
        if self.past.empty:
            past_n = pd.Series(dtype="int64")
        else:
            past_n = self.past.groupby("unique_id").size()
        past_slots_by_uid = (
            rfa["unique_id"].map(past_n).fillna(0).astype(int) * rfa["t_visit_slots"]
        )
        past_by_prog = (
            pd.Series(past_slots_by_uid, index=rfa.index)
            .groupby(rfa["program_code"])
            .sum()
        )
        self.programs["past_slots"] = (
            past_by_prog.reindex(self.programs.index).fillna(0).astype(int)
        )

    @cached_property
    def _boost_by_uid(self):
        """dict[unique_id -> boost factor], or ``None`` if no boost was passed."""
        if self.boost is None:
            return None
        return dict(
            zip(self.boost["unique_id"].astype(str), self.boost["boost"].astype(float))
        )

    def _program_slot_value(self):
        """dict[program_code -> float] of scheduled slot-time at the
        *current* Gurobi solution. Not cached -- ``.X`` changes every round."""
        return {
            p: sum(
                self.Yrds[k].X * n for k, n in zip(g["rds"], g["t_visit_slots"])
            )
            for p, g in self.request_slots.groupby("program_code")
        }

    # ---- throttling & bonus round ----

    def constraint_throttle(self, throttle_grace=1.0):
        """
        Not described in Lubin et al. 2025.

        Ensure that no program is scheduled for more time than they bring to
        the queue (within a grace amount). Past usage is counted over ALL
        request rows (active and inactive) via ``self.programs["past_slots"]``,
        while only active targets contribute schedulable slots (inactive
        targets have no ``Yrds`` variables).
        """
        logs.info("Constraint: Throttling over-requested programs.")
        awarded_slots_grace_by_program = (
            self.programs["awarded_slots"] * throttle_grace
        ).astype(int)

        # Past budget: ALL rows (active + inactive).
        past_used_slots_by_program = self.programs["past_slots"]

        # Schedulable budget: only ACTIVE targets get Yrds variables
        # (request_slots is active-only).
        slot_expr_by_program = {
            p: gp.quicksum(
                self.Yrds[k] * n for k, n in zip(g["rds"], g["t_visit_slots"])
            )
            for p, g in self.request_slots.groupby("program_code")
        }

        clamped = []
        for program, awarded_slots_grace in awarded_slots_grace_by_program.items():
            awarded_slots_grace = int(awarded_slots_grace)
            schedulable_slots = slot_expr_by_program.get(program, 0)
            past_used = int(past_used_slots_by_program.get(program, 0))
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

    def constraint_fix_global_shortfall(self, slack_factor):
        """Upcoming-night round: cap weighted shortfall at Round-1 optimum * slack."""
        round1_theta = self.round_state["round1_weighted_theta"]
        if round1_theta is None:
            raise RuntimeError(
                "Round 1 must be solved before fixing global shortfall."
            )
        cap = round1_theta * slack_factor
        logs.info(
            "Constraint: weighted shortfall <= Round-1 optimum * %g "
            "(cap=%.3f from Round-1 weighted shortfall=%.3f)",
            slack_factor,
            cap,
            round1_theta,
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
        t_visit = self.requests_active.set_index("unique_id")["t_visit_slots"]
        return gp.quicksum(
            self.theta[uid] * t_visit[uid] for uid in self.schedulable_uids
        )

    def _eval_weighted_theta(self):
        """Evaluate weighted shortfall at the current Gurobi solution."""
        t_visit = self.requests_active.set_index("unique_id")["t_visit_slots"]
        return sum(
            self.theta[uid].X * t_visit[uid] for uid in self.schedulable_uids
        )

    def set_objective_minimize_theta_time_normalized(self):
        """See Equation 1 in Lubin et al. 2025."""
        theta_obj = self._weighted_theta_expr()
        if self.boost is not None:
            boost_by_uid = self._boost_by_uid
            d_today = self.today_starting_night
            boost_terms = [
                boost_by_uid[uid] * self.Yrds[uid, d, s]
                for uid, d, s in self.yrds_tuples
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

    def set_objective_maximize_slots_used_tonight(self):
        """Upcoming-night round: maximize filled slots on ``current_day``."""
        d_today = self.today_starting_night
        current_day = self.config.get("global", "current_day")
        logs.info(
            "Objective: Maximize slot usage on upcoming night current_day=%s (d=%d).",
            current_day,
            d_today,
        )
        t_visit = self.requests_active.set_index("unique_id")["t_visit_slots"]
        self.model.setObjective(
            gp.quicksum(
                t_visit[uid] * self.Yrds[uid, d, s]
                for uid, d, s in self.yrds_tuples
                if d == d_today
            ),
            GRB.MAXIMIZE,
        )

    def set_objective_minimize_empty_slots(self):
        """Bonus round: minimize empty slots."""
        logs.info("Objective: Minimize the number of empty slots.")
        t_visit = self.requests_active.set_index("unique_id")["t_visit_slots"]
        total_slots = self.semester_length * self.access_obj.nslots
        self.model.setObjective(
            (total_slots - gp.quicksum(
                t_visit[uid] * self.Yrds[uid, d, s]
                for uid, d, s in self.yrds_tuples
            )),
            GRB.MINIMIZE,
        )

    def build_model_round2(self):
        """Round 2: maximize the summed program fill factors, holding each
        program at or above its Round-1 fill. Not in Lubin et al. 2025."""
        program_slot_value = self._program_slot_value()

        fill_factor = {}
        for p, g in self.request_slots.groupby("program_code"):
            awarded_slots = self.programs["awarded_slots"].get(p)
            if awarded_slots is None or awarded_slots <= 0:
                logs.warning(
                    "Program %s missing or has non-positive awarded slots; "
                    "skipping fill factor.",
                    p,
                )
                continue
            awarded_slots = float(awarded_slots)
            slot_expr = gp.quicksum(
                self.Yrds[k] * n for k, n in zip(g["rds"], g["t_visit_slots"])
            )
            fill_factor[p] = slot_expr / awarded_slots

            slots_used_val = program_slot_value.get(p, 0.0)
            logs.info(
                "Round2 priority: program %s awarded_slots=%.0f "
                "slots_used=%.0f fill_factor=%.3f",
                p,
                awarded_slots,
                slots_used_val,
                slots_used_val / awarded_slots,
            )
            self.model.addConstr(
                slots_used_val / awarded_slots <= fill_factor[p],
                "maintain_fill_factor_for_program_" + p,
            )

        self.model.setObjective(
            gp.quicksum(fill_factor.values()), GRB.MAXIMIZE
        )

    def constraint_hold_program_fill_factors(self, alpha=0.0):
        """Constrain each program's scheduled slot-time to stay within
        ``alpha`` (relative tolerance) of its value at the current solution.
        Applied by Round 3 to hold the Round-2 fills."""
        self.round_state["hold_fill_alpha"] = float(alpha)
        logs.info("Constraint: Holding program fill factors.")

        program_slot_value = self._program_slot_value()

        for p, g in self.request_slots.groupby("program_code"):
            r2_slots = program_slot_value.get(p, 0.0)
            slot_expr = gp.quicksum(
                self.Yrds[k] * n for k, n in zip(g["rds"], g["t_visit_slots"])
            )
            logs.info(
                f"Holding program {p} to at least {alpha * 100:.1f}% less than "
                f"Round 2 scheduled slots: {r2_slots:.0f}"
            )
            self.model.addConstr(
                r2_slots * (1 - alpha) <= slot_expr,
                "hold_program_fill_factors_lower_" + p,
            )

    def build_model_round3(self):
        """Round 3: hold Round-2 program fills, optimize intra-program
        priorities (``splan_weight``). Not in Lubin et al. 2025."""
        self.constraint_hold_program_fill_factors()
        self.set_objective_priorities_intra()

    def set_objective_priorities_intra(self):
        """Maximize slot-time weighted by inverse ``splan_weight`` (lower
        weight = higher priority) within each program."""
        logs.info("Objective: Intra-program priorities.")

        rf_by_uid = self.requests_active.set_index("unique_id")
        weight_by_id = rf_by_uid["splan_weight"]
        t_visit = rf_by_uid["t_visit_slots"]
        uids_by_program = (
            self.requests_active.groupby("program_code")["unique_id"]
            .apply(set)
            .to_dict()
        )

        self.model.setObjective(
            gp.quicksum(
                gp.quicksum(
                    (1.0 / weight_by_id.loc[r])
                    * self.Yrds[r, d, s]
                    * t_visit[r]
                    for r, d, s in self.yrds_tuples
                    if r in uids_by_program[p]
                )
                for p in uids_by_program
            ),
            GRB.MAXIMIZE,
        )

    def remove_constraint_throttle(self):
        """Remove throttle constraints from a prior round."""
        logs.info("Constraint: Removing previous throttle constraints.")
        for program in self.programs.index:
            rm_const = self.model.getConstrByName(f"throttle_program_{program}")
            if rm_const is not None:
                self.model.remove(rm_const)

    def constraint_fix_previous_completion_rates(self):
        """Keep each target's shortfall no worse than the prior round's
        solved value (lower ``theta`` is better)."""
        logs.info(
            "Constraint: Per-target shortfall must stay at or below prior round."
        )
        for uid in self.schedulable_uids:
            self.model.addConstr(
                self.theta[uid] <= float(self.theta[uid].X),
                f"theta_le_prior_{uid}",
            )

    def build_model_round4(self):
        """Round 4: minimize empty slots with a widened throttle, holding
        per-target shortfalls at their Round-3 values."""
        self.constraint_fix_previous_completion_rates()
        self.remove_constraint_throttle()
        self.constraint_throttle(throttle_grace=_ROUND4_THROTTLE_GRACE)
        self.round_state["throttle_grace"] = _ROUND4_THROTTLE_GRACE
        self.set_objective_minimize_empty_slots()

    # ==================================================================
    # Model orchestration.
    # ==================================================================

    def build_model_round1(self):
        """Round 1 objective per Lubin et al. 2025 (structural constraints
        already live on the model from :meth:`build_model`)."""
        self.set_objective_minimize_theta_time_normalized()

    def build_model_upcoming_night(self):
        """Upcoming-night round: cap global shortfall, fill ``current_day``."""
        slack = self.config.getfloat(
            "semester", "global_shortfall_slack", fallback=_DEFAULT_GLOBAL_SHORTFALL_SLACK
        )
        self.constraint_fix_global_shortfall(slack_factor=slack)
        self.set_objective_maximize_slots_used_tonight()

    def _resolve_mode(self):
        """Resolve ``[semester] mode`` to a round-label sequence.

        Two modes are supported: ``round1`` (default) and ``full`` (Rounds
        1-5). The pre-refactor boolean keys are rejected with a migration
        hint rather than silently ignored.
        """
        for legacy in ("run_bonus_round", "run_upcoming_night_round"):
            if self.config.has_option("semester", legacy):
                raise ValueError(
                    f"[semester] {legacy} was replaced by "
                    f"'mode = round1 | full'; update the config."
                )
        mode = self.config.get("semester", "mode", fallback="round1")
        mode = mode.strip().lower()
        if mode not in MODE_SEQUENCES:
            valid = " | ".join(MODE_SEQUENCES)
            raise ValueError(f"[semester] mode={mode!r} invalid; expected {valid}")
        return mode, MODE_SEQUENCES[mode]

    def _begin_round(self, round_label, *, round_num, round_total):
        """Log a visible banner at the start of each scheduling round."""
        desc, _ = _ROUND_SPECS[round_label]
        logs.info(
            "===== Semester Round %d/%d: %s =====", round_num, round_total, desc
        )

    def _log_solver_config_once(self, n_rounds):
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
            n_rounds,
        )
        return method, norel_budget, milp_time

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
        """Solve the configured round sequence over the structural model.

        The sequence comes from ``[semester] mode`` (:meth:`_resolve_mode`);
        values later rounds need from earlier solves are captured in
        ``self.round_state`` right after each solve.
        """
        mode, sequence = self._resolve_mode()
        self.round_state = _initial_round_state()
        self._solver_method, self._solver_norel_budget, self._solver_milp_time = (
            self._log_solver_config_once(len(sequence))
        )
        logs.info("Semester scheduling mode: %s", mode)

        for round_num, round_label in enumerate(sequence, start=1):
            self._begin_round(
                round_label, round_num=round_num, round_total=len(sequence)
            )
            _, build_method = _ROUND_SPECS[round_label]
            t_build = time.time()
            getattr(self, build_method)()
            self.optimize_model(round_label, build_secs=time.time() - t_build)
            if round_label == "Round1":
                self.round_state["round1_weighted_theta"] = (
                    self._eval_weighted_theta()
                )
            self._finalize_round(round_label)

        logs.info("Scheduling complete, clear skies!")

    def _finalize_round(self, round_label):
        """Build schedule, log report, persist per-night handoff + snapshot."""
        self.build_schedule()
        if round_label == "Round2":
            self.round_state["round2_slots_by_program"] = (
                self._program_slot_value()
            )
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
            self.requests_active[["unique_id", "target"]],
            on="unique_id",
            how="left",
        )
        sparse["target"] = sparse["target"].fillna("NO MATCHING NAME")
        self._ensure_output_dir()
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
        t_visit_slots = self.requests_active.set_index("unique_id")["t_visit_slots"]
        slots_per_visit_future = sched_future["unique_id"].map(t_visit_slots).fillna(1)
        slots_per_visit_today = sched_today["unique_id"].map(t_visit_slots).fillna(1)
        future_reserved = int(slots_per_visit_future.sum())
        today_reserved = int(slots_per_visit_today.sum())

        summary = pd.Series(
            {
                "Total requests": len(self.requests),
                "Total requests (active)": len(self.requests_active),
                "Total allocated slots": allocated,
                "Total slots requested": slot_demand_slots(self.requests),
                "Total slots requested (active)": slot_demand_slots(
                    self.requests_active
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
        awarded = self.programs["awarded_hours"]
        awarded_slots = self.programs["awarded_slots"]
        past_slots_by_prog = self.programs["past_slots"]
        past_by_prog = past_slots_by_prog.astype("float64") / slots_per_hour

        rf = self.requests_active.copy()
        rf["requested_h"] = (
            rf["t_visit_slots"] * rf["n_intra_max"] * rf["n_inter_max"]
        ) / slots_per_hour
        requested_by_prog = rf.groupby("program_code")["requested_h"].sum()

        sched_with_prog = sched.merge(
            self.requests_active[["unique_id", "program_code", "t_visit_slots"]],
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

        # Grace / hold values reflect what is actually enforced on the model
        # this round (captured in round_state as rounds solve).
        throttle_grace = self.round_state["throttle_grace"]
        round2_slots = self.round_state["round2_slots_by_program"]
        hold_alpha = self.round_state["hold_fill_alpha"]
        hold_active = round_label in ("Round3", "Round4", "UpcomingNight")

        miff_pct = []
        maff_pct = []
        for prog, row in table.iterrows():
            aw_slots = float(awarded_slots.get(prog, 0.0))
            past_slots = float(past_slots_by_prog.get(prog, 0))

            if aw_slots > 0:
                grace_slots = int(aw_slots * throttle_grace)
                if grace_slots < past_slots:
                    grace_slots = int(past_slots)
                maff_pct.append(100.0 * grace_slots / aw_slots)
            else:
                maff_pct.append(0.0)

            if hold_active and round2_slots is not None:
                r2_slots = float(round2_slots.get(prog, 0.0))
                min_slots = past_slots + r2_slots * (1.0 - hold_alpha)
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
        round1_theta = self.round_state["round1_weighted_theta"]
        if round_label == "UpcomingNight" and round1_theta is not None:
            slack = self.config.getfloat(
                "semester", "global_shortfall_slack",
                fallback=_DEFAULT_GLOBAL_SHORTFALL_SLACK,
            )
            logs.info(
                "UpcomingNight: Round-1 weighted shortfall=%.3f cap=%.3f "
                "post-round weighted shortfall=%.3f",
                round1_theta,
                round1_theta * slack,
                self._eval_weighted_theta(),
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
        selected_df = self.requests_active[
            self.requests_active["unique_id"].isin(selected)
        ].copy()
        selected_df["nplan_weight"] = 1.0
        self._ensure_output_dir()
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
        self._ensure_output_dir()
        tmp_path = hdf5_path + ".tmp"
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

        self.requests.to_hdf(
            tmp_path, key="requests", mode="a", format="table"
        )
        past_fmt = "fixed" if self.past.empty else "table"
        self.past.to_hdf(tmp_path, key="past", mode="a", format=past_fmt)
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
        consumers (plot.py, nplan.py) read requests, schedule,
        access_record, past, and queue, all of which are restored.
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

        requests = pd.read_hdf(hdf5_path, key="requests")
        try:
            past = pd.read_hdf(hdf5_path, key="past")
        except KeyError:
            past = pd.DataFrame(columns=PAST_COLS)
        try:
            schedule = pd.read_hdf(hdf5_path, key="schedule")
        except KeyError:
            schedule = None

        instance = cls.__new__(cls)
        instance._config_ini_text = config_ini_text
        instance.config = ConfigParser()
        instance.config.read_string(config_ini_text)
        instance.queue = astroq.queue.from_config(instance.config)

        instance._attach_slot_columns(requests)
        instance.requests = requests

        instance.past = past
        instance._attach_past_columns()
        instance.programs = instance._load_frame("programs")
        instance._attach_program_columns()
        # Downstream consumers only use the rehydrated access_obj for
        # coordinate-based queries (accessible_at, slotmidpoints); the
        # allocation/custom cubes live in the persisted access_record.
        instance.allocation = None
        instance.custom = None
        instance.access_obj = ac.Access.from_planner(instance)
        instance.access_record = access_record
        instance.schedule = schedule
        instance.round_state = _initial_round_state()

        logs.info(f"SemesterPlanner loaded from HDF5: {hdf5_path}")
        return instance