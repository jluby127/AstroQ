"""
Module that defines the SemesterPlanner class. This class is responsible for defining,
building, and solving the Gurobi model for semester-level observation planning. It is
nearly completely agnostic to all astronomy knowledge.
"""

import logging
import os
from configparser import ConfigParser
from functools import cached_property
from pathlib import Path

import gurobipy as gp
import h5py
import numpy as np
import pandas as pd
from astropy.time import Time
from gurobipy import GRB
import astroq.access as ac
import astroq.io
import astroq.queue

logs = logging.getLogger(__name__)

# Schema for h5 serialization bump when the on-disk layout changes
SEMESTER_PLANNER_H5_SCHEMA = 6

# Request columns denormalized onto the observability grid to form
# ``request_slots`` -- the single relational table the model is defined over.
STRATEGY_COLS = [
    "r",
    "target",
    "program_code",
    "n_intra_min",
    "n_intra_max",
    "n_inter_max",
    "tau_inter",
    "t_visit_slots",
    "tau_intra_slots",
]

REQUEST_SLOT_COLS = ("t_visit_slots", "tau_intra_slots")
REQUEST_PAST_COLS = (
    "past_nights_observed",
    "past_n_exposures",
    "past_date_last_observed",
    "desired_max_obs",
)
REQUEST_DERIVED_COLS = REQUEST_SLOT_COLS + REQUEST_PAST_COLS
PROGRAM_DERIVED_COLS = ("awarded_hours", "awarded_slots", "past_slots")

PROGRAM_LEDGER_COLS = (
    "past_hours",
    "sched_hours",
    "proj_hours",
    "requested_hours",
    "fill_proj",
    "fill_min",
    "fill_max",
)

TIMELINE_VALUE_COLS = (
    "unique_id",
    "program_code",
    "t_visit_slots",
)



_DEFAULT_GLOBAL_SHORTFALL_SLACK = 1.1

# Pipelines keyed by ``[semester] mode`` (mode is required in config; no default).
_MODE_PIPELINES = {
    "shortfall": "run_model_shortfall",
    "shortfall,balance,prioritize,fill-empty,fill-current-day": (
        "run_model_shortfall_balance_prioritize_fillempty_fillcurrentday"
    ),
}
_ALLOWED_MODES = list(_MODE_PIPELINES)

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
        - ``programs.csv`` -- awarded hours per program (drives throttling).

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

    def __init__(self, cf):
        """See class docstring."""

        # Read config as text so we can persist it verbatim and recreate the
        # parser on from_hdf5.
        self._config_ini_text = Path(cf).read_text()
        self.config = ConfigParser()
        self.config.optionxform = str
        self.config.read_string(self._config_ini_text)
        self.queue = astroq.queue.from_config(self.config)
        self.schedule = None

        # Load input data
        self.requests = self._load_frame("request")
        self.past = self._load_past()
        self.allocation = self._load_frame("allocation")
        self.custom = self._load_frame("custom")
        self.programs = self._load_frame("programs")

        # Add additional columns needed for model on the fly (not saved)
        self._add_request_columns()
        self._add_program_columns()
        self._validate_program_coverage()

        # Observability cube (single source of truth for which slots are valid).
        self.access_obj = ac.Access.from_planner(self)
        self.access_record = self.access_obj.build_access()

        # Build the request_slots table (sparse r,d,s table) 
        request_slots = (
            self.access_obj.observability(self.access_record.is_observable)
            .rename(columns={"unique_id": "r"})
            .merge(self.requests_active[STRATEGY_COLS], on="r")
        )
        request_slots["rds"] = request_slots[["r", "d", "s"]].apply(tuple, axis=1)
        self.request_slots = request_slots
        self.build_model()

    def _load_frame(self, kind):
        """Load a validated CSV frame. ``kind`` maps to ``{kind}_file`` in config."""
        key = f"{kind}_file"
        raw = self.config.get("data", key)
        path = raw if os.path.isabs(raw) else os.path.join(
            self.config.get("global", "workdir"), raw
        )
        return astroq.io.read_csv(path, kind)

    def _load_past(self):
        """Load past.csv and validate timestamps against the semester window."""
        past = self._load_frame("past")
        astroq.io.validate_past_in_semester(
            past,
            self.config.get("global", "semester_start_day"),
            self.config.get("global", "semester_end_day"),
        )
        return past

    # ------------------------------------------------------------------
    # Properties (paths derived from config).
    # ------------------------------------------------------------------

    @property
    def output_directory(self):
        return os.path.join(self.config.get("global", "workdir"), "outputs")

    @property
    def requests_active(self):
        return self.requests[~self.requests["inactive"]]

    @cached_property
    def semester_length(self):
        """Inclusive semester span in nights (computed once)."""
        start = self.config.get("global", "semester_start_day")
        end = self.config.get("global", "semester_end_day")
        start = Time(start, format="iso", scale="utc")
        end = Time(end, format="iso", scale="utc")
        return int(round(end.jd - start.jd)) + 1

    @cached_property
    def slot_size(self):
        """Minutes per scheduling slot from ``[semester] slot_size``."""
        return self.config.getfloat("semester", "slot_size")

    @property
    def slots_per_hour(self):
        """Slots per awarded/executed hour (inverse of slot duration)."""
        return 60 / self.slot_size

    # ------------------------------------------------------------------
    # Construction helpers.
    # ------------------------------------------------------------------

    def _add_request_columns(self):
        """Add REQUEST_DERIVED_COLS to ``self.requests``.

        Slot cols from :meth:`Queue.visit_seconds`; past cols from
        ``self.past``. Mutates in place (idempotent).
        """
        rf = self.requests
        rf["r"] = rf["unique_id"]
        self.past["r"] = self.past["unique_id"]
        visit_s = self.queue.visit_seconds(rf["exptime"], rf["n_exp"])
        rf["t_visit_slots"] = (
            (visit_s / (self.slot_size * 60.0)).round().clip(lower=1).astype(int)
        )
        rf["tau_intra_slots"] = (
            (rf["tau_intra"] * self.slots_per_hour).round().astype(int)
        )

        rs = rf["r"]
        night = self.past["timestamp"].str[:10]
        g = self.past.assign(_night=night).groupby("r")
        agg = pd.DataFrame({
            "nights": g["_night"].nunique(),
            "n_exp": g.size(),
            "last": g["_night"].max(),
        }).reindex(rs).fillna({"nights": 0, "n_exp": 0, "last": ""})

        rf["past_nights_observed"] = agg["nights"].astype(int).to_numpy()
        rf["past_n_exposures"] = agg["n_exp"].astype(int).to_numpy()
        rf["past_date_last_observed"] = agg["last"].astype(str).to_numpy()

        n_max = rf["n_inter_max"].to_numpy()
        past = rf["past_nights_observed"].to_numpy()
        over = past > n_max
        rf["desired_max_obs"] = np.where(over, past, n_max - past).astype(int)

    def _add_program_columns(self):
        """Add PROGRAM_DERIVED_COLS to ``self.programs``.

        ``past_slots`` counts ALL request rows (active + inactive), matching
        throttle semantics. Requires ``t_visit_slots`` on ``self.requests``.
        """
        if "min_fillfactor" not in self.programs.columns:
            self.programs["min_fillfactor"] = 0.0
        if "max_fillfactor" not in self.programs.columns:
            self.programs["max_fillfactor"] = 1.25

        self.programs["awarded_hours"] = self.programs["hours"]
        self.programs["awarded_slots"] = (
            self.programs["hours"] * self.slots_per_hour
        )

        req_cols = ["r", "program_code", "t_visit_slots"]
        past_by_prog = (
            self.past.merge(self.requests[req_cols], on="r", how="inner")
            .groupby("program_code")["t_visit_slots"]
            .sum()
        )
        self.programs["past_slots"] = (
            past_by_prog.reindex(self.programs.index).fillna(0).astype(int)
        )

    def _timeline_past(self):
        """One row per past.csv visit; indexed by night ``d``."""
        req_cols = ["r", "program_code", "t_visit_slots"]
        past = self.past.merge(self.requests[req_cols], on="r", how="inner")
        past["unique_id"] = past["unique_id"].astype(str)
        date_to_d = {
            d: i for i, d in enumerate(self.access_obj.all_dates_array)
        }
        past["d"] = past["timestamp"].astype(str).str[:10].map(date_to_d)
        return past.set_index("d")[list(TIMELINE_VALUE_COLS)]

    def _timeline_future(self):
        """One row per scheduled visit; indexed by night ``d``."""
        req_cols = ["r", "program_code", "t_visit_slots"]
        future = self.schedule.merge(
            self.requests[req_cols],
            left_on="unique_id",
            right_on="r",
            how="left",
        )
        future["unique_id"] = self.schedule["unique_id"].astype(str).values
        future["d"] = future["d"].astype(int)
        future["t_visit_slots"] = future["t_visit_slots"].fillna(1).astype(int)
        return future.set_index("d")[list(TIMELINE_VALUE_COLS)]

    def _invalidate_timeline(self):
        self.__dict__.pop("timeline", None)

    @cached_property
    def timeline(self):
        """Executed past visits plus scheduled future visits (one row per visit)."""
        if self.schedule is None:
            raise RuntimeError("call build_schedule() before accessing timeline")
        return pd.concat([self._timeline_past(), self._timeline_future()])

    @property
    def programs_ledger(self):
        """``programs`` enriched with timeline charged hours and fill factors.
        Requires :meth:`build_schedule` so ``timeline`` is defined. Does not
        mutate ``self.programs``.
        """
        if self.schedule is None:
            raise RuntimeError("No schedule created")
        sph = self.slots_per_hour
        idx = self.programs.index
        ps = self.timeline
        today = self.access_obj.current_night_index
        def hours(mask):
            return (
                ps.loc[mask, "t_visit_slots"]
                .groupby(ps.loc[mask, "program_code"])
                .sum()
                .reindex(idx, fill_value=0)
                / sph
            )
        rf = self.requests_active
        req_h = (
            rf["t_visit_slots"] * rf["n_intra_max"] * rf["n_inter_max"] / sph
        ).groupby(rf["program_code"]).sum().reindex(idx, fill_value=0)
        out = self.programs.assign(
            past_hours=hours(ps.index < today),
            sched_hours=hours(ps.index >= today),
            requested_hours=req_h,
        )
        out["proj_hours"] = out["past_hours"] + out["sched_hours"]
        for col, attr in (
            ("fill_proj", "X"),
            ("fill_min", "LB"),
            ("fill_max", "UB"),
        ):
            out[col] = idx.map(pd.Series(self.model.getAttr(attr, self.F)))
        return out

    def _validate_program_coverage(self):
        """Every request program_code must appear in programs.csv."""
        codes = set(self.requests["program_code"].dropna().astype(str))
        missing = sorted(codes - set(self.programs.index.astype(str)))
        if missing:
            raise ValueError(
                f"program_code(s) missing from programs.csv: {missing}"
            )

    def build_model(self):
        """Gurobi variables plus every structural constraint, inline.

        Assumes ``self.request_slots`` (built in ``__init__``) -- one row per
        observable ``(r, d, s)`` with an ``rds`` column holding each row's key
        into ``self.Yrds``. Everything here is required by every scheduling
        mode; pipeline steps layer objectives and step-specific constraints on top.

        Paper map (Lubin et al. 2025 -> code): ``Y_{r,d,s}`` -> ``Yrds``,
        ``W_{r,d}`` -> ``Wrd``, shortfall -> ``theta``; constraint numbers
        below refer to that paper. Set-building is relational (merges /
        groupbys over ``request_slots``); Python loops only emit ``addConstr``.
        """
        logs.debug("Building the SemesterPlanner.")
        logs.debug("Initializing complete.")
        self.model = gp.Model("splan")

        rs = self.request_slots

        # ---- variables ----
        self.Yrds = self.model.addVars(rs["rds"], vtype=GRB.BINARY, name="Yrds")
        self.theta = self.model.addVars(self.requests_active["r"], name="Theta")
        wrd_keys = rs.loc[rs.n_intra_max > 1].groupby(["r", "d"], sort=False).groups
        if wrd_keys:
            self.Wrd = self.model.addVars(wrd_keys.keys(), vtype=GRB.BINARY, name="Wrd")

        # ---- Eq. 3: shortfall definition. theta_r >= remaining nights owed
        # minus scheduled nights (visits / n_intra_max); lb=0 from addVars. ----
        logs.info("Constraint: Build theta variable")
        rf_indexed = self.requests_active.set_index("r")
        for r, grp_keys in rs.groupby("r", sort=False)["rds"]:
            row = rf_indexed.loc[r]
            self.model.addConstr(
                self.theta[r]
                >= row["n_inter_max"]
                - row["past_nights_observed"]
                - gp.quicksum(self.Yrds[k] for k in grp_keys) / row["n_intra_max"],
                f"greater_than_nobs_shortfall_{r}",
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
            rds_off = (
                rs_multislot["rds"].loc[[ds_on]] if ds_on in rs_multislot.index else []
            )
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
            .groupby("r", sort=False)["rds"]
            .agg(list)
        )
        for r, rds_keys in r_maxnights_single.items():
            self.model.addConstr(
                gp.quicksum(self.Yrds[rds] for rds in rds_keys) <= desired_max_obs.loc[r],
                "max_desired_unique_nights_for_request_{0}".format(r),
            )
        r_maxnights_multi = (
            rs.query("n_intra_max > 1")
            .drop_duplicates(["r", "d"])
            .groupby("r", sort=False)["d"]
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
        # (r, d) -> the rds that are forbidden when the request is observed on
        # night d.
        logs.info("Constraint: Enforce inter-night cadence.")
        rd_internight = (
            rs.query("tau_inter > 1")[["r", "d", "tau_inter"]]
            .drop_duplicates(["r", "d"])
            .merge(rs[["r", "d", "rds"]], suffixes=["", "_future"], on="r")
            .query("d < d_future < d + tau_inter")
            .groupby(["r", "d"], sort=False)["rds"]
            .agg(list)
        )
        for rd_on, grp in rs.groupby(["r", "d"], sort=False):
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
                rs.query("n_intra_max > 1")[["r", "d", "s", "tau_intra_slots"]],
                rs.query("n_intra_max > 1")[["r", "d", "s", "rds"]],
                suffixes=["", "_future"],   # left s -> s, right s -> s_future
                on=["r", "d"],
            )
            .query("s < s_future < s + tau_intra_slots")
            .groupby(["r", "d", "s"], sort=False)["rds"]
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
            .groupby(["r", "d"], sort=False)
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
            .groupby(["r", "d"], sort=False)["rds"]
            .agg(list)
        )
        for rd_on, rds_on in rd_single.items():
            self.model.addConstr(
                gp.quicksum(self.Yrds[rds] for rds in rds_on) <= 1,
                "enforce_max_visits_{0}_{1}d".format(*rd_on),
            )

        # ---- Per-program fill factor F[p] and linking constraints ----
        logs.info("Constraint: Program fill factors (F[p]).")
        self.sched_slots_by_program = {
            p: gp.quicksum(
                self.Yrds[k] * n for k, n in zip(g["rds"], g["t_visit_slots"])
            )
            for p, g in rs.groupby("program_code")
        }
        self.program_keys = self.programs.index[
            self.programs["awarded_slots"] > 0
        ]
        self.F = self.model.addVars(self.program_keys, name="F")
        awarded = self.programs["awarded_slots"]
        past = self.programs["past_slots"]
        for p in self.program_keys:
            self.model.addConstr(
                self.F[p] * float(awarded[p])
                == float(past[p]) + self.sched_slots_by_program.get(p, 0),
                f"f_link_{p}",
            )
        self._constraint_fillfactor()

    def _constraint_fillfactor(self, *, min_fillfactor=None, max_fillfactor=None):
        """Set F[p] lower/upper bounds from self.programs, with optional overrides.

        F[p] tracks projected fill factor (past + scheduled) / awarded. Because
        sched_slots >= 0, F[p] is structurally floored at past_ff = past/awarded;
        a ceiling below that floor would be infeasible, so we clamp UB up to
        past_ff (which pins F=past_ff and forces sched=0 -- the graceful
        "no new scheduling for over-budget programs" behavior of the old throttle).

        Args:
            min_fillfactor: optional floor override, scalar (all programs) or
                pd.Series keyed by program. Applied as lb = max(csv_min, override).
            max_fillfactor: optional ceiling override, scalar or pd.Series.
                Replaces csv_max for this call (ub = override).

        Bounds not overridden reset to their csv values on every call, so a
        step must re-pass any floor/ceiling it wants to preserve.
        """
        def _at(val, p):
            if isinstance(val, pd.Series):
                return float(val[p])
            return float(val)

        awarded = self.programs["awarded_slots"]
        past = self.programs["past_slots"]

        for p in self.F:
            lb = float(self.programs.at[p, "min_fillfactor"])
            ub = float(self.programs.at[p, "max_fillfactor"])

            if min_fillfactor is not None:
                lb = max(lb, _at(min_fillfactor, p))
            if max_fillfactor is not None:
                ub = _at(max_fillfactor, p)

            past_ff = float(past[p]) / float(awarded[p])
            if past_ff > ub:
                logs.warning(
                    "Program %s over ceiling from past alone "
                    "(past_ff=%.3f > max_fillfactor=%.3f); pinning F, sched=0.",
                    p,
                    past_ff,
                    ub,
                )
                ub = past_ff

            lb = min(lb, ub)
            self.F[p].LB = lb
            self.F[p].UB = ub

        if getattr(self, "model", None) is not None:
            self.model.update()

    # ==================================================================
    # Objectives.
    # ==================================================================

    def _objective_weighted_theta(self):
        """Time-weighted global shortfall (Lubin et al. Eq. 1)."""
        t_visit = self.requests_active.set_index("r")["t_visit_slots"]
        return gp.quicksum(
            self.theta[r] * t_visit[r] for r in self.theta
        )

    def _objective_slots_used_tonight(self):
        """Filled slot-time on ``current_day``."""
        d_today = self.access_obj.current_night_index
        t_visit = self.requests_active.set_index("r")["t_visit_slots"]
        return gp.quicksum(
            t_visit[r] * self.Yrds[r, d, s]
            for r, d, s in self.request_slots["rds"]
            if d == d_today
        )

    def _objective_prioritize_intra(self):
        """Weight scheduled slot-time by inverse ``splan_weight`` per target."""
        logs.info("Objective: Intra-program priorities (inverse splan_weight).")
        splan_weight = self.requests_active.set_index("r")["splan_weight"]
        return gp.quicksum(
            (1.0 / splan_weight.loc[r])
            * self.Yrds[k]
            * n
            for k, n, r in zip(
                self.request_slots["rds"],
                self.request_slots["t_visit_slots"],
                self.request_slots["r"],
            )
        )

    def _objective_minimize_empty_slots(self):
        """Minimize unscheduled allocated slot-time."""
        logs.info("Objective: Minimize empty allocated slots.")
        total_allocated = int(self.access_record["is_allocated"][0].sum())
        scheduled = gp.quicksum(
            self.Yrds[k] * n
            for k, n in zip(
                self.request_slots["rds"],
                self.request_slots["t_visit_slots"],
            )
        )
        return total_allocated - scheduled

    # ==================================================================
    # Model orchestration.
    # ==================================================================

    def optimize_model(self, step):
        """Apply per-step Gurobi params from config and solve."""
        params = self.model.Params
        for section in ("semester.default.gurobi", f"semester.{step}.gurobi"):
            if not self.config.has_section(section):
                continue
            for key in self.config.options(section):
                if not hasattr(params, key):
                    logs.warning("Ignoring unknown Gurobi param %s", key)
                    continue
                template = getattr(params, key)
                if isinstance(template, bool):
                    val = self.config.getboolean(section, key)
                else:
                    raw = self.config.get(section, key)
                    val = (
                        int(raw)
                        if isinstance(template, int) and "." not in raw
                        else float(raw)
                    )
                setattr(params, key, val)
        self.model.update()
        self.model.optimize()

        if self.model.Status == GRB.INFEASIBLE:
            raise RuntimeError(
                f"{step} solve infeasible (status={self.model.Status})"
            )

    def run_model(self):
        """Dispatch to the semester scheduling pipeline named in ``[semester] mode``."""
        mode = self.config.get("semester", "mode").strip().lower()
        if mode not in _MODE_PIPELINES:
            raise ValueError(
                f"[semester] mode={mode!r} invalid; expected one of {_ALLOWED_MODES!r}"
            )
        getattr(self, _MODE_PIPELINES[mode])()
        logs.info("Scheduling complete, clear skies!")

    def run_model_shortfall(self):
        """Shortfall-only pipeline: minimize weighted theta, write outputs."""
        self._constraint_fillfactor(max_fillfactor=1.0)
        self.model.setObjective(self._objective_weighted_theta(), GRB.MINIMIZE)
        self.optimize_model("shortfall")
        self.build_schedule()
        self.log_report("shortfall")
        self.write_request_selected()
        self.to_hdf5()

    def run_model_shortfall_balance_prioritize_fillempty_fillcurrentday(self):
        """Full pipeline: shortfall through fill-current-day."""
        # ===== shortfall =====
        self._constraint_fillfactor(max_fillfactor=1.0)
        self.model.setObjective(self._objective_weighted_theta(), GRB.MINIMIZE)
        self.optimize_model("shortfall")
        objective_shortfall_min = self.model.ObjVal
        self.build_schedule()
        self.log_report("shortfall")

        # ===== balance =====
        self._constraint_fillfactor(
            min_fillfactor=pd.Series(self.model.getAttr("X", self.F))
        )
        self.model.setObjective(self.F.sum(), GRB.MAXIMIZE)
        self.optimize_model("balance")
        self.build_schedule()
        self.log_report("balance")

        # ===== prioritize =====
        hold_alpha = self.config.getfloat(
            "semester.prioritize",
            "hold_fill_alpha",
            fallback=0.0,
        )
        hold_scale = 1.0 - hold_alpha
        self._constraint_fillfactor(
            min_fillfactor=pd.Series(self.model.getAttr("X", self.F)) * hold_scale,
        )
        self.model.setObjective(self._objective_prioritize_intra(), GRB.MAXIMIZE)
        self.optimize_model("prioritize")
        self.build_schedule()
        self.log_report("prioritize", hold_fill_alpha=hold_alpha)

        # ===== fill-empty =====
        logs.info("Constraint: Per-target shortfall frozen at prior values.")
        for r in self.theta:
            self.model.addConstr(
                self.theta[r] <= float(self.theta[r].X),
                f"theta_le_prior_{r}",
            )
        self.model.setObjective(self._objective_minimize_empty_slots(), GRB.MINIMIZE)
        self.optimize_model("fill-empty")
        self.build_schedule()
        self.log_report("fill-empty", hold_fill_alpha=hold_alpha)

        # ===== fill-current-day =====
        slack = self.config.getfloat(
            "semester.fill-current-day",
            "global_shortfall_slack",
            fallback=_DEFAULT_GLOBAL_SHORTFALL_SLACK,
        )
        cap = objective_shortfall_min * slack
        logs.info(
            "Constraint: weighted shortfall <= shortfall optimum * %g "
            "(cap=%.3f from objective_shortfall_min=%.3f)",
            slack,
            cap,
            objective_shortfall_min,
        )
        self.model.addConstr(
            self._objective_weighted_theta() <= cap,
            "fix_global_shortfall_upcoming_night",
        )
        d_today = self.access_obj.current_night_index
        current_day = self.config.get("global", "current_day")
        logs.info(
            "Objective: Maximize slot usage on upcoming night current_day=%s (d=%d).",
            current_day,
            d_today,
        )
        self.model.setObjective(self._objective_slots_used_tonight(), GRB.MAXIMIZE)
        self.optimize_model("fill-current-day")

        self.build_schedule()
        self.log_report(
            "fill-current-day",
            objective_shortfall_min=objective_shortfall_min,
            hold_fill_alpha=hold_alpha,
        )
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
        df = pd.DataFrame(self.Yrds.keys(), columns=["r", "d", "s"])
        df["value"] = [self.Yrds[k].x for k in self.Yrds.keys()]
        sparse = df.query("value > 0").drop(columns=["value"]).copy()
        sparse = sparse.merge(
            self.requests_active[["r", "unique_id", "target"]],
            on="r",
            how="left",
        )
        sparse = sparse.drop(columns=["r"])
        sparse["target"] = sparse["target"].fillna("NO MATCHING NAME")
        os.makedirs(self.output_directory, exist_ok=True)
        sparse.to_csv(
            os.path.join(self.output_directory, "semester_plan.csv"),
            index=False,
            na_rep="",
        )
        self.schedule = sparse
        self._invalidate_timeline()

    def to_string_summary(
        self,
        *,
        header="Semester Planner Statistics",
    ):
        """Global run-report summary (request/slot counts and fill factors).

        Requires that :meth:`build_schedule` has been called so
        ``self.schedule`` is set.
        """
        if self.schedule is None:
            raise RuntimeError("call build_schedule() before to_string_summary()")

        def slot_demand_slots(frame):
            return int(
                (
                    frame["t_visit_slots"]
                    * frame["n_intra_max"]
                    * frame["n_inter_max"]
                ).sum()
            )

        today_idx = self.access_obj.current_night_index
        active_with_future_slots = (
            self.request_slots.loc[self.request_slots["d"] >= today_idx, "r"]
            .unique()
        )
        n_active_future_slots = sum(
            r in active_with_future_slots for r in self.requests_active["r"]
        )
        is_alloc_2d = self.access_record["is_allocated"][0]
        allocated = int(is_alloc_2d.sum())
        allocated_future = int(is_alloc_2d[today_idx:].sum())
        allocated_today = int(is_alloc_2d[today_idx].sum())

        sched = self.schedule
        sched_future = sched[sched["d"] >= today_idx]
        sched_today = sched[sched["d"] == today_idx]
        t_visit_slots = self.requests_active.set_index("r")["t_visit_slots"]
        slots_per_visit_future = sched_future["unique_id"].map(t_visit_slots).fillna(1)
        slots_per_visit_today = sched_today["unique_id"].map(t_visit_slots).fillna(1)
        future_reserved = int(slots_per_visit_future.sum())
        today_reserved = int(slots_per_visit_today.sum())

        summary = pd.Series(
            {
                "Total requests": len(self.requests),
                "Total requests (active)": len(self.requests_active),
                "Total requests (active, future slots > 0)": n_active_future_slots,
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

        divider = "-" * 54
        parts = [
            header,
            divider,
            summary.to_string(float_format=lambda x: f"{int(round(x))}"),
        ]
        return "\n".join(parts) + "\n"

    def to_string_programs(self):
        """Per-program hours table and fill-factor statistics.

        Requires that :meth:`build_schedule` has been called so
        ``self.schedule`` is set.
        """
        if self.schedule is None:
            raise RuntimeError("call build_schedule() before to_string_programs()")

        hour_cols = ("aw", "req", "past", "proj")
        pct_cols = ("miff%", "maff%", "past%", "proj%")
        ledger = self.programs_ledger.sort_index()
        aw = ledger["awarded_hours"]
        table = ledger.rename(
            columns={
                "awarded_hours": "aw",
                "requested_hours": "req",
                "past_hours": "past",
                "proj_hours": "proj",
            }
        )[(*hour_cols,)]
        table["proj%"] = 100.0 * ledger["fill_proj"]
        table["miff%"] = 100.0 * ledger["fill_min"]
        table["maff%"] = 100.0 * ledger["fill_max"]
        table["past%"] = np.where(aw > 0, 100.0 * ledger["past_hours"] / aw, 0.0)
        table = table[[*hour_cols, *pct_cols]]
        table[["proj%", "miff%", "maff%"]] = table[["proj%", "miff%", "maff%"]].fillna(
            0.0
        )
        for col in hour_cols:
            table[col] = table[col].map("{:.1f}".format)
        for col in pct_cols:
            table[col] = table[col].map(lambda x: f"{int(round(x))}%")

        stats_divider = "-" * 19 + " Program Statistics " + "-" * 19
        return "\n".join(
            [stats_divider, table.to_string(), "", _PROGRAM_STATS_KEY, ""]
        ) + "\n"

    def log_report(self, step, **report_ctx):
        """Emit the run-report text to stdout (no log prefix on table lines)."""
        objective_shortfall_min = report_ctx.get("objective_shortfall_min")
        if step == "fill-current-day" and objective_shortfall_min is not None:
            slack = self.config.getfloat(
                "semester.fill-current-day",
                "global_shortfall_slack",
                fallback=_DEFAULT_GLOBAL_SHORTFALL_SLACK,
            )
            logs.info(
                "fill-current-day: shortfall objective=%.3f cap=%.3f "
                "post-step shortfall objective=%.3f",
                objective_shortfall_min,
                objective_shortfall_min * slack,
                self._objective_weighted_theta().getValue(),
            )
        logs.info("Run report (%s):", step)
        print(self.to_string_summary().rstrip(), flush=True)
        print()
        print(self.to_string_programs().rstrip(), flush=True)

    def write_request_selected(self):
        """Write ``request_selected.csv`` -- the handoff to ``NightPlanner``."""
        today_idx = self.access_obj.current_night_index
        selected = {
            k[0] for k, v in self.Yrds.items() if v.x > 0 and k[1] == today_idx
        }
        selected_df = self.requests_active[
            self.requests_active["r"].isin(selected)
        ].copy()
        selected_df["nplan_weight"] = 1.0
        os.makedirs(self.output_directory, exist_ok=True)
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
        os.makedirs(self.output_directory, exist_ok=True)
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
        past = pd.read_hdf(hdf5_path, key="past")
        try:
            schedule = pd.read_hdf(hdf5_path, key="schedule")
        except KeyError:
            schedule = None

        instance = cls.__new__(cls)
        instance._config_ini_text = config_ini_text
        instance.config = ConfigParser()
        instance.config.optionxform = str
        instance.config.read_string(config_ini_text)
        instance.queue = astroq.queue.from_config(instance.config)

        instance.requests = requests
        instance.past = past
        instance.programs = instance._load_frame("programs")
        instance._add_request_columns()
        instance._add_program_columns()
        # Downstream consumers only use the rehydrated access_obj for
        # coordinate-based queries (accessible_at, slotmidpoints); the
        # allocation/custom cubes live in the persisted access_record.
        instance.allocation = None
        instance.custom = None
        instance.access_obj = ac.Access.from_planner(instance)
        instance.access_record = access_record
        instance.schedule = schedule

        logs.info(f"SemesterPlanner loaded from HDF5: {hdf5_path}")
        return instance