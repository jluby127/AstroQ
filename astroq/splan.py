"""
Module that defines the SemesterPlanner class. This class is responsible for defining,
building, and solving the Gurobi model for semester-level observation planning. It is
nearly completely agnostic to all astronomy knowledge.
"""

import logging
import os
import time
from configparser import ConfigParser
from datetime import datetime
from pathlib import Path

import gurobipy as gp
import h5py
import numpy as np
import pandas as pd
from gurobipy import GRB
import astroq.access as ac
import astroq.queue

logs = logging.getLogger(__name__)

# Schema for h5 serialization bump when the on-disk layout changes
SEMESTER_PLANNER_H5_SCHEMA = 4

# Canonical past.csv column schema. ``junk`` is optional.
PAST_COLS = ["unique_id", "target", "timestamp", "exposure_time"]


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
        """Read ``past.csv``, drop junk-flagged visits, return a DataFrame.

        Empty/missing files yield an empty frame with the canonical schema.
        Junk filter: drop a visit (``unique_id, timestamp`` group) when at
        least half of its rows are flagged ``junk=True``.
        """
        path = self._resolve_path("past_file")
        if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
            return pd.DataFrame(columns=PAST_COLS)
        try:
            df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=PAST_COLS)
        if "junk" in df.columns:
            df["junk"] = df["junk"].fillna(False).astype(bool)
            keep = df.groupby(["unique_id", "timestamp"])["junk"].transform(
                lambda s: s.sum() < len(s) / 2
            )
            df = df.loc[keep]
        return df.reset_index(drop=True)

    # ------------------------------------------------------------------
    # Properties (date-derived; path attrs are set in __init__).
    # ------------------------------------------------------------------

    @property
    def semester_length(self):
        start = datetime.strptime(
            self.config.get("global", "semester_start_day"), "%Y-%m-%d"
        )
        end = datetime.strptime(
            self.config.get("global", "semester_end_day"), "%Y-%m-%d"
        )
        return int((end - start).days + 1)

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
        """Read request.csv, clean, validate. Returns ``(all_frame, active_frame)``.

        Cleaning rules (applied once at CSV ingest; not repeated on HDF5
        rehydrate): tolerate "None" strings left over from the early-2025B
        webform, ensure ``comments`` column exists, normalize ``unique_id``
        and ``target`` to strings, and fail on duplicate active unique_id
        (which would otherwise produce a cryptic Gurobi error later).

        Also appends two derived slot-unit columns to the active frame:

        - ``t_visit_slots`` -- full per-visit duration in slots, computed
          from :meth:`astroq.queue.base.Queue.visit_seconds` (includes
          inter-shot readouts and slew overhead), rounded and clipped to
          >= 1. This is the slot reservation charged by every Gurobi
          consumer (constraint_reserve_multislot_exposures, objectives,
          throttle, Access multi-slot windowing).
        - ``tau_intra_slots`` -- minimum intra-night spacing between
          visits, in slots.

        Original units of ``exptime`` (seconds) and ``tau_intra`` (hours)
        are left untouched.
        """
        request_file = self._resolve_path("request_file")
        if not os.path.exists(request_file):
            raise FileNotFoundError(f"Requests file not found: {request_file}")

        rfa = pd.read_csv(request_file)
        if "comments" not in rfa.columns:
            rfa["comments"] = ""
        rfa["inactive"] = rfa["inactive"].fillna(False).astype(bool)
        mask = ~rfa["inactive"]
        logs.warning(
            f"There are {len(rfa[~mask])} inactive of {len(rfa)} requests."
        )
        rf = rfa[mask].reset_index(drop=True).copy()

        for col, default in (("n_intra_max", 1), ("n_intra_min", 1), ("tau_intra", 0)):
            rf[col] = rf[col].replace("None", np.nan).fillna(default)
        for band in (1, 2, 3):
            col = f"weather_band_{band}"
            if col in rf.columns:
                rf[col] = rf[col].replace("None", np.nan).fillna(False)
        rf["unique_id"] = rf["unique_id"].astype(str)
        rf["target"] = rf["target"].astype(str)

        dup_mask = rf["unique_id"].duplicated(keep=False)
        if dup_mask.any():
            dup_ids = sorted(rf.loc[dup_mask, "unique_id"].unique())
            raise ValueError(
                f"Duplicate unique_id among active requests in {request_file!r}: "
                f"{dup_ids}. Remove or merge duplicate rows so each active "
                f"request has one row."
            )

        self._attach_slot_columns(rf)
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
        slot_size = self.config.getint("semester", "slot_size")
        visit_s = self.queue.visit_seconds(
            rf["exptime"].astype(float),
            rf["n_exp"].astype(int),
            rf["n_intra_max"].astype(int),
        )
        rf["t_visit_slots"] = (
            (visit_s / (slot_size * 60.0)).round().clip(lower=1).astype(int)
        )
        rf["tau_intra_slots"] = (
            (rf["tau_intra"].astype(float) * 60 / slot_size).round().astype(int)
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
        self.build_observation_chains()
        self.build_yrds_tuples()

    def build_observation_chains(self):
        """Group single-shot targets; write ``ss_chain`` / ``ss_chain_idx``."""
        rf = self.requests_frame.copy()
        rf["ss_chain"] = pd.Series(pd.NA, index=rf.index, dtype="Int64")
        rf["ss_chain_idx"] = pd.Series(pd.NA, index=rf.index, dtype="Int64")

        if not self.config.getboolean("semester", "use_observation_chains"):
            self.requests_frame = rf
            return

        max_chain_len = self.config.getint("semester", "chain_max_length")
        min_overlap = self.config.getint("semester", "chain_min_overlap")
        max_chain_slots = int(
            self.config.getfloat("semester", "chain_max_minutes")
            / self.config.getint("semester", "slot_size")
        )

        pool = rf[
            (rf["n_inter_max"].astype(int) == 1)
            & (rf["n_intra_max"].astype(int) == 1)
        ]
        pool_uids = pool["unique_id"].tolist()
        uid_to_tvisit = pool.set_index("unique_id")["t_visit_slots"].astype(int)

        obs = self.observability[["unique_id", "d", "s"]]
        rng = np.random.default_rng(self.config.getint("semester", "random_seed"))

        chain_id = 0
        n_chained = 0
        remaining = set(pool_uids)
        uid_order = list(pool_uids)
        rng.shuffle(uid_order)

        while remaining:
            seed_uid = next(u for u in uid_order if u in remaining)
            remaining.remove(seed_uid)

            chain = [seed_uid]
            chain_slots_used = int(uid_to_tvisit.at[seed_uid])
            offset_next = chain_slots_used

            while len(chain) < max_chain_len and remaining:
                anchor_uid = chain[0]
                offset = offset_next

                anchor = obs.query("unique_id == @anchor_uid")[["d", "s"]].rename(
                    columns={"s": "s_anchor"},
                )
                cands = obs.loc[
                    obs["unique_id"].isin(remaining), ["unique_id", "d", "s"]
                ]
                aligned = anchor.merge(cands, on="d").query("s == s_anchor + @offset")
                scores = aligned.groupby("unique_id").size()
                scores = scores[
                    scores.index.map(
                        lambda uid: chain_slots_used + int(uid_to_tvisit.at[uid])
                        <= max_chain_slots
                    )
                ]
                if scores.empty or scores.max() < min_overlap:
                    break

                best_uid = scores.idxmax()
                remaining.remove(best_uid)
                chain.append(best_uid)
                chain_slots_used += int(uid_to_tvisit.at[best_uid])
                offset_next = chain_slots_used

            if len(chain) < 2:
                continue

            for chain_idx, uid in enumerate(chain):
                mask = rf["unique_id"] == uid
                rf.loc[mask, "ss_chain"] = chain_id
                rf.loc[mask, "ss_chain_idx"] = chain_idx
            chain_id += 1
            n_chained += len(chain)

        self.requests_frame = rf
        logs.info(
            "Built %d observation chains covering %d / %d pool targets",
            chain_id,
            n_chained,
            len(pool_uids),
        )

    def build_yrds_tuples(self):
        """``yrds_tuples``: full observability minus orphan slots for chained followers."""
        if not self.config.getboolean("semester", "use_observation_chains"):
            self.yrds_tuples = list(self.observability_tuples)
            self._index_yrds_tuples()
            return

        keep = set(self.observability_tuples)
        obs = self.observability
        rf = self.requests_frame
        chained = rf[rf["ss_chain"].notna()]

        if not chained.empty:
            for _chain_id, grp in chained.groupby("ss_chain"):
                members = grp.sort_values("ss_chain_idx")
                anchor = members.loc[
                    members["ss_chain_idx"] == 0, "unique_id"
                ].iloc[0]
                offsets = (
                    members.set_index("ss_chain_idx")["t_visit_slots"]
                    .astype(int)
                    .cumsum()
                    .shift(1, fill_value=0)
                )
                anchor_starts = obs.loc[
                    obs["unique_id"] == anchor, ["d", "s"]
                ].rename(columns={"s": "s_anchor"})

                for _, row in members[members["ss_chain_idx"] > 0].iterrows():
                    follower = row["unique_id"]
                    offset = int(offsets.loc[row["ss_chain_idx"]])
                    follower_at = obs.loc[
                        obs["unique_id"] == follower, ["d", "s"]
                    ]
                    on_chain = anchor_starts.merge(follower_at, on="d").query(
                        "s == s_anchor + @offset"
                    )[["d", "s"]]
                    on_chain_set = {
                        (int(d), int(s))
                        for d, s in on_chain.itertuples(index=False, name=None)
                    }

                    for d, s in obs.loc[
                        obs["unique_id"] == follower, ["d", "s"]
                    ].itertuples(index=False, name=None):
                        if (int(d), int(s)) not in on_chain_set:
                            keep.discard((follower, int(d), int(s)))

        self.yrds_tuples = list(keep)
        yrds_df = pd.DataFrame(self.yrds_tuples, columns=["unique_id", "d", "s"])
        self.all_valid_ds_for_request = yrds_df.groupby("unique_id")[["d", "s"]].agg(
            list
        )
        self._index_yrds_tuples()
        logs.info(
            "Yrds tuples: %d (observability %d)",
            len(self.yrds_tuples),
            len(self.observability_tuples),
        )

    def _index_yrds_tuples(self):
        """Build ``_yrds_keys`` and ``_yrds_by_ds`` from ``yrds_tuples``."""
        self._yrds_keys = set(self.yrds_tuples)
        by_ds = {}
        for uid, d, s in self.yrds_tuples:
            by_ds.setdefault((int(d), int(s)), set()).add(uid)
        self._yrds_by_ds = by_ds

    def _build_chain_link_rows(self):
        """Rows for batched chain equalities: anchor, follower, d, s, offset."""
        obs = self.observability
        rf = self.requests_frame
        parts = []

        for _chain_id, grp in rf[rf["ss_chain"].notna()].groupby("ss_chain"):
            members = grp.sort_values("ss_chain_idx")
            anchor = members.loc[members["ss_chain_idx"] == 0, "unique_id"].iloc[0]
            offsets = (
                members.set_index("ss_chain_idx")["t_visit_slots"]
                .astype(int)
                .cumsum()
                .shift(1, fill_value=0)
            )
            anchor_starts = obs.loc[
                obs["unique_id"] == anchor, ["d", "s"]
            ].rename(columns={"s": "s_anchor"})

            for _, row in members[members["ss_chain_idx"] > 0].iterrows():
                follower = row["unique_id"]
                offset = int(offsets.loc[row["ss_chain_idx"]])
                follower_at = obs.loc[
                    obs["unique_id"] == follower, ["d", "s"]
                ].rename(columns={"s": "s_follower"})
                links = anchor_starts.merge(follower_at, on="d").query(
                    "s_follower == s_anchor + @offset"
                )
                parts.append(
                    links.assign(anchor=anchor, follower=follower, offset=offset)
                    .rename(columns={"s_anchor": "s"})[
                        ["anchor", "follower", "d", "s", "offset"]
                    ]
                )

        if not parts:
            return pd.DataFrame(columns=["anchor", "follower", "d", "s", "offset"])
        return pd.concat(parts, ignore_index=True)

    def _log_boost_current_day_slots(self):
        """Report observable slot counts on current_day for each boosted target."""
        if self.boost is None:
            return
        current_day = self.config.get("global", "current_day")
        d_today = self.today_starting_night
        boost_by_uid = dict(
            zip(
                self.boost["unique_id"].astype(str),
                self.boost["boost"].astype(float),
            )
        )
        uid_to_target = dict(
            zip(
                self.requests_frame_all["unique_id"].astype(str),
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
        joiner_uids = self.joiner["unique_id"].astype(str)
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
            night = self.past_df["timestamp"].astype(str).str[:10]
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
        n_max = rf["n_inter_max"].astype(int).to_numpy()
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

    def constraint_link_observation_chains(self):
        """Tie follower Yrds to anchor Yrds via one batched addConstrs call."""
        if not self.config.getboolean("semester", "use_observation_chains"):
            return
        if self.requests_frame["ss_chain"].notna().sum() == 0:
            return

        chain_links = self._build_chain_link_rows()
        if chain_links.empty:
            return

        logs.info(
            "Constraint: Link observation chains (%d equalities).",
            len(chain_links),
        )
        self.model.addConstrs(
            (
                self.Yrds[r, d, s + o] == self.Yrds[a, d, s]
                for a, r, d, s, o in chain_links.itertuples(index=False)
            ),
            name="chain",
        )

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

    @property
    def multi_visit_uids(self):
        """uids that may receive >1 visit per night (Wrd is defined for these)."""
        return set(
            self.joiner.loc[self.joiner["n_intra_max"] > 1, "unique_id"].unique()
        )

    # ---- throttling & bonus round ----

    def constraint_throttle(self):
        """
        Not described in Lubin et al. 2025.

        Ensure that no program is scheduled for more time than they bring to
        the queue (within a grace amount).
        """
        logs.info("Constraint: Throttling over-requested programs.")
        program_frame = pd.read_csv(self.programs_file).set_index("program")
        slot_size = self.config.getint("semester", "slot_size")
        hours_per_night = self.config.getfloat("semester", "hours_per_night")
        throttle_grace = self.config.getfloat("semester", "throttle_grace")

        program_frame["awarded_slots"] = (
            program_frame["nights"] * hours_per_night * 60 / slot_size
        )
        program_frame["awarded_slots_grace"] = (
            program_frame["awarded_slots"] * throttle_grace
        ).astype(int)

        rf = self.requests_frame
        t_visit_slots = dict(zip(rf["unique_id"], rf["t_visit_slots"]))
        merged_df = rf[
            ["unique_id", "program_code", "past_n_exposures", "t_visit_slots"]
        ].copy()
        merged_df["past_slots_used"] = (
            merged_df["past_n_exposures"].astype(int)
            * merged_df["t_visit_slots"].astype(int)
        )
        merged_df = merged_df.drop(columns=["past_n_exposures", "t_visit_slots"])
        merged_df = merged_df.merge(
            program_frame[["nights"]],
            left_on="program_code",
            right_index=True,
            how="inner",
        )

        past_used_slots_by_program = (
            merged_df.groupby("program_code")["past_slots_used"].sum().to_dict()
        )
        program_requests_map = (
            merged_df.groupby("program_code")["unique_id"].apply(set).to_dict()
        )

        for program, row in program_frame.iterrows():
            awarded_slots_grace = int(row["awarded_slots_grace"])
            uids_for_program = program_requests_map.get(program, set())
            schedulable_slots = gp.quicksum(
                self.Yrds[r, d, s] * t_visit_slots[r]
                for r, d, s in self.yrds_tuples
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
            gp.quicksum(self.theta[uid] for uid in self.requests_frame["unique_id"])
            <= self.model.objval + epsilon,
            "fix_previous_objective",
        )

    # ==================================================================
    # Objectives.
    # ==================================================================

    def set_objective_minimize_theta_time_normalized(self):
        """See Equation 1 in Lubin et al. 2025."""
        schedulable_uids = list(self.joiner["unique_id"].unique())
        t_visit_slots = dict(
            zip(self.requests_frame["unique_id"], self.requests_frame["t_visit_slots"])
        )
        theta_obj = gp.quicksum(
            self.theta[uid] * t_visit_slots[uid] for uid in schedulable_uids
        )
        if self.boost is not None:
            boost_by_uid = dict(
                zip(
                    self.boost["unique_id"].astype(str),
                    self.boost["boost"].astype(float),
                )
            )
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
        t_visit_slots = dict(
            zip(self.requests_frame["unique_id"], self.requests_frame["t_visit_slots"])
        )
        self.model.setObjective(
            gp.quicksum(
                t_visit_slots[uid] * self.Yrds[uid, d, s]
                for uid, d, s in self.yrds_tuples
            ),
            GRB.MAXIMIZE,
        )

    # ==================================================================
    # Model orchestration.
    # ==================================================================

    def build_model_round1(self):
        """Round 1 constraints + objective per Lubin et al. 2025."""
        t1 = time.time()
        self.constraint_link_observation_chains()
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
        self.model.params.MIPFocus = 1 # 0 means balance objective and feasibility, 1 means feasibility, 2 means optimality
        # -1 means default degeneracy moves. helps find feasible solutions in highly degenerate cases.
        self.model.params.DegenMoves = -1
        norel_heur_time = self.config.getfloat("semester", "norel_heur_time", fallback=0.0)
        if norel_heur_time > 0:
            self.model.params.NoRelHeurTime = norel_heur_time
            logs.info("Gurobi NoRelHeurTime = %g s", norel_heur_time)
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
        logs.info(f"Time to finish solver: {time.time() - t1:.3f}")

    def run_model(self):
        """Construct and solve the Gurobi model (with optional bonus round)."""
        self.build_model_round1()
        self.optimize_model()
        self._finalize_round("Round1")
        if self.config.getboolean("semester", "run_bonus_round"):
            self.build_model_round2()
            self.optimize_model()
            self._finalize_round("Round2")
        logs.info("Scheduling complete, clear skies!")

    def _finalize_round(self, round_label):
        """Build schedule, log report, persist per-night handoff + snapshot."""
        self.build_schedule()
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

    def to_string(self, *, header="Stats for Round1"):
        """Run report: top-level summary Series + per-program hours DataFrame.

        Requires that :meth:`build_schedule` has been called so
        ``self.schedule`` is set.
        """
        if self.schedule is None:
            raise RuntimeError("call build_schedule() before to_string()")

        slot_size = self.config.getint("semester", "slot_size")
        hours_per_night = self.config.getfloat("semester", "hours_per_night")
        slots_per_hour = 60 / slot_size

        # ---- top-level summary as a Series ----
        is_alloc_2d = self.access_record["is_allocated"][0]
        sched = self.schedule
        t_visit_slots = self.requests_frame.set_index("unique_id")["t_visit_slots"]
        slots_per_visit = sched["unique_id"].map(t_visit_slots).fillna(1)
        scheduled_starting = len(sched)
        reserved = int((slots_per_visit - 1).clip(lower=0).sum())
        total_scheduled = scheduled_starting + reserved
        allocated = int(is_alloc_2d.sum())
        rf_slots = self.requests_frame["t_visit_slots"]
        total_requested = int(
            (
                rf_slots
                * self.requests_frame["n_intra_max"]
                * self.requests_frame["n_inter_max"]
            ).sum()
        )

        summary = pd.Series(
            {
                "N slots in semester": is_alloc_2d.size,
                "N available slots": allocated,
                "N starting slots scheduled": scheduled_starting,
                "N reserved slots": reserved,
                "N total slots scheduled": total_scheduled,
                "N slots left empty": allocated - total_scheduled,
                "N slots requested (total)": total_requested,
                "Utilization (% of available slots)": (
                    100 * total_scheduled / allocated if allocated else 0.0
                ),
                "Utilization (% of requested slots)": (
                    100 * total_scheduled / total_requested
                    if total_requested
                    else 0.0
                ),
            }
        )

        # ---- per-program table (hours only) ----
        progs = (
            pd.read_csv(self.programs_file)
            .rename(columns={"program": "program_code", "nights": "awarded_nights"})
            .set_index("program_code")
        )
        awarded = progs["awarded_nights"] * hours_per_night

        rf = self.requests_frame.copy()
        rf["requested_h"] = (
            rf["t_visit_slots"] * rf["n_intra_max"] * rf["n_inter_max"]
        ) / slots_per_hour
        rf["past_h"] = (rf["past_n_exposures"] * rf["t_visit_slots"]) / slots_per_hour
        by_prog = rf.groupby("program_code")[["requested_h", "past_h"]].sum()

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
            pd.DataFrame({"awarded": awarded})
            .join(by_prog, how="left")
            .join(scheduled_h, how="left")
            .fillna(0.0)
            .rename(
                columns={
                    "requested_h": "requested",
                    "past_h": "past",
                    "scheduled_h": "scheduled",
                }
            )
        )
        done = table["past"] + table["scheduled"]
        table["req/aw%"] = np.where(
            table["awarded"] > 0, 100 * table["requested"] / table["awarded"], 0.0
        )
        table["done/req%"] = np.where(
            table["requested"] > 0, 100 * done / table["requested"], 0.0
        )
        table["done/aw%"] = np.where(
            table["awarded"] > 0, 100 * done / table["awarded"], 0.0
        )
        table = table.sort_index()

        divider = "-" * 54
        parts = [
            header,
            divider,
            summary.to_string(float_format=lambda x: f"{x:.2f}"),
            "",
            "Program Statistics (hours):",
            divider,
            table.to_string(float_format=lambda x: f"{x:.2f}"),
            "",
        ]
        return "\n".join(parts) + "\n"

    def log_report(self, round_label):
        """Log the run-report text (the same content that used to land in runReport.txt)."""
        report = self.to_string(header=f"Stats for {round_label}")
        for line in report.splitlines():
            logs.info(line)

    def write_request_selected(self):
        """Write ``request_selected.csv`` -- the handoff to ``NightPlanner``."""
        today_idx = self.all_dates_dict[self.config.get("global", "current_day")]
        selected = {
            k[0] for k, v in self.Yrds.items() if v.x > 0 and k[1] == today_idx
        }
        self.requests_frame[
            self.requests_frame["unique_id"].isin(selected)
        ].to_csv(
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

        # Cleaning was done at original CSV ingest; on rehydrate we just split
        # active vs. all, then re-derive the slot/past columns (they're pure
        # functions of the persisted data so we don't ship them on disk).
        instance.requests_frame_all = requests_frame_all
        instance.requests_frame = (
            requests_frame_all[~requests_frame_all["inactive"].astype(bool)]
            .reset_index(drop=True)
            .copy()
        )
        instance._attach_slot_columns(instance.requests_frame)

        instance.past_df = past_df
        instance._attach_past_columns()
        instance.access_obj = ac.Access.from_planner(instance)
        instance.access_record = access_record
        instance.schedule = schedule

        logs.info(f"SemesterPlanner loaded from HDF5: {hdf5_path}")
        return instance
