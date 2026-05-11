"""
In-house Gurobi ILP for night-level scheduling.

Replaces the external TTP wrapper (`astroq.nplan.NightPlanner`) by reusing
the splan packing/cadence formulation restricted to a single night, at fine
slot granularity, with a slew floor baked into the cover constraint.

Plumbing reuses:
  - `astroq.access.Access` for the alt/az + moon + custom + allocation cube
    (instantiated for ONE night, fine slot_size, restricted to splan's
    selected targets).
  - `astroq.splan.SemesterPlanner` loaded from h5 for `past_history`, etc.
  - `astroq.nplan.get_nightly_times_from_allocation` for the night window.

Variables:
  Y[i, s] in {0, 1}      (i, s) in A
  W[i]    in {0, 1}      i in multi-visit

C1a (cover form: at most one visit occupies any slot)
  For every slot s:
      Sum_{(i, s') : s' <= s < s' + t_visit[i] + slew_floor}  Y[i, s']  <=  1
  This is the aggregated splan-style "reserve multislot exposures" form
  (Lubin et al. 2025 Constraint 1), applied to one night. `slew_floor` is a
  per-target slot inflation (default 1 slot) that bakes in the minimum slew
  between any two consecutive visits.

  Aggregating identical clique inequalities into one row per slot (rather
  than one row per (i, s)) is the textbook LP-tightening for time-indexed
  scheduling formulations and is the reason splan's LP relaxation is
  integer-tight at the semester scale.

C1b (pair-specific slew refinement)
  C1a uses a constant `slew_floor` for every target pair. The true slew
  between consecutive visits j and i, evaluated at the slot s where i
  starts, can exceed that floor. For every (i, s) we add the residual
  forbidden range only:
      Y[i, s] + Σ_{j != i}
                 Σ_{s' in [s - t_visit[j] - slew(j, i, s) + 1,
                           s - t_visit[j] - slew_floor + 1)} Y[j, s']
      <=  1
  These rows are sparse (only the slew-excess slots), so the cover form
  still does the heavy LP lifting while exact pair-wise slew is honored.

C4 (intra-night cadence)
  For each (i, s) in A with n_intra_max[i] > 1 and tau_intra_slots[i] > 1:
      Y[i, s] + Sum_{s' in (s, s + tau_intra_slots[i]) : (i, s') in A} Y[i, s']
      <= W[i].

C5 (per-night visit count)
  Multi-visit: n_intra_min[i] * W[i] <= Sum_s Y[i, s] <= n_intra_max[i] * W[i].
  Single-visit:                          Sum_s Y[i, s] <= 1.

Objective
  maximize Sum_{(i, s) in A}  Y[i, s]      (count of exposures fit)

References
----------
- Lubin+ 2025 (arXiv:2506.08195): semester MILP we mirror.
- Handley+ 2024 (arXiv:2310.18497): TTP, replaced by this module.
"""

import logging
import os
import time
from configparser import ConfigParser

import gurobipy as gp
import numpy as np
import pandas as pd
from astropy.time import Time
from gurobipy import GRB

from astroq.access import Access
from astroq.nplan import get_nightly_times_from_allocation
from astroq.slew import slew_seconds_at
from astroq.splan import SemesterPlanner

logs = logging.getLogger(__name__)


class NightPlannerMILP:
    """In-house Gurobi ILP night planner.

    Reuses `astroq.access.Access` for the access cube; this class only owns
    the ILP itself (Y/W variables, slew-aware C1, C4, C5, count objective).
    """

    def __init__(self, config_file):
        config = ConfigParser()
        config.read(config_file)
        self.config = config

        workdir = config.get("global", "workdir")
        self.workdir = workdir
        self.current_day = config.get("global", "current_day")
        self.observatory_string = config.get("global", "observatory")
        self.output_directory = os.path.join(workdir, "outputs")
        os.makedirs(self.output_directory, exist_ok=True)

        def _resolve(section, key):
            val = config.get(section, key)
            return val if os.path.isabs(val) else os.path.join(workdir, val)

        self.allocation_file = _resolve("data", "allocation_file")
        self.custom_file = _resolve("data", "custom_file")

        def _get(section, key, fallback, kind=str):
            if not config.has_option(section, key):
                return fallback
            if kind is int:
                return config.getint(section, key)
            if kind is float:
                return config.getfloat(section, key)
            if kind is bool:
                return config.getboolean(section, key)
            return config.get(section, key)

        # 30 s default. On the 2026-05-09 full-band1 reference night the
        # splan-style aggregated cover + pair-specific slew tail closes to
        # OPTIMAL at every slot resolution tested, and total scheduled slew
        # is minimized at 30 s (51.1 min, vs 57.7 min at 60 s and 56.8 min
        # at 15 s); 15 s costs ~4x the solve time (74 s vs 18 s) without
        # cutting slew further. See plan's Empirical record.
        self.slot_size_seconds = _get("night", "slot_size_seconds", 30, int)
        self.max_solve_time = _get("night", "max_solve_time", 300, int)
        self.max_solve_gap = _get("night", "max_solve_gap", 0.005, float)
        self.show_gurobi_output = _get("night", "show_gurobi_output", True, bool)

        sp_path = os.path.join(self.output_directory, "semester_planner.h5")
        if not os.path.exists(sp_path):
            raise FileNotFoundError(
                f"semester_planner.h5 not found at {sp_path}. "
                "Run plan-semester first."
            )
        self.semester_planner = SemesterPlanner.from_hdf5(sp_path)

        selected_path = os.path.join(self.output_directory, "request_selected.csv")
        if not os.path.exists(selected_path):
            raise FileNotFoundError(
                f"request_selected.csv not found at {selected_path}."
            )
        self.selected_df = self._load_selected(selected_path)
        self.N = len(self.selected_df)

        self._warm_start_schedule = None
        self._warm_start_ratio = 1

    def set_warm_start_from_coarser(self, schedule_df, ratio):
        """Provide a coarser-resolution schedule to use as MIP start.

        Args:
            schedule_df: schedule DataFrame from a previous run at slot size
                `ratio * self.slot_size_seconds`.
            ratio: integer ratio (e.g. 2 means previous slot was 2x coarser).
        """
        self._warm_start_schedule = schedule_df
        self._warm_start_ratio = int(ratio)

    @staticmethod
    def _load_selected(path):
        df = pd.read_csv(path)
        df["unique_id"] = df["unique_id"].astype(str)
        # Defensive: request_selected.csv comes from splan post-cleaning so
        # these should already be numeric / non-"None", but the upstream
        # webform has flaked before. Cf. splan.SemesterPlanner.__init__.
        for col, default in [
            ("n_intra_max", 1),
            ("n_intra_min", 1),
            ("tau_intra", 0.0),
            ("tau_inter", 1.0),
            ("minimum_elevation", 33.0),
            ("minimum_moon_separation", 30.0),
        ]:
            if col in df.columns:
                df[col] = df[col].replace("None", np.nan).fillna(default)
        for col in [
            "ra", "dec", "exptime", "n_exp", "n_intra_max", "n_intra_min",
            "tau_intra", "tau_inter", "minimum_elevation",
            "minimum_moon_separation",
        ]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        df.reset_index(drop=True, inplace=True)
        return df

    def _build_t_visit_dict(self):
        """Slots needed per visit (exposure + readout only) at fine slot_size.

        Differs from `splan._build_slots_required_dictionary` by dropping the
        `60 * n_intra_max` slew-budget term. Splan folds an implicit per-visit
        slew allowance into t_visit because it has no explicit slew model;
        this module models slew explicitly via the pair-specific C1
        constraint, so including it here would double-count.

        Formula:
            t_visit_seconds = exptime * n_exp + 45 * (n_exp - 1)
        """
        out = {}
        for _, row in self.selected_df.iterrows():
            n_exp = int(row["n_exp"])
            exptime = float(row["exptime"])
            total_s = exptime * n_exp + 45.0 * (n_exp - 1)
            slots = int(np.ceil(total_s / self.slot_size_seconds))
            out[row["unique_id"]] = max(slots, 1)
        return out

    def build_access(self):
        """Use astroq.access.Access for the one-night access cube."""
        self.slots_needed_for_exposure_dict = self._build_t_visit_dict()
        self.t_visit = np.array(
            [
                self.slots_needed_for_exposure_dict[uid]
                for uid in self.selected_df["unique_id"]
            ]
        )

        # slot_size in Access is "minutes" but is used as a multiplier with
        # u.min, so floats work fine -- 0.25 min = 15 s.
        slot_size_min = self.slot_size_seconds / 60.0

        self.access_obj = Access(
            semester_start_date=self.current_day,
            semester_length=1,
            n_nights_in_semester=1,
            today_starting_night=0,
            current_day=self.current_day,
            all_dates_dict={self.current_day: 0},
            all_dates_array=[self.current_day],
            slot_size=slot_size_min,
            slots_needed_for_exposure_dict=self.slots_needed_for_exposure_dict,
            custom_file=self.custom_file,
            allocation_file=self.allocation_file,
            past_history={},  # splan already accounted for tau_inter on selection
            output_directory=self.output_directory,
            run_weather_loss=False,
            run_band3=False,
            observatory_string=self.observatory_string,
            request_frame=self.selected_df,
        )
        access_record = self.access_obj.produce_ultimate_map()
        self.access_record = access_record
        # Single-night arrays: shape (N, n_slots).
        self.is_observable = access_record.is_observable[:, 0, :]
        self.is_alloc = access_record.is_alloc[:, 0, :]
        self.n_slots = self.is_observable.shape[1]
        # Slot midpoints for slew time computation.
        self.slot_midpoints = self.access_obj.slotmidpoints[0]
        # First / last allocated slot indices (assumed contiguous for now).
        alloc_any = self.is_alloc.any(axis=0)
        if not alloc_any.any():
            raise RuntimeError(
                "No allocated slots tonight per the access cube; check allocation.csv."
            )
        self.alloc_start_slot = int(np.argmax(alloc_any))
        self.alloc_stop_slot = int(self.n_slots - np.argmax(alloc_any[::-1]))

    def compute_altaz(self):
        """Per-target alt/az at every slot midpoint."""
        altaz = self.access_obj.observatory.altaz(
            self.slot_midpoints,
            self.access_obj.targets,
            grid_times_targets=True,
        )
        self.alts = altaz.alt.deg  # (N, n_slots)
        self.azs = altaz.az.deg

    def slew_slots_at(self, s):
        """Pair-specific slew time in slots, evaluated at slot s."""
        slew_s = slew_seconds_at(self.alts[:, s], self.azs[:, s])
        slots = np.ceil(slew_s / self.slot_size_seconds).astype(int)
        np.fill_diagonal(slots, 0)
        return slots

    def build_model(self):
        """Construct the Gurobi ILP."""
        m = gp.Model("NightPlannerMILP")
        m.Params.OutputFlag = 1 if self.show_gurobi_output else 0
        m.Params.TimeLimit = self.max_solve_time
        m.Params.MIPGap = self.max_solve_gap
        # Barrier root LP and aggressive heuristics are both important at
        # fine slot resolutions. With the cover-form C1 the presolve picture
        # is benign: empirically Presolve=1 (conservative) removes ~180 rows
        # in ~1.4 s and trims ~2 s overall, while Presolve=2 (aggressive)
        # burns the full TimeLimit chasing marginal row removals and never
        # gets to the root LP. Keep it overridable so we can re-sweep.
        m.Params.Method = 2
        m.Params.Presolve = getattr(self, "presolve_override", 1)
        m.Params.MIPFocus = 1
        m.Params.Heuristics = 0.5

        rows, slots = np.where(self.is_observable)
        if len(rows) == 0:
            raise RuntimeError(
                "is_observable is empty everywhere; nothing to schedule."
            )
        obs_pairs = list(zip(rows.tolist(), slots.tolist()))
        self.obs_pairs = obs_pairs

        Y = m.addVars(obs_pairs, vtype=GRB.BINARY, name="Y")
        self.Y = Y

        n_intra_max = self.selected_df["n_intra_max"].astype(int).values
        n_intra_min = self.selected_df["n_intra_min"].astype(int).values
        tau_intra_hr = self.selected_df["tau_intra"].astype(float).values
        # tau_intra is documented as hours in request.csv.
        tau_intra_slots = np.ceil(
            tau_intra_hr * 3600.0 / self.slot_size_seconds
        ).astype(int)

        multi_visit_idx = [i for i in range(self.N) if n_intra_max[i] > 1]
        if multi_visit_idx:
            W = m.addVars(multi_visit_idx, vtype=GRB.BINARY, name="W")
        else:
            W = {}
        self.W = W
        self.multi_visit_idx = set(multi_visit_idx)

        # Per-slot index of who can start where, plus per-target slot lists.
        starts_at = [[] for _ in range(self.n_slots)]
        starts_for = [[] for _ in range(self.N)]
        for i, s in obs_pairs:
            starts_at[s].append(i)
            starts_for[i].append(s)
        # Each starts_for[i] is naturally sorted ascending (np.where output).

        # ----- C1a: splan-style aggregated cover constraint -----
        # For each slot s, at most one visit may be occupying s:
        #
        #   Σ_{(i, s') : s - t_block[i] + 1 ≤ s' ≤ s} Y[i, s'] ≤ 1
        #
        # with t_block[i] = t_visit[i] + slew_floor. This is the splan
        # "reserve multislot exposures" form (Lubin et al. 2025 Constraint 1)
        # adapted to one night: aggregating identical clique inequalities
        # into one row per slot is what gives splan its integer-tight LP
        # relaxation, and the same effect carries over here.
        #
        # The cover form alone treats slew as a constant per-target floor of
        # `slew_floor` slots regardless of the actual transit between the
        # two consecutive visits. C1b below adds the *pair-specific* tail.
        slew_floor = 1
        t_block = self.t_visit + slew_floor
        max_t_block = int(t_block.max())

        for s in range(self.n_slots):
            cover_terms = [(i, s) for i in starts_at[s]]
            for delta in range(1, max_t_block):
                s_prev = s - delta
                if s_prev < 0:
                    break
                for i in starts_at[s_prev]:
                    if t_block[i] >= delta + 1:
                        cover_terms.append((i, s_prev))
            if len(cover_terms) < 2:
                continue
            m.addConstr(
                gp.quicksum(Y[i, s_prev] for (i, s_prev) in cover_terms) <= 1,
                name=f"cover_{s}",
            )

        # ----- C1b: pair-specific slew refinement -----
        # The cover constraint above forbids any earlier start (j, s_prev)
        # with s_prev ≥ s - t_visit[j] - slew_floor + 1 from coexisting
        # with Y[i, s]. The *true* pair-specific feasibility requires
        #
        #     s - s_prev ≥ t_visit[j] + slew_slots(j, i, evaluated at s)
        #
        # which, when slew_slots(j, i, s) > slew_floor, forbids an
        # additional `slew_slots(j, i, s) - slew_floor` slots of earlier
        # starts. We add only the residual forbidden range here so the
        # cover form still does the heavy LP lifting:
        #
        #     Y[i, s] + Σ_j Σ_{s_prev ∈ R_j(i, s)} Y[j, s_prev]  ≤  1
        #
        # where R_j(i, s) = [s - t_visit[j] - slew(j, i, s) + 1,
        #                    s - t_visit[j] - slew_floor + 1).
        #
        # Constraints are sparse (only the slew-excess slots, not the full
        # t_visit window), so total nonzero count stays comparable to the
        # cover form.
        t_visit = self.t_visit
        for s in range(self.n_slots):
            if not starts_at[s]:
                continue
            slew_at_s = self.slew_slots_at(s)
            for i in starts_at[s]:
                extra_terms = []
                for j in range(self.N):
                    if j == i:
                        continue
                    sl = int(slew_at_s[j, i])
                    if sl <= slew_floor:
                        continue
                    lo = max(0, s - int(t_visit[j]) - sl + 1)
                    hi = s - int(t_visit[j]) - slew_floor + 1
                    if lo >= hi:
                        continue
                    sj = starts_for[j]
                    a = np.searchsorted(sj, lo, side="left")
                    b = np.searchsorted(sj, hi, side="left")
                    for s_prev in sj[a:b]:
                        extra_terms.append(Y[j, s_prev])
                if extra_terms:
                    m.addConstr(
                        Y[i, s] + gp.quicksum(extra_terms) <= 1,
                        name=f"pair_slew_{i}_{s}",
                    )

        # ----- C4 + C5: intra-night cadence and per-night visit count -----
        for i in range(self.N):
            slots_i = starts_for[i]
            if not slots_i:
                continue
            if i in self.multi_visit_idx:
                w_i = W[i]
                m.addConstr(
                    gp.quicksum(Y[i, s] for s in slots_i)
                    <= int(n_intra_max[i]) * w_i,
                    name=f"max_visits_{i}",
                )
                m.addConstr(
                    gp.quicksum(Y[i, s] for s in slots_i)
                    >= int(n_intra_min[i]) * w_i,
                    name=f"min_visits_{i}",
                )
                tau = int(tau_intra_slots[i])
                if tau > 1:
                    slots_arr = np.array(slots_i)
                    for s in slots_i:
                        s_within = slots_arr[
                            (slots_arr > s) & (slots_arr < s + tau)
                        ]
                        if len(s_within) == 0:
                            continue
                        m.addConstr(
                            Y[i, s]
                            + gp.quicksum(Y[i, int(s2)] for s2 in s_within)
                            <= w_i,
                            name=f"intra_cad_{i}_{s}",
                        )
            else:
                m.addConstr(
                    gp.quicksum(Y[i, s] for s in slots_i) <= 1,
                    name=f"single_visit_{i}",
                )

        # ----- Objective: maximize visits scheduled -----
        # splan already chose which targets are valuable; the night planner
        # should fit as many as possible. Mirrors Handley+ 2024 eq. 10.
        m.setObjective(
            gp.quicksum(Y[i, s] for (i, s) in obs_pairs), GRB.MAXIMIZE
        )

        # ----- Optional MIP warm start from an outer schedule -----
        # The fine-grained model (30 s slots) is dominated by the time
        # Gurobi spends in root-node cutting before it ever branches; a
        # known-feasible schedule from a coarser run lets it skip the
        # search effort it would otherwise need to discover the same
        # incumbent. Lifted schedules from a 2x-coarser solve are always
        # feasible because slew_slots increases by at most a factor of 2
        # under slot halving, and t_visit_slots scales the same way.
        if self._warm_start_schedule is not None:
            uid_to_idx = {
                uid: i for i, uid in enumerate(self.selected_df["unique_id"])
            }
            ratio = self._warm_start_ratio
            # Build the full assignment: every Y is 0 unless explicitly the
            # lifted start. Partial starts let Gurobi try to fill in the
            # blanks with an LP, which often picks integer values that
            # violate C1 against our enforced ones.
            start_set = set()
            for _, row in self._warm_start_schedule.iterrows():
                uid = row["unique_id"]
                if uid not in uid_to_idx:
                    continue
                i = uid_to_idx[uid]
                s_fine = int(row["slot"]) * ratio
                if (i, s_fine) in Y:
                    start_set.add((i, s_fine))
            n_loaded = n_zero = 0
            for key, var in Y.items():
                if key in start_set:
                    var.Start = 1.0
                    n_loaded += 1
                else:
                    var.Start = 0.0
                    n_zero += 1
            # W variables (if any) consistent with the lifted starts.
            for i in self.multi_visit_idx:
                W[i].Start = 1.0 if any(k[0] == i for k in start_set) else 0.0
            logs.info(
                f"warm start: loaded={n_loaded} zeroed={n_zero} "
                f"(ratio={ratio})"
            )
        elif os.environ.get("NPLAN_MILP_WARM_START", "0") == "1":
            warm = self._greedy_warm_start(starts_for, starts_at)
            for (i, s), val in warm.items():
                Y[i, s].Start = val

        self.model = m
        return m

    def _greedy_warm_start(self, starts_for, starts_at):
        """Return a feasible {(i, s): 1.0} dict for a greedy schedule.

        Time-ordered greedy: at each slot, place the unscheduled target with
        the smallest t_visit that is observable now AND whose slew from every
        prior placed visit (evaluated at slot s) keeps the prior end before
        s. The "every prior" check is essential because slot-and-pair-
        specific slew can make a 2-back visit's slew to i exceed the 1-back
        visit's slew to i, so a myopic 1-back check can produce schedules
        that violate C1.
        """
        N = self.N
        t_visit = self.t_visit
        n_intra_max = self.selected_df["n_intra_max"].astype(int).values
        scheduled_count = np.zeros(N, dtype=int)
        # (prev_idx, prev_end_slot) for every previously placed visit.
        placed = []
        starts = {}
        s = self.alloc_start_slot
        while s < self.alloc_stop_slot:
            candidates = [
                i for i in starts_at[s]
                if scheduled_count[i] < n_intra_max[i]
            ]
            if not candidates:
                s += 1
                continue
            if placed:
                slew_at_s = self.slew_slots_at(s)
                ok = []
                for i in candidates:
                    # All previously placed visits must clear i's slew window.
                    feasible = True
                    for (j, end_j) in placed:
                        if end_j + int(slew_at_s[j, i]) > s:
                            feasible = False
                            break
                    if feasible:
                        ok.append(i)
                if not ok:
                    s += 1
                    continue
                candidates = ok
            i = min(candidates, key=lambda k: t_visit[k])
            end_slot = s + int(t_visit[i])
            if end_slot > self.alloc_stop_slot:
                s += 1
                continue
            starts[(i, s)] = 1.0
            scheduled_count[i] += 1
            placed.append((i, end_slot))
            s = end_slot
        return starts

    def solve(self):
        t0 = time.time()
        self.model.optimize()
        wall = time.time() - t0
        return {
            "status": int(self.model.Status),
            "wall_seconds": wall,
            "obj": float(self.model.ObjVal) if self.model.SolCount > 0 else None,
            "bound": float(self.model.ObjBound) if self.model.SolCount > 0 else None,
            "mip_gap": float(self.model.MIPGap) if self.model.SolCount > 0 else None,
            "node_count": int(self.model.NodeCount),
            "n_vars": int(self.model.NumVars),
            "n_constrs": int(self.model.NumConstrs),
            "n_nonzeros": int(self.model.NumNZs),
        }

    def extract_schedule(self):
        """Read out (i, s) pairs where Y == 1, sorted by slot."""
        if self.model.SolCount == 0:
            return pd.DataFrame(
                columns=[
                    "unique_id", "starname", "slot", "start_utc",
                    "t_visit_slots", "exptime", "ra", "dec",
                ]
            )
        rows = []
        for (i, s), var in self.Y.items():
            if var.X > 0.5:
                rows.append(
                    {
                        "unique_id": self.selected_df["unique_id"].iloc[i],
                        "starname": self.selected_df["starname"].iloc[i],
                        "slot": s,
                        "start_utc": str(self.slot_midpoints[s].iso),
                        "t_visit_slots": int(self.t_visit[i]),
                        "exptime": float(self.selected_df["exptime"].iloc[i]),
                        "ra": float(self.selected_df["ra"].iloc[i]),
                        "dec": float(self.selected_df["dec"].iloc[i]),
                    }
                )
        return pd.DataFrame(rows).sort_values("slot").reset_index(drop=True)

    def compute_slew_summary(self, schedule_df):
        """Total slew time in seconds for the actual scheduled order."""
        if len(schedule_df) < 2:
            return {"n_transitions": 0, "total_slew_seconds": 0.0,
                    "total_slew_minutes": 0.0}
        uid_to_idx = {uid: i for i, uid in enumerate(self.selected_df["unique_id"])}
        slots = schedule_df["slot"].tolist()
        idxs = [uid_to_idx[u] for u in schedule_df["unique_id"]]
        total = 0.0
        for k in range(len(idxs) - 1):
            j, i = idxs[k], idxs[k + 1]
            # Use slew at the *next* exposure's slot (matches C1's convention).
            s_next = slots[k + 1]
            slew_s = slew_seconds_at(self.alts[:, s_next], self.azs[:, s_next])
            total += float(slew_s[j, i])
        return {
            "n_transitions": len(idxs) - 1,
            "total_slew_seconds": total,
            "total_slew_minutes": total / 60.0,
        }

    def run(self, log_prefix=""):
        """End-to-end build + solve. Returns (summary, schedule_df, slew_stats)."""
        t0 = time.time()
        self.build_access()
        n_obs = int(self.is_observable.sum())
        logs.info(
            f"{log_prefix}access: N={self.N} n_slots={self.n_slots} "
            f"alloc_slots=[{self.alloc_start_slot}, {self.alloc_stop_slot}) "
            f"valid_(r,s)={n_obs} ({time.time()-t0:.2f}s)"
        )
        self.compute_altaz()
        logs.info(f"{log_prefix}altaz computed ({time.time()-t0:.2f}s)")
        self.build_model()
        logs.info(
            f"{log_prefix}model: vars={self.model.NumVars} "
            f"constrs={self.model.NumConstrs} nzs={self.model.NumNZs} "
            f"({time.time()-t0:.2f}s build)"
        )
        summary = self.solve()
        logs.info(
            f"{log_prefix}solve: status={summary['status']} "
            f"obj={summary['obj']} gap={summary['mip_gap']} "
            f"nodes={summary['node_count']} "
            f"wall={summary['wall_seconds']:.2f}s"
        )
        schedule_df = self.extract_schedule()
        slew_stats = self.compute_slew_summary(schedule_df)
        logs.info(
            f"{log_prefix}schedule: {len(schedule_df)} visits, "
            f"total slew {slew_stats['total_slew_minutes']:.2f} min"
        )
        return summary, schedule_df, slew_stats
