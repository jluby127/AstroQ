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
  B[i, s] in [0, 1]      (i, s) with s in [s_lo[i], s_hi[i]]  (continuous)

C1a (delta-form occupancy cover: at most one visit occupies any slot)
  Equivalence-preserving rewrite of the aggregated splan-style "reserve
  multislot exposures" cover (Lubin et al. 2025 Constraint 1), expressed as
  a sliding-sum recurrence on a continuous occupancy indicator:
      B[i, s] = Σ_{s' = s - t_block[i] + 1}^{s}  Y[i, s']
  where t_block[i] = t_visit[i] + slew_floor and slew_floor (default 1 slot)
  bakes in the minimum slew between any two consecutive visits.

  The cover is then three sparse blocks:
      B_init :  B[i, s_lo[i]] = Y[i, s_lo[i]]                  (1 row per i)
      B_rec  :  B[i, s] - B[i, s-1] - Y[i, s] + Y[i, s - t_block[i]] = 0
                                                                (1 row per (i, s))
      B_cap  :  Σ_i B[i, s] <= 1                               (1 row per s)
  Each recurrence row has at most 4 nonzeros; the cap row has at most one
  entry per target with active B support at slot s. Treating B as continuous
  is exact: the recurrence forces B integral whenever Y is integral, and
  projecting onto Y reproduces the original aggregated cover polytope (same
  LP relaxation). The delta form sidesteps the sliding-shift density that
  Gurobi's Sparsify otherwise has to discover at presolve time.

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
from astropy.time import Time, TimeDelta
from gurobipy import GRB

from astroq.access import Access
from astroq.nplan import get_nightly_times_from_allocation
from astroq.slew import slew_seconds_at
from astroq.splan import SemesterPlanner

logs = logging.getLogger(__name__)


class NightPlannerILP:
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
        # Optional in some configs but expected by save_night_planner_hdf5.
        if config.has_option("data", "filler_file"):
            self.filler_file = _resolve("data", "filler_file")
        else:
            self.filler_file = ""

        # Aliases / extras the legacy NightPlanner sets and that the HDF5
        # writer + downstream consumers (webapp, astroq plot) expect.
        self.upstream_path = workdir
        self.semester_directory = workdir
        self.reports_directory = self.output_directory

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
        self.past_history = self.semester_planner.past_history

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
        m = gp.Model("NightPlannerILP")
        m.Params.OutputFlag = 1 if self.show_gurobi_output else 0
        m.Params.TimeLimit = self.max_solve_time
        m.Params.MIPGap = self.max_solve_gap
        gurobi_log_file = getattr(self, "gurobi_log_file", None)
        if gurobi_log_file:
            m.Params.LogFile = gurobi_log_file
        # Barrier root LP and aggressive heuristics are both important at
        # fine slot resolutions. With the cover-form C1 the presolve picture
        # is benign: empirically Presolve=1 (conservative) removes ~180 rows
        # in ~1.4 s and trims ~2 s overall, while Presolve=2 (aggressive)
        # burns the full TimeLimit chasing marginal row removals and never
        # gets to the root LP. Keep it overridable so we can re-sweep.
        m.Params.Method = 2
        m.Params.Presolve = 2
        m.Params.Symmetry = 2
        m.Params.Aggregate = 2
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

        # ----- C1a: delta-form occupancy cover -----
        # Equivalence-preserving rewrite of the splan-style aggregated
        # cover constraint. Introduces a continuous occupancy indicator
        #
        #     B[i, s] = Σ_{s' = s - t_block[i] + 1}^{s} Y[i, s']
        #
        # i.e. B[i, s] = 1 iff target i is mid-visit at slot s, with the
        # per-target footprint length t_block[i] = t_visit[i] + slew_floor
        # baking in a constant slew floor. C1b below adds the pair-specific
        # slew tail.
        #
        # Three sparse blocks replace the dense sliding-window cover:
        #
        #   B_init : B[i, s_lo[i]] = Y[i, s_lo[i]]           (one row per i)
        #   B_rec  : B[i, s] - B[i, s-1] - Y[i, s] + Y[i, s - t_block[i]] = 0
        #                                                     (one row per (i, s))
        #   B_cap  : Σ_i B[i, s] <= 1                        (one row per s)
        #
        # Each recurrence row carries at most 4 nonzeros; the cap row
        # carries one entry per target with active B support at slot s.
        # Treating B as continuous in [0, 1] does not relax the integer
        # program: B is forced integral by the recurrence whenever Y is
        # integral. Projecting onto Y reproduces today's cover polytope
        # exactly, so the root LP bound is unchanged.
        slew_floor = 1
        t_block = self.t_visit + slew_floor

        obs_pairs_set = set(obs_pairs)

        # Per-target B support [s_lo[i], s_hi[i]]. Targets with no
        # observable start (starts_for[i] empty) contribute no B rows.
        s_lo = np.full(self.N, -1, dtype=int)
        s_hi = np.full(self.N, -1, dtype=int)
        for i in range(self.N):
            if not starts_for[i]:
                continue
            s_lo[i] = int(starts_for[i][0])
            s_hi[i] = min(
                self.n_slots - 1,
                int(starts_for[i][-1]) + int(t_block[i]) - 1,
            )

        B_pairs = [
            (i, s)
            for i in range(self.N)
            if s_lo[i] >= 0
            for s in range(int(s_lo[i]), int(s_hi[i]) + 1)
        ]
        B = m.addVars(
            B_pairs, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="B"
        )
        self.B = B

        # Per-slot covering set for the cap row.
        covers_slot = [[] for _ in range(self.n_slots)]
        for (i, s) in B_pairs:
            covers_slot[s].append(i)

        # B_init: anchor B at the first slot of each target's support.
        for i in range(self.N):
            if s_lo[i] < 0:
                continue
            s0 = int(s_lo[i])
            # By construction s0 = starts_for[i][0], so (i, s0) ∈ obs_pairs.
            m.addConstr(B[i, s0] == Y[i, s0], name=f"B_init_{i}")

        # B_rec: telescoping recurrence over each target's B support.
        for i in range(self.N):
            if s_lo[i] < 0:
                continue
            tb = int(t_block[i])
            for s in range(int(s_lo[i]) + 1, int(s_hi[i]) + 1):
                y_in = Y[i, s] if (i, s) in obs_pairs_set else 0
                s_out = s - tb
                y_out = (
                    Y[i, s_out]
                    if s_out >= 0 and (i, s_out) in obs_pairs_set
                    else 0
                )
                m.addConstr(
                    B[i, s] - B[i, s - 1] - y_in + y_out == 0,
                    name=f"B_rec_{i}_{s}",
                )

        # B_cap: at most one target mid-visit at any given slot.
        for s in range(self.n_slots):
            if not covers_slot[s]:
                continue
            m.addConstr(
                gp.quicksum(B[i, s] for i in covers_slot[s]) <= 1,
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
        elif os.environ.get("NPLAN_ILP_WARM_START", "0") == "1":
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
        """End-to-end build + solve. Returns (summary, schedule_df, slew_stats).

        Also exposes a TTP-shaped `self.solution = [obj]` so downstream consumers
        (webapp, `astroq plot`, `hirescps.write_starlist`, HDF5 round-trip)
        keep working without code changes.
        """
        from astroq.nplan import get_nightly_times_from_allocation
        try:
            self.observation_start_time, self.observation_stop_time = (
                get_nightly_times_from_allocation(self.allocation_file, self.current_day)
            )
        except ValueError:
            logs.warning(
                f"{log_prefix}No allocation for {self.current_day}; "
                "ILP cannot run, returning empty result."
            )
            self.solution = None
            return (
                {"status": -1, "wall_seconds": 0.0, "obj": None,
                 "bound": None, "mip_gap": None, "node_count": 0,
                 "n_vars": 0, "n_constrs": 0, "n_nonzeros": 0},
                pd.DataFrame(),
                {"n_transitions": 0, "total_slew_seconds": 0.0,
                 "total_slew_minutes": 0.0},
            )

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

        self.schedule_df = schedule_df
        self.summary = summary
        self.slew_stats = slew_stats
        self.solution = [self._build_solution(schedule_df, summary, slew_stats)]
        return summary, schedule_df, slew_stats

    # ------------------------------------------------------------------
    # TTP-shaped solution shim + on-disk artifacts
    # ------------------------------------------------------------------

    @staticmethod
    def _keck1_observatory():  # exposed via helper below for clarity
        from ttp import telescope
        return telescope.Keck1()

    def _get_first_last_minutes(self, target_idx):
        """First/last observable slot (relative to obs start time) in minutes.

        Mirrors the legacy `NightPlanner.get_first_last_indices` (which used a
        coarser slot grid via `semester_planner.access_record`) but operates on
        the fine ILP access cube and returns numeric minutes-from-night-start
        to match the TTP `solution.plotly['First/Last Available']` convention.
        """
        observable = self.is_observable[target_idx, :]
        true_idx = np.where(observable)[0]
        slot_size_min = self.slot_size_seconds / 60.0
        if len(true_idx) == 0:
            # Sentinel: clamp to night end. TTP returns the literal allocation
            # window in this case; matching that is fine since downstream code
            # only renders these for display.
            return 0.0, 0.0
        first_slot = int(true_idx[0])
        last_slot = int(true_idx[-1])
        first_min = first_slot * slot_size_min - self._slot0_offset_min()
        last_min = last_slot * slot_size_min - self._slot0_offset_min()
        return first_min, last_min

    def _slot0_offset_min(self):
        """Minutes between Access's slot 0 (00:00 UTC of current_day) and obs start."""
        return (
            (self.observation_start_time
             - self.access_obj.daily_start).sec / 60.0
        )

    def _build_solution(self, schedule_df, summary, slew_stats):
        """Build a TTP-shaped `solution` object from the ILP schedule.

        Downstream consumers (write_starlist, get_script_plan, get_ladder,
        plot_path_2D_interactive, get_slew_animation_plotly, save_night_planner_hdf5)
        treat `solution.plotly` as a dict-of-lists keyed by Starname / Start Exposure
        / Stop Exposure / First Available / Last Available / Total Exp Time (min) /
        Minutes the from Start of the Night / human_starname / UTC Start Time.
        Times are minutes from `nightstarts`.
        """
        from types import SimpleNamespace
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        slot_size_min = self.slot_size_seconds / 60.0
        slot0_offset_min = self._slot0_offset_min()
        obs_start = self.observation_start_time
        obs_stop = self.observation_stop_time
        dur_min = float((obs_stop.jd - obs_start.jd) * 24.0 * 60.0)

        # Build a per-row id_to_idx into self.selected_df.
        uid_to_idx = {uid: i for i, uid in enumerate(self.selected_df["unique_id"])}

        starnames = []
        human = []
        starts = []
        stops = []
        firsts = []
        lasts = []
        totals = []
        exposure_mins = []
        n_shots_list = []
        utc_strs = []
        stars = []
        times_list = []
        az_path = []
        alt_path = []

        for _, row in schedule_df.iterrows():
            uid = row["unique_id"]
            i = uid_to_idx[uid]
            s = int(row["slot"])
            t_visit_slots = int(row["t_visit_slots"])
            start_min = s * slot_size_min - slot0_offset_min
            total_min = t_visit_slots * slot_size_min
            stop_min = start_min + total_min
            first_min, last_min = self._get_first_last_minutes(i)

            starnames.append(uid)
            human.append(str(self.selected_df["starname"].iloc[i]))
            starts.append(round(start_min, 2))
            stops.append(round(stop_min, 2))
            firsts.append(round(first_min, 2))
            lasts.append(round(last_min, 2))
            totals.append(round(total_min, 2))
            exposure_mins.append(round(float(self.selected_df["exptime"].iloc[i]) / 60.0, 2))
            n_shots_list.append(int(self.selected_df["n_exp"].iloc[i]))

            start_abs = obs_start + TimeDelta(start_min * 60.0, format="sec")
            utc_strs.append(str(start_abs)[11:16])
            times_list.append(start_abs)

            star = SimpleNamespace(
                name=uid,
                target=SkyCoord(
                    ra=float(self.selected_df["ra"].iloc[i]) * u.deg,
                    dec=float(self.selected_df["dec"].iloc[i]) * u.deg,
                ),
            )
            stars.append(star)
            az_path.append(float(self.azs[i, s]))
            alt_path.append(float(self.alts[i, s]))

        plotly = {
            "Starname": starnames,
            "human_starname": human,
            "Start Exposure": starts,
            "Stop Exposure": stops,
            "First Available": firsts,
            "Last Available": lasts,
            "Total Exp Time (min)": totals,
            "Exposure Time (min)": exposure_mins,
            "N_shots": n_shots_list,
            "Minutes the from Start of the Night": starts,  # alias TTP also sets
            "UTC Start Time": utc_strs,
            # Fixed priority 10 mirrors the legacy `to_ttp` "Priority": 10 used by
            # the TTP path; downstream ladder plotting colors by this.
            "Priority": [10] * len(starnames),
        }

        # Empty `extras` (no unscheduled-targets section per spec).
        extras = {k: [] for k in plotly}

        # `schedule` is consumed by get_slew_animation_plotly which needs
        # parallel `Starname` and `Time` (JD) arrays.
        schedule = {
            "Starname": list(starnames),
            "Time": [t.jd for t in times_list],
        }

        time_exposing = sum(totals)
        time_slewing = float(slew_stats["total_slew_minutes"])
        time_idle = max(0.0, dur_min - time_exposing - time_slewing)
        num_scheduled = len(schedule_df)

        solution = SimpleNamespace(
            plotly=plotly,
            extras=extras,
            schedule=schedule,
            stars=stars,
            times=times_list,
            nightstarts=obs_start,
            nightends=obs_stop,
            az_path=np.array(az_path),
            alt_path=np.array(alt_path),
            dur=dur_min,
            time_exposing=time_exposing,
            time_slewing=time_slewing,
            time_idle=time_idle,
            num_scheduled=num_scheduled,
            # TTP exposes N = real_targets + 2 anchors; mirror that so legacy
            # text like "Observations Requested: N - 2" still works.
            N=num_scheduled + 2,
            solve_time=float(summary.get("wall_seconds", 0.0)),
            observatory=self._keck1_observatory(),
        )
        return solution

    def write_outputs(self):
        """Write all on-disk artifacts the legacy TTP path produced.

        Files written into `self.output_directory`:
          - ttp_prepared.csv        (debug parity with legacy)
          - ObserveOrder_<date>.txt
          - script_<date>_nominal.txt
          - TTPstatistics.txt
        """
        from astroq.queue import hirescps

        if self.solution is None or len(self.solution[0].plotly["Starname"]) == 0:
            logs.warning(
                "ILP produced no scheduled visits; not writing TTP artifacts."
            )
            return False

        solution = self.solution[0]
        plotly = solution.plotly
        observers_path = self.output_directory
        os.makedirs(observers_path, exist_ok=True)

        # ---- ttp_prepared.csv (debug parity) ----
        to_ttp = pd.DataFrame({
            "Starname": self.selected_df["unique_id"],
            "RA": self.selected_df["ra"],
            "Dec": self.selected_df["dec"],
            "Exposure Time": self.selected_df["exptime"],
            "Exposures Per Visit": self.selected_df["n_exp"],
            "Visits In Night": self.selected_df["n_intra_max"],
            "Intra_Night_Cadence": self.selected_df["tau_intra"],
            "Priority": 10,
            "First Available": "",
            "Last Available": "",
        })
        to_ttp.to_csv(os.path.join(observers_path, "ttp_prepared.csv"), index=False)

        # ---- ObserveOrder_<date>.txt ----
        observe_order_file = os.path.join(
            observers_path, f"ObserveOrder_{self.current_day}.txt"
        )
        order_rows = []
        for i, uid in enumerate(plotly["Starname"]):
            order_rows.append({
                "unique_id": str(uid),
                "Target": plotly["human_starname"][i],
                "StartExposure": plotly["UTC Start Time"][i],
            })
        pd.DataFrame(order_rows).to_csv(observe_order_file, index=False)

        # ---- script_<date>_nominal.txt ----
        hirescps.write_starlist(
            self.selected_df,
            plotly,
            self.observation_start_time,
            solution.extras,
            [],
            str(self.current_day),
            observers_path,
            all_active_requests=self.semester_planner.requests_frame,
            past_history=self.past_history,
        )

        # ---- TTPstatistics.txt + matching stdout for check_night_plans.py ----
        stats_lines = [
            "Stats for ILP Solution",
            "------------------------------------",
            f"    Model ran for {solution.solve_time:.2f} seconds",
            f"     Observations Requested: {self.N}",
            f"     Observations Scheduled: {solution.num_scheduled}",
            "------------------------------------",
            f"   Observing Duration (min): {solution.dur:.2f}",
            f"  Time Spent Exposing (min): {solution.time_exposing:.2f}",
            f"      Time Spent Idle (min): {solution.time_idle:.2f}",
            f"   Time Spent Slewing (min): {solution.time_slewing:.2f}",
            "------------------------------------",
        ]
        ttp_stats_path = os.path.join(observers_path, "TTPstatistics.txt")
        with open(ttp_stats_path, "w") as f:
            f.write("\n".join(stats_lines) + "\n")
        # Emit to stdout so check_night_plans.py (which greps astroq.log) finds them.
        print("\n" + "\n".join(stats_lines))
        return True

    def to_hdf5(self, hdf5_path=None):
        """Write night_planner.h5 using the same format as the legacy engine."""
        from astroq.nplan import save_night_planner_hdf5
        if self.solution is None:
            logs.warning("ILP has no solution; skipping night_planner.h5 write.")
            return None
        return save_night_planner_hdf5(self, hdf5_path)
