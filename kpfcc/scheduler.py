"""
Module that defines the Scheduler class. This class is responsible for defining, building, and solving the
Gurobi model. It is nearly completely agnostic to all astronomy knowledge.

"""
import sys
import time
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB

# NOTE TO SELF: this is the last obstacle to totally divorcing from all astronomy knowledge
import kpfcc.access as ac

class Scheduler(object):
    """A Scheduler object, from which we can define a Gurobi model, build constraints, and solve."""

    def __init__(self, request_set, manager):

        print("Building the Scheduler.")
        self.manager = manager
        self.request_set = request_set

        if self.manager.run_optimal_allocation:
            if self.manager.semester_grid == [] or self.manager.quarters_grid == []:
                print("WARNING: Missing necessary parameters to run Optimal Allocation.")
                print("---- A semester grid and a quarters grid is required to run optimal allocation.")
            if self.manager.max_quarters == 0 or self.manager.max_unique_nights == 0:
                print("WARNING: Missing necessary parameters to run Optimal Allocation.")
                print("---- Neither Max Quarters and Max Unique Nights parameters cannot be zero.")

        self.model = gp.Model('Semester_Scheduler')
        # Yrs is technically a 1D matrix indexed by tuples.
        # But in practice best think of it as a 2D square matrix of requests r and slots s, with gaps.
        # Slot s for request r will be 1 to indicate starting an exposure for that request in that slot
        # self.Yrds = self.model.addVars(self.request_set.Aset, vtype = GRB.BINARY, name = 'Requests_Slots')
        tuple_array = [(row.id, row.d, row.s) for row in self.request_set.Aset.itertuples(index=False)]
        self.Yrds = self.model.addVars(tuple_array, vtype = GRB.BINARY, name = 'Requests_Slots')
        # theta is the "shortfall" variable, continous in natural numbers.
        self.theta = self.model.addVars(self.request_set.schedulable_requests, name = 'Shortfall')

        if self.manager.run_optimal_allocation:
            # Anq is a 2D matrix of N_nights_in_semester by N_quarters_in_night
            # element will be 1 if that night/quarter is allocated to KPF and 0 otherwise
            self.Anq = self.model.addVars(self.manager.semester_grid, self.manager.quarters_grid, vtype = GRB.BINARY, name = 'Allocation')

            # Un is a 1D matrix of N_nights_in_semester
            # element will be 1 if at least one quarter in that night is allocated
            self.Un = self.model.addVars(self.manager.semester_grid, vtype = GRB.BINARY, name = 'On-Sky')
            self.model.update()

        if self.request_set.Wset != []:
            # Wrd is technically a 1D matrix indexed by tuples.
            # But in practice best think of it as a 2D square matrix of requests r and nights d, with gaps.
            # Night d for request r will be 1 to indicate at least one exposure is scheduled for this night.
            # Note that Wrd is only valid for requests r which have at least 2 visits requested in the night.
            self.Wrd = self.model.addVars(self.request_set.Wset, vtype = GRB.BINARY, name = 'Nightly_Observation_Tracker')

        desired_max_obs_allowed_dict = {}
        absolute_max_obs_allowed_dict = {}
        past_nights_observed_dict = {}
        for name in self.request_set.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['Starname']==name][0]
            if self.manager.database_info_dict == {}:
                past_nights_observed = 0
            else:
                past_nights_observed = len(self.manager.database_info_dict[name][1])

            # Safety valve for if the target is over-observed for any reason
            if past_nights_observed > self.manager.requests_frame['# of Nights Per Semester'][idx] + \
                        int(self.manager.requests_frame['# of Nights Per Semester'][idx]*self.manager.max_bonus):
                desired_max_obs = past_nights_observed
            else:
                desired_max_obs = (self.manager.requests_frame['# of Nights Per Semester'][idx] - past_nights_observed)
                absolute_max_obs = (self.manager.requests_frame['# of Nights Per Semester'][idx] - past_nights_observed) \
                        + int(self.manager.requests_frame['# of Nights Per Semester'][idx]*self.manager.max_bonus)
                # second safety valve
                if past_nights_observed > absolute_max_obs:
                    absolute_max_obs = past_nights_observed
            past_nights_observed_dict[name] = past_nights_observed
            desired_max_obs_allowed_dict[name] = desired_max_obs
            absolute_max_obs_allowed_dict[name] = absolute_max_obs
            self.desired_max_obs_allowed_dict = desired_max_obs_allowed_dict
            self.absolute_max_obs_allowed_dict = absolute_max_obs_allowed_dict
            self.past_nights_observed_dict = past_nights_observed_dict

    def constraint_build_theta(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 0: Build theta variable")
        aframe_slots_for_request = self.request_set.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.request_set.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['Starname']==name][0]
            self.model.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.model.addConstr(self.theta[name] >= ((self.manager.requests_frame['# of Nights Per Semester'][idx] - \
                        self.past_nights_observed_dict[name]) - gp.quicksum(self.Yrds[name, d, s] for d,s in available)), \
                        'greater_than_nobs_shortfall_' + str(name))

    # Note: Likely can delete this constraint
    def constraint_build_theta_time_normalized(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 0: Build theta variable")
        aframe_slots_for_request = self.request_set.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.request_set.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['Starname']==name][0]
            self.model.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.model.addConstr(self.theta[name] >= ((self.manager.requests_frame['# of Nights Per Semester'][idx] - \
                        self.past_nights_observed_dict[name]) - gp.quicksum(self.Yrds[name, d, s] for d,s in available))*self.manager.slots_needed_for_exposure_dict[name], \
                        'greater_than_nobs_shortfall_' + str(name))

    # Note: Likely can delete this constraint
    def constraint_build_theta_program_normalized(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 0: Build theta variable")
        aframe_slots_for_request = self.request_set.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.request_set.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['Starname']==name][0]
            self.model.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.model.addConstr(self.theta[name] >= ((self.manager.requests_frame['# of Nights Per Semester'][idx] - \
                        self.past_nights_observed_dict[name]) - gp.quicksum(self.Yrds[name, d, s] for d,s in available))/self.manager.requests_frame['# of Nights Per Semester'][idx], \
                        'greater_than_nobs_shortfall_' + str(name))

    def constraint_one_request_per_slot(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 1: Enforce one request per slot.")
        self.request_set.Aset.rename(columns={'id': 'r'}, inplace=True) #temporary!
        for key in self.Yrds.keys():
            print(f"(name, d, s) = {key}")
        # Get all requests which are valid in slot (d, s)
        Aframe_ds = pd.merge(
            self.request_set.Aset.drop_duplicates(['d', 's']),
            self.request_set.Aset[['r', 'd', 's']],
            suffixes=['', '2'],
            on=['d', 's'])
        ds_requests = Aframe_ds.groupby(['d','s'])[['r2']].agg(list)
        Aset_ds_no_duplicates = self.request_set.Aset.copy()
        Aset_ds_no_duplicates = self.request_set.Aset.drop_duplicates(subset=['d', 's'])
        # Construct the constraint
        for i, row in Aset_ds_no_duplicates.iterrows():
            requests_valid_in_slot = list(ds_requests.loc[(row.d, row.s)])[0]
            self.model.addConstr((gp.quicksum(self.Yrds[name,row.d,row.s] for name in requests_valid_in_slot) <= 1),
                            'one_request_per_slot_' + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_reserve_multislot_exposures(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 2: Reserve slots for for multi-slot exposures.")
        # Get all requests that are valid in (d,s+e) pair for a given (d,s,1..e)
        print(self.request_set.Aframe)
        requires_multislot = self.request_set.Aframe.t_visit > 1
        # Note: earlier we had a full query and merge of an outer join via pandas but found that this quickly
        # ballooned in size in terms of memory required to complete the merge. This is an equal shortcut.
        frame_holder = []
        for day in range(self.manager.n_nights_in_semester):
            today_mask = self.request_set.Aset.d==day
            frame_holder.append(pd.merge(self.request_set.Aframe[today_mask&requires_multislot]['r d s e'.split()] \
                ,self.request_set.Aset['r d s'.split()],on=['d'],suffixes=['','2']) \
                .query('s < s2 < s + e').groupby('r d s'.split()).agg(list))
        requests_valid_in_reserved_slots = pd.concat(frame_holder)
        # If request requires only 1 slot to complete, then no constraint on reserving additional slots
        Aframe_multislots = self.request_set.Aset[requires_multislot]
        for i, row in Aframe_multislots.iterrows():
            # construct list of (r,d,s) indices to be constrained. These are all requests that are
            # valid in slots (d, s+1) through (d, s + e)
            # Ensuring slot (d, s) is not double filled already taken care of in Constraint 1.
            if (row.r, row.d, row.s) in requests_valid_in_reserved_slots.index:
                allr = list(requests_valid_in_reserved_slots.loc[row.r, row.d, row.s]['r2'])
                alls = list(requests_valid_in_reserved_slots.loc[row.r, row.d, row.s]['s2'])
                all_reserved_slots = list(zip(allr, [row.d]*len(allr), alls))
                self.model.addConstr((row.e*(1 - self.Yrds[row.r,row.d,row.s])) >= gp.quicksum(self.Yrds[c]
                                                                for c in all_reserved_slots),
                                'reserve_multislot_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_max_visits_per_night(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 3: Schedule request's maximum observations per night.")
        # Get all valid slots s for request r on day d
        aframe_slots_on_day = pd.merge(
            self.request_set.Aset.drop_duplicates(['r','d',]),
            self.request_set.Aset[['r','d','s']],
            suffixes=['','2'],on=['r']
        ).query('d == d2')
        slots_on_day = aframe_slots_on_day.groupby(['r','d'])[['s2']].agg(list)
        self.slots_on_day = slots_on_day
        unique_request_day_pairs = self.request_set.Aset.drop_duplicates(['r','d'])
        for i, row in unique_request_day_pairs.iterrows():
            constrained_slots_tonight = np.array(slots_on_day.loc[(row.r, row.d)][0])
            self.model.addConstr((gp.quicksum(self.Yrds[row.r,row.d,ss] for ss in constrained_slots_tonight) <= row.maxv),
                    'max_observations_per_night_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_enforce_internight_cadence(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 4: Enforce inter-night cadence.")
        aframe_slots_on_day = pd.merge(
            self.request_set.Aset.drop_duplicates(['r','d',]),
            self.request_set.Aset[['r','d','s']],
            suffixes=['','2'],on=['r']
        ).query('d == d2')
        slots_on_day = aframe_slots_on_day.groupby(['r','d'])[['s2']].agg(list)
        self.slots_on_day = slots_on_day
        aframe_intercadence = pd.merge(
            self.request_set.Aset.drop_duplicates(['r','d',]),
            self.request_set.Aset[['r','d','s']], #
            suffixes=['','2'],on=['r']
        ).query('d + 0 < d2 < d + ter') # When +1, it excludes the first day of the semester always
        intercadence = aframe_intercadence.groupby(['r','d'])[['d2','s2']].agg(list)
        # When inter-night cadence is 1, there will be no keys to constrain so skip
        # While the if/else statement would catch these, by shrinking the list here we do fewer
        # total steps in the loop.
        mask_inter_1 = self.request_set.Aframe['tau_inter'] > 1
        Aset_inter = self.request_set.Aframe[mask_inter_1]
        # We don't want duplicate slots on day d because we only need this constraint once per day
        # With duplicates, the same constraint would be applied to (r, d, s) and (r, d, s+1) which
        # is superfluous since we are summing over tonight's slots
        Aset_inter_no_duplicates = Aset_inter.copy()
        Aset_inter_no_duplicates = Aset_inter_no_duplicates.drop_duplicates(subset=['r', 'd'])
        for i, row in Aset_inter_no_duplicates.iterrows():
            constrained_slots_tonight = np.array(self.slots_on_day.loc[(row.r, row.d)][0])
            # Get all slots for pair (r, d) where valid
            if (row.r, row.d) in intercadence.index:
                slots_to_constrain_future = intercadence.loc[(row.r, row.d)]
                ds_pairs = list(zip(slots_to_constrain_future.d2, slots_to_constrain_future.s2))
                self.model.addConstr((gp.quicksum(self.Yrds[row.r,row.d,s2] for s2 in constrained_slots_tonight)/row.maxv \
                     <= (1 - (gp.quicksum(self.Yrds[row.r,d3,s3] for d3, s3 in ds_pairs)))), \
                    'enforce_internight_cadence_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")
            # else:
            #     # For request r, there are no (d,s) pairs that are within "inter cadence" days of the
            #     # given day d, therefore nothing to constrain. If I can find a way to filter out these
            #     # rows as a "mask_inter_2", then the if/else won't be needed
            #     continue

    def constraint_set_max_quarters_allocated(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: setting max number of quarters that can be allocated.")
        self.model.addConstr(gp.quicksum(self.Anq[d,q] for d in range(self.manager.n_nights_in_semester) \
                        for q in range(self.manager.n_quarters_in_night))
                        <= self.manager.max_quarters, "maximumQuartersAllocated")

    def constraint_set_max_onsky_allocated(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: setting max number of unique nights that can be allocated.")
        self.model.addConstr(gp.quicksum(self.Un[d] for d in range(self.manager.n_nights_in_semester))
                            <= self.manager.max_unique_nights, "maximumNightsAllocated")

    def constraint_relate_allocation_and_onsky(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: relating allocation map and unique night allocation map.")
        for d in range(self.manager.n_nights_in_semester):
            for q in range(self.manager.n_quarters_in_night):
                self.model.addConstr(self.Un[d] >= self.Anq[d,q], "relatedUnique_andNonUnique_lowerbound_" + str(d) + "d_" + str(q) + "q")
            self.model.addConstr(self.Un[d] <= gp.quicksum(self.Anq[d,q] for q in range(self.manager.n_quarters_in_night)), "relatedUnique_andNonUnique_upperbound_" + str(d) + "d")

    def constraint_all_portions_of_night_represented(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: setting minimum number of times each quarter to be allocated.")
        self.model.addConstr(gp.quicksum(self.Anq[d,0] for d in range(self.manager.n_nights_in_semester)) >= self.manager.min_represented, "minQuarterSelection_0q")
        self.model.addConstr(gp.quicksum(self.Anq[d,1] for d in range(self.manager.n_nights_in_semester)) >= self.manager.min_represented, "minQuarterSelection_1q")
        self.model.addConstr(gp.quicksum(self.Anq[d,2] for d in range(self.manager.n_nights_in_semester)) >= self.manager.min_represented, "minQuarterSelection_2q")
        self.model.addConstr(gp.quicksum(self.Anq[d,3] for d in range(self.manager.n_nights_in_semester)) >= self.manager.min_represented, "minQuarterSelection_3q")

    def constraint_forbidden_quarter_patterns(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: forbid certain patterns of quarter night allocations within night.")
        for d in range(self.manager.n_nights_in_semester):
            # Cannot have 1st and 3rd quarter allocated without also allocating 2nd quarter (no gap), regardless of if 4th quarter is allocated or not
            self.model.addConstr(self.Anq[d,0] + (self.Un[d]-self.Anq[d,1]) + self.Anq[d,2] <= 2*self.Un[d], "NoGap2_" + str(d) + "d")
            # Cannot have 2nd and 4th quarter allocated without also allocating 3rd quarter (no gap), regardless of if 1st quarter is allocated or not
            self.model.addConstr(self.Anq[d,1] + (self.Un[d]-self.Anq[d,2]) + self.Anq[d,3] <= 2*self.Un[d], "NoGap3_" + str(d) + "d")
            # Cannot have only 2nd and 3rd quarters allocated (no middle half)
            self.model.addConstr((self.Un[d]-self.Anq[d,0]) + self.Anq[d,1] + self.Anq[d,2] + (self.Un[d]-self.Anq[d,3]) <= 3*self.Un[d], "NoMiddleHalf_" + str(d) + "d")
            if self.manager.allow_single_quarters == False:
                # Cannot have only 1st and 4th quarters allocated (no end-cap half)
                self.model.addConstr(self.Anq[d,0] + (self.Un[d]-self.Anq[d,1]) + (self.Un[d]-self.Anq[d,2]) + self.Anq[d,3] <= 3*self.Un[d], "NoEndCapHalf_" + str(d) + "d")
                # Cannot choose single quarter allocations
                self.model.addConstr(self.Anq[d,0] + (self.Un[d]-self.Anq[d,1]) + (self.Un[d]-self.Anq[d,2]) + (self.Un[d]-self.Anq[d,3]) <= 3*self.Un[d], "No1stQOnly_" + str(d) + "d")
                self.model.addConstr((self.Un[d]-self.Anq[d,0]) + self.Anq[d,1] + (self.Un[d]-self.Anq[d,2]) + (self.Un[d]-self.Anq[d,3]) <= 3*self.Un[d], "No2ndQOnly_" + str(d) + "d")
                self.model.addConstr((self.Un[d]-self.Anq[d,0]) + (self.Un[d]-self.Anq[d,1]) + self.Anq[d,2] + (self.Un[d]-self.Anq[d,3]) <= 3*self.Un[d], "No3rdQOnly_" + str(d) + "d")
                self.model.addConstr((self.Un[d]-self.Anq[d,0]) + (self.Un[d]-self.Anq[d,1]) + (self.Un[d]-self.Anq[d,2]) + self.Anq[d,3] <= 3*self.Un[d], "No4thQOnly_" + str(d) + "d")
                # Cannot choose 3/4 allocations
                self.model.addConstr(self.Anq[d,0] + self.Anq[d,1] + self.Anq[d,2] + (self.Un[d]-self.Anq[d,3]) <= 3*self.Un[d], "No3/4Q_v1_" + str(d) + "d")
                self.model.addConstr((self.Un[d]-self.Anq[d,0]) + self.Anq[d,1] + self.Anq[d,2] + self.Anq[d,3] <= 3*self.Un[d], "No3/4Q_v2_" + str(d) + "d")

    def constraint_cannot_observe_if_not_allocated(self, twilight_map_remaining_2D):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: cannot observe if night/quarter is not allocated.")
        # if quarter is not allocated, all slots in quarter must be zero
        # note that the twilight times at the front and end of the night have to be respected
        for r, d, s in self.request_set.Aset:
            split_1st2nd, split_2nd3rd, split_3rd4th = ac.convert_slot_to_quarter(twilight_map_remaining_2D[d])
            if s < split_1st2nd:
                q = 0
            elif s >= split_1st2nd and s < split_2nd3rd:
                q = 1
            elif s >= split_2nd3rd and s < split_3rd4th:
                q = 2
            elif s >= split_3rd4th:
                q = 3
            else:
                q = 100
                print("Houston, we've had a problem.")
            self.model.addConstr(self.Yrds[r, d, s] <= self.Anq[d, q], "dontSched_ifNot_Allocated_"+ str(d) + "d_" + str(q) + "q_" + str(s) + "s_" + r)

    def constraint_max_consecutive_onsky(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: setting maximum number of consecutive unique nights allocated.")
        for d in range(self.manager.n_nights_in_semester - self.manager.max_consecutive):
            self.model.addConstr(gp.quicksum(self.Un[d + t] for t in range(self.manager.max_consecutive)) <= self.manager.max_consecutive - 1, "consecutiveNightsMax_" + str(d) + "d")

    def constraint_minimum_consecutive_offsky(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: setting minimum gap in days between allocated unique nights.")
        # Enforce at least one day allocated every X days (no large gaps)
        # Note you cannot run this when observatory/instrument has an extended shutdown
        for d in range(self.manager.n_nights_in_semester - self.manager.min_consecutive):
            self.model.addConstr(gp.quicksum(self.Un[d + t] for t in range(self.manager.min_consecutive)) >= 2, "noLargeGaps_" + str(d) + "d")

    def constraint_enforce_restricted_nights(self, limit):
        """
        According to Eq X in Lubin et al. 2025.

        Args:
            limit (int): either 0 or 1 to enforce blackout and whiteout, respectively
        """
        if limit not in [0, 1]:
            print("Limit must be an integer, either 0 (blackout) or 1 (whiteout).")
        elif limit == 0:
            filename = self.manager.blackout_file
        else:
            filename = self.manager.whiteout_file
        print("Constraint: enforcing quarters that cannot be chosen.")
        selections = pd.read_csv(filename)
        for s in range(len(selections)):
            night = self.manager.all_dates_dict[selections['Date'][s]]
            quarter = selections['Quarter'][s]
            self.model.addConstr(self.Anq[night,quarter] == limit, "enforced_" + str(limit) + "_" + str(night) + "d_" + str(quarter) + 'q')

    def constraint_maximize_baseline(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: maximize the baseline of unique nights allocated.")
        self.model.addConstr(gp.quicksum(self.Un[0 + t] for t in range(self.manager.max_baseline)) >= 1, "maxBase_early")
        self.model.addConstr(gp.quicksum(self.Un[self.manager.n_nights_in_semester - self.manager.max_baseline + t] for t in range(self.manager.max_baseline)) >= 1, "maxBase_late")

    def constraint_fix_previous_objective(self, epsilon=5):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: Fixing the previous solution's objective value.")
        self.model.addConstr(gp.quicksum(self.theta[name] for name in self.manager.requests_frame['Starname']) <= \
                        self.model.objval + epsilon)

    def set_objective_maximize_slots_used(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Objective: Maximize the number of slots used.")
        self.model.setObjective(gp.quicksum(self.manager.slots_needed_for_exposure_dict[r]*self.Yrds[r,d,s]
                            for r, d, s in self.request_set.Aset), GRB.MAXIMIZE)

    # This objective can likely be deleted.
    def set_objective_minimize_theta(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        self.model.setObjective(gp.quicksum(self.theta[name] for name in self.request_set.schedulable_requests), GRB.MINIMIZE)

    def set_objective_minimize_theta_time_normalized(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        self.model.setObjective(gp.quicksum(self.theta[name]*self.manager.slots_needed_for_exposure_dict[name] for name in self.request_set.schedulable_requests), GRB.MINIMIZE)

    # This objective can likely be deleted.
    def set_objective_minimize_theta_prog_norm(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        self.model.setObjective(gp.quicksum(self.theta[name]/self.manager.requests_frame.loc[self.manager.requests_frame['Starname'] == name, '# of Nights Per Semester'] for name in self.request_set.schedulable_requests), GRB.MINIMIZE)

    def constraint_connect_Wrd_and_Yrds(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint -1: Connect W and Y for all requests.")
        # Get all slots s that are valid for a given r and d
        grouped_s = self.request_set.Aframe.groupby(['r', 'd'])['s'].unique().reset_index()
        grouped_s.set_index(['r', 'd'], inplace=True)
        for i, row in self.request_set.Aframe.iterrows():
            all_valid_slots_tonight = list(grouped_s.loc[(row.r, row.d)]['s'])
            self.model.addConstr(gp.quicksum(self.Yrds[row.r,row.d,s3] for s3 in all_valid_slots_tonight) <= \
                row.maxv*self.Wrd[row.r, row.d],
                'connect_W_and_Y' + row.r + "_" + str(row.d) + "d")

    def constraint_build_theta_multivisit(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 0: Build theta variable")
        aframe_slots_for_request = self.request_set.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.request_set.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['Starname']==name][0]
            self.model.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            all_d = list(set(list(aframe_slots_for_request.loc[name].d)))
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.model.addConstr(self.theta[name] >= ((self.manager.requests_frame['# of Nights Per Semester'][idx] - \
                        self.past_nights_observed_dict[name]) - (gp.quicksum(self.Yrds[name, d, s] for d, s in available))/self.manager.requests_frame['Desired Visits per Night'][idx]), \
                        'greater_than_nobs_shortfall_' + str(name))

    def constraint_set_max_desired_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraining desired maximum observations.")
        aframe_slots_for_request = self.request_set.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.request_set.schedulable_requests:
            all_d = list(set(list(aframe_slots_for_request.loc[name].d)))
            self.model.addConstr(gp.quicksum(self.Wrd[name, d] for d in all_d) <=
                        self.desired_max_obs_allowed_dict[name],
                        'max_desired_unique_nights_for_request_' + str(name))

    def remove_constraint_set_max_desired_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Removing previous maximum observations.")
        for name in self.request_set.schedulable_requests:
            rm_const = self.model.getConstrByName("max_desired_unique_nights_for_request_" + str(name))
            self.model.remove(rm_const)

    def constraint_set_max_absolute_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraining absolute maximum observations.")
        aframe_slots_for_request = self.request_set.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.request_set.schedulable_requests:
            all_d = list(set(list(aframe_slots_for_request.loc[name].d)))
            self.model.addConstr(gp.quicksum(self.Wrd[name, d] for d in all_d) <=
                    self.absolute_max_obs_allowed_dict[name],
                    'max_absolute_unique_nights_for_request_' + str(name))

    # This constraint can likely be deleted.
    def constraint_set_max_desired_unique_nights_Yrds(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        aframe_slots_for_request = self.request_set.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.request_set.schedulable_requests:
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.model.addConstr(gp.quicksum(self.Yrds[name, d, s] for d,s in available) <=
                    self.desired_max_obs_allowed_dict[name],
                    'max_desired_unique_nights_for_request_' + str(name))

    # This constraint can likely be deleted.
    def constraint_set_max_absolute_unique_nights_Yrds(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        aframe_slots_for_request = self.request_set.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.request_set.schedulable_requests:
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.model.addConstr(gp.quicksum(self.Yrds[name, d, s] for d,s in available) <=
                    self.absolute_max_obs_allowed_dict[name],
                    'max_absolute_unique_nights_for_request_' + str(name))

    def constraint_build_enforce_intranight_cadence(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        # # When intra-night cadence is 0, there will be no keys to constrain so skip
        mask_intra_1 = self.request_set.Aframe['tra'] >= 1
        aframe_intra = self.request_set.Aframe[mask_intra_1]
        print("Constraint 6: Enforce intra-night cadence.")
        # get all combos of slots that must be constrained if given slot is scheduled
        aframe_intracadence = pd.merge(
            aframe_intra.drop_duplicates(['r','d','s']),
            aframe_intra[['r','d','s']],
            suffixes=['','2'],on=['r', 'd']
            ).query('s + 0 < s2 <= s + tra')
        intracadence = aframe_intracadence.groupby(['r','d','s'])[['s2']].agg(list)
        for i, row in aframe_intra.iterrows():
            if (row.r, row.d, row.s) in intracadence.index:
                # Get all slots tonight which are too soon after given slot for another visit
                slots_to_constrain_tonight_intra = list(intracadence.loc[(row.r, row.d, row.s)][0])
                self.model.addConstr((self.Yrds[row.r,row.d,row.s] <= (self.Wrd[row.r, row.d] - (gp.quicksum(self.Yrds[row.r,row.d,s3] \
                    for s3 in slots_to_constrain_tonight_intra)))), \
                    'enforce_intranight_cadence_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_set_min_max_visits_per_night(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 7: Allow minimum and maximum visits.")
        aframe_intra = self.request_set.Aframe.copy()
        aframe_intra_no_duplicates = aframe_intra.drop_duplicates(subset=['r', 'd'])
        grouped_s = self.request_set.Aframe.groupby(['r', 'd'])['s'].unique().reset_index()
        grouped_s.set_index(['r', 'd'], inplace=True)
        for i, row in aframe_intra_no_duplicates.iterrows():
            all_valid_slots_tonight = list(grouped_s.loc[(row.r, row.d)]['s'])
            self.model.addConstr((((gp.quicksum(self.Yrds[row.r,row.d,s3] for s3 in all_valid_slots_tonight)))) <= \
                row.maxv*self.Wrd[row.r, row.d], 'enforce_max_visits1_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")
            self.model.addConstr((((gp.quicksum(self.Yrds[row.r,row.d,s3] for s3 in all_valid_slots_tonight)))) >= \
                row.minv*self.Wrd[row.r, row.d], 'enforce_min_visits_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def optimize_model(self):
        print("Begin model solve.")

        self.model.params.TimeLimit = self.manager.solve_time_limit
        self.model.Params.OutputFlag = self.manager.gurobi_output
        # Allow stop at 5% gap to prevent from spending lots of time on marginally better solution
        self.model.params.MIPGap = self.manager.solve_max_gap
        # More aggressive presolve gives better solution in shorter time
        self.model.params.Presolve = 2
        self.model.update()
        self.model.optimize()

        if self.model.Status == GRB.INFEASIBLE:
            print('Model remains infeasible. Searching for invalid constraints')
            search = self.model.computeIIS()
            print("Printing bad constraints:")
            for c in self.model.getConstrs():
                if c.IISConstr:
                    print('%s' % c.ConstrName)
            for c in m.getGenConstrs():
                if c.IISGenConstr:
                    print('%s' % c.GenConstrName)
        else:
            print("Model Successfully Solved.")

    def run_model(self):
        self.round_info = 'Round1'
        self.build_model_round1()
        self.solve_model()
        self.serialize_results_csv()
        self.round_info = 'Round2'
        # np.savetxt(manager.output_directory + 'Round2_Requests.txt', [], delimiter=',', fmt="%s")
        if self.manager.run_round_two:
            self.build_model_round2()
            self.solve_model()
        self.serialize_results_csv()
        print("Scheduling complete, clear skies!")

    def build_model_round1(self):
        t1 = time.time()
        # self.model = cf.GorubiModel(manager, request_set)

        self.constraint_one_request_per_slot()
        self.constraint_reserve_multislot_exposures()
        self.constraint_enforce_internight_cadence()

        self.constraint_set_max_desired_unique_nights_Wrd()
        self.constraint_build_enforce_intranight_cadence()
        self.constraint_set_min_max_visits_per_night()
        self.constraint_build_theta_multivisit()

        if self.manager.run_optimal_allocation:
            self.constraint_set_max_quarters_allocated()
            self.constraint_set_max_onsky_allocated()
            self.constraint_relate_allocation_and_onsky()
            self.constraint_all_portions_of_night_represented()
            self.constraint_forbidden_quarter_patterns()
            self.constraint_cannot_observe_if_not_allocated(self.manager.twilight_map_remaining_2D)
            if os.path.exists(self.manager.blackout_file):
                self.constraint_enforce_restricted_nights(limit=0)
            if os.path.exists(self.manager.whiteout_file):
                self.constraint_enforce_restricted_nights(limit=1)
            if manager.include_aesthetic:
                self.constraint_max_consecutive_onsky()
                self.constraint_minimum_consecutive_offsky()
                self.constraint_maximize_baseline()

        self.set_objective_minimize_theta_time_normalized()
        print("Time to build constraints: ", np.round(time.time()-t1,3))

    def build_model_round2(self):
        t1 = time.time()
        self.remove_constraint_set_max_desired_unique_nights_Wrd()
        self.constraint_set_max_absolute_unique_nights_Wrd()
        self.constraint_fix_previous_objective()
        self.set_objective_maximize_slots_used()
        print("Time to build constraints: ", np.round(time.time()-t1,3))

    def solve_model(self):
        t1 = time.time()
        self.optimize_model()
        print("Time to finish solver: ", np.round(time.time()-t1,3))

    def serialize_results_csv(self):
        print("Building human readable schedule.")
        if self.manager.run_optimal_allocation:
            self.retrieve_ois_solution()
        self.human_read_available = io.write_available_human_readable(self.manager)
        self.human_read_schedule = io.write_stars_schedule_human_readable(self.human_read_available, self.model.Yrds, self.manager, self.round_info)
        io.build_fullness_report(self.human_read_schedule, self.manager, self.round_info)
        io.write_out_results(self.manager, self.model, self.round_info)

    def retrieve_ois_solution(self):
        print("Retrieving results of Optimal Instrument Allocation set of nights.")
        allocation_schedule_1d = []
        for v in self.model.Anq.values():
            if np.round(v.X,0) == 1:
                allocation_schedule_1d.append(1)
            else:
                allocation_schedule_1d.append(0)
        allocation_schedule = np.reshape(allocation_schedule_1d, (self.manager.n_nights_in_semester, self.manager.n_quarters_in_night))
        manager.allocation_map_2D_NQ = allocation_schedule
        weather_holder = np.zeros(np.shape(allocation_schedule))
        allocation_map_1D, allocation_map_2D, weathered_map = mp.build_allocation_map(self.manager, allocation_schedule, weather_holder)
        mp.convert_allocation_array_to_binary(self.manager)
        self.manager.allocation_all = allocation_map_1D
        self.manager.allocation_1D = allocation_map_1D
        self.manager.allocation_map_2D = allocation_map_2D
        self.manager.weathered_map = weathered_map
