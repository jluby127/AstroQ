"""
Module for building constraints

This module constructs constraints for the gurobi model.

Example usage:
    import constraints_functions as cf
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

from kpfcc import DATADIR
import kpfcc.access as ac

class GorubiModel(object):
    """A Gorubi Model object, from which we can add constraints and share variables easily."""

    def __init__(self, manager, Aset, Aframe, schedulable_requests, Wset):

        print("Defining model.")
        self.manager = manager
        self.Aset = Aset
        self.Aframe = Aframe
        self.schedulable_requests = schedulable_requests

        run_direct_multi_visit = True
        if run_direct_multi_visit:
            self.Wset = Wset
        else:
            self.Wset = []

        if manager.run_optimal_allocation:
            if self.manager.semester_grid == [] or self.manager.quarters_grid == []:
                print("WARNING: Missing necessary parameters to run Optimal Allocation.")
                print("---- A semester grid and a quarters grid is required to run optimal allocation.")
            if self.manager.max_quarters == 0 or self.manager.max_unique_nights == 0:
                print("WARNING: Missing necessary parameters to run Optimal Allocation.")
                print("---- Neither Max Quarters and Max Unique Nights parameters cannot be zero.")

        self.m = gp.Model('Semester_Scheduler')
        # Yrs is technically a 1D matrix indexed by tuples.
        # But in practice best think of it as a 2D square matrix of requests r and slots s, with gaps.
        # Slot s for request r will be 1 to indicate starting an exposure for that request in that slot
        self.Yrds = self.m.addVars(self.Aset, vtype = GRB.BINARY, name = 'Requests_Slots')
        # theta is the "shortfall" variable, continous in natural numbers.
        self.theta = self.m.addVars(self.schedulable_requests, name = 'Shortfall')

        if self.manager.run_optimal_allocation:
            # Anq is a 2D matrix of N_nights_in_semester by N_quarters_in_night
            # element will be 1 if that night/quarter is allocated to KPF and 0 otherwise
            self.Anq = self.m.addVars(self.manager.semester_grid, self.manager.quarters_grid, vtype = GRB.BINARY, name = 'Allocation')

            # Un is a 1D matrix of N_nights_in_semester
            # element will be 1 if at least one quarter in that night is allocated
            self.Un = self.m.addVars(self.manager.semester_grid, vtype = GRB.BINARY, name = 'On-Sky')
            self.m.update()

        if self.Wset != []:
            # Wrd is technically a 1D matrix indexed by tuples.
            # But in practice best think of it as a 2D square matrix of requests r and nights d, with gaps.
            # Night d for request r will be 1 to indicate at least one exposure is scheduled for this night.
            # Note that Wrd is only valid for requests r which have at least 2 visits requested in the night.
            self.Wrd = self.m.addVars(Wset, vtype = GRB.BINARY, name = 'Nightly_Observation_Tracker')

        desired_max_obs_allowed_dict = {}
        absolute_max_obs_allowed_dict = {}
        past_nights_observed_dict = {}
        for name in self.schedulable_requests:
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
        aframe_slots_for_request = self.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['Starname']==name][0]
            self.m.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.m.addConstr(self.theta[name] >= ((self.manager.requests_frame['# of Nights Per Semester'][idx] - \
                        self.past_nights_observed_dict[name]) - gp.quicksum(self.Yrds[name, d, s] for d,s in available)), \
                        'greater_than_nobs_shortfall_' + str(name))

    def constraint_build_theta_time_normalized(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 0: Build theta variable")
        aframe_slots_for_request = self.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['Starname']==name][0]
            self.m.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.m.addConstr(self.theta[name] >= ((self.manager.requests_frame['# of Nights Per Semester'][idx] - \
                        self.past_nights_observed_dict[name]) - gp.quicksum(self.Yrds[name, d, s] for d,s in available))*self.manager.slots_needed_for_exposure_dict[name], \
                        'greater_than_nobs_shortfall_' + str(name))

    def constraint_build_theta_program_normalized(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 0: Build theta variable")
        aframe_slots_for_request = self.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['Starname']==name][0]
            self.m.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.m.addConstr(self.theta[name] >= ((self.manager.requests_frame['# of Nights Per Semester'][idx] - \
                        self.past_nights_observed_dict[name]) - gp.quicksum(self.Yrds[name, d, s] for d,s in available))/self.manager.requests_frame['# of Nights Per Semester'][idx], \
                        'greater_than_nobs_shortfall_' + str(name))

    def constraint_one_request_per_slot(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 1: Enforce one request per slot.")
        # Get all requests which are valid in slot (d, s)
        Aframe_ds = pd.merge(
            self.Aframe.drop_duplicates(['d', 's']),
            self.Aframe[['r', 'd', 's']],
            suffixes=['', '2'],
            on=['d', 's'])
        ds_requests = Aframe_ds.groupby(['d','s'])[['r2']].agg(list)
        Aset_ds_no_duplicates = self.Aframe.copy()
        Aset_ds_no_duplicates = self.Aframe.drop_duplicates(subset=['d', 's'])
        # Construct the constraint
        for i, row in Aset_ds_no_duplicates.iterrows():
            requests_valid_in_slot = list(ds_requests.loc[(row.d, row.s)])[0]
            self.m.addConstr((gp.quicksum(self.Yrds[name,row.d,row.s] for name in requests_valid_in_slot) <= 1),
                            'one_request_per_slot_' + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_reserve_multislot_exposures(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 2: Reserve slots for for multi-slot exposures.")
        n_multi_slot_requests = np.sum(self.Aframe.e > 1)
        if n_multi_slot_requests > 0:
            # Get all requests that are  valid in (d,s+e) pair for a given (d,s,1..e)
            requests_valid_in_reserved_slots = pd.merge(self.Aframe.query('e > 1 ')['r d s e'.split()] \
                ,self.Aframe['r d s'.split()],on=['d'],suffixes=['','2']) \
                .query('s < s2 < s + e').groupby('r d s'.split()).agg(list)
            # If request requires only 1 slot to complete, then no constraint on reserving additional slots
            Aframe_multislots = self.Aframe[self.Aframe.e > 1]
            print('LENGTH OF AFRAME_MULTISLOTS: ', len(Aframe_multislots))
            for i, row in Aframe_multislots.iterrows():
                # construct list of (r,d,s) indices to be constrained. These are all requests that are
                # valid in slots (d, s+1) through (d, s + e)
                # Ensuring slot (d, s) is not double filled already taken care of in Constraint 1.
                if (row.r, row.d, row.s) in requests_valid_in_reserved_slots.index:
                    allr = list(requests_valid_in_reserved_slots.loc[row.r, row.d, row.s]['r2'])
                    alls = list(requests_valid_in_reserved_slots.loc[row.r, row.d, row.s]['s2'])
                    all_reserved_slots = list(zip(allr, [row.d]*len(allr), alls))
                    self.m.addConstr((row.e*(1 - self.Yrds[row.r,row.d,row.s])) >= gp.quicksum(self.Yrds[c]
                                                                    for c in all_reserved_slots),
                                    'reserve_multislot_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")
        else:
            print("No multi slot exposures!")

    def constraint_max_visits_per_night(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 3: Schedule request's maximum observations per night.")
        # Get all valid slots s for request r on day d
        aframe_slots_on_day = pd.merge(
            self.Aframe.drop_duplicates(['r','d',]),
            self.Aframe[['r','d','s']],
            suffixes=['','2'],on=['r']
        ).query('d == d2')
        slots_on_day = aframe_slots_on_day.groupby(['r','d'])[['s2']].agg(list)
        self.slots_on_day = slots_on_day
        unique_request_day_pairs = self.Aframe.drop_duplicates(['r','d'])
        # Build the constraint
        for i, row in unique_request_day_pairs.iterrows():
            constrained_slots_tonight = np.array(slots_on_day.loc[(row.r, row.d)][0])
            self.m.addConstr((gp.quicksum(self.Yrds[row.r,row.d,ss] for ss in constrained_slots_tonight) <= row.maxv),
                    'max_observations_per_night_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_enforce_internight_cadence(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 4: Enforce inter-night cadence.")
        aframe_intercadence = pd.merge(
            self.Aframe.drop_duplicates(['r','d',]),
            self.Aframe[['r','d','s']], #
            suffixes=['','2'],on=['r']
        ).query('d + 0 < d2 < d + ter') # When +1, it excludes the first day of the semester always
        intercadence = aframe_intercadence.groupby(['r','d'])[['d2','s2']].agg(list)
        # When inter-night cadence is 1, there will be no keys to constrain so skip
        # While the if/else statement would catch these, by shrinking the list here we do fewer
        # total steps in the loop.
        mask_inter_1 = self.Aframe['ter'] > 1
        Aset_inter = self.Aframe[mask_inter_1]
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
                self.m.addConstr((gp.quicksum(self.Yrds[row.r,row.d,s2] for s2 in constrained_slots_tonight)/row.v \
                     <= (1 - (gp.quicksum(self.Yrds[row.r,d3,s3] for d3, s3 in ds_pairs)))), \
                    'enforce_internight_cadence_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")
            else:
                # For request r, there are no (d,s) pairs that are within "inter cadence" days of the
                # given day d, therefore nothing to constrain. If I can find a way to filter out these
                # rows as a "mask_inter_2", then the if/else won't be needed
                continue

    def constraint_set_max_quarters_allocated(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: setting max number of quarters allocated.")
        # No more than a maximum number of quarters can be allocated
        self.m.addConstr(gp.quicksum(self.Anq[d,q] for d in range(self.manager.n_nights_in_semester) \
                        for q in range(self.manager.n_quarters_in_night))
                        <= self.manager.max_quarters, "maximumQuartersAllocated")

    def constraint_set_max_onsky_allocated(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: setting max number of unique nights allocated.")
        # No more than a maximum number of unique nights can be allocated
        self.m.addConstr(gp.quicksum(self.Un[d] for d in range(self.manager.n_nights_in_semester))
                            <= self.manager.max_unique_nights, "maximumNightsAllocated")

    def constraint_relate_allocation_and_onsky(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: relating allocation map and unique night allocation map.")
        # relate unique_allocation and allocation
        # if any one of the q in allocation[given date, q] is 1, then unique_allocation[given date] must be 1, zero otherwise
        for d in range(self.manager.n_nights_in_semester):
            for q in range(self.manager.n_quarters_in_night):
                self.m.addConstr(self.Un[d] >= self.Anq[d,q], "relatedUnique_andNonUnique_lowerbound_" + str(d) + "d_" + str(q) + "q")
            self.m.addConstr(self.Un[d] <= gp.quicksum(self.Anq[d,q] for q in range(self.manager.n_quarters_in_night)), "relatedUnique_andNonUnique_upperbound_" + str(d) + "d")

    def constraint_all_portions_of_night_represented(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: setting min number each quarter to be allocated.")
        # Minimum number of each quarter must be allocated
        self.m.addConstr(gp.quicksum(self.Anq[d,0] for d in range(self.manager.n_nights_in_semester)) >= self.manager.min_represented, "minQuarterSelection_0q")
        self.m.addConstr(gp.quicksum(self.Anq[d,1] for d in range(self.manager.n_nights_in_semester)) >= self.manager.min_represented, "minQuarterSelection_1q")
        self.m.addConstr(gp.quicksum(self.Anq[d,2] for d in range(self.manager.n_nights_in_semester)) >= self.manager.min_represented, "minQuarterSelection_2q")
        self.m.addConstr(gp.quicksum(self.Anq[d,3] for d in range(self.manager.n_nights_in_semester)) >= self.manager.min_represented, "minQuarterSelection_3q")

    def constraint_forbidden_quarter_patterns(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: forbid certain patterns of quarter night allocations within night.")
        # Disallow certain patterns of quarters selected within same night
        for d in range(self.manager.n_nights_in_semester):
            # Cannot have 1st and 3rd quarter allocated without also allocating 2nd quarter (no gap), regardless of if 4th quarter is allocated or not
            self.m.addConstr(self.Anq[d,0] + (self.Un[d]-self.Anq[d,1]) + self.Anq[d,2] <= 2*self.Un[d], "NoGap2_" + str(d) + "d")
            # Cannot have 2nd and 4th quarter allocated without also allocating 3rd quarter (no gap), regardless of if 1st quarter is allocated or not
            self.m.addConstr(self.Anq[d,1] + (self.Un[d]-self.Anq[d,2]) + self.Anq[d,3] <= 2*self.Un[d], "NoGap3_" + str(d) + "d")
            # Cannot have only 2nd and 3rd quarters allocated (no middle half)
            self.m.addConstr((self.Un[d]-self.Anq[d,0]) + self.Anq[d,1] + self.Anq[d,2] + (self.Un[d]-self.Anq[d,3]) <= 3*self.Un[d], "NoMiddleHalf_" + str(d) + "d")
            if self.manager.allow_single_quarters == False:
                # Cannot have only 1st and 4th quarters allocated (no end-cap half)
                self.m.addConstr(self.Anq[d,0] + (self.Un[d]-self.Anq[d,1]) + (self.Un[d]-self.Anq[d,2]) + self.Anq[d,3] <= 3*self.Un[d], "NoEndCapHalf_" + str(d) + "d")
                # Cannot choose single quarter allocations
                self.m.addConstr(self.Anq[d,0] + (self.Un[d]-self.Anq[d,1]) + (self.Un[d]-self.Anq[d,2]) + (self.Un[d]-self.Anq[d,3]) <= 3*self.Un[d], "No1stQOnly_" + str(d) + "d")
                self.m.addConstr((self.Un[d]-self.Anq[d,0]) + self.Anq[d,1] + (self.Un[d]-self.Anq[d,2]) + (self.Un[d]-self.Anq[d,3]) <= 3*self.Un[d], "No2ndQOnly_" + str(d) + "d")
                self.m.addConstr((self.Un[d]-self.Anq[d,0]) + (self.Un[d]-self.Anq[d,1]) + self.Anq[d,2] + (self.Un[d]-self.Anq[d,3]) <= 3*self.Un[d], "No3rdQOnly_" + str(d) + "d")
                self.m.addConstr((self.Un[d]-self.Anq[d,0]) + (self.Un[d]-self.Anq[d,1]) + (self.Un[d]-self.Anq[d,2]) + self.Anq[d,3] <= 3*self.Un[d], "No4thQOnly_" + str(d) + "d")
                # Cannot choose 3/4 allocations
                self.m.addConstr(self.Anq[d,0] + self.Anq[d,1] + self.Anq[d,2] + (self.Un[d]-self.Anq[d,3]) <= 3*self.Un[d], "No3/4Q_v1_" + str(d) + "d")
                self.m.addConstr((self.Un[d]-self.Anq[d,0]) + self.Anq[d,1] + self.Anq[d,2] + self.Anq[d,3] <= 3*self.Un[d], "No3/4Q_v2_" + str(d) + "d")

    def constraint_cannot_observe_if_not_allocated(self, twilight_map_remaining_2D):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: cannot observe if night/quarter is not allocated.")
        # if quarter is not allocated, all slots in quarter must be zero
        # note that the twilight times at the front and end of the night have to be respected
        for r, d, s in self.Aset:
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
            self.m.addConstr(self.Yrds[r, d, s] <= self.Anq[d, q], "dontSched_ifNot_Allocated_"+ str(d) + "d_" + str(q) + "q_" + str(s) + "s_" + r)

    def constraint_max_consecutive_onsky(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: setting max number of consecutive unique nights allocated.")
        # Don't allocate more than X consecutive nights
        for d in range(self.manager.n_nights_in_semester - self.manager.max_consecutive):
            self.m.addConstr(gp.quicksum(self.Un[d + t] for t in range(self.manager.max_consecutive)) <= self.manager.max_consecutive - 1, "consecutiveNightsMax_" + str(d) + "d")

    def constraint_minimum_consecutive_offsky(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: setting min gap in days between allocated unique nights.")
        # Enforce at least one day allocated every X days (no large gaps)
        # Note you cannot run this when observatory/instrument has an extended shutdown
        for d in range(self.manager.n_nights_in_semester - self.manager.min_consecutive):
            self.m.addConstr(gp.quicksum(self.Un[d + t] for t in range(self.manager.min_consecutive)) >= 2, "noLargeGaps_" + str(d) + "d")

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
            self.m.addConstr(self.Anq[night,quarter] == limit, "enforced_" + str(limit) + "_" + str(night) + "d_" + str(quarter) + 'q')

    def constraint_maximize_baseline(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: maximize the baseline of unique nights allocated.")
        # Enforce a night to be allocated within the first X nights and the last X nights of the semester (max baseline)
        self.m.addConstr(gp.quicksum(self.Un[0 + t] for t in range(self.manager.max_baseline)) >= 1, "maxBase_early")
        self.m.addConstr(gp.quicksum(self.Un[self.manager.n_nights_in_semester - self.manager.max_baseline + t] for t in range(self.manager.max_baseline)) >= 1, "maxBase_late")

    def constraint_fix_previous_objective(self, epsilon=5):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint: Fixing the previous solution's objective value.")
        self.m.addConstr(gp.quicksum(self.theta[name] for name in self.manager.requests_frame['Starname']) <= \
                        self.m.objval + epsilon)

    def set_objective_maximize_slots_used(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Objective: Maximize the number of slots used.")
        self.m.setObjective(gp.quicksum(self.manager.slots_needed_for_exposure_dict[r]*self.Yrds[r,d,s]
                            for r, d, s in self.Aset), GRB.MAXIMIZE)

    def set_objective_minimize_theta(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        self.m.setObjective(gp.quicksum(self.theta[name] for name in self.schedulable_requests), GRB.MINIMIZE)

    def set_objective_minimize_theta_time_norm(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        self.m.setObjective(gp.quicksum(self.theta[name]*self.manager.slots_needed_for_exposure_dict[name] for name in self.schedulable_requests), GRB.MINIMIZE)

    def set_objective_minimize_theta_prog_norm(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        self.m.setObjective(gp.quicksum(self.theta[name]/self.manager.requests_frame.loc[self.manager.requests_frame['Starname'] == name, '# of Nights Per Semester'] for name in self.schedulable_requests), GRB.MINIMIZE)

    def solve_model(self):
        print("Begin model solve.")

        self.m.params.TimeLimit = self.manager.solve_time_limit
        self.m.Params.OutputFlag = self.manager.gurobi_output
        # Allow stop at 5% gap to prevent from spending lots of time on marginally better solution
        self.m.params.MIPGap = self.manager.solve_max_gap
        # More aggressive presolve gives better solution in shorter time
        self.m.params.Presolve = 2
        self.m.update()
        self.m.optimize()

        if self.m.Status == GRB.INFEASIBLE:
            print('Model remains infeasible. Searching for invalid constraints')
            search = self.m.computeIIS()
            print("Printing bad constraints:")
            for c in self.m.getConstrs():
                if c.IISConstr:
                    print('%s' % c.ConstrName)
            for c in m.getGenConstrs():
                if c.IISGenConstr:
                    print('%s' % c.GenConstrName)
        else:
            print("Model Successfully Solved.")

    def constraint_connect_Wrd_and_Yrds(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint -1: Connect W and Y for all requests.")
        # Get all slots s that are valid for a given r and d
        grouped_s = self.Aframe.groupby(['r', 'd'])['s'].unique().reset_index()
        grouped_s.set_index(['r', 'd'], inplace=True)
        for i, row in self.Aframe.iterrows():
            all_valid_slots_tonight = list(grouped_s.loc[(row.r, row.d)]['s'])
            self.m.addConstr(gp.quicksum(self.Yrds[row.r,row.d,s3] for s3 in all_valid_slots_tonight) <= \
                row.maxv*self.Wrd[row.r, row.d],
                'connect_W_and_Y' + row.r + "_" + str(row.d) + "d")

    def constraint_build_theta_direct_multivisit_v1(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 0: Build theta variable")
        aframe_slots_for_request = self.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['Starname']==name][0]

            self.m.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            all_d = list(set(list(aframe_slots_for_request.loc[name].d)))
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.m.addConstr(self.theta[name] >= ((self.manager.requests_frame['# of Nights Per Semester'][idx] - \
                        self.past_nights_observed_dict[name]) - (gp.quicksum(self.Wrd[name, d] for d in all_d))), \
                        'greater_than_nobs_shortfall_' + str(name))

    def constraint_build_theta_direct_multivisit_v2(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 0: Build theta variable")
        aframe_slots_for_request = self.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['Starname']==name][0]
            self.m.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            all_d = list(set(list(aframe_slots_for_request.loc[name].d)))
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.m.addConstr(self.theta[name] >= ((self.manager.requests_frame['# of Nights Per Semester'][idx] - \
                        self.past_nights_observed_dict[name]) - (gp.quicksum(self.Yrds[name, d, s] for d, s in available))/self.manager.requests_frame['Desired Visits per Night'][idx]), \
                        'greater_than_nobs_shortfall_' + str(name))

    def constraint_set_max_desired_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraining desired maximum observations.")
        aframe_slots_for_request = self.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.schedulable_requests:
            all_d = list(set(list(aframe_slots_for_request.loc[name].d)))
            self.m.addConstr(gp.quicksum(self.Wrd[name, d] for d in all_d) <=
                        self.desired_max_obs_allowed_dict[name],
                        'max_desired_unique_nights_for_request_' + str(name))

    def remove_constraint_set_max_desired_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Removing previous maximum observations.")
        for name in self.schedulable_requests:
            rm_const = self.m.getConstrByName("max_desired_unique_nights_for_request_" + str(name))
            self.m.remove(rm_const)

    def constraint_set_max_absolute_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraining absolute maximum observations.")
        aframe_slots_for_request = self.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.schedulable_requests:
            all_d = list(set(list(aframe_slots_for_request.loc[name].d)))
            self.m.addConstr(gp.quicksum(self.Wrd[name, d] for d in all_d) <=
                    self.absolute_max_obs_allowed_dict[name],
                    'max_absolute_unique_nights_for_request_' + str(name))

    def constraint_set_max_desired_unique_nights_Yrds(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        aframe_slots_for_request = self.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.schedulable_requests:
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.m.addConstr(gp.quicksum(self.Yrds[name, d, s] for d,s in available) <=
                    self.desired_max_obs_allowed_dict[name],
                    'max_desired_unique_nights_for_request_' + str(name))

    def constraint_set_max_absolute_unique_nights_Yrds(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        aframe_slots_for_request = self.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.schedulable_requests:
            available = list(zip(list(aframe_slots_for_request.loc[name].d), list(aframe_slots_for_request.loc[name].s)))
            self.m.addConstr(gp.quicksum(self.Yrds[name, d, s] for d,s in available) <=
                    self.absolute_max_obs_allowed_dict[name],
                    'max_absolute_unique_nights_for_request_' + str(name))

    def constraint_build_enforce_intranight_cadence(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        # # When intra-night cadence is 0, there will be no keys to constrain so skip
        mask_intra_1 = self.Aframe['tra'] >= 1
        aframe_intra = self.Aframe[mask_intra_1]
        print("Constraint 6: Enforce intra-night cadence.")
        # get all combos of slots that must be constrained if given slot is scheduled
        aframe_intracadence = pd.merge(
            aframe_intra.drop_duplicates(['r','d','s']),
            aframe_intra[['r','d','s']],
            suffixes=['','2'],on=['r', 'd']
            ).query('s + 0 < s2 <= s + tra')
        intracadence = aframe_intracadence.groupby(['r','d','s'])[['s2']].agg(list)
        count8 = 0
        for i, row in aframe_intra.iterrows():
            if (row.r, row.d, row.s) in intracadence.index:
                # Get all slots tonight which are too soon after given slot for another visit
                slots_to_constrain_tonight_intra = list(intracadence.loc[(row.r, row.d, row.s)][0])
                self.m.addConstr((self.Yrds[row.r,row.d,row.s] <= (self.Wrd[row.r, row.d] - (gp.quicksum(self.Yrds[row.r,row.d,s3] \
                    for s3 in slots_to_constrain_tonight_intra)))), \
                    'enforce_intranight_cadence_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_set_min_max_visits_per_night(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 7: Allow minimum and maximum visits.")
        grouped_s = self.Aframe.groupby(['r', 'd'])['s'].unique().reset_index()
        grouped_s.set_index(['r', 'd'], inplace=True)
        mask_intra_1 = self.Aframe['tra'] >= 1
        aframe_intra = self.Aframe[mask_intra_1]
        aframe_intra_no_duplicates = aframe_intra.drop_duplicates(subset=['r', 'd'])
        aframe_intracadence = pd.merge(
            aframe_intra.drop_duplicates(['r','d','s']),
            aframe_intra[['r','d','s']],
            suffixes=['','2'],on=['r', 'd']
            ).query('s + 0 < s2 <= s + tra')
        intracadence = aframe_intracadence.groupby(['r','d','s'])[['s2']].agg(list)
        for i, row in aframe_intra_no_duplicates.iterrows():
            if (row.r, row.d, row.s) in intracadence.index:
                all_valid_slots_tonight = list(grouped_s.loc[(row.r, row.d)]['s'])
                self.m.addConstr((((gp.quicksum(self.Yrds[row.r,row.d,s3] for s3 in all_valid_slots_tonight)))) <= \
                    row.maxv*self.Wrd[row.r, row.d], 'enforce_max_visits1_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")
                self.m.addConstr((((gp.quicksum(self.Yrds[row.r,row.d,s3] for s3 in all_valid_slots_tonight)))) >= \
                    row.minv*self.Wrd[row.r, row.d], 'enforce_min_visits_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")


def define_slot_index_frame(manager, available_indices_for_request):
    """
    Using the dictionary of indices where each request is available, define a dataframe for which
    we will use to cut/filter/merge r,d,s tuples

    Args:
        requests_frame (dataframe): the pandas dataframe containing the request information
        slots_needed_for_exposure_dict (dictionary): keys are the request names, values are the
                        number of slots required for one exposure of that request
        available_indices_for_request (dictionary): keys are the starnames and values are a 1D array
                                                  the indices where available_slots_for_request is 1.
    Returns:
        Aframe (dataframe): a pandas dataframe containing one row for each valid r,d,s tuple
        Aset (array): an array of [r, d, s] tuples that make up the indices of Yrds
        schedulable_requests (list): the requests that are valid for this schedule

    """
    # Define the tuples of request and available slot for each request.
    # This becomes the grid over which the Gurobi variables are defined.
    # Now, slots that were never possible for scheduling are not included in the model.
    Aset = []
    Wset = []
    Aframe_keys = []
    for n,row in manager.requests_frame.iterrows():
        name = row['Starname']
        max_n_visits = int(row['Desired Visits per Night'])
        min_n_visits = int(row['Accepted Visits per Night'])
        intra = int(row['Minimum Intra-Night Cadence'])
        inter = int(row['Minimum Inter-Night Cadence'])
        slots_needed = manager.slots_needed_for_exposure_dict[name]
        for d in range(len(available_indices_for_request[name])):
            Wset.append((name, d))
            for s in available_indices_for_request[name][d]:
                Aset.append((name, d, s))
                Aframe_keys.append([name, d, s, slots_needed, max_n_visits, min_n_visits, intra, inter])

    Aframe = pd.DataFrame(Aframe_keys, columns =['r', 'd', 's', 'e', 'maxv', 'minv', 'tra', 'ter'])
    schedulable_requests = list(Aframe['r'].unique())
    for name in list(manager.requests_frame['Starname']):
        if name not in schedulable_requests:
            print("WARNING: Target " + name + " has no valid day/slot pairs and therefore is effectively removed from the model.")
    # duplicate columns for easy indexing later
    Aframe['rr'] = Aframe['r']
    Aframe['dd'] = Aframe['d']
    Aframe['ss'] = Aframe['s']

    return Aframe, Aset, schedulable_requests, Wset
