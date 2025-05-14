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

import kpfcc.io as io
import kpfcc.management as mn
import kpfcc.request as rq
import kpfcc.maps as mp

class Scheduler(object):
    """A Scheduler object, from which we can define a Gurobi model, build constraints, and solve."""

    def __init__(self, request_set, cf):
        print("Building the Scheduler.")
        self.start_the_clock = time.time()

        manager = mn.data_admin(cf)
        if manager.current_day != request_set.meta['d1_date']:
            print("Mismatch on 'current_date' between config file and request_set file! Using the request_set file's value.")
            manager.current_day = request_set.meta['d1_date']
        manager.run_admin()

        self.manager = manager
        request_set = rq.cull_from_weather(request_set, self.manager.weathered_days)
        self.request_set = request_set
        self.observability_tuples = list(request_set.observability.itertuples(index=False, name=None))

        self.joiner = pd.merge(self.request_set.strategy, self.request_set.observability, on=['id'])
        # add dummy columns for easier joins
        self.joiner['id2'] = self.joiner['id']
        self.joiner['d2'] = self.joiner['d']
        self.joiner['s2'] = self.joiner['s']
        self.joiner['tau_intra'] *= int(60/self.manager.slot_size) # convert hours to slots
        self.joiner['tau_intra'] += self.joiner['t_visit'] # start the minimum intracadence time from the end of the previous exposure, not the beginning

        # Prepare information by construction observability_nights (Wset) and schedulable_requests
        self.observability_nights = self.joiner[self.joiner['n_intra_max'] > 1][['id', 'd']].drop_duplicates().copy()
        self.multi_visit_requests = list(self.observability_nights['id'].unique())

        self.all_requests = list(manager.requests_frame['starname'])
        self.schedulable_requests =  list(self.joiner['id'].unique())
        self.single_visit_requests = [item for item in self.schedulable_requests if item not in self.multi_visit_requests]
        for name in list(manager.requests_frame['starname']):
            if name not in self.schedulable_requests:
                print("WARNING: Target " + name + " has no valid day/slot pairs and therefore is effectively removed from the model.")

        # Construct a few useful joins
        # Get each request's full list of valid d/s pairs
        self.all_valid_ds_for_request = self.joiner.groupby(['id'])[['d', 's']].agg(list)
        # Get all requests which are valid in slot (d, s)
        requests_valid_for_ds = pd.merge(
            self.joiner.drop_duplicates(['d', 's']),
            self.joiner[['id', 'd', 's']],
            suffixes=['', '3'],
            on=['d', 's'])
        self.requests_valid_for_ds = requests_valid_for_ds.groupby(['d','s'])[['id3']].agg(list)

        # Get all valid d/s pairs
        self.valid_ds_pairs = self.joiner.copy()
        # self.valid_ds_pairs = self.joiner.drop_duplicates(subset=['d', 's'])
        # Get all requests that require multiple slots to complete one observation
        self.multislot_mask = self.joiner.t_visit > 1
        self.multi_slot_frame = self.joiner[self.multislot_mask]
        # Get all valid slots s for request r on day d
        valid_s_for_rd = pd.merge(
            self.joiner.drop_duplicates(['id','d',]),
            self.joiner[['id','d','s']],
            suffixes=['','3'],on=['id']
        ).query('d == d3')
        self.slots_on_day_for_r = valid_s_for_rd.groupby(['id','d'])[['s3']].agg(list)
        # Get all request id's that are valid on a given day
        self.unique_request_on_day_pairs = self.joiner.copy().drop_duplicates(['id','d'])

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
        observability_array = list(self.request_set.observability.itertuples(index=False, name=None))
        self.Yrds = self.model.addVars(observability_array, vtype = GRB.BINARY, name = 'Requests_Slots')

        if len(self.observability_nights) != 0:
            # Wrd is technically a 1D matrix indexed by tuples.
            # But in practice best think of it as a 2D square matrix of requests r and nights d, with gaps.
            # Night d for request r will be 1 to indicate at least one exposure is scheduled for this night.
            # Note that Wrd is only valid for requests r which have at least 2 visits requested in the night.
            observability_array_onsky = list(self.observability_nights.itertuples(index=False, name=None))
            self.Wrd = self.model.addVars(observability_array_onsky, vtype = GRB.BINARY, name = 'OnSky')

        # theta is the "shortfall" variable, continous in natural numbers.
        self.theta = self.model.addVars(self.all_requests, name = 'Shortfall', vtype = GRB.INTEGER)

        if self.manager.run_optimal_allocation:
            # Anq is a 2D matrix of N_nights_in_semester by N_quarters_in_night
            # element will be 1 if that night/quarter is allocated to KPF and 0 otherwise
            self.Anq = self.model.addVars(self.manager.semester_grid, self.manager.quarters_grid, vtype = GRB.BINARY, name = 'Allocation')

            # Un is a 1D matrix of N_nights_in_semester
            # element will be 1 if at least one quarter in that night is allocated
            self.Un = self.model.addVars(self.manager.semester_grid, vtype = GRB.BINARY, name = 'On-Sky')
            self.model.update()

        desired_max_obs_allowed_dict = {}
        absolute_max_obs_allowed_dict = {}
        past_nights_observed_dict = {}
        for name in self.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['starname']==name][0]
            if self.manager.database_info_dict == {}:
                past_nights_observed = 0
            else:
                past_nights_observed = len(self.manager.database_info_dict[name][1])

            # Safety valve for if the target is over-observed for any reason
            if past_nights_observed > self.manager.requests_frame['n_inter_max'][idx] + \
                        int(self.manager.requests_frame['n_inter_max'][idx]*self.manager.max_bonus):
                desired_max_obs = past_nights_observed
            else:
                desired_max_obs = (self.manager.requests_frame['n_inter_max'][idx] - past_nights_observed)
                absolute_max_obs = (self.manager.requests_frame['n_inter_max'][idx] - past_nights_observed) \
                        + int(self.manager.requests_frame['n_inter_max'][idx]*self.manager.max_bonus)
                # second safety valve
                if past_nights_observed > absolute_max_obs:
                    absolute_max_obs = past_nights_observed
            past_nights_observed_dict[name] = past_nights_observed
            desired_max_obs_allowed_dict[name] = desired_max_obs
            absolute_max_obs_allowed_dict[name] = absolute_max_obs
            self.desired_max_obs_allowed_dict = desired_max_obs_allowed_dict
            self.absolute_max_obs_allowed_dict = absolute_max_obs_allowed_dict
            self.past_nights_observed_dict = past_nights_observed_dict
        print("Initializing complete.")

    def constraint_one_request_per_slot(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 1: Enforce one request per slot.")
        # Construct the constraint
        tmpjoiner = self.joiner.drop_duplicates(subset=['d', 's'])
        for i, row in tmpjoiner.iterrows():
            requests_valid_in_slot = list(self.requests_valid_for_ds.loc[(row.d, row.s)])[0]
            self.model.addConstr((gp.quicksum(self.Yrds[name,row.d,row.s] for name in requests_valid_in_slot) <= 1),
                            'one_request_per_slot_' + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_reserve_multislot_exposures(self):
        """
        According to Eq X in Lubin et al. 2025.

        Get requests that are valid in (d,s+t_visit) pair for a given (d,s,1...t_visit)
        Note: earlier we had a full query and merge of an outer join via pandas but found that this quickly
        # ballooned in size in terms of memory required to complete the merge. This is an equal shortcut.

        # The Yr,d,s matrix represents the day/slot pairs where an exposure r will begin. Therefore, when a re-
        # quest's exposure time exceeds the length of time repre-
        # sented by a single slot, then multiple consecutive slots
        # must be reserved for the request. In this case, no re-
        # quest, including the one being scheduled, may be sched-
        # uled into the required consecutive slots
        """
        print("Constraint 2: Reserve slots for for multi-slot exposures.")
        obs = self.request_set.observability
        ms = self.multi_slot_frame
        # If request requires only 1 slot to complete, then no constraint on reserving additional slots

        for d in ms.d.drop_duplicates():
            obs_day = obs[obs.d == d]
            ms_day = ms[ms.d == d]

            # Precompute all requests that are valid on day d, and slot s.
            valid_requests = obs_day.groupby('s')['id'].apply(set).to_dict()

            # Loop over all multi-slot requests on day d and enforce no other request in
            # slots within t_visit from slot starting slot
            ids = ms_day.id.values
            ss = ms_day.s.values
            t_visits = ms_day.t_visit.values
            for i in range(len(ids)):
                id = ids[i]
                s = ss[i]
                t_visit = t_visits[i]
                rhs = []
                for delta in range(1,t_visit):
                    s_shift = s + delta
                    if s_shift in valid_requests:
                        lhs = 1 - self.Yrds[id,d,s]
                        rhs = gp.quicksum(self.Yrds[r,d,s_shift] for r in valid_requests[s_shift])
                        name = ""
                        constr = (lhs >= rhs)
                        self.model.addConstr(constr, name)

    # this function can likely be deleted - Jack 4/28/25
    def constraint_max_visits_per_night(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 3: Schedule request's maximum observations per night.")
        for i, row in self.unique_request_on_day_pairs.iterrows():
            constrained_slots_tonight = np.array(self.slots_on_day_for_r.loc[(row.id, row.d)][0])
            print(row.id, row.n_intra_max)
            self.model.addConstr((gp.quicksum(self.Yrds[row.id,row.d,ss] for ss in constrained_slots_tonight) <= row.n_intra_max),
                    'max_observations_per_night_' + row.id + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_enforce_internight_cadence(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 4: Enforce inter-night cadence.")
        # Get all (d',s') pairs for a request that must be zero if a (d,s) pair is selected
        intercadence = pd.merge(
            self.joiner.drop_duplicates(['id','d',]),
            self.joiner[['id','d','s']],
            suffixes=['','3'],on=['id']
        ).query('d + 0 < d3 < d + tau_inter')
        self.intercadence_tracker = intercadence.groupby(['id','d'])[['d3','s3']].agg(list)
        # When inter-night cadence is 1, there will be no keys to constrain so skip
        # While the if/else statement would catch these, by shrinking the list here we do fewer
        # total steps in the loop.
        intercadence_valid_tuples = self.joiner.copy()[self.joiner['tau_inter'] > 1]
        # We don't want duplicate slots on day d because we only need this constraint once per day
        # With duplicates, the same constraint would be applied to (r, d, s) and (r, d, s+1) which
        # is superfluous since we are summing over tonight's slots
        intercadence_valid_tuples = intercadence_valid_tuples.drop_duplicates(subset=['id', 'd'])
        for i, row in intercadence_valid_tuples.iterrows():
            constrained_slots_tonight = np.array(self.slots_on_day_for_r.loc[(row.id2, row.d2)][0])
            # Get all slots for pair (r, d) where valid
            # print(row.id, row.d)
            if (row.id, row.d) in self.intercadence_tracker.index:
                slots_to_constrain_future = self.intercadence_tracker.loc[(row.id2, row.d2)]
                ds_pairs = zip(list(np.array(slots_to_constrain_future.d3).flatten()), list(np.array(slots_to_constrain_future.s3).flatten()))
                self.model.addConstr((gp.quicksum(self.Yrds[row.id,row.d,s2] for s2 in constrained_slots_tonight)/row.n_intra_max \
                     <= (1 - (gp.quicksum(self.Yrds[row.id,d3,s3] for d3, s3 in ds_pairs)))), \
                    'enforce_internight_cadence_' + row.id + "_" + str(row.d) + "d_" + str(row.s) + "s")
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
        for id, d, s in zip(self.request_set.observability['id'], self.request_set.observability['d'], self.request_set.observability['s']):
            q = rq.convert_slot_to_quarter(d, s, twilight_map_remaining_2D[d])
            self.model.addConstr(self.Yrds[id, d, s] <= self.Anq[d, q], "dontSched_ifNot_Allocated_"+ str(d) + "d_" + str(q) + "q_" + str(s) + "s_" + id, d, id)

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
        self.model.addConstr(gp.quicksum(self.theta[name] for name in self.manager.requests_frame['starname']) <= \
                        self.model.objval + epsilon)

    def set_objective_maximize_slots_used(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Objective: Maximize the number of slots used.")
        self.model.setObjective(gp.quicksum(self.manager.slots_needed_for_exposure_dict[id]*self.Yrds[id,d,s]
                            for id, d, s in self.observability_tuples), GRB.MAXIMIZE)
                            # for id, d, s in self.joiner), GRB.MAXIMIZE)

    def set_objective_minimize_theta_time_normalized(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        self.model.setObjective(gp.quicksum(self.theta[name]*self.manager.slots_needed_for_exposure_dict[name] for name in self.schedulable_requests), GRB.MINIMIZE)

    def constraint_build_theta_multivisit(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 0: Build theta variable")
        for name in self.schedulable_requests:
            idx = self.manager.requests_frame.index[self.manager.requests_frame['starname']==name][0]
            self.model.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            all_d = list(set(list(self.all_valid_ds_for_request.loc[name].d)))
            available = list(zip(list(self.all_valid_ds_for_request.loc[name].d), list(self.all_valid_ds_for_request.loc[name].s)))
            self.model.addConstr(self.theta[name] >= ((self.manager.requests_frame['n_inter_max'][idx] - \
                        self.past_nights_observed_dict[name]) - (gp.quicksum(self.Yrds[name, d, s] for d, s in available))/self.manager.requests_frame['n_intra_max'][idx]), \
                        'greater_than_nobs_shortfall_' + str(name))

    def constraint_set_max_desired_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraining desired maximum observations.")
        for name in self.multi_visit_requests:
            all_d = list(set(list(self.all_valid_ds_for_request.loc[name].d)))
            self.model.addConstr(gp.quicksum(self.Wrd[name, d] for d in all_d) <=
                        self.desired_max_obs_allowed_dict[name],
                        'max_desired_unique_nights_for_request_' + str(name))
        for name in self.single_visit_requests:
            available = list(zip(list(self.all_valid_ds_for_request.loc[name].d), list(self.all_valid_ds_for_request.loc[name].s)))
            self.model.addConstr(gp.quicksum(self.Yrds[name, d, s] for d, s in available) <=
                        self.desired_max_obs_allowed_dict[name],
                        'max_desired_unique_nights_for_request_' + str(name))

    def remove_constraint_set_max_desired_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Removing previous maximum observations.")
        for name in self.multi_visit_requests:
            rm_const = self.model.getConstrByName("max_desired_unique_nights_for_request_" + str(name))
            self.model.remove(rm_const)

    def constraint_set_max_absolute_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraining absolute maximum observations.")
        for name in self.multi_visit_requests:
            all_d = list(set(list(self.all_valid_ds_for_request.loc[name].d)))
            self.model.addConstr(gp.quicksum(self.Wrd[name, d] for d in all_d) <=
                    self.absolute_max_obs_allowed_dict[name],
                    'max_absolute_unique_nights_for_request_' + str(name))

    def constraint_build_enforce_intranight_cadence(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 6: Enforce intra-night cadence.")
        # get all combos of slots that must be constrained if given slot is scheduled
        # # When intra-night cadence is 0, there will be no keys to constrain so skip
        intracadence_valid_tuples = self.joiner.copy()[self.joiner['n_intra_max'] > 1]
        intracadence_frame = pd.merge(
            intracadence_valid_tuples.drop_duplicates(['id','d','s']),
            intracadence_valid_tuples[['id','d','s']],
            suffixes=['','3'],on=['id', 'd']
            ).query('s + 0 < s3 <= s + tau_intra')
        intracadence_frame = intracadence_frame.groupby(['id','d','s'])[['s3']].agg(list)
        for i, row in intracadence_valid_tuples.iterrows():
            if (row.id, row.d, row.s) in intracadence_frame.index:
                # Get all slots tonight which are too soon after given slot for another visit
                slots_to_constrain_tonight_intra = list(intracadence_frame.loc[(row.id, row.d, row.s)][0])
                self.model.addConstr((self.Yrds[row.id,row.d,row.s] <= (self.Wrd[row.id, row.d] - (gp.quicksum(self.Yrds[row.id,row.d,s3] \
                    for s3 in slots_to_constrain_tonight_intra)))), \
                    'enforce_intranight_cadence_' + row.id + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_set_min_max_visits_per_night(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 7: Allow minimum and maximum visits.")
        intracadence_frame_on_day = self.joiner.copy().drop_duplicates(subset=['id', 'd'])
        grouped_s = self.joiner.copy().groupby(['id', 'd'])['s'].unique().reset_index()
        grouped_s.set_index(['id', 'd'], inplace=True)
        for i, row in intracadence_frame_on_day.iterrows():
            all_valid_slots_tonight = list(grouped_s.loc[(row.id, row.d)]['s'])
            if row.id in self.multi_visit_requests:
                self.model.addConstr((((gp.quicksum(self.Yrds[row.id, row.d,s3] for s3 in all_valid_slots_tonight)))) <= \
                    row.n_intra_max*self.Wrd[row.id, row.d], 'enforce_max_visits1_' + row.id + "_" + str(row.d) + "d_" + str(row.s) + "s")
                self.model.addConstr((((gp.quicksum(self.Yrds[row.id,row.d,s3] for s3 in all_valid_slots_tonight)))) >= \
                    row.n_intra_min*self.Wrd[row.id, row.d], 'enforce_min_visits_' + row.id + "_" + str(row.d) + "d_" + str(row.s) + "s")
            else:
                self.model.addConstr((((gp.quicksum(self.Yrds[row.id,row.d,s3] for s3 in all_valid_slots_tonight)))) <= \
                    row.n_intra_max, 'enforce_min_visits_' + row.id + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def optimize_model(self):
        print("Begin model solve.")
        t1 = time.time()
        self.model.params.TimeLimit = self.manager.solve_time_limit
        self.model.params.TimeLimit = 900
        self.model.Params.OutputFlag = self.manager.gurobi_output
        # Allow stop at 5% gap to prevent from spending lots of time on marginally better solution
        self.model.params.MIPGap = self.manager.solve_max_gap
        # More aggressive presolve gives better solution in shorter time
        self.model.params.Presolve = 0
#        self.model.params.Symmetry = 2
        #self.model.params.Presolve = 0
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
        print("Time to finish solver: ", np.round(time.time()-t1,3))


    def run_model(self):
        self.round_info = 'Round1'
        self.build_model_round1()
        self.optimize_model()
        self.serialize_results_csv()
        self.round_info = 'Round2'
        # np.savetxt(manager.output_directory + 'Round2_Requests.txt', [], delimiter=',', fmt="%s")
        if self.manager.run_round_two:
            self.build_model_round2()
            self.optimize_model()
        self.serialize_results_csv()
        print("Scheduling complete, clear skies!")

    def build_model_round1(self):
        t1 = time.time()

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
            if self.manager.include_aesthetic:
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

    def serialize_results_csv(self):
        print("Building human readable schedule.")
        if self.manager.run_optimal_allocation:
            self.retrieve_ois_solution()
        self.human_read_available = io.write_available_human_readable(self.manager)
        self.human_read_schedule = io.write_stars_schedule_human_readable(self.human_read_available, self.Yrds, self.manager, self.round_info)
        io.build_fullness_report(self.human_read_schedule, self.manager, self.round_info)
        io.write_out_results(self.manager, self.theta, self.round_info, self.start_the_clock)
        mn.get_gap_filler_targets(self.manager)
        io.serialize_schedule(self.Yrds, self.manager,)

    def retrieve_ois_solution(self):
        print("Retrieving results of Optimal Instrument Allocation set of nights.")
        allocation_schedule_1d = []
        for v in self.Anq.values():
            if np.round(v.X,0) == 1:
                allocation_schedule_1d.append(1)
            else:
                allocation_schedule_1d.append(0)
        allocation_schedule = np.reshape(allocation_schedule_1d, (self.manager.n_nights_in_semester, self.manager.n_quarters_in_night))
        self.manager.allocation_map_2D_NQ = allocation_schedule
        print("Printing quicklook optimal instrument allocation statistics.")
        print(self.manager.allocation_map_2D_NQ)
        print(np.sum(self.manager.allocation_map_2D_NQ))
        print(np.sum(self.manager.allocation_map_2D_NQ, axis=1))
        print(np.sum(self.manager.allocation_map_2D_NQ, axis=0))
        weather_holder = np.zeros(np.shape(allocation_schedule))
        allocation_map_1D, allocation_map_2D, weathered_map = mp.build_allocation_map(self.manager, allocation_schedule, weather_holder)
        mp.convert_allocation_array_to_binary(self.manager)
        self.manager.allocation_all = allocation_map_1D
        self.manager.allocation_1D = allocation_map_1D
        self.manager.allocation_map_2D = allocation_map_2D
        self.manager.weathered_map = weathered_map
