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
import kpfcc.helper_functions as hf
import kpfcc.twilight_functions as tw
import kpfcc.reporting_functions as rf
import kpfcc.processing_functions as pf
import kpfcc.mapping_functions as mf

class GorubiModel(object):
    """A Gorubi Model object, from which we can add constraints and share variables easily."""

    def __init__(self, Aset, Aframe, schedulable_requests, requests_frame, database_info_dict,
                 solve_time_limit, gurobi_output, max_gap=0.05, max_bonus_observations_pct = 0.5,
                 run_optimal_allocation=False, semester_grid = [], quarters_grid = []):

        print("Defining model.")
        self.Aset = Aset
        self.Aframe = Aframe
        self.schedulable_requests = schedulable_requests
        self.requests_frame = requests_frame
        self.database_info_dict = database_info_dict
        self.solve_time_limit = solve_time_limit
        self.gurobi_output = gurobi_output
        self.max_gap = max_gap
        self.max_bonus_observations_pct = max_bonus_observations_pct

        self.run_optimal_allocation = run_optimal_allocation
        self.semester_grid = semester_grid
        self.quarters_grid = quarters_grid
        if self.run_optimal_allocation and (self.semester_grid == [] or self.quarters_grid == []):
            print("A semester grid and a quarters grid is required to run optimal allocation.")

        self.m = gp.Model('Semester_Scheduler')
        # Yrs is technically a 1D matrix indexed by tuples.
        # But in practice best think of it as a 2D square matrix of requests r and slots s, with gaps.
        # Slot s for request r will be 1 to indicate starting an exposure for that request in that slot
        self.Yrds = self.m.addVars(self.Aset, vtype = GRB.BINARY, name = 'Requests_Slots')
        # theta is the "shortfall" variable, continous in natural numbers.
        self.theta = self.m.addVars(self.schedulable_requests, name = 'Shortfall')

        if self.run_optimal_allocation:
            # Anq is a 2D matrix of N_nights_in_semester by N_quarters_in_night
            # element will be 1 if that night/quarter is allocated to KPF and 0 otherwise
            self.Anq = self.m.addVars(self.semester_grid, self.quarters_grid, vtype = GRB.BINARY, name = 'Allocation')

            # Un is a 1D matrix of N_nights_in_semester
            # element will be 1 if at least one quarter in that night is allocated
            self.Un = self.m.addVars(self.semester_grid, vtype = GRB.BINARY, name = 'On-Sky')
            self.m.update()

    def add_constraint_build_theta(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 0: Build theta variable")
        aframe_slots_for_request = self.Aframe.groupby(['r'])[['d', 's']].agg(list)
        for name in self.schedulable_requests:
            idx = self.requests_frame.index[self.requests_frame['Starname']==name][0]

            if self.database_info_dict == {}:
                past_nights_observed = 0
            else:
                past_nights_observed = len(self.database_info_dict[name][1])

            # Safety valve for if the target is over-observed for any reason
            # Example: a cadence target that is also an RM target will have many more past
            # observations than is requested.
            # When we move over to parsing the Keck database for a unique request ID, instead of
            # parsing the Jump database on non-unique star name, this can be removed.
            if past_nights_observed > self.requests_frame['# of Nights Per Semester'][idx] + \
                        int(self.requests_frame['# of Nights Per Semester'][idx]*self.max_bonus_observations_pct):
                true_max_obs = past_nights_observed
            else:
                true_max_obs = (self.requests_frame['# of Nights Per Semester'][idx] - past_nights_observed)\
                       + int(self.requests_frame['# of Nights Per Semester'][idx]*self.max_bonus_observations_pct)

            self.m.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            available = list(zip(list(aframe_slots_for_request.loc[name].d), \
                                                            list(aframe_slots_for_request.loc[name].s)))
            self.m.addConstr(self.theta[name] >= ((self.requests_frame['# of Nights Per Semester'][idx] - \
                        past_nights_observed) - gp.quicksum(self.Yrds[name, d, s] for d,s in available)), \
                        'greater_than_nobs_shortfall_' + str(name))
            self.m.addConstr(gp.quicksum(self.Yrds[name, d, s] for d,s in available) <=
                        true_max_obs, 'max_unique_nights_for_request_' + str(name))

    def add_constraint_one_request_per_slot(self):
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
        # Aset_inter_no_duplicates = Aframe.copy()
        Aset_inter_no_duplicates = self.Aframe.drop_duplicates(subset=['d', 's'])
        # Construct the constraint
        for i, row in Aset_ds_no_duplicates.iterrows():
            requests_valid_in_slot = list(ds_requests.loc[(row.d, row.s)])[0]
            self.m.addConstr((gp.quicksum(self.Yrds[name,row.d,row.s] for name in requests_valid_in_slot) <= 1),
                            'one_request_per_slot_' + str(row.d) + "d_" + str(row.s) + "s")

    def add_constraint_reserve_multislot_exposures(self):
        """
        According to Eq X in Lubin et al. 2025.
        """
        print("Constraint 2: Reserve slots for for multi-slot exposures.")
        # Get all requests that are  valid in (d,s+e) pair for a given (d,s,1..e)
        requests_valid_in_reserved_slots = pd.merge(self.Aframe.query('e > 1 ')['r d s e'.split()] \
            ,self.Aframe['r d s'.split()],on=['d'],suffixes=['','2']) \
            .query('s < s2 < s + e').groupby('r d s'.split()).agg(list)
        # If request requires only 1 slot to complete, then no constraint on reserving additional slots
        Aframe_multislots = self.Aframe[self.Aframe.e > 1]
        for i, row in Aframe_multislots.iterrows():
            # construct list of (r,d,s) indices to be constrained. These are all requests that are
            # valid in slots (d, s+1) through (d, s + e)
            # Ensuring slot (d, s) is not double filled already taken care of in Constraint 1.
            allr = list(requests_valid_in_reserved_slots.loc[row.r, row.d, row.s]['r2'])
            alls = list(requests_valid_in_reserved_slots.loc[row.r, row.d, row.s]['s2'])
            all_reserved_slots = list(zip(allr, [row.d]*len(allr), alls))
            self.m.addConstr((row.e*(1 - self.Yrds[row.r,row.d,row.s])) >= gp.quicksum(self.Yrds[c]
                                                            for c in all_reserved_slots),
                            'reserve_multislot_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def add_constraint_max_visits_per_night(self):
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
        unique_request_day_pairs = self.Aframe.drop_duplicates(['r','d'])
        # Build the constraint
        for i, row in unique_request_day_pairs.iterrows():
            constrained_slots_tonight = np.array(slots_on_day.loc[(row.r, row.d)][0])
            self.m.addConstr((gp.quicksum(self.Yrds[row.r,row.d,ss] for ss in constrained_slots_tonight) <= row.v),
                    'max_observations_per_night_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def add_constraint_enforce_internight_cadence(self):
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
            constrained_slots_tonight = np.array(slots_on_day.loc[(row.r, row.d)][0])
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





    def solve_model(self):
        print("Begin model solve.")
        self.m.setObjective(gp.quicksum(self.theta[name] for name in self.schedulable_requests), GRB.MINIMIZE)

        self.m.params.TimeLimit = self.solve_time_limit
        self.m.Params.OutputFlag = self.gurobi_output
        # Allow stop at 5% gap to prevent from spending lots of time on marginally better solution
        self.m.params.MIPGap = self.max_gap
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
            print("Round 1 Model Solved.")
