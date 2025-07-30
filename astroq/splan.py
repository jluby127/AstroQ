"""
Module that defines the SemesterPlanner class. This class is responsible for defining, building, and solving the
Gurobi model for semester-level observation planning. It is nearly completely agnostic to all astronomy knowledge.

"""
import sys
import time
import os
import math
import warnings
warnings.filterwarnings('ignore')
import logging
from configparser import ConfigParser
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB
from astropy.time import Time, TimeDelta

import astroq.io as io
import astroq.access as ac
import astroq.history as hs

logs = logging.getLogger(__name__)

class SemesterPlanner(object):
    """A SemesterPlanner object, from which we can define a Gurobi model, build constraints, and solve semester-level observation schedules."""

    def __init__(self, cf):
        logs.debug("Building the SemesterPlanner.")
        self.start_the_clock = time.time()

        # Read config file directly
        config = ConfigParser()
        config.read(cf)
        
        # Extract configuration parameters directly
        upstream_path = eval(config.get('required', 'folder'), {"os": os})
        self.current_day = str(config.get('required', 'current_day'))
        self.observatory = config.get('required', 'observatory')
        self.slot_size = int(config.get('other', 'slot_size'))
        self.n_quarters_in_night = int(config.get('other', 'quarters_in_night'))
        self.n_hours_in_night = int(config.get('other', 'hours_in_night'))
        self.daily_starting_time = str(config.get('other', 'daily_starting_time'))
        
        semester_directory = upstream_path
        
        # Solver configuration
        self.solve_time_limit = int(config.get('gurobi', 'max_solve_time'))
        self.gurobi_output = config.get('gurobi', 'show_gurobi_output').strip().lower() == "true"
        self.solve_max_gap = float(config.get('gurobi', 'max_solve_gap'))
        
        # Other parameters
        self.max_bonus = float(config.get('other', 'maximum_bonus_size'))
        
        # Output directory
        self.output_directory = upstream_path + "outputs/"
        check = os.path.isdir(self.output_directory)
        if not check:
            os.makedirs(self.output_directory)
            file = open(self.output_directory + "runReport.txt", "w")
            file.close()
        
        # Set up file paths
        self.semester_directory = upstream_path
        self.past_file = os.path.join(self.semester_directory, "inputs/past.csv")
        self.custom_file = os.path.join(self.semester_directory, "inputs/custom.csv")
        
        # Resolve allocation file path
        allocation_file_config = str(config.get('options', 'allocation_file'))
        if os.path.isabs(allocation_file_config):
            self.allocation_file = allocation_file_config
        else:
            self.allocation_file = os.path.join(self.semester_directory, allocation_file_config)
        
        # Load data files
        self.requests_frame = pd.read_csv(os.path.join(self.semester_directory, "inputs/requests.csv"))
        self.strategy = self.requests_frame[['starname','n_intra_min','n_intra_max','tau_intra','n_inter_max','tau_inter']]
        self.strategy = self.strategy.rename(columns={'starname':'id'})
        self.strategy['t_visit'] = (self.requests_frame['exptime'] / 60 / self.slot_size).clip(lower=1).round().astype(int) 


        # Load past history
        self.past_history = hs.process_star_history(self.past_file)
        
        # Build slots needed dictionary
        self.slots_needed_for_exposure_dict = self._build_slots_required_dictionary()
        
        # Calculate semester info
        self._calculate_semester_info()
        
        # Build date dictionary
        self.all_dates_dict, self.all_dates_array = self._build_date_dictionary()
        
        # Calculate slot info
        self._calculate_slot_info()
        
        # Build strategy and observability data
        self.observability = self._build_observability()
        # Build meta data
        daily_starting_time = str(config.get('other', 'daily_starting_time'))
        current_day = str(config.get('required', 'current_day'))
        slot_size = int(config.get('other', 'slot_size'))
        self.meta = {"s1_time":daily_starting_time, "d1_date":current_day, "slot_duration":slot_size}

        self.observability_tuples = list(self.observability.itertuples(index=False, name=None))

        self.joiner = pd.merge(self.strategy, self.observability, on=['id'])
        # add dummy columns for easier joins
        self.joiner['id2'] = self.joiner['id']
        self.joiner['d2'] = self.joiner['d']
        self.joiner['s2'] = self.joiner['s']
        self.joiner['tau_intra'] *= int(60/self.slot_size) # convert hours to slots
        self.joiner['tau_intra'] += self.joiner['t_visit'] # start the minimum intracadence time from the end of the previous exposure, not the beginning

        # Prepare information by construction observability_nights (Wset) and schedulable_requests
        self.observability_nights = self.joiner[self.joiner['n_intra_max'] > 1][['id', 'd']].drop_duplicates().copy()
        self.multi_visit_requests = list(self.observability_nights['id'].unique())

        self.all_requests = list(self.requests_frame['starname'])
        self.schedulable_requests =  list(self.joiner['id'].unique())
        self.single_visit_requests = [item for item in self.schedulable_requests if item not in self.multi_visit_requests]
        for name in list(self.requests_frame['starname']):
            if name not in self.schedulable_requests:
                logs.warning("Target " + name + " has no valid day/slot pairs and therefore is effectively removed from the model.")

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

        self.model = gp.Model('Semester_Scheduler')
        # Yrs is technically a 1D matrix indexed by tuples.
        # But in practice best think of it as a 2D square matrix of requests r and slots s, with gaps.
        # Slot s for request r will be 1 to indicate starting an exposure for that request in that slot
        observability_array = list(self.observability.itertuples(index=False, name=None))
        self.Yrds = self.model.addVars(observability_array, vtype = GRB.BINARY, name = 'Requests_Slots')

        if len(self.observability_nights) != 0:
            # Wrd is technically a 1D matrix indexed by tuples.
            # But in practice best think of it as a 2D square matrix of requests r and nights d, with gaps.
            # Night d for request r will be 1 to indicate at least one exposure is scheduled for this night.
            # Note that Wrd is only valid for requests r which have at least 2 visits requested in the night.
            observability_array_onsky = list(self.observability_nights.itertuples(index=False, name=None))
            self.Wrd = self.model.addVars(observability_array_onsky, vtype = GRB.BINARY, name = 'OnSky')

        # theta is the "shortfall" variable, continous in natural numbers.
        self.theta = self.model.addVars(self.all_requests, name = 'Shortfall')

        desired_max_obs_allowed_dict = {}
        absolute_max_obs_allowed_dict = {}
        past_nights_observed_dict = {}
        for name in self.all_requests:
            idx = self.requests_frame.index[self.requests_frame['starname']==name][0]
            if name in list(self.past_history.keys()):
                past_nights_observed = self.past_history[name].total_n_unique_nights
            else:
                past_nights_observed = 0

            # Safety valve for if the target is over-observed for any reason
            if past_nights_observed > self.requests_frame['n_inter_max'][idx] + \
                        int(self.requests_frame['n_inter_max'][idx]*self.max_bonus):
                desired_max_obs = past_nights_observed
            else:
                desired_max_obs = (self.requests_frame['n_inter_max'][idx] - past_nights_observed)
                absolute_max_obs = (self.requests_frame['n_inter_max'][idx] - past_nights_observed) \
                        + int(self.requests_frame['n_inter_max'][idx]*self.max_bonus)
                # second safety valve
                if past_nights_observed > absolute_max_obs:
                    absolute_max_obs = past_nights_observed
            past_nights_observed_dict[name] = past_nights_observed
            desired_max_obs_allowed_dict[name] = desired_max_obs
            absolute_max_obs_allowed_dict[name] = absolute_max_obs
            self.desired_max_obs_allowed_dict = desired_max_obs_allowed_dict
            self.absolute_max_obs_allowed_dict = absolute_max_obs_allowed_dict
            self.past_nights_observed_dict = past_nights_observed_dict
        logs.debug("Initializing complete.")

    def _calculate_semester_info(self):
        """Calculate semester information based on current day."""
        current_date = datetime.strptime(self.current_day, '%Y-%m-%d')
        
        # Determine semester based on date
        if current_date.month in [8, 9, 10, 11, 12]:
            semester_letter = 'A'
            semester_year = current_date.year
        else:
            semester_letter = 'B'
            semester_year = current_date.year - 1
        
        # Set semester boundaries
        if semester_letter == 'A':
            semester_start_date = f"{semester_year}-08-01"
            semester_end_date = f"{semester_year + 1}-01-31"
        else:
            semester_start_date = f"{semester_year + 1}-02-01"
            semester_end_date = f"{semester_year + 1}-07-31"
        
        # Calculate semester length
        start_date = datetime.strptime(semester_start_date, '%Y-%m-%d')
        end_date = datetime.strptime(semester_end_date, '%Y-%m-%d')
        semester_length = (end_date - start_date).days + 1
        
        self.semester_start_date = semester_start_date
        self.semester_length = semester_length
        self.semester_letter = semester_letter

    def _build_date_dictionary(self):
        """Build date dictionary for the semester."""
        start_date = datetime.strptime(self.semester_start_date, '%Y-%m-%d')
        end_date = datetime.strptime(self.semester_start_date, '%Y-%m-%d') + timedelta(days=self.semester_length - 1)
        
        all_dates_dict = {}
        all_dates_array = []
        
        current_date = start_date
        day_index = 0
        
        while current_date <= end_date:
            date_str = current_date.strftime('%Y-%m-%d')
            all_dates_dict[date_str] = day_index
            all_dates_array.append(date_str)
            current_date += timedelta(days=1)
            day_index += 1
        
        return all_dates_dict, all_dates_array

    def _calculate_slot_info(self):
        """Calculate slot-related information."""
        # Calculate slots per quarter and night
        self.n_slots_in_quarter = int(((self.n_hours_in_night * 60) / self.n_quarters_in_night) / self.slot_size)
        self.n_slots_in_night = self.n_slots_in_quarter * self.n_quarters_in_night
        
        # Calculate remaining semester info
        self.n_nights_in_semester = len(self.all_dates_dict) - self.all_dates_dict[self.current_day]
        self.n_slots_in_semester = self.n_slots_in_night * self.n_nights_in_semester
        
        # Calculate today's starting positions
        self.today_starting_slot = self.all_dates_dict[self.current_day] * self.n_slots_in_night
        self.today_starting_night = self.all_dates_dict[self.current_day]

    def _build_slots_required_dictionary(self, always_round_up_flag=False):
        """Build dictionary mapping star names to required slots."""
        slots_needed_for_exposure_dict = {}
        for n, row in self.requests_frame.iterrows():
            name = row['starname']
            exposure_time = float(row['exptime'])
            
            if always_round_up_flag:
                slots_needed = int(np.ceil(exposure_time / (self.slot_size * 60.0)))
            else:
                slots_needed = int(np.round(exposure_time / (self.slot_size * 60.0)))
            
            slots_needed_for_exposure_dict[name] = slots_needed
        
        return slots_needed_for_exposure_dict

    def _build_observability(self):
        """
        Build strategy and observability dataframes directly from config file.
        This replaces the need for a RequestSet object.
        """
        # Create Access object with parameters from config
        access_obj = ac.Access(
            semester_start_date=self.semester_start_date,
            semester_length=self.semester_length,
            slot_size=self.slot_size,
            observatory=self.observatory,
            current_day=self.current_day,
            all_dates_dict=self.all_dates_dict,
            custom_file=self.custom_file,
            allocation_file=self.allocation_file,
            past_history=self.past_history,
            today_starting_night=self.today_starting_night,
            slots_needed_for_exposure_dict=self.slots_needed_for_exposure_dict
        )
        observability = access_obj.observability(self.requests_frame)

        return observability

    def _compute_slots_required_for_exposure(self, exposure_time, slot_size, always_round_up_flag):
        """
        Compute the number of slots required for a given exposure time.
        
        Args:
            exposure_time: Exposure time in minutes
            slot_size: Slot size in minutes
            always_round_up_flag: If True, always round up to the next slot
            
        Returns:
            Number of slots required
        """
        time_per_slot = slot_size / 60.0  # Convert slot_size to hours
        
        if always_round_up_flag:
            return math.ceil(exposure_time / time_per_slot)
        else:
            return round(exposure_time / time_per_slot)

    def constraint_reserve_multislot_exposures(self):
        """
        According to Eq X in Lubin et al. 2025.

        Round 1:

        Reserve multiple time slots for exposures that require
        more than one time slot to complete, and ensure that
        no other observations are scheduled during these slots.
        """
        logs.info("Constraint 2: Reserve slots for for multi-slot exposures.")
        max_t_visit = self.strategy.t_visit.max() # longest exposure time
        R_ds = self.observability.groupby(['d','s'])['id'].apply(set).to_dict()
        R_geq_t_visit = {} # dictionary of requests where t_visit is greater than or equal to t_visit
        strategy = self.strategy
        for t_visit in range(1,max_t_visit+1):
            R_geq_t_visit[t_visit] = set(strategy[strategy.t_visit >= t_visit]['id'])

        for d,s in self.observability.drop_duplicates(['d','s'])[['d','s']].itertuples(index=False, name=None):
            rhs = []
            for delta in range(1,max_t_visit):
                s_shift = s - delta
                if (d,s_shift) in R_ds:
                    rhs.extend(self.Yrds[r,d,s_shift] for r in R_ds[d,s_shift] & R_geq_t_visit[delta+1])

            lhs = 1 - gp.quicksum(self.Yrds[r,d,s] for r in R_ds[d,s])
            rhs = gp.quicksum(rhs)
            self.model.addConstr(lhs >= rhs, f'reserve_multislot_{d}d_{s}s')

    def constraint_enforce_internight_cadence(self):
        """
        According to Eq X in Lubin et al. 2025.

        Round 1:

        Ensure that the minimum number of days pass between
        consecutive observations of a given target. If a
        target is scheduled for observation on a given date,
        prevent it from being scheduled again until the
        minimum number of days have passed.
        """
        logs.info("Constraint 4: Enforce inter-night cadence.")
        # Get all (d',s') pairs for a request that must be zero if a (d,s) pair is selected
        # Ensure tau_inter is numeric before the query
        joiner_for_intercadence = self.joiner.copy()
        joiner_for_intercadence['tau_inter'] = pd.to_numeric(joiner_for_intercadence['tau_inter'], errors='coerce')
        
        intercadence = pd.merge(
            joiner_for_intercadence.drop_duplicates(['id','d',]),
            joiner_for_intercadence[['id','d','s']],
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

    def constraint_fix_previous_objective(self, epsilon=0.03):
        """
        According to Eq X in Lubin et al. 2025.

        Round 2:

        This constraint ensures that the
        objective function value calculated during
        Round 2 be within a given tolerance of the
        Round 1 value. This constraint ensures that
        Round 2 result in only small changes to the
        optimal solution found in Round 1.
        """
        logs.info("Constraint: Fixing the previous solution's objective value.")
        self.model.addConstr(gp.quicksum(self.theta[name] for name in self.requests_frame['starname']) <= \
                    self.model.objval + epsilon)

    def set_objective_maximize_slots_used(self):
        """
        According to Eq X in Lubin et al. 2025.

        In Round 2, maximize the number of filled slots,
        i.e., slots during which an exposure occurs.
        """
        logs.info("Objective: Maximize the number of slots used.")
        self.model.setObjective(gp.quicksum(self.slots_needed_for_exposure_dict[id]*self.Yrds[id,d,s]
                        for id, d, s in self.observability_tuples), GRB.MAXIMIZE)
                        # for id, d, s in self.joiner), GRB.MAXIMIZE)

    def set_objective_minimize_theta_time_normalized(self):
        """
        According to Eq X in Lubin et al. 2025.

        Minimize the total shortfall for the number
        of targets that receive their requested number
        of observations, weighted by the time needed
        to complete one observation.
        """
        self.model.setObjective(gp.quicksum(self.theta[name]*self.slots_needed_for_exposure_dict[name] for name in self.schedulable_requests), GRB.MINIMIZE)

    def constraint_build_theta_multivisit(self):
        """
        According to Eq X in Lubin et al. 2025.

        Definition of the "shortfall" matrix, Theta.
        The shortfall is defined for each target,
        giving for each target the difference between
        the number of requested nights for that target
        and the sum of the past and future scheduled
        observations of that target.
        """
        logs.info("Constraint 0: Build theta variable")
        for name in self.schedulable_requests:
            idx = self.requests_frame.index[self.requests_frame['starname']==name][0]
            self.model.addConstr(self.theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
            # Get all (d,s) pairs for which this request is valid.
            all_d = list(set(list(self.all_valid_ds_for_request.loc[name].d)))
            available = list(zip(list(self.all_valid_ds_for_request.loc[name].d), list(self.all_valid_ds_for_request.loc[name].s)))
            self.model.addConstr(self.theta[name] >= ((self.requests_frame['n_inter_max'][idx] - \
                    self.past_nights_observed_dict[name]) - (gp.quicksum(self.Yrds[name, d, s] for d, s in available))/self.requests_frame['n_intra_max'][idx]), \
                    'greater_than_nobs_shortfall_' + str(name))

    def constraint_set_max_desired_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.

        Round 1:

        Limit the number of observations scheduled for a given
        target to the maximum value provided by the PI. This
        constraint may later be relaxed if Round 2 of scheduling
        is invoked.
        """
        logs.info("Constraining desired maximum observations.")
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

        Round 2:

        Remove the maximum number of observations set by
        constraints_set_max_desired_unique_nights_Wrd.
        """
        logs.info("Removing previous maximum observations.")
        for name in self.multi_visit_requests:
            rm_const = self.model.getConstrByName("max_desired_unique_nights_for_request_" + str(name))
            self.model.remove(rm_const)

    def constraint_set_max_absolute_unique_nights_Wrd(self):
        """
        According to Eq X in Lubin et al. 2025.

        Round 2:

        Set the maximum number of observations for a target to
        150% of the original requested number.
        """
        logs.info("Constraining absolute maximum observations.")
        for name in self.multi_visit_requests:
            all_d = list(set(list(self.all_valid_ds_for_request.loc[name].d)))
            self.model.addConstr(gp.quicksum(self.Wrd[name, d] for d in all_d) <=
                    self.absolute_max_obs_allowed_dict[name],
                    'max_absolute_unique_nights_for_request_' + str(name))

    def constraint_build_enforce_intranight_cadence(self):
        """
        According to Eq X in Lubin et al. 2025.

        Ensure that the minimum number of hours pass between
        consecutive observations of a given target on the same
        night. If a target is scheduled for observation at
        a given time, prevent it from being scheduled again
        until the minimum number of hours have passed.
        """
        logs.info("Constraint 6: Enforce intra-night cadence.")
        # get all combos of slots that must be constrained if given slot is scheduled
        # # When intra-night cadence is 0, there will be no keys to constrain so skip
        intracadence_valid_tuples = self.joiner.copy()[self.joiner['n_intra_max'] > 1]
        intracadence_frame = pd.merge(
            intracadence_valid_tuples.drop_duplicates(['id','d','s']),
            intracadence_valid_tuples[['id','d','s']],
            suffixes=['','3'],on=['id', 'd']
            ).query('s + 0 < s3 < s + tau_intra')
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

        Require that the number of scheduled visits to a target
        in a given night falls between the minimum and maximum
        values supplied by the PI.
        """
        logs.info("Constraint 7: Allow minimum and maximum visits.")
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
        logs.debug("Begin model solve.")
        t1 = time.time()
        self.model.params.TimeLimit = self.solve_time_limit
        self.model.Params.OutputFlag = self.gurobi_output
        # Allow stop at 5% gap to prevent from spending lots of time on marginally better solution
        self.model.params.MIPGap = self.solve_max_gap
        # More aggressive presolve gives better solution in shorter time
        self.model.params.Presolve = 2
        #self.model.params.Presolve = 0
        self.model.update()
        self.model.optimize()

        if self.model.Status == GRB.INFEASIBLE:
            logs.critical('Model remains infeasible. Searching for invalid constraints.')
            search = self.model.computeIIS()
            logs.critical("Printing bad constraints:")
            for c in self.model.getConstrs():
                if c.IISConstr:
                    logs.critical('%s' % c.ConstrName)
            for c in self.model.getGenConstrs():
                if c.IISGenConstr:
                    logs.critical('%s' % c.GenConstrName)
        else:
            logs.debug("Model Successfully Solved.")
        logs.info("Time to finish solver: {:.3f}".format(time.time()-t1))


    def run_model(self):
        self.round_info = 'Round1'
        self.build_model_round1()
        self.optimize_model()
        self.serialize_results_csv()
        logs.info("Scheduling complete, clear skies!")

    def build_model_round1(self):
        t1 = time.time()

        #self.constraint_one_request_per_slot()
        self.constraint_reserve_multislot_exposures()
        self.constraint_enforce_internight_cadence()

        self.constraint_set_max_desired_unique_nights_Wrd()
        self.constraint_build_enforce_intranight_cadence()
        self.constraint_set_min_max_visits_per_night()
        self.constraint_build_theta_multivisit()

        self.set_objective_minimize_theta_time_normalized()
        logs.info("Time to build constraints: ", np.round(time.time()-t1,3))

    def build_model_round2(self):
        t1 = time.time()
        self.remove_constraint_set_max_desired_unique_nights_Wrd()
        self.constraint_set_max_absolute_unique_nights_Wrd()
        self.constraint_fix_previous_objective()
        self.set_objective_maximize_slots_used()
        logs.info("Time to build constraints: ", np.round(time.time()-t1,3))

    def serialize_results_csv(self):
        logs.debug("Building human readable schedule.")
        
        io.serialize_schedule(self.Yrds, self)

        # --- New: Output selected requests for current day ---
        today_idx = self.all_dates_dict[self.current_day]
        selected = [k[0] for k, v in self.Yrds.items() if v.x > 0 and k[1] == today_idx]
        selected = list(set(selected))
        selected_df = self.requests_frame[self.requests_frame['starname'].isin(selected)].copy()
        # Save to CSV with new name
        selected_df.to_csv(os.path.join(self.output_directory, 'request_selected.csv'), index=False)
