"""
Module that defines the SemesterPlanner class. This class is responsible for defining, building, and solving the
Gurobi model for semester-level observation planning. It is nearly completely agnostic to all astronomy knowledge.

"""

# Standard library imports
import logging
import os
import time
import warnings
from configparser import ConfigParser
from datetime import datetime, timedelta
import pickle
import json
import h5py

# Third-party imports
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
import astroplan as apl

# Local imports
import astroq.access as ac
import astroq.history as hs
import astroq.io as io

# Suppress warnings
warnings.filterwarnings('ignore')

logs = logging.getLogger(__name__)

class SemesterPlanner(object):
    """
    Define the SemesterPlanner object. This is the heart of AstroQ. 
    This object: 
        - manages and holds the parameters defined in the config.ini 
        - constructs additional metadata for easy sharing/storage across functions
        - builds the Gurobi model
        - defines the constraints
        - sets the objective function
        - kicks off the model solver
        - serializes the results to a csv file
        - saves the object to an hdf5 file for use later by the nplan and plot modules
     """

    def __init__(self, cf, run_band3):
        """
        Initialize the SemesterPlanner object.

        Args:
            cf (str): the path to the config.ini file
            run_band3 (bool): whether to run the band 3 weather loss model (this will be unnecessary in the 2026A semester)

        Returns:
            None
        """

        logs.debug("Building the SemesterPlanner.")
        self.start_the_clock = time.time()

        # Read config file directly
        config = ConfigParser()
        config.read(cf)
        self.run_band3 = run_band3
        self.config = config 

        # Extract configuration parameters from new format
        workdir = str(config.get('global', 'workdir'))
        self.semester_directory = workdir
        self.current_day = str(config.get('global', 'current_day'))
        self.observatory = config.get('global', 'observatory')
        
        # Get semester parameters from semester section
        self.slot_size = config.getint('semester', 'slot_size')
        self.run_weather_loss = config.getboolean('semester', 'run_weather_loss')
        self.solve_time_limit = config.getint('semester', 'max_solve_time')
        self.gurobi_output = config.getboolean('semester', 'show_gurobi_output')
        self.solve_max_gap = config.getfloat('semester', 'max_solve_gap')
        self.max_bonus = config.getfloat('semester', 'maximum_bonus_size')
        self.run_bonus_round = config.getboolean('semester', 'run_bonus_round')
        
        # Output directory
        self.output_directory = workdir + "outputs/"
        check = os.path.isdir(self.output_directory)
        if not check:
            os.makedirs(self.output_directory)
        
        # Set up file paths from data section
        allocation_file_config = str(config.get('data', 'allocation_file'))
        if os.path.isabs(allocation_file_config):
            self.allocation_file = allocation_file_config
        else:
            self.allocation_file = os.path.join(self.semester_directory, allocation_file_config)

        # Define the input files 
        if self.run_band3:
            request_file_config = str(config.get('data', 'filler_file'))
            self.add_twilights()
        else:
            request_file_config = str(config.get('data', 'request_file'))
        if os.path.isabs(request_file_config):
            self.request_file = request_file_config
        else:
            self.request_file = os.path.join(self.semester_directory, request_file_config)
        
        past_file_config = str(config.get('data', 'past_file'))
        if os.path.isabs(past_file_config):
            self.past_file = past_file_config
        else:
            self.past_file = os.path.join(self.semester_directory, past_file_config)
        
        custom_file_config = str(config.get('data', 'custom_file'))
        if os.path.isabs(custom_file_config):
            self.custom_file = custom_file_config
        else:
            self.custom_file = os.path.join(self.semester_directory, custom_file_config)
        
        if not os.path.exists(self.request_file):
            raise FileNotFoundError(f"Requests file not found: {self.request_file}")
        self.requests_frame = pd.read_csv(self.request_file)

        # Data cleaning
        # Fill NaN values with defaults --- for now in early 2025B since we had issues with the webform.
        # Replace "None" strings with NaN first, then fill with defaults to make sure we get them all 
        self.requests_frame['n_intra_max'] = self.requests_frame['n_intra_max'].replace('None', np.nan).fillna(1)
        self.requests_frame['n_intra_min'] = self.requests_frame['n_intra_min'].replace('None', np.nan).fillna(1)
        self.requests_frame['tau_intra'] = self.requests_frame['tau_intra'].replace('None', np.nan).fillna(0)
        # Handle weather band columns - process each band column individually
        for band_num in [1, 2, 3]:
            weather_band_col = f'weather_band_{band_num}'
            if weather_band_col in self.requests_frame.columns:
                self.requests_frame[weather_band_col] = self.requests_frame[weather_band_col].replace('None', np.nan).fillna(False)
        self.requests_frame['unique_id'] = self.requests_frame['unique_id'].astype(str)
        self.requests_frame['starname'] = self.requests_frame['starname'].astype(str)

        # Build the "strategy" dataframe. Note exptime is in minutes and tau_intra is in hours they are both converted to slots here
        strategy = self.requests_frame[['starname', 'unique_id', 'n_intra_min','n_intra_max','n_inter_max','tau_inter']]
        strategy['t_visit'] = (self.requests_frame['exptime'] / 60 / self.slot_size).clip(lower=1).round().astype(int) 
        strategy['tau_intra'] = (self.requests_frame['tau_intra'] * 60 / self.slot_size).round().astype(int) 
        self.strategy = strategy

        # Compile additional data and metadata 
        self.past_history = hs.process_star_history(self.past_file)
        self.slots_needed_for_exposure_dict = self._build_slots_required_dictionary()
        self._calculate_semester_info()
        self.all_dates_dict, self.all_dates_array = self._build_date_dictionary()
        self._calculate_slot_info()
        
        # Observability represents the indices of the slots where targets are observable
        self.observability = self._build_observability()
        self.observability_tuples = list(self.observability.itertuples(index=False, name=None))

        # Joiner combines strategy and observability
        self.joiner = pd.merge(self.strategy, self.observability, on=['unique_id'])
        # add dummy columns for easier joins
        self.joiner['unique_id2'] = self.joiner['unique_id']
        self.joiner['d2'] = self.joiner['d']
        self.joiner['s2'] = self.joiner['s']

        # Determine the nights where multi-visit requests are observable and the list of multi-visit requests
        self.observability_nights = self.joiner[self.joiner['n_intra_max'] > 1][['unique_id', 'd']].drop_duplicates().copy()
        self.multi_visit_requests = list(self.observability_nights['unique_id'].unique())

        # Define subsets of requests 
        self.all_requests = list(self.requests_frame['unique_id'])
        self.schedulable_requests =  list(self.joiner['unique_id'].unique())
        self.single_visit_requests = [item for item in self.schedulable_requests if item not in self.multi_visit_requests]
        for starid in list(self.requests_frame['unique_id']):
            if starid not in self.schedulable_requests:
                starname = self.requests_frame[self.requests_frame['unique_id']==starid]['starname'].values[0]
                logs.warning("Target " + starname + " with unique id " + starid +  " has no valid day/slot pairs and therefore is effectively removed from the model.")

        # Get each request's full list of valid d/s pairs
        self.all_valid_ds_for_request = self.joiner.groupby(['unique_id'])[['d', 's']].agg(list)
        # Get all requests which are valid in slot (d, s)
        requests_valid_for_ds = pd.merge(
            self.joiner.drop_duplicates(['d', 's']),
            self.joiner[['unique_id', 'd', 's']],
            suffixes=['', '3'],
            on=['d', 's'])
        self.requests_valid_for_ds = requests_valid_for_ds.groupby(['d','s'])[['unique_id3']].agg(list)

        # Get all valid d/s pairs
        self.valid_ds_pairs = self.joiner.copy()
        # Get all requests that require multiple slots to complete one observation
        self.multislot_mask = self.joiner.t_visit > 1
        self.multi_slot_frame = self.joiner[self.multislot_mask]
        # Get all valid slots s for request r on day d
        valid_s_for_rd = pd.merge(
            self.joiner.drop_duplicates(['unique_id','d',]),
            self.joiner[['unique_id','d','s']],
            suffixes=['','3'],on=['unique_id']
        ).query('d == d3')
        self.slots_on_day_for_r = valid_s_for_rd.groupby(['unique_id','d'])[['s3']].agg(list)
        # Get all request id's that are valid on a given day
        self.unique_request_on_day_pairs = self.joiner.copy().drop_duplicates(['unique_id','d'])

        # Define the Gurobi model
        self.model = gp.Model('Semester_Scheduler')
        # Yrds is technically a 1D matrix indexed by tuples.
        # But in practice best think of it as a 3D ragged matrix of requests r, nights d, and slots s, with gaps.
        # Day d / Slot s for request r will be 1 to indicate starting an exposure for that request in that day/slot
        observability_array = list(self.observability.itertuples(index=False, name=None))
        self.Yrds = self.model.addVars(observability_array, vtype = GRB.BINARY, name = 'Requests_Slots')

        if len(self.observability_nights) != 0:
            # Wrd is technically a 1D matrix indexed by tuples.
            # But in practice best think of it as a 2D ragged matrix of requests r and nights d, with gaps.
            # Night d for request r will be 1 to indicate at least one exposure is scheduled for this night.
            # Note that Wrd is only valid for requests r which have at least 2 visits requested in the night.
            observability_array_onsky = list(self.observability_nights.itertuples(index=False, name=None))
            self.Wrd = self.model.addVars(observability_array_onsky, vtype = GRB.BINARY, name = 'OnSky')

        # theta is the "shortfall" variable, continous in natural numbers.
        self.theta = self.model.addVars(self.all_requests, name = 'Shortfall')

        # Set the max allowed number of observations for each request based on the past history and the requested number of observations
        desired_max_obs_allowed_dict = {}
        absolute_max_obs_allowed_dict = {}
        past_nights_observed_dict = {}
        for name in self.all_requests:
            idx = self.requests_frame.index[self.requests_frame['unique_id']==name][0]
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
        """
        Compile parameters defining the semester based on the current day.

        Returns:
            None
        """

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
        """
        Construct useful data structures that are used throughout the semester planner.

        Returns:
            all_dates_dict (dict): a dictionary where keys are the dates in the semester and values are the day index
            all_dates_array (list): a list of the dates in the semester
        """
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
        """
        Compute important numbers relating to the quantity of slots.

        Returns:
            None
        """
        # Calculate slots per quarter and night
        self.n_slots_in_night = int(24 * 60 / self.slot_size)
        
        # Calculate remaining semester info
        self.n_nights_in_semester = len(self.all_dates_dict) - self.all_dates_dict[self.current_day]
        self.n_slots_in_semester = self.n_slots_in_night * self.n_nights_in_semester
        
        # Calculate today's starting positions
        self.today_starting_slot = self.all_dates_dict[self.current_day] * self.n_slots_in_night
        self.today_starting_night = self.all_dates_dict[self.current_day]

    def _build_slots_required_dictionary(self, always_round_up_flag=False):
        """
        Determine the number of slots required to complete for each visit of a given request.

        Returns:
            slots_needed_for_exposure_dict (dict): a dictionary where keys are the star names and values are the number of slots required for each exposure
        """
        slots_needed_for_exposure_dict = {}
        for n, row in self.requests_frame.iterrows():
            starid = row['unique_id']
            exposure_time = float(row['exptime'])
            overhead = 45*float(row['n_exp'] - 1) #+ 180*float(row['n_intra_max'])
            
            if always_round_up_flag:
                slots_needed = int(np.ceil((exposure_time + overhead) / (self.slot_size * 60.0)))
            else:
                slots_needed = int(np.round((exposure_time + overhead) / (self.slot_size * 60.0)))
            if slots_needed < 1:
                slots_needed = 1
            
            slots_needed_for_exposure_dict[starid] = slots_needed
        
        return slots_needed_for_exposure_dict

    def _build_observability(self):
        """
        Determine the indices of the slots where targets are observable using the Access object.

        Returns:
            observability (dict): a dictionary where keys are the star names and values are the indices of the slots where the target is observable
        """
        # Create Access object with parameters from config
        self.access_obj = ac.Access(
            semester_start_date=self.semester_start_date,
            semester_length=self.semester_length,
            slot_size=self.slot_size,
            observatory=self.observatory,
            current_day=self.current_day,
            all_dates_dict=self.all_dates_dict,
            all_dates_array=self.all_dates_array,
            n_nights_in_semester=self.n_nights_in_semester,
            custom_file=self.custom_file,
            allocation_file=self.allocation_file,
            past_history=self.past_history,
            today_starting_night=self.today_starting_night,
            slots_needed_for_exposure_dict=self.slots_needed_for_exposure_dict,
            run_weather_loss=self.run_weather_loss,
            output_directory=self.output_directory,
            run_band3=self.run_band3
        )
        # Store the full access record array for later use
        self.access_record = self.access_obj.produce_ultimate_map(self.requests_frame)
        observability = self.access_obj.observability(self.requests_frame, access=self.access_record)

        return observability

    def _compute_slots_required_for_exposure(self, exposure_time, slot_size, always_round_up_flag):
        """
        Compute the number of slots required for a given exposure time.
        
        Args:
            exposure_time: Exposure time in minutes
            slot_size: Slot size in minutes
            always_round_up_flag: If True, always round up to the next slot
            
        Returns:
            slots_needed (int): the number of slots required for the given exposure time
        """
        time_per_slot = slot_size / 60.0  # Convert slot_size to hours
        
        if always_round_up_flag:
            return math.ceil(exposure_time / time_per_slot)
        else:
            return round(exposure_time / time_per_slot)

    def constraint_reserve_multislot_exposures(self):
        """
        See Constraint 1 in Lubin et al. 2025.

        Reserve multiple time slots for exposures that require
        more than one time slot to complete, and ensure that
        no other observations are scheduled during these slots.
        """
        logs.info("Constraint: Reserve slots for for multi-slot exposures.")
        max_t_visit = self.strategy.t_visit.max() # longest exposure time
        R_ds = self.observability.groupby(['d','s'])['unique_id'].apply(set).to_dict()
        R_geq_t_visit = {} # dictionary of requests where t_visit is greater than or equal to t_visit
        strategy = self.strategy
        for t_visit in range(1,max_t_visit+1):
            R_geq_t_visit[t_visit] = set(strategy[strategy.t_visit >= t_visit]['unique_id'])

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
        See Constraint 3 in Lubin et al. 2025.

        Ensure that the minimum number of days pass between
        consecutive observations of a given target. If a
        target is scheduled for observation on a given date,
        prevent it from being scheduled again until the
        minimum number of days have passed.
        """
        logs.info("Constraint: Enforce inter-night cadence.")
        # Get all (d',s') pairs for a request that must be zero if a (d,s) pair is selected
        # Ensure tau_inter is numeric before the query
        joiner_for_intercadence = self.joiner.copy()
        joiner_for_intercadence['tau_inter'] = pd.to_numeric(joiner_for_intercadence['tau_inter'], errors='coerce')
        
        intercadence = pd.merge(
            joiner_for_intercadence.drop_duplicates(['unique_id','d',]),
            joiner_for_intercadence[['unique_id','d','s']],
            suffixes=['','3'],on=['unique_id']
        ).query('d + 0 < d3 < d + tau_inter')
        self.intercadence_tracker = intercadence.groupby(['unique_id','d'])[['d3','s3']].agg(list)
        # When inter-night cadence is 1, there will be no keys to constrain so skip
        # While the if/else statement would catch these, by shrinking the list here we do fewer
        # total steps in the loop.
        intercadence_valid_tuples = self.joiner.copy()[self.joiner['tau_inter'] > 1]
        # We don't want duplicate slots on day d because we only need this constraint once per day
        # With duplicates, the same constraint would be applied to (r, d, s) and (r, d, s+1) which
        # is superfluous since we are summing over tonight's slots
        intercadence_valid_tuples = intercadence_valid_tuples.drop_duplicates(subset=['unique_id', 'd'])
        for i, row in intercadence_valid_tuples.iterrows():
            constrained_slots_tonight = np.array(self.slots_on_day_for_r.loc[(row.unique_id2, row.d2)][0])
            # Get all slots for pair (r, d) where valid
            if (row.unique_id, row.d) in self.intercadence_tracker.index:
                slots_to_constrain_future = self.intercadence_tracker.loc[(row.unique_id2, row.d2)]
                ds_pairs = zip(list(np.array(slots_to_constrain_future.d3).flatten()), list(np.array(slots_to_constrain_future.s3).flatten()))

                lhs = gp.quicksum(self.Yrds[row.unique_id,row.d,s2] for s2 in constrained_slots_tonight)/row.n_intra_max 
                rhs = 1 - (gp.quicksum(self.Yrds[row.unique_id,d3,s3] for d3, s3 in ds_pairs))
                self.model.addConstr(lhs <= rhs, 'enforce_internight_cadence_' + row.unique_id + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_fix_previous_objective(self, epsilon=0.03):
        """
        Bonus round constraint: not featured in Lubin et al. 2025.

        This constraint ensures that the
        objective function value calculated during
        Round 2 be within a given tolerance of the
        Round 1 value. This constraint ensures that
        Round 2 result in only small changes to the
        optimal solution found in Round 1.
        """
        logs.info("Constraint: Fixing the previous solution's objective value.")
        lhs = gp.quicksum(self.theta[name] for name in self.requests_frame['unique_id'])
        rhs = self.model.objval + epsilon
        self.model.addConstr(lhs <= rhs, 'fix_previous_objective')

    def set_objective_maximize_slots_used(self):
        """
        Bonus round constraint: not featured in Lubin et al. 2025.

        In Round 2, maximize the number of filled slots,
        i.e., slots during which an exposure occurs.
        """
        logs.info("Objective: Maximize the number of slots used.")
        self.model.setObjective(gp.quicksum(self.slots_needed_for_exposure_dict[id]*self.Yrds[id,d,s]
                        for id, d, s in self.observability_tuples), GRB.MAXIMIZE)

    def set_objective_minimize_theta_time_normalized(self):
        """
        See Equation 1 in Lubin et al. 2025.

        Minimize the total shortfall for the number
        of targets that receive their requested number
        of observations, weighted by the time needed
        to complete one observation.
        """
        self.model.setObjective(gp.quicksum(self.theta[name]*self.slots_needed_for_exposure_dict[name] for name in self.schedulable_requests), GRB.MINIMIZE)

    def constraint_build_theta_multivisit(self):
        """
        See Equation 3 in Lubin et al. 2025.

        Definition of the "shortfall" matrix, Theta.
        The shortfall is defined for each target,
        giving for each target the difference between
        the number of requested nights for that target
        and the sum of the past and future scheduled
        observations of that target.
        """
        logs.info("Constraint 0: Build theta variable")
        for starid in self.schedulable_requests:
            idx = self.requests_frame.index[self.requests_frame['unique_id']==starid][0]
            lhs1 = self.theta[starid]
            rhs1 = 0
            self.model.addConstr(lhs1 >= rhs1, 'greater_than_zero_shortfall_' + str(starid))
            
            # Get all (d,s) pairs for which this request is valid.
            all_d = list(set(list(self.all_valid_ds_for_request.loc[starid].d)))
            available = list(zip(list(self.all_valid_ds_for_request.loc[starid].d), list(self.all_valid_ds_for_request.loc[starid].s)))
            lhs2 = self.theta[starid]
            rhs2 = self.requests_frame['n_inter_max'][idx] - self.past_nights_observed_dict[starid] - (gp.quicksum(self.Yrds[starid, d, s] for d, s in available))/self.requests_frame['n_intra_max'][idx]
            self.model.addConstr(lhs2 >= rhs2, 'greater_than_nobs_shortfall_' + str(starid))

    def constraint_set_max_desired_unique_nights_Wrd(self):
        """
        See Constraint 2 in Lubin et al. 2025.

        Limit the number of observations scheduled for a given
        target to the maximum value provided by the PI. This
        constraint may later be relaxed if Round 2 of scheduling
        is invoked.
        """
        logs.info("Constraint: Set desired maximum observations.")
        for name in self.multi_visit_requests:
            all_d = list(set(list(self.all_valid_ds_for_request.loc[name].d)))
            lhs1 = gp.quicksum(self.Wrd[name, d] for d in all_d)
            rhs1 = self.desired_max_obs_allowed_dict[name]
            self.model.addConstr(lhs1 <= rhs1, 'max_desired_unique_nights_for_request_' + str(name))
        
        for name in self.single_visit_requests:
            available = list(zip(list(self.all_valid_ds_for_request.loc[name].d), list(self.all_valid_ds_for_request.loc[name].s)))
            lhs2 = gp.quicksum(self.Yrds[name, d, s] for d, s in available)
            rhs2 = self.desired_max_obs_allowed_dict[name]
            self.model.addConstr(lhs2 <= rhs2, 'max_desired_unique_nights_for_request_' + str(name))

    def remove_constraint_set_max_desired_unique_nights_Wrd(self):
        """
        Bonus round constraint: not featured in Lubin et al. 2025.

        Remove the maximum number of observations set by
        constraints_set_max_desired_unique_nights_Wrd.
        """
        logs.info("Constraint: Removing previous maximum observations constraint.")
        for name in self.multi_visit_requests:
            rm_const = self.model.getConstrByName("max_desired_unique_nights_for_request_" + str(name))
            self.model.remove(rm_const)

    def constraint_set_max_absolute_unique_nights_Wrd(self):
        """
        Bonus round constraint: not featured in Lubin et al. 2025.

        Set the maximum number of observations for a target to
        150% of the original requested number.
        """
        logs.info("Constraint: Set absolute maximum observations.")
        for name in self.multi_visit_requests:
            all_d = list(set(list(self.all_valid_ds_for_request.loc[name].d)))
            lhs = gp.quicksum(self.Wrd[name, d] for d in all_d)
            rhs = self.absolute_max_obs_allowed_dict[name]
            self.model.addConstr(lhs <= rhs, 'max_absolute_unique_nights_for_request_' + str(name))

    def constraint_build_enforce_intranight_cadence(self):
        """
        Constraint 4 in Lubin et al. 2025.

        Ensure that the minimum number of hours pass between
        consecutive observations of a given target on the same
        night. If a target is scheduled for observation at
        a given time, prevent it from being scheduled again
        until the minimum number of hours have passed.
        """
        logs.info("Constraint: Enforce intra-night cadence.")
        # get all combos of slots that must be constrained if given slot is scheduled
        # # When intra-night cadence is 0, there will be no keys to constrain so skip
        intracadence_valid_tuples = self.joiner.copy()[self.joiner['n_intra_max'] > 1]
        intracadence_frame = pd.merge(
            intracadence_valid_tuples.drop_duplicates(['unique_id','d','s']),
            intracadence_valid_tuples[['unique_id','d','s']],
            suffixes=['','3'],on=['unique_id', 'd']
            ).query('s + 0 < s3 < s + tau_intra')
        intracadence_frame = intracadence_frame.groupby(['unique_id','d','s'])[['s3']].agg(list)
        for i, row in intracadence_valid_tuples.iterrows():
            if (row.unique_id, row.d, row.s) in intracadence_frame.index:
                # Get all slots tonight which are too soon after given slot for another visit
                slots_to_constrain_tonight_intra = list(intracadence_frame.loc[(row.unique_id, row.d, row.s)][0])
                lhs = self.Yrds[row.unique_id,row.d,row.s]
                rhs = self.Wrd[row.unique_id, row.d] - (gp.quicksum(self.Yrds[row.unique_id,row.d,s3] for s3 in slots_to_constrain_tonight_intra))
                self.model.addConstr(lhs <= rhs, 'enforce_intranight_cadence_' + row.unique_id + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def constraint_set_min_max_visits_per_night(self):
        """
        See Constraint 5 in Lubin et al. 2025.

        Require that the number of scheduled visits to a target
        in a given night falls between the minimum and maximum
        values supplied by the PI.
        """
        logs.info("Constraint: Bound minimum and maximum visits per night.")
        intracadence_frame_on_day = self.joiner.copy().drop_duplicates(subset=['unique_id', 'd'])
        grouped_s = self.joiner.copy().groupby(['unique_id', 'd'])['s'].unique().reset_index()
        grouped_s.set_index(['unique_id', 'd'], inplace=True)
        for i, row in intracadence_frame_on_day.iterrows():
            all_valid_slots_tonight = list(grouped_s.loc[(row.unique_id, row.d)]['s'])
            if row.unique_id in self.multi_visit_requests:
                lhs1 = gp.quicksum(self.Yrds[row.unique_id, row.d,s3] for s3 in all_valid_slots_tonight)
                rhs1 = row.n_intra_max*self.Wrd[row.unique_id, row.d]
                self.model.addConstr(lhs1 <= rhs1, 'enforce_max_visits1_' + row.unique_id + "_" + str(row.d) + "d_" + str(row.s) + "s")
                
                lhs2 = gp.quicksum(self.Yrds[row.unique_id,row.d,s3] for s3 in all_valid_slots_tonight)
                rhs2 = row.n_intra_min*self.Wrd[row.unique_id, row.d]
                self.model.addConstr(lhs2 >= rhs2, 'enforce_min_visits_' + row.unique_id + "_" + str(row.d) + "d_" + str(row.s) + "s")
            else:
                lhs3 = gp.quicksum(self.Yrds[row.unique_id,row.d,s3] for s3 in all_valid_slots_tonight)
                rhs3 = row.n_intra_max
                self.model.addConstr(lhs3 <= rhs3, 'enforce_min_visits_' + row.unique_id + "_" + str(row.d) + "d_" + str(row.s) + "s")

    def optimize_model(self):
        """
        Solve the Gurobi model.

        Returns:
            None
        """

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
        """
        Construct and solve the Gurobi model.

        Returns:
            None
        """
        self.round_info = 'Round1'
        self.build_model_round1()
        self.optimize_model()
        self.serialize_results_csv()
        logs.info("Round 1 complete.")
        if self.run_bonus_round:
            self.round_info = 'Round2'
            self.build_model_round2()
            self.optimize_model()
            self.serialize_results_csv()
            logs.info("Round 2 complete.")
        logs.info("Scheduling complete, clear skies!")

    def build_model_round1(self):
        """
        Implement the constraints and objective function for Round 1 as described in Lubin et al. 2025.

        Returns:
            None
        """
        t1 = time.time()
        self.constraint_reserve_multislot_exposures()
        self.constraint_enforce_internight_cadence()
        self.constraint_set_max_desired_unique_nights_Wrd()
        self.constraint_build_enforce_intranight_cadence()
        self.constraint_set_min_max_visits_per_night()
        self.constraint_build_theta_multivisit()
        self.set_objective_minimize_theta_time_normalized()
        logs.info(f"Time to build constraints: {np.round(time.time()-t1,3):.3f}")

    def build_model_round2(self):
        """
        Implement the constraints and objective function for Round 2. Not described in Lubin et al. 2025.

        Returns:
            None
        """
        t1 = time.time()
        self.remove_constraint_set_max_desired_unique_nights_Wrd()
        self.constraint_set_max_absolute_unique_nights_Wrd()
        self.constraint_fix_previous_objective()
        self.set_objective_maximize_slots_used()
        logs.info(f"Time to build constraints: {np.round(time.time()-t1,3):.3f}")

    def serialize_results_csv(self):
        """
        Serialize the results to a CSV file.

        Returns:
            None
        """

        logs.debug("Building human readable schedule.")
        serialized_schedule = io.serialize_schedule(self.Yrds, self)
        self.serialized_schedule = serialized_schedule
        today_idx = self.all_dates_dict[self.current_day]
        selected = [k[0] for k, v in self.Yrds.items() if v.x > 0 and k[1] == today_idx]
        selected = list(set(selected))
        selected_df = self.requests_frame[self.requests_frame['unique_id'].isin(selected)].copy()
        selected_df.to_csv(os.path.join(self.output_directory, 'request_selected.csv'), index=False)
        self.to_hdf5()

    def to_hdf5(self, hdf5_path=None):
        """
        Save the SemesterPlanner object to an HDF5 file.
        
        Args:
            hdf5_path (str, optional): Path to save the HDF5 file. 
                                      If None, saves to output_directory/semester_planner.h5
        """
        if hdf5_path is None:
            hdf5_path = os.path.join(self.output_directory, 'semester_planner.h5')
        
        # Remove existing file if it exists
        if os.path.exists(hdf5_path):
            os.remove(hdf5_path)
        
        # Save DataFrames using pandas HDF5 support
        dataframes_to_save = {
            'requests_frame': self.requests_frame,
            'serialized_schedule': self.serialized_schedule,
        }
        
        for key, df in dataframes_to_save.items():
            if df is not None:
                df.to_hdf(hdf5_path, key=key, mode='a', format='table')
        
        # Save numpy arrays and other data using h5py
        with h5py.File(hdf5_path, 'a') as f:
            # Save numpy arrays
            if hasattr(self, 'access_record') and self.access_record is not None:
                # access_record is a structured array, save its fields
                for field_name in self.access_record.dtype.names:
                    f.create_dataset(f'access_record/{field_name}', 
                                   data=self.access_record[field_name], 
                                   compression='gzip')
            
            # Save scalar attributes
            scalars = {
                'current_day': self.current_day,
                'semester_start_date': self.semester_start_date,
                'semester_length': self.semester_length,
                'semester_letter': self.semester_letter,
                'slot_size': self.slot_size,
                'n_slots_in_night': self.n_slots_in_night,
                'n_nights_in_semester': self.n_nights_in_semester,
                'n_slots_in_semester': self.n_slots_in_semester,
                'today_starting_slot': self.today_starting_slot,
                'today_starting_night': self.today_starting_night,
                'run_band3': self.run_band3,
                'observatory': self.observatory,
                'output_directory': self.output_directory,
                'run_weather_loss': self.run_weather_loss,
                'solve_time_limit': self.solve_time_limit,
                'gurobi_output': self.gurobi_output,
                'solve_max_gap': self.solve_max_gap,
                'max_bonus': self.max_bonus,
                'run_bonus_round': self.run_bonus_round,
                'semester_directory': self.semester_directory,
                'custom_file': getattr(self, 'custom_file', None),
                'allocation_file': getattr(self, 'allocation_file', None),
            }
            
            for key, value in scalars.items():
                if value is not None:
                    f.attrs[key] = value
            
            # Save lists
            if self.all_dates_array is not None:
                f.create_dataset('all_dates_array', data=np.array(self.all_dates_array, dtype='S'))

            # Save dictionaries as JSON strings
            dicts_to_save = {
                'all_dates_dict': getattr(self, 'all_dates_dict', None),
                'slots_needed_for_exposure_dict': getattr(self, 'slots_needed_for_exposure_dict', None),
                'past_nights_observed_dict': getattr(self, 'past_nights_observed_dict', None),
            }
            
            for key, value in dicts_to_save.items():
                if value is not None:
                    f.attrs[f'{key}_json'] = json.dumps(value)
            
            # Save past_history dictionary (convert StarHistory namedtuples to dict)
            if hasattr(self, 'past_history') and self.past_history is not None:
                past_history_serialized = {}
                for star_id, star_hist in self.past_history.items():
                    past_history_serialized[star_id] = {
                        'name': star_hist.name if hasattr(star_hist, 'name') else star_id,
                        'date_last_observed': star_hist.date_last_observed if hasattr(star_hist, 'date_last_observed') else None,
                        'total_n_exposures': star_hist.total_n_exposures if hasattr(star_hist, 'total_n_exposures') else 0,
                        'total_n_visits': star_hist.total_n_visits if hasattr(star_hist, 'total_n_visits') else 0,
                        'total_n_unique_nights': star_hist.total_n_unique_nights if hasattr(star_hist, 'total_n_unique_nights') else 0,
                        'total_open_shutter_time': star_hist.total_open_shutter_time if hasattr(star_hist, 'total_open_shutter_time') else 0,
                        'n_obs_on_nights': star_hist.n_obs_on_nights if hasattr(star_hist, 'n_obs_on_nights') else [],
                        'n_visits_on_nights': star_hist.n_visits_on_nights if hasattr(star_hist, 'n_visits_on_nights') else [],
                    }
                f.attrs['past_history_json'] = json.dumps(past_history_serialized)
        
        logs.info(f"SemesterPlanner saved to HDF5: {hdf5_path}")
        return hdf5_path

    @classmethod
    def from_hdf5(cls, hdf5_path):
        """
        Load a SemesterPlanner object from an HDF5 file.
        
        Args:
            hdf5_path (str): Path to the HDF5 file
            
        Returns:
            SemesterPlanner: Reconstructed SemesterPlanner object
        """        
        # Create a new instance without calling __init__
        instance = cls.__new__(cls)
        
        # Load DataFrames
        try:
            instance.requests_frame = pd.read_hdf(hdf5_path, key='requests_frame')
        except KeyError:
            instance.requests_frame = None
        
        try:
            instance.serialized_schedule = pd.read_hdf(hdf5_path, key='serialized_schedule')
        except KeyError:
            instance.serialized_schedule = None
        
        # Load other data from HDF5
        with h5py.File(hdf5_path, 'r') as f:
            # Load scalar attributes
            for key in ['current_day', 'semester_start_date', 'semester_length', 'semester_letter',
                       'slot_size', 'n_slots_in_night', 'n_nights_in_semester', 'n_slots_in_semester',
                       'today_starting_slot', 'today_starting_night', 'run_band3', 'observatory',
                       'output_directory', 'run_weather_loss', 'solve_time_limit', 'gurobi_output',
                       'solve_max_gap', 'max_bonus', 'run_bonus_round', 'semester_directory',
                       'custom_file', 'allocation_file']:
                if key in f.attrs:
                    setattr(instance, key, f.attrs[key])
            
            # Load access_record (structured array)
            if 'access_record' in f:
                # Reconstruct structured array from saved fields
                field_names = list(f['access_record'].keys())
                if field_names:
                    # Load all field data first
                    field_data = {}
                    for field_name in field_names:
                        field_data[field_name] = f[f'access_record/{field_name}'][:]
                    
                    # Determine the number of records (first dimension of first field)
                    first_field = field_data[field_names[0]]
                    n_records = first_field.shape[0]
                    
                    # Create dtype list with proper shapes for multidimensional fields
                    dtype_list = []
                    for field_name in field_names:
                        data = field_data[field_name]
                        if data.ndim == 1:
                            # 1D field: just use the dtype
                            dtype_list.append((field_name, data.dtype))
                        else:
                            # Multidimensional field: include shape (excluding first dimension)
                            dtype_list.append((field_name, data.dtype, data.shape[1:]))
                    
                    # Create structured array and populate it
                    struct_array = np.zeros(n_records, dtype=dtype_list)
                    for field_name in field_names:
                        struct_array[field_name] = field_data[field_name]
                    
                    # Convert to recarray so we can use dot notation (e.g., access_record.is_observable)
                    instance.access_record = struct_array.view(np.recarray)
                else:
                    instance.access_record = None
            else:
                instance.access_record = None
            
            # Load lists
            if 'all_dates_array' in f:
                instance.all_dates_array = [d.decode('utf-8') if isinstance(d, bytes) else d 
                                           for d in f['all_dates_array'][:]]
            
            # Load dictionaries from JSON
            for key in ['all_dates_dict', 'slots_needed_for_exposure_dict', 
                       'past_nights_observed_dict']:
                json_key = f'{key}_json'
                if json_key in f.attrs:
                    setattr(instance, key, json.loads(f.attrs[json_key]))
            
            # Load past_history (reconstruct StarHistory namedtuples from dict)
            if 'past_history_json' in f.attrs:
                past_history_data = json.loads(f.attrs['past_history_json'])
                # Import StarHistory namedtuple
                from astroq.history import StarHistory
                instance.past_history = {}
                for star_id, hist_data in past_history_data.items():
                    # Create StarHistory namedtuple with all required fields
                    star_hist = StarHistory(
                        name=hist_data.get('name', star_id),
                        date_last_observed=hist_data.get('date_last_observed'),
                        total_n_exposures=hist_data.get('total_n_exposures', 0),
                        total_n_visits=hist_data.get('total_n_visits', 0),
                        total_n_unique_nights=hist_data.get('total_n_unique_nights', 0),
                        total_open_shutter_time=hist_data.get('total_open_shutter_time', 0),
                        n_obs_on_nights=hist_data.get('n_obs_on_nights', []),
                        n_visits_on_nights=hist_data.get('n_visits_on_nights', [])
                    )
                    instance.past_history[star_id] = star_hist
            else:
                instance.past_history = {}
        
        # Recreate access_obj using the loaded parameters
        # Only recreate if we have the necessary attributes and requests_frame
        if (hasattr(instance, 'requests_frame') and instance.requests_frame is not None and
            hasattr(instance, 'semester_start_date') and hasattr(instance, 'all_dates_dict')):
            try:
                instance.access_obj = ac.Access(
                    semester_start_date=instance.semester_start_date,
                    semester_length=instance.semester_length,
                    slot_size=instance.slot_size,
                    observatory=instance.observatory,
                    current_day=instance.current_day,
                    all_dates_dict=getattr(instance, 'all_dates_dict', {}),
                    all_dates_array=getattr(instance, 'all_dates_array', []),
                    n_nights_in_semester=instance.n_nights_in_semester,
                    custom_file=getattr(instance, 'custom_file', None),
                    allocation_file=getattr(instance, 'allocation_file', None),
                    past_history=getattr(instance, 'past_history', {}),
                    today_starting_night=instance.today_starting_night,
                    slots_needed_for_exposure_dict=getattr(instance, 'slots_needed_for_exposure_dict', {}),
                    run_weather_loss=getattr(instance, 'run_weather_loss', False),
                    output_directory=instance.output_directory,
                    run_band3=getattr(instance, 'run_band3', False)
                )
                logs.info("access_obj recreated from HDF5")
            except Exception as e:
                logs.warning(f"Could not recreate access_obj: {e}")
                instance.access_obj = None
        else:
            instance.access_obj = None
        
        logs.info(f"SemesterPlanner loaded from HDF5: {hdf5_path}")
        return instance

    def add_twilights(self):
        """Add 20-minute buffer to allocation times that match 12-degree twilight."""
        observatory = self.config.get('global', 'observatory')
        keck = apl.Observer.at_site(observatory)
        allocation_df = pd.read_csv(self.allocation_file)
        
        for idx, row in allocation_df.iterrows():
            # Get date from start time (first 10 chars = YYYY-MM-DD)
            date_str = str(row['start'])[:10]
            day = Time(date_str, format='iso', scale='utc')
            
            # Get 12-degree twilight times
            evening_12 = keck.twilight_evening_nautical(day, which='next')
            morning_12 = keck.twilight_morning_nautical(day, which='next')
            
            # Check if start time matches evening twilight (within 10 min)
            start_time = Time(row['start'])
            if abs(start_time - evening_12) <= TimeDelta(10, format='jd') / 1440:  # 10 minutes
                adjusted_start = start_time - TimeDelta(20, format='jd') / 1440
                allocation_df.loc[idx, 'start'] = adjusted_start.strftime('%Y-%m-%dT%H:%M')
                logs.info(f"Adjusted start time for {date_str}: subtracted 20 min")
            
            # Check if stop time matches morning twilight (within 10 min)
            stop_time = Time(row['stop'])
            if abs(stop_time - morning_12) <= TimeDelta(10, format='jd') / 1440:  # 10 minutes
                adjusted_stop = stop_time + TimeDelta(20, format='jd') / 1440
                allocation_df.loc[idx, 'stop'] = adjusted_stop.strftime('%Y-%m-%dT%H:%M')
                logs.info(f"Adjusted stop time for {date_str}: added 20 min")
        
        # Save updated allocation file
        allocation_df.to_csv(self.allocation_file, index=False)
        logs.info("Allocation file updated with twilight adjustments")
