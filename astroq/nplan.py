"""
Night Planning Module (nplan.py)

Module for night-level observation planning and optimization.
Uses the Target & Time Planner (TTP) to optimize nightly observation sequences.

Main Functions:
- run_ttp(config_file): Optimize nightly observation sequences
- produce_bright_backups(config_file): Create backup target lists for poor weather
- prepare_for_ttp(...): Prepare data for the TTP system

See https://github.com/lukehandley/ttp/tree/main for more info about the TTP
"""

# Standard library imports
import os
import pickle
from configparser import ConfigParser

# Third-party imports
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
        
# Local imports
import astroq.access as ac
import astroq.io as io

# TTP imports (assuming TTP is installed separately)
import ttp.formatting as formatting
import ttp.telescope as telescope
import ttp.plotting as plotting
import ttp.model as model

class NightPlanner(object):
    """
    Night Planner for optimizing observation sequences within a single night.
    
    Uses the Target & Time Planner (TTP) to create optimal observation schedules
    and backup target lists for poor weather conditions.
    
    Reads configuration directly from config file and uses SemesterPlanner
    for semester calculations and data management.
    """
    
    def __init__(self, config_file):
        """
        Initialize the Night Planner with a config file.
        
        Args:
            config_file: Path to configuration file
        """
        
        # Parse config file directly for paths (following SemesterPlanner pattern)
        from configparser import ConfigParser
        config = ConfigParser()
        config.read(config_file)
        
        # Get workdir from global section
        workdir = str(config.get('global', 'workdir'))
        self.upstream_path = workdir
        self.semester_directory = self.upstream_path
        self.current_day = str(config.get('global', 'current_day'))
        self.output_directory = self.upstream_path + "outputs/"
        self.reports_directory = self.upstream_path + "outputs/"

        # Get night plan specific parameters
        self.max_solve_gap = config.getfloat('night', 'max_solve_gap')
        self.max_solve_time = config.getint('night', 'max_solve_time')
        self.show_gurobi_output = config.getboolean('night', 'show_gurobi_output')
        
        # Set up allocation file path from data section
        allocation_file_config = str(config.get('data', 'allocation_file'))
        if os.path.isabs(allocation_file_config):
            self.allocation_file = allocation_file_config
        else:
            self.allocation_file = os.path.join(self.semester_directory, allocation_file_config)
            
        # Set up backup file path
        # DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),'data')
        # self.backup_file = os.path.join(DATADIR, "bright_backups_frame.csv")
        filler_file_config = str(config.get('data', 'filler_file'))
        if os.path.isabs(filler_file_config):
            self.filler_file = filler_file_config
        else:
            self.filler_file = os.path.join(self.semester_directory, filler_file_config)

        
        # Set up custom file path from data section
        custom_file_config = str(config.get('data', 'custom_file'))
        if os.path.isabs(custom_file_config):
            self.custom_file = custom_file_config
        else:
            self.custom_file = os.path.join(self.semester_directory, custom_file_config)
        
        # Load SemesterPlanner from pickle file instead of creating new one
        config = ConfigParser()
        config.read(config_file)
        workdir = str(config.get('global', 'workdir')) + "/outputs/"
        semester_planner_pkl = os.path.join(workdir, 'semester_planner.pkl')
        
        with open(semester_planner_pkl, 'rb') as f:
            self.semester_planner = pickle.load(f)
        
        # Pull properties from SemesterPlanner for consistency
        self.semester_start_date = self.semester_planner.semester_start_date
        self.semester_length = self.semester_planner.semester_length
        self.all_dates_dict = self.semester_planner.all_dates_dict
        self.all_dates_array = self.semester_planner.all_dates_array
        self.today_starting_night = self.semester_planner.today_starting_night
        self.past_history = self.semester_planner.past_history
        self.slots_needed_for_exposure_dict = self.semester_planner.slots_needed_for_exposure_dict
        self.run_weather_loss = self.semester_planner.run_weather_loss
        
    def run_ttp(self):
        """
        Produce the TTP solution given the results of the autoscheduler.
        Optimizes the nightly observation sequence for scheduled targets.
        """
        observers_path = self.semester_directory + 'outputs/'
        check1 = os.path.isdir(observers_path)
        if not check1:
            os.makedirs(observers_path)

        observatory = telescope.Keck1()
        # Get start/stop times from allocation file
        observation_start_time, observation_stop_time = get_nightly_times_from_allocation(self.allocation_file, self.current_day)
        total_time = np.round((observation_stop_time.jd-observation_start_time.jd)*24,3)
        print("Time in Night for Observations: " + str(total_time) + " hours.")

        # Use only request_selected.csv as the source of scheduled targets
        selected_path = os.path.join(self.output_directory, 'request_selected.csv')
        if not os.path.exists(selected_path):
            raise FileNotFoundError(f"{selected_path} not found. Please run the scheduler first.")
        selected_df = pd.read_csv(selected_path)
        # Fill NaN values with defaults --- for now in early 2025B since we had issues with the webform.c
        # Replace "None" strings with NaN first, then fill with defaults
        selected_df['n_intra_max'] = selected_df['n_intra_max'].replace('None', np.nan).fillna(1)
        selected_df['n_intra_min'] = selected_df['n_intra_min'].replace('None', np.nan).fillna(1)
        selected_df['tau_intra'] = selected_df['tau_intra'].replace('None', np.nan).fillna(0.0)
        selected_df['jmag'] = selected_df['jmag'].replace('None', np.nan).fillna(0.0)
        selected_df['gmag'] = selected_df['gmag'].replace('None', np.nan).fillna(0.0)
        selected_df['pmra'] = selected_df['pmra'].replace('None', np.nan).fillna(0.0)
        selected_df['pmdec'] = selected_df['pmdec'].replace('None', np.nan).fillna(0.0)
        selected_df['epoch'] = selected_df['epoch'].replace('None', np.nan).fillna(0.0)

        # Prepare the TTP input DataFrame (matching the old prepare_for_ttp output)
        to_ttp = pd.DataFrame({
            "Starname": selected_df["unique_id"],
            "RA": selected_df["ra"],
            "Dec": selected_df["dec"],
            "Exposure Time": selected_df["exptime"],
            "Exposures Per Visit": selected_df["n_exp"],
            "Visits In Night": selected_df["n_intra_max"],
            "Intra_Night_Cadence": selected_df["tau_intra"],
            "Priority": 10  # Default priority, or you can add logic if needed
        })

        filename = os.path.join(self.output_directory, 'request_selected.txt')
        to_ttp.to_csv(filename, index=False)
        target_list = formatting.theTTP(filename)

        solution = model.TTPModel(observation_start_time, observation_stop_time, target_list,
                                    observatory, observers_path, runtime=self.max_solve_time, optgap=self.max_solve_gap, useHighEl=False)

        gurobi_model_backup = solution.gurobi_model  # backup the attribute, probably don't need this
        del solution.gurobi_model                   # remove attribute so pickle works
    
        # add human readable starname to the solution so that it can be used in the plotting functions
        id_to_name = dict(zip(selected_df['unique_id'], selected_df['starname']))
        solution.plotly['human_starname'] = [
            id_to_name.get(uid, "NO MATCHING NAME") for uid in solution.plotly['Starname']
        ]
        self.solution = [solution]

        observe_order_file = os.path.join(observers_path, f"ObserveOrder_{self.current_day}.txt")
        # Convert solution.plotly to a DataFrame for easier handling
        plotly_df = pd.DataFrame(solution.plotly)
        use_starnames = []
        use_star_ids = []
        use_start_exposures = []
        for i in range(len(plotly_df)):
            adjusted_timestamp = TimeDelta(plotly_df['Start Exposure'].iloc[i]*60,format='sec') + observation_start_time
            use_start_exposures.append(str(adjusted_timestamp)[11:16])            
            use_starnames.append(selected_df[selected_df['unique_id'] == plotly_df['Starname'].iloc[i]]['starname'].iloc[0])
            use_star_ids.append(str(plotly_df['Starname'].iloc[i]))
        # Convert solution.extras to a DataFrame for consistency
        extras_df = pd.DataFrame(solution.extras)

        for j in range(len(extras_df)):
            use_start_exposures.append('24:00')
            use_star_ids.append(str(extras_df['Starname'].iloc[j]))
            use_starnames.append(selected_df[selected_df['unique_id'] == extras_df['Starname'].iloc[j]]['starname'].iloc[0])
        use_frame = pd.DataFrame({'unique_id': use_star_ids, 'Target': use_starnames, 'StartExposure': use_start_exposures})
        use_frame.to_csv(observe_order_file, index=False)

        # Check if the file was created before trying to read it
        if os.path.exists(observe_order_file):
            obs_and_times = pd.read_csv(observe_order_file)
        else:
            print(f"Warning: {observe_order_file} was not created by writeStarList")
            obs_and_times = pd.DataFrame()  # Create empty DataFrame as fallback
        io.write_starlist(selected_df, solution.plotly, observation_start_time, solution.extras,
                            [], str(self.current_day), observers_path)
        print("The optimal path through the sky for the selected stars is found. Clear skies!")

        return #save_data

def get_nightly_times_from_allocation(allocation_file, current_day):
    """
    Extract start and stop times for a specific date from allocation.csv.
    
    Args:
        allocation_file (str): path to the allocation file
        current_day (str): the date to look for in YYYY-MM-DD format
        
    Returns:
        tuple: (start_time, stop_time) as Time objects
    """
    allocated_times_frame = pd.read_csv(allocation_file)
    allocated_times_frame['start'] = allocated_times_frame['start'].apply(Time)
    allocated_times_frame['stop'] = allocated_times_frame['stop'].apply(Time)
    
    # Filter for the current day
    current_day_str = str(current_day)
    day_allocations = []
    for _, row in allocated_times_frame.iterrows():
        start_datetime = str(row['start'])[:10]  # Extract date part (YYYY-MM-DD)
        if start_datetime == current_day_str:
            day_allocations.append(row)
    
    if not day_allocations:
        raise ValueError(f"No allocation found for date {current_day_str}")
    
    # For multiple allocations on the same day, use the earliest start and latest stop
    start_times = [row['start'] for row in day_allocations]
    stop_times = [row['stop'] for row in day_allocations]
    
    earliest_start = min(start_times)
    latest_stop = max(stop_times)
    
    return earliest_start, latest_stop
        