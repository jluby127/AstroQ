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
        
        # Set up allocation file path from data section
        allocation_file_config = str(config.get('data', 'allocation_file'))
        if os.path.isabs(allocation_file_config):
            self.allocation_file = allocation_file_config
        else:
            self.allocation_file = os.path.join(self.semester_directory, allocation_file_config)
            
        # Set up backup file path
        DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),'data')
        self.backup_file = os.path.join(DATADIR, "bright_backups_frame.csv")
        
        # Set up custom file path from data section
        custom_file_config = str(config.get('data', 'custom_file'))
        if os.path.isabs(custom_file_config):
            self.custom_file = custom_file_config
        else:
            self.custom_file = os.path.join(self.semester_directory, custom_file_config)
        
        # Create SemesterPlanner to get all the helper methods and properties
        # This avoids duplicating semester calculation logic
        import astroq.splan as splan
        self.semester_planner = splan.SemesterPlanner(config_file)
        
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
            "Starname": selected_df["starname"],
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
                                    observatory, observers_path, runtime=10, optgap=0.01, useHighEl=False)

        gurobi_model_backup = solution.gurobi_model  # backup the attribute, probably don't need this
        del solution.gurobi_model                   # remove attribute so pickle works
        save_data = [solution]
        with open(self.reports_directory + 'ttp_data.pkl', 'wb') as f:
            pickle.dump(save_data, f)

        observe_order_file = os.path.join(observers_path,'night_plan.csv')
        observe_order_txt = os.path.join(observers_path)#, f"ObserveOrder_{self.current_day}.txt")
        plotting.writeStarList(solution.plotly, observation_start_time, self.current_day,observe_order_txt)
        # temporarily open the ObserveOrder file and write in the id numbers.
        obsord = pd.read_csv(observe_order_txt + f"ObserveOrder_{self.current_day}.txt")
        obsord['id'] = obsord['Target'].apply(lambda x: selected_df[selected_df['starname'] == x]['unique_id'].iloc[0] if x in selected_df['starname'].values else x)
        obsord.to_csv(observe_order_txt + f"ObserveOrder_{self.current_day}.txt", index=False)

        # plotting.plot_path_2D(solution,outputdir=observers_path)
        # plotting.nightPlan(solution.plotly, self.current_day, outputdir=observers_path)
        
        # Check if the file was created before trying to read it
        if os.path.exists(observe_order_file):
            obs_and_times = pd.read_csv(observe_order_file)
        else:
            print(f"Warning: {observe_order_file} was not created by writeStarList")
            obs_and_times = pd.DataFrame()  # Create empty DataFrame as fallback
        io.write_starlist(selected_df, solution.plotly, observation_start_time, solution.extras,
                            [], str(self.current_day), observers_path)
        print("The optimal path through the sky for the selected stars is found. Clear skies!")
        return save_data

    def produce_bright_backups(self, nstars_max=100):
        """
        Produce backup target lists for poor weather conditions.
        
        Args:
            nstars_max (int): Maximum number of backup stars to include
        """
        backups_path = self.semester_directory + 'outputs/'
        check = os.path.isdir(backups_path)
        if not check:
            os.makedirs(backups_path)

        # Get start/stop times from allocation file
        observation_start_time, observation_stop_time = get_nightly_times_from_allocation(self.allocation_file, self.current_day)
        diff_minutes = int(abs((observation_stop_time - observation_start_time).to('min').value))
        print("Minutes on sky: ", diff_minutes)

        backup_starlist = pd.read_csv(self.backup_file)
        self.requests_frame = backup_starlist
        # Create Access object with required parameters
        access_obj = ac.Access(
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
            output_directory=self.output_directory
        )
        available_indices = access_obj.produce_ultimate_map(self.requests_frame, running_backup_stars=True)
        slots_available_tonight_for_star = {k: len(v[0]) for k, v in available_indices.items()}
        stars_with_sufficient_availability_tonight = [k for k, v in slots_available_tonight_for_star.items() if v > int(0.25*int(diff_minutes/5))]

        self.requests_frame = backup_starlist
        isTonight = backup_starlist['starname'].isin(stars_with_sufficient_availability_tonight)
        hasDR3name = backup_starlist['gaia_id'].str.startswith('Gaia DR2')
        pool_tonight = self.requests_frame[isTonight&hasDR3name]
        pool_tonight = pool_tonight.sample(frac=1).reset_index(drop=True)
        pool_tonight = pool_tonight[:nstars_max]

        ready_for_ttp = self.prepare_for_ttp(pool_tonight, list(pool_tonight['starname']), [])
        ready_for_ttp.to_csv(backups_path + "selected_stars.csv", index=False)
        target_list = formatting.theTTP(backups_path + "selected_stars.csv")

        observatory = telescope.Keck1()
        solution_b = model.TTPModel(observation_start_time, observation_stop_time,
                                    target_list, observatory, backups_path,
                                    runtime=10, optgap=0.05)

        backup_order_file = os.path.join(backups_path, 'night_plan.csv')
        backup_order_txt = os.path.join(backups_path, f"ObserveOrder_{self.current_day}.txt")
        plotting.writeStarList(solution_b.plotly, observation_start_time, self.current_day,
                                    backup_order_txt)
        plotting.plot_path_2D(solution_b, outputdir = backups_path)
        plotting.nightPlan(solution_b.plotly, self.current_day, outputdir = backups_path)
        obs_and_times_b = pd.read_csv(backups_path + 'night_plan.csv')
        io.write_starlist(pool_tonight, solution_b.plotly, observation_start_time,
                            solution_b.extras, [], self.current_day, backups_path, "backups")
        print("Bright backups script created.")

    def prepare_for_ttp(self, request_frame, night_plan, round_two_targets):
        """
        Prepare tonight's scheduled stars for their run through the TTP.

        Args:
            request_frame (dataframe): the pandas dataframe of PI requests
            night_plan (array): the n'th row of combined_semester_schedule array
            round_two_targets (array): a 1D list of the stars that were added in the bonus round

        Returns:
            to_ttp (dataframe): the data on the stars to be observed tonight, formatted for TTP
        """
        ignore = ['*', 'W', '', '*X', 'X']
        selected_stars = []
        for i, item in enumerate(night_plan):
            if night_plan[i] not in ignore and night_plan[i][:4] != "RM___":
                selected_stars.append(night_plan[i])
                ignore.append(night_plan[i])

        starnames = []
        ras = []
        decs = []
        exposure_times = []
        exposures_per_visit = []
        visits_in_night = []
        cadences = []
        priorities = []
        for j, item in enumerate(selected_stars):
            idx = request_frame.index[request_frame['starname']==str(selected_stars[j])][0]
            starnames.append(str(request_frame['starname'][idx]))
            ras.append(request_frame['ra'][idx])
            decs.append(request_frame['dec'][idx])
            exposure_times.append(int(request_frame['exptime'][idx]))
            exposures_per_visit.append(int(request_frame['n_exp'][idx]))
            visits_in_night.append(int(request_frame['n_intra_max'][idx]))
            cadences.append(int(request_frame['tau_intra'][idx]))
            # higher numbers are higher priorities, filler targets get low priority
            if str(selected_stars[j]) in round_two_targets:
                prior = 1
            else:
                prior = 10
            priorities.append(prior)
        to_ttp = pd.DataFrame({"Starname":starnames,"RA":ras,"Dec":decs,
                              "Exposure Time":exposure_times,
                              "Exposures Per Visit":exposures_per_visit,
                              "Visits In Night":visits_in_night, "Intra_Night_Cadence":cadences,
                              "Priority":priorities})
        return to_ttp




# Legacy function wrappers for backwards compatibility
def run_ttp(config_file):
    """
    Legacy wrapper function for backwards compatibility.
    Use NightPlanner class for new code.
    
    Args:
        config_file: Path to configuration file
    """
    planner = NightPlanner(config_file)
    planner.run_ttp()


def produce_bright_backups(config_file, nstars_max=100):
    """
    Legacy wrapper function for backwards compatibility.
    Use NightPlanner class for new code.
    
    Args:
        config_file: Path to configuration file
        nstars_max: Maximum number of backup stars to include
    """
    planner = NightPlanner(config_file)
    planner.produce_bright_backups(nstars_max)


def prepare_for_ttp(request_frame, night_plan, round_two_targets):
    """
    Legacy wrapper function for backwards compatibility.
    Use NightPlanner class for new code.
    """
    # This function is stateless, so we can just call it directly
    # Create a temporary planner instance
    temp_planner = NightPlanner(None)
    return temp_planner.prepare_for_ttp(request_frame, night_plan, round_two_targets)

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
        