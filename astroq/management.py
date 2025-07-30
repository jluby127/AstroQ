"""
Module defining the admin object, which stores, manipulates, and passes information across
all requests easily.

Example usage:
    import admin_functions as af
"""
import os
import math

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta
from configparser import ConfigParser
import logging
logs = logging.getLogger(__name__)

# from kpfcc import DATADIR
DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),'data')
EXAMPLEDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),'examples')

import astroq.history as hs

class data_admin(object):
    """A Data Admin object, from which we can easily pass around information.

    Args:
        - config_path (str) = the path and file name to the config.ini file.
        - current_day (str) = the calendar date of the night to produce a script. Format: YYYY-MM-DD.
    Returns:
        None
    """

    def __init__(self, config_path):

        config = ConfigParser()
        config.read(config_path)

        self.upstream_path = eval(config.get('required', 'folder'), {"os": os})

        self.current_day = str(config.get('required', 'current_day'))
        self.semester_directory = self.upstream_path
        self.reports_directory = self.upstream_path + 'outputs/'
        self.observatory = config.get('required', 'observatory')

        self.slot_size = int(config.get('other', 'slot_size'))
        self.n_quarters_in_night = int(config.get('other', 'quarters_in_night'))
        self.n_hours_in_night = int(config.get('other', 'hours_in_night'))
        self.daily_starting_time = str(config.get('other', 'daily_starting_time'))
        self.daily_ending_time  = f"{(int(self.daily_starting_time.split(':')[0]) + self.n_hours_in_night) % 24:02d}:{int(self.daily_starting_time.split(':')[1]):02d}"

        # Resolve allocation file path relative to semester directory
        allocation_file_config = str(config.get('options', 'allocation_file'))
        if os.path.isabs(allocation_file_config):
            self.allocation_file = allocation_file_config
        else:
            self.allocation_file = os.path.join(self.semester_directory, allocation_file_config)
        self.past_file = os.path.join(self.semester_directory, "inputs/past.csv")
        # Removed special_map_file, zero_out_file, nonqueue_map_file - these will be handled differently
        self.custom_file = os.path.join(self.semester_directory, "inputs/custom.csv")
        self.nonqueue_file = os.path.join(self.semester_directory, "inputs/NonQueueMap.csv")
        self.future_forecast = os.path.join(self.semester_directory, "outputs/semester_plan.csv")

        self.random_seed = int(config.get('options', 'random_seed'))
        np.random.seed(self.random_seed)
        self.run_weather_loss = eval(config.get('options', 'run_weather_loss'))#.strip().lower() == "true"

        self.DATADIR = DATADIR
        self.gurobi_output = config.get('gurobi', 'show_gurobi_output').strip().lower() == "true"
        self.plot_results = config.get('options', 'run_plots').strip().lower() == "true"
        self.run_plots = config.get('options', 'run_plots').strip().lower() == "true"
        self.solve_time_limit = int(config.get('gurobi', 'max_solve_time'))
        self.solve_max_gap = float(config.get('gurobi', 'max_solve_gap'))

        self.run_scheduler = config.get('options', 'run_scheduler').strip().lower() == "true"
        self.run_ttp = config.get('options', 'run_ttp').strip().lower() == "true"
        self.run_round_two = config.get('options', 'run_bonus_round').strip().lower() == "true"
        self.max_bonus = float(config.get('other', 'maximum_bonus_size'))

        self.folder_forecasts = os.path.join(self.semester_directory, "/data/first_forecasts/")
        self.run_backup_scripts = config.get('other', 'generate_backup_script').strip().lower() == "true"
        self.backup_file = os.path.join(DATADIR,"bright_backups_frame.csv")
        self.backup_observability_file = os.path.join(DATADIR,"bright_backup_observability.csv")

    def run_admin(self):
        """
        Given today's date, collate all important information about the semester.
        """
        # build out some paths here, so that if current_day changes due to request_set.json file, it is reflected properly
        self.requests_frame = pd.read_csv(os.path.join(self.semester_directory, "inputs/requests.csv"))

        self.output_directory = self.upstream_path + "outputs/"
        self.reports_directory = self.upstream_path + 'outputs/'
        self.folder_cadences = os.path.join(self.semester_directory, "/outputs/" + str(self.current_day) + "/cadences/")
        # Suggest your output directory be something so that it doesn't autosave
        # to the same directory as the run files and crowds up the GitHub repo.
        check = os.path.isdir(self.output_directory)
        if not check:
            os.makedirs(self.output_directory)
            file = open(self.output_directory + "runReport.txt", "w") # later append to this file
            file.close()

        # Get semester parameters and define important quantities
        self.past_history = hs.process_star_history(self.past_file)
        self.slots_needed_for_exposure_dict = self.build_slots_required_dictionary()

        # Calculate semester info based on current_day
        from datetime import datetime
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
        
        all_dates_dict, all_dates_array = build_date_dictionary(semester_start_date, semester_length)
        n_nights_in_semester = len(all_dates_dict) - all_dates_dict[self.current_day]
        logs.debug("Total semester length: ", semester_length)
        logs.debug("There are " + str(n_nights_in_semester) + " calendar nights remaining in the semester.")

        # Compute slot grid sizes
        n_slots_in_quarter = int(((self.n_hours_in_night*60)/self.n_quarters_in_night)/self.slot_size)
        n_slots_in_night = n_slots_in_quarter*self.n_quarters_in_night
        n_slots_in_semester = n_slots_in_night*n_nights_in_semester

        # Define the slot and night represents the today's date
        today_starting_slot = all_dates_dict[self.current_day]*n_slots_in_night
        today_starting_night =  all_dates_dict[self.current_day]
        logs.debug("There are " + str(n_slots_in_semester) + " slots remaining in the semester.")

        self.semester_start_date = semester_start_date
        self.semester_length = semester_length
        self.semester_letter = semester_letter
        self.all_dates_dict = all_dates_dict
        self.all_dates_array = all_dates_array
        self.n_slots_in_night = n_slots_in_night
        self.n_nights_in_semester = n_nights_in_semester
        self.n_slots_in_semester = n_slots_in_semester
        self.today_starting_slot = today_starting_slot
        self.today_starting_night = today_starting_night
        self.semester_grid = np.arange(0, self.n_nights_in_semester, 1)
        self.quarters_grid = np.arange(0, self.n_quarters_in_night, 1)
    def build_slots_required_dictionary(self, always_round_up_flag=False):
        """
        Computes the slots needed for a given exposure for all requests.
        When always_round_up_flag is false, we can round up or down.
        Example: with 5 minute slot sizes, we want a 6 minute exposure to only require one slot
                (round down) as opposed to 2 slots (round up)

        Args:
            manager (obj): a data_admin object
            always_round_up_flag (boolean): if true, slots needed is always larger than exposure_time
        """
        logs.info("Determining slots needed for exposures.")
        # schedule multi-shots and multi-visits as if a single, long exposure.
        # When n_shots and n_visits are both 1, this reduces down to just the stated exposure time.
        slots_needed_for_exposure_dict = {}
        for n,row in self.requests_frame.iterrows():
            name = row['starname']
            exposure_time = row['exptime']*row['n_exp'] + 45*(row['n_exp'] - 1)
            slots_needed_for_exposure_dict[name] = compute_slots_required_for_exposure(exposure_time, self.slot_size, always_round_up_flag)
        return slots_needed_for_exposure_dict



def build_date_dictionary(semester_start_date, semester_length):
    """
    Builds a dictionary where keys are the calendar dates within the semester and values are the
    corresponding day numbers of the semester.

    Args:
        semester_start_date (str): the first day of the semester, format YYYY-MM-DD
        semester_length (int): the number of days in the full semester
    """
    all_dates = {}
    date_formal = Time(semester_start_date, format='iso',scale='utc')
    date = str(date_formal)[:10]
    all_dates[date] = 0
    for i in range(1, semester_length):
        date_formal += TimeDelta(1,format='jd')
        date = str(date_formal)[:10]
        all_dates[date] = i
    return all_dates, list(all_dates.keys())

def compute_slots_required_for_exposure(exposure_time, slot_size, always_round_up_flag):
    slot_size = slot_size * 60 # converting to seconds
    if always_round_up_flag:
        slots_needed_for_exposure = math.ceil(exposure_time/slot_size)
    else:
        if exposure_time > slot_size:
            slots_needed_for_exposure = int(round(exposure_time/slot_size))
        else:
            slots_needed_for_exposure = 1
    return slots_needed_for_exposure

def get_gap_filler_targets(manager):
    """
    Using the results of the two rounds of scheduling, determine what is different between them,
    i.e. which targets were added in Round 2

    Args:
        manager (obj): a data_admin object

    Returns:
        None
    """
    round1_results = manager.output_directory + 'raw_combined_semester_schedule_Round1.txt'
    round2_results = manager.output_directory + 'raw_combined_semester_schedule_Round2.txt'
    if os.path.exists(round1_results) and os.path.exists(round2_results):
        scheduleR1 = np.loadtxt(round1_results, delimiter=',', dtype=str)
        scheduleR2 = np.loadtxt(round2_results, delimiter=',', dtype=str)
        new = scheduleR2[manager.all_dates_dict[manager.current_day]]
        old = scheduleR1[manager.all_dates_dict[manager.current_day]]
        gap_fillers = [x for x in new if x not in old]
    else:
        gap_fillers = []
    np.savetxt(manager.output_directory + 'Round2_Requests.txt', gap_fillers, delimiter=',', fmt="%s")

