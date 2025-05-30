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
from types import SimpleNamespace

# from kpfcc import DATADIR
DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),'data')
import kpfcc.history as hs
import kpfcc.access as ac
import kpfcc.maps as mp

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
        self.reports_directory = self.upstream_path + 'reports/'
        self.observatory = config.get('required', 'observatory')

        self.slot_size = int(config.get('other', 'slot_size'))
        self.n_quarters_in_night = int(config.get('other', 'quarters_in_night'))
        self.n_hours_in_night = int(config.get('other', 'hours_in_night'))
        self.daily_starting_time = str(config.get('other', 'daily_starting_time'))
        self.daily_ending_time  = f"{(int(self.daily_starting_time.split(':')[0]) + self.n_hours_in_night) % 24:02d}:{int(self.daily_starting_time.split(':')[1]):02d}"

        # self.allocation_file = os.path.join(self.semester_directory, "inputs/allocation_schedule.txt")
        self.allocation_file = str(config.get('oia', 'allocation_file'))
        print("Using allocation map as defined in: ", self.allocation_file)
        self.requests_frame = pd.read_csv(os.path.join(self.semester_directory, "inputs/Requests.csv"))
        self.twilight_frame = pd.read_csv(os.path.join(self.semester_directory, "inputs/twilight_times.csv"), parse_dates=True)
        self.past_database_file = os.path.join(self.semester_directory, "inputs/queryJumpDatabase.csv")
        self.accessibilities_file = os.path.join(self.semester_directory, "inputs/accessibilities_" + str(self.slot_size) + "minSlots.json")
        self.special_map_file = os.path.join(self.semester_directory, "inputs/specialMaps_" + str(self.slot_size) + "minSlots.txt")
        self.zero_out_file = os.path.join(self.semester_directory, "inputs/zero_out.csv")
        self.nonqueue_map_file = os.path.join(self.semester_directory, "inputs//NonQueueMap"  + str(self.slot_size) + ".txt")
        self.nonqueue_file = os.path.join(self.semester_directory, "inputs/NonQueueMap.csv")

        self.random_seed = int(config.get('options', 'random_seed'))
        np.random.seed(self.random_seed)
        self.run_weather_loss = eval(config.get('options', 'run_weather_loss'))#.strip().lower() == "true"

        self.run_optimal_allocation = config.get('oia', 'run_optimal_allocation').strip().lower() == "true"
        self.include_aesthetic = config.get('oia', 'run_with_aesthetics').strip().lower() == "true"
        self.max_quarters = int(config.get('oia', 'maximum_allocated_quarters'))
        self.max_unique_nights = int(config.get('oia', 'maximum_allocated_nights'))
        self.min_represented = 5
        self.whiteout_file = os.path.join(self.semester_directory, "inputs/whiteout_dates.csv")
        self.blackout_file = os.path.join(self.semester_directory, "inputs/blackout_dates.csv")
        self.allow_single_quarters = config.get('oia', 'allow_single_quarter_allocations').strip().lower() == "true"
        self.max_consecutive = int(config.get('oia', 'maximum_consecutive_onsky'))
        self.min_consecutive = int(config.get('oia', 'minimum_consecutive_offsky'))
        self.max_baseline = int(config.get('oia', 'maximum_baseline'))

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
        self.turn_on_off_file = os.path.join(self.semester_directory, "inputs/turnOnOffDates.csv")
        self.starmap_template_filename = os.path.join(self.semester_directory, "inputs/cadenceTemplateFile.csv")
        self.build_starmaps = config.get('other', 'build_starmaps').strip().lower() == "true"

        self.nightly_start_stop_times_file = os.path.join(self.semester_directory, "inputs/nightly_start_stop_times.csv")
        self.run_backup_scripts = config.get('other', 'generate_backup_script').strip().lower() == "true"
        self.backup_file = os.path.join(DATADIR,"bright_backups_frame.csv")
        self.backup_observability_file = os.path.join(DATADIR,"bright_backup_observability.csv")

    def run_admin(self):
        """
        Given today's date, collate all important information about the semester.
        """
        # build out some paths here, so that if current_day changes due to request_set.json file, it is reflected properly
        self.output_directory = self.upstream_path  + "outputs/" + str(self.current_day) + "/"
        self.folder_cadences = os.path.join(self.semester_directory, "/outputs/" + str(self.current_day) + "/cadences/")
        self.future_forecast = os.path.join(self.semester_directory, "outputs/" + str(self.current_day) + "/raw_combined_semester_schedule_Round2.txt")
        # Suggest your output directory be something so that it doesn't autosave
        # to the same directory as the run files and crowds up the GitHub repo.
        check = os.path.isdir(self.output_directory)
        if not check:
            os.makedirs(self.output_directory)
        file = open(self.output_directory + "runReport.txt", "w") # later append to this file
        file.close()

        # Get semester parameters and define important quantities
        self.database_info_dict = hs.build_past_history(self.past_database_file, self.requests_frame, self.twilight_frame)
        self.slots_needed_for_exposure_dict = self.build_slots_required_dictionary()

        semester_start_date, semester_end_date, semester_length, semester_year, semester_letter = \
            get_semester_info(self.current_day)
        all_dates_dict, all_dates_array = build_date_dictionary(semester_start_date, semester_length)
        n_nights_in_semester = len(all_dates_dict) - all_dates_dict[self.current_day]
        print("Total semester length: ", semester_length)
        print("There are " + str(n_nights_in_semester) + " calendar nights remaining in the semester.")

        # Compute slot grid sizes
        n_slots_in_quarter = int(((self.n_hours_in_night*60)/self.n_quarters_in_night)/self.slot_size)
        n_slots_in_night = n_slots_in_quarter*self.n_quarters_in_night
        n_slots_in_semester = n_slots_in_night*n_nights_in_semester

        # Define the slot and night represents the today's date
        today_starting_slot = all_dates_dict[self.current_day]*n_slots_in_night
        today_starting_night =  all_dates_dict[self.current_day]
        print("There are " + str(n_slots_in_semester) + " slots remaining in the semester.")

        self.semester_start_date = semester_start_date
        self.semester_length = semester_length
        self.all_dates_dict = all_dates_dict
        self.all_dates_array = all_dates_array
        self.n_slots_in_night = n_slots_in_night
        self.n_nights_in_semester = n_nights_in_semester
        self.n_slots_in_semester = n_slots_in_semester
        self.today_starting_slot = today_starting_slot
        self.today_starting_night = today_starting_night
        self.semester_grid = np.arange(0, self.n_nights_in_semester, 1)
        self.quarters_grid = np.arange(0, self.n_quarters_in_night, 1)

        twilight_map_remaining_2D, available_slots_in_each_night = ac.construct_twilight_map(self)
        self.twilight_map_remaining_2D = twilight_map_remaining_2D
        self.available_slots_in_each_night = available_slots_in_each_night
        # When running a normal schedule, include the observatory's allocation map
        if self.run_optimal_allocation == False:
                weather_diff_remaining, allocation_map_1D, allocation_map_2D, weathered_map = \
                                    mp.prepare_allocation_map(self)
        else:
            # When running the optimal allocation, all dates are possible except for those specifically blacked out
            # Weather arrays are zeros since no weather losses are modeled
            weather_diff_remaining = np.zeros(self.n_nights_in_semester, dtype='int')
            weathered_map = np.zeros((self.n_nights_in_semester, self.n_slots_in_night), dtype='int')

            # allocation maps are ones because all nights are possible to be allocated
            allocation_map_1D = np.ones(self.n_slots_in_semester, dtype='int')
            allocation_map_2D = np.ones((self.n_nights_in_semester, self.n_slots_in_night), dtype='int')

        self.available_slots_in_each_night = available_slots_in_each_night
        self.weather_diff_remaining = weather_diff_remaining
        self.allocation_map_1D = allocation_map_1D
        self.allocation_map_2D = allocation_map_2D
        self.weathered_map = weathered_map
        if np.sum(self.weather_diff_remaining) == 0:
            self.weathered_days = []
        else:
            self.weathered_days = np.where(np.any(self.weather_diff_remaining == 1, axis=1))[0]

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
        print("Determining slots needed for exposures.")
        # schedule multi-shots and multi-visits as if a single, long exposure.
        # When n_shots and n_visits are both 1, this reduces down to just the stated exposure time.
        slots_needed_for_exposure_dict = {}
        for n,row in self.requests_frame.iterrows():
            name = row['starname']
            exposure_time = row['exptime']*row['n_exp'] + \
                45*(row['n_exp'] - 1)
            slots_needed_for_exposure_dict[name] = compute_slots_required_for_exposure(exposure_time, self.slot_size, always_round_up_flag)
        return slots_needed_for_exposure_dict

def get_semester_info(current_day):
    """
    Given today's date, return information about the semester we are currently in.

    Args:
        manager (obj): a data_admin object
    """
    year_flag = False
    # "A" semester runs from Feb 01 through July 31
    if current_day[5:7] in ['02', '03', '04', '05', '06', '07']:
        semester_letter = 'A'
    # "B" semester runs from Aug 1 through Jan 01
    elif current_day[5:7] in ['08', '09', '10', '11', '12', '01']:
        semester_letter = 'B'
        if current_day[5:7] == '01':
            year_flag = True
    else:
        print("Invalid date. Exiting.")
        return None
    semester_year = current_day[:4]

    if semester_letter == 'A':
        semester_start_date = semester_year + '-02-01'
        semester_end_date = semester_year + '-07-31'
        # check if this is a leap year
        this_year = 2024
        year_limit = 2074
        # Note from Jack Lubin in the year 2024: The year 2074 is arbitrary. In this year,
        # you, the current queue manager, will have to update this line for another 50 years.
        # I could have extended this thousands of years in the future, but thought it would be more
        # fun if one day this line breaks, and every 50 years and someone manually updates it.
        # If/when you need to update it, please send me an email because I'd like to know that Keck
        # is still using this software! Also when you update the line, please sign the list of
        # queue managers below:
        #
        # --------------------------
        # KPF-CC Queue Managers:
        # --------------------------
        # 2024 - Jack Lubin
        # 2074 - your_name_here
        # --------------------------
        if int(semester_year) > year_limit:
            print("Time to update the leap year array!!! See line 255 in management.py!!!")
            semester_length = 0
        elif int(semester_year) in np.arange(this_year, year_limit, 4):
            semester_length = 182
        else:
            semester_length = 181
    elif semester_letter == 'B':
        if year_flag:
            semester_start_date = str(int(semester_year) - 1) + '-08-01'
            semester_end_date = semester_year + '-01-31'
        else:
            semester_start_date = semester_year + '-08-01'
            semester_end_date = str(int(semester_year) + 1) + '-01-31'
        semester_length = 184
    else:
        print("Unrecognized semester letter designation!")
        semester_start_date = semester_year + '-01-01'
        semester_end_date = str(semester_year)+ '-01-01'
        semester_length = 0
    return semester_start_date, semester_end_date, semester_length, semester_year, semester_letter

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
    # manager.all_dates_dict = all_dates
    # manager.all_dates_array = list(all_dates.keys())

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


def prepare_new_semester(config_path):

    config = ConfigParser()
    config.read(config_path)

    little_manager = SimpleNamespace()


    # Set up important variables
    # -----------------------------------------------------------------------------------------
    little_manager.current_day = str(config.get('required', 'current_day'))
    little_manager.run_optimal_allocation = config.get('oia', 'run_optimal_allocation').strip().lower() == "true"

    little_manager.upstream_path = eval(config.get('required', 'folder'), {"os": os})
    little_manager.observatory = config.get('required', 'observatory')

    little_manager.slot_size = int(config.get('other', 'slot_size'))
    little_manager.n_quarters_in_night = int(config.get('other', 'quarters_in_night'))
    little_manager.n_hours_in_night = int(config.get('other', 'hours_in_night'))
    little_manager.daily_starting_time = str(config.get('other', 'daily_starting_time'))
    little_manager.daily_ending_time  = f"{(int(little_manager.daily_starting_time.split(':')[0]) + little_manager.n_hours_in_night) % 24:02d}:{int(little_manager.daily_starting_time.split(':')[1]):02d}"

    print(little_manager.upstream_path)
    print(os.path.join(little_manager.upstream_path, "inputs/Requests.csv"))
    little_manager.requests_frame = pd.read_csv(os.path.join(little_manager.upstream_path, "inputs/Requests.csv"))
    try:
        little_manager.nonqueue_frame = pd.read_csv(os.path.join(little_manager.upstream_path, "inputs/NonQueueMap"  + str(little_manager.slot_size) + ".csv"))
    except:
        print("There are no times reserved for non-queue observations.")
        little_manager.nonqueue_frame = None

    little_manager.semester_start_date, little_manager.semester_end_date, little_manager.semester_length, little_manager.semester_year, little_manager.semester_letter = get_semester_info(little_manager.current_day)
    little_manager.all_dates_dict, little_manager.all_dates_array = build_date_dictionary(little_manager.semester_start_date, little_manager.semester_length)
    little_manager.n_nights_in_semester = len(little_manager.all_dates_dict) - little_manager.all_dates_dict[little_manager.current_day]

    # Create the template file for the cadence plots including true weather map
    # -----------------------------------------------------------------------------------------
    print("Generate the cadence plot template file.")
    dateslist = []
    quarterlist = []
    for d in range(len(little_manager.all_dates_array)):
        for q in range(4):
            dateslist.append(little_manager.all_dates_array[d])
            quarterlist.append(0.5+q)
    falselist = [False]*len(dateslist)
    template_frame = pd.DataFrame({'Date':dateslist, 'Quarter':quarterlist,'Allocated':falselist,'Weathered':falselist})
    template_frame.to_csv(little_manager.upstream_path + 'inputs/cadenceTemplateFile.csv', index=False)

    # Create the template file for the cadence plots including true weather map
    # -----------------------------------------------------------------------------------------
    print("Computing twilight times for the semester.")
    twilight_frame = ac.generate_twilight_times(little_manager.all_dates_array)
    twilight_frame.to_csv(little_manager.upstream_path + 'inputs/twilight_times.csv', index=False)

    # Create the json file containing the accessiblity for each target (elevation + moon safe)
    # -----------------------------------------------------------------------------------------
    print("Computing access maps for all stars.")
    print("This could take some time.")
    little_manager.accessibilities_file = os.path.join(little_manager.upstream_path, "inputs/accessibilities_" + str(little_manager.slot_size) + "minSlots.json")
    default_access_maps = ac.construct_access_dict(little_manager)

    # Create the csv file containing the turn on and turn off dates for each target in each quarter
    # -----------------------------------------------------------------------------------------
    ########### code goes here, use the separation slots between quarters and sum access maps

    if little_manager.run_optimal_allocation == False:

        little_manager.allocation_file = os.path.join(little_manager.upstream_path, "inputs/Observatory_Allocation.csv")
        if os.path.exists(little_manager.allocation_file):
            allocation = mp.format_keck_allocation_info(little_manager.allocation_file)
            allocation_binary = mp.convert_allocation_info_to_binary(little_manager, allocation)
        else:
            print("No Keck Observatory instruments schedule found. See https://www2.keck.hawaii.edu/observing/keckSchedule/queryForm.php to download and then try again.")
