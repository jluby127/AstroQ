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

from kpfcc import DATADIR
import kpfcc.history as hs
import kpfcc.access as ac
import kpfcc.maps as mp

class data_admin(object):
    """A Data Admin object, from which we can easily pass around information.

    Args:
        - current_day (str) = the calendar date of the night to produce a script.
                              Sets the "first" day of the semester from which to compute the
                            semester schedule solution from this day forward. Format: YYYY-MM-DD.
        - requests_file (str) = the path and file name to the CSV with all the PI requests.
                                Confirm that column names are correct.
        - allocation_file (str) = the path and file name to the binary map of allocated nights.
        - accessibilities_file (str) = the path and file name to the pickle file containing a
                                       dictionary of target names and associated pre-computed 1D
                                       accessibility maps of length equal to n_slots_in_semester.
        - twilight_file (str) = the path and file name to the CSV with precomputed twilight times.
        - output_directory (str) = the path where all outputs of this function should be saved.
                                    It is recommended that the path be outside the git repo.
        - slot_size (int) = the time, in minutes, for a single slot.
        - run_round_two (boolean) = when True, run the bonus round.
                                    When False, do not run the bonus round.
        - past_observations_file (str) = the path and file name of the CSV containing information
                                         on all previous observations in the semester. If file
                                         does not exist, then we are ignoring prior observations.
        - semester_template_file (str) = the path and file name of the CSV containing the visual
                                         template of the semester. For plotting purposes only.
        - turn_off_on_file (str) = the path and file name of the CSV containing the pre-computed
                                   first and last day of accessiblity for each target.
                                   For plotting purposes only.
        - nonqueue_map_file (str) = the path and file name of the CSV containining a grid of
                                    n_nights_in_semester by n_slots_in_night elements where slots
                                    reserved for non-queue observations are filled with target name.
        - special_map_file (str) = the path and file name of the CSV containining a grid of
                                    n_nights_in_semester by n_slots_in_night elements which contains
                                    information on the custom set of slots a request can be
                                    scheduled into for various reasons of the PI
        - zero_out_file (str) = the path and file name of list of stars that cannot be scheduled
                                tonight for any reason. Often this is empty.
        - run_weather_loss (boolean) = if False, then no nights are lost to weather.
        - run_optimal_allocation (boolean) = if True, then run in optimal allocation mode.
        - include_aesthetic (boolean) = if True, then run optimal allocation with aesthetic constraints.
        - max_quarters (int) = the maximum number of quarters that can be allocated in optimal allocation
        - max_nights (int) = the maximum number of nights on sky that can be allocated in optimal allocation
        - whiteout_file (str) = the path and filename to the list of whiteout dates.
        - blackout_file (str) = the path and filename to the list of blackout dates.
        - gurobi_output (boolean) = a flag to turn off or on the feature of Gurobi printing
                                    to the terminal as it solves the model.
        - plot_results (boolean) = a flag to turn off or on the plotting outputs.
        - solve_time_limit (int) = the maximum time, in seconds, to allow Gurobi to solve the model.
    Returns:
        None
    """

    def __init__(self,
                output_directory,
                current_day,
                observatory,
                slot_size,
                requests_file,
                twilight_file,
                past_observations_file,
                allocation_file,
                accessibilities_file,
                special_map_file,
                zero_out_file,
                nonqueue_map_file,
                nonqueue_file,
                run_weather_loss,
                run_optimal_allocation,
                include_aesthetic,
                max_quarters,
                max_nights,
                min_represented,
                allow_single_quarters,
                max_consecutive,
                min_consecutive,
                max_baseline,
                whiteout_file,
                blackout_file,
                gurobi_output,
                plot_results,
                solve_time_limit,
                solve_max_gap,
                run_round_two,
                max_bonus,
                folder_forecasts,
                folder_cadences,
                turn_on_off_file,
                starmap_template_filename,
                future_forecast,
                build_starmaps,
                nightly_start_stop_times,
                backup_file,
                backup_observability_file,
                ):

        self.observatory = observatory
        self.current_day = current_day
        self.semester_directory = output_directory
        self.output_directory = output_directory  + "outputs/" + str(self.current_day) + "/"
        self.reports_directory = output_directory + 'reports/'

        # Suggest your output directory be something so that it doesn't autosave
        # to the same directory as the run files and crowds up the GitHub repo.
        check = os.path.isdir(self.output_directory)
        if not check:
            os.makedirs(self.output_directory)
        file = open(self.output_directory + "runReport.txt", "w") # later append to this file
        file.close()

        self.slot_size = slot_size
        self.n_quarters_in_night = 4
        self.n_hours_in_night = 14

        self.requests_frame = pd.read_csv(requests_file)
        self.twilight_frame = pd.read_csv(twilight_file, parse_dates=True)
        self.database_info_dict = hs.build_past_history(past_observations_file, self.requests_frame, self.twilight_frame)
        self.slots_needed_for_exposure_dict = self.build_slots_required_dictionary()

        self.past_database_file = past_observations_file
        self.allocation_file = allocation_file
        self.accessibilities_file = accessibilities_file
        self.special_map_file = special_map_file
        self.zero_out_file = zero_out_file
        self.nonqueue_map_file = nonqueue_map_file
        self.nonqueue_file = nonqueue_file

        self.run_weather_loss = run_weather_loss

        self.run_optimal_allocation = run_optimal_allocation
        self.include_aesthetic = include_aesthetic
        self.max_quarters = max_quarters
        self.max_unique_nights = max_nights
        self.min_represented = min_represented
        self.whiteout_file = whiteout_file
        self.blackout_file = blackout_file
        self.allow_single_quarters = allow_single_quarters
        self.max_consecutive = max_consecutive
        self.min_consecutive = min_consecutive
        self.max_baseline = max_baseline

        self.DATADIR = DATADIR
        self.gurobi_output = gurobi_output
        self.plot_results = plot_results
        self.solve_time_limit = solve_time_limit
        self.solve_max_gap = solve_max_gap

        self.run_round_two = run_round_two
        self.max_bonus = max_bonus

        self.folder_forecasts = folder_forecasts
        self.folder_cadences = folder_cadences
        self.turn_on_off_file = turn_on_off_file
        self.starmap_template_filename = starmap_template_filename
        self.future_forecast = future_forecast
        self.build_starmaps = build_starmaps

        self.nightly_start_stop_times = nightly_start_stop_times
        self.backup_file = backup_file
        self.backup_observability_file = backup_observability_file

    def run_admin(self):
        """
        Given today's date, collate all important information about the semester.
        """
        # Get semester parameters and define important quantities
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
            name = row['Starname']
            exposure_time = row['Nominal Exposure Time [s]']*row['# of Exposures per Visit'] + \
                45*(row['Desired Visits per Night'] - 1)
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
    scheduleR1 = np.loadtxt(manager.output_directory + 'raw_combined_semester_schedule_Round1.txt',
        delimiter=',', dtype=str)
    scheduleR2 = np.loadtxt(manager.output_directory + 'raw_combined_semester_schedule_Round2.txt',
        delimiter=',', dtype=str)
    new = scheduleR2[manager.all_dates_dict[manager.current_day]]
    old = scheduleR1[manager.all_dates_dict[manager.current_day]]
    gap_fillers = [x for x in new if x not in old]
    np.savetxt(manager.output_directory + 'Round2_Requests.txt', gap_fillers, delimiter=',', fmt="%s")
