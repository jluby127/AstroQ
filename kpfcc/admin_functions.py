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
import kpfcc.helper_functions as hf
import kpfcc.processing_functions as pf

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
        self.database_info_dict = pf.build_past_history(past_observations_file, self.requests_frame, self.twilight_frame)
        self.slots_needed_for_exposure_dict = hf.build_slots_required_dictionary(self.requests_frame, self.slot_size)

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
            hf.get_semester_info(self.current_day)
        all_dates_dict = hf.build_date_dictionary(semester_start_date, semester_length)
        all_dates_array = list(all_dates_dict.keys())
        n_nights_in_semester = hf.current_day_tracker(self.current_day, all_dates_dict)
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
