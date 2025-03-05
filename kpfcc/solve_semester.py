"""
Module for running the semester solver.

This module organizes, builds, and solves the gurobi model. It also writes and organizes the
outputs of the schedule into both machine human readable forms. It is designed to produce info
that is in correct formatting for running the TTP. Further designed to be only run as a function
call from the generateScript.py script.

Example usage:
    import solve_semester as ssm
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
import kpfcc.constraint_functions as cf
import kpfcc.admin_functions as af

# import line_profiler
# @profile
def run_kpfcc(current_day,
               requests_file,
               allocation_file,
               accessibilities_file,
               twilight_file,
               output_directory,
               slot_size,
               run_round_two,
               past_observations_file,
               semester_template_file,
               turn_off_on_file,
               nonqueue_map_file,
               special_map_file,
               zero_out_file,
               run_weather_loss,
               optional_allo = False,
               include_aesthetic = False,
               gurobi_output = True,
               plot_results = True,
               solve_time_limit = 300,
               whiteout_file = 'nofilename.txt',
               blackout_file = 'nofilename.txt'):

    """run_kpfcc
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
        - optional_allo (boolean) = if True, then run in optimal allocation mode. Default is False.
        - gurobi_output (boolean) = a flag to turn off or on the feature of Gurobi printing
                                    to the terminal as it solves the model.
        - plot_results (boolean) = a flag to turn off or on the plotting outputs.
        - solve_time_limit (int) = the maximum time, in seconds, to allow Gurobi to solve the model.
    Returns:
        None
    """
    start_the_clock = time.time()

    check = os.path.isdir(output_directory)
    if not check:
        os.makedirs(output_directory)

    # admin = af.data_admin(output_directory, current_day[0])

    # Set up logistics.
    # Only today is important to the semester solver.
    # Any additional dates in the array are for the TTP.
    current_day = current_day[0]
    semester_start_date, semester_length, all_dates_dict, all_dates_array, n_slots_in_night, \
            n_nights_in_semester, today_starting_slot, today_starting_night = af.run_admin(current_day, slot_size)
    n_slots_in_semester = n_nights_in_semester*n_slots_in_night

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Read in files and prep targets
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    print("Reading inputs and prepping files.")
    requests_frame = pd.read_csv(requests_file)
    twilight_frame = pd.read_csv(twilight_file, parse_dates=True)
    database_info_dict = pf.build_past_history(past_observations_file, requests_frame, twilight_frame)
    slots_needed_for_exposure_dict = hf.build_slots_required_dictionary(requests_frame, slot_size)
    default_access_maps = mf.construct_access_dict(accessibilities_file, requests_frame)
    custom_access_maps = mf.construct_custom_map_dict(special_map_file)
    zero_out_names = mf.construct_zero_out_arr(zero_out_file)
    twilight_map_remaining_flat, twilight_map_remaining_2D, available_slots_in_each_night \
                                = tw.construct_twilight_map(current_day, \
                                twilight_frame, slot_size, all_dates_dict, \
                                n_slots_in_night, n_nights_in_semester)

    nonqueue_map_file_slots_ints = mf.construct_nonqueue_arr(nonqueue_map_file, today_starting_slot)

    # When running a normal schedule, include the observatory's allocation map
    if optional_allo == False:
        weather_diff_remaining, allocation_map_1D, allocation_map_2D, weathered_map = \
                    mf.prepare_allocation_map(allocation_file, current_day, semester_length, DATADIR, \
                    all_dates_dict, all_dates_array, run_weather_loss, n_slots_in_night,
                    available_slots_in_each_night, today_starting_night, output_directory)
        semester_grid = []
        quarters_grid = []
    # When running the optimal allocation, all dates are possible except for those specifically blacked out
    else:
        semester_grid = np.arange(0, n_nights_in_semester, 1)
        n_quarters_in_night = 4
        quarters_grid = np.arange(0, n_quarters_in_night, 1)

        # Weather arrays are zeros since no weather losses are modeled
        weather_diff_remaining = np.zeros(n_nights_in_semester, dtype='int')
        weathered_map = np.zeros((n_nights_in_semester, n_slots_in_night), dtype='int')

        # allocation maps are ones because all nights are possible to be allocated
        allocation_map_1D = np.ones(n_slots_in_semester, dtype='int')
        allocation_map_2D = np.ones((n_nights_in_semester, n_slots_in_night), dtype='int')

        # note here add blackout restrictions to cut down parameter space

    available_indices_for_request = mf.produce_ultimate_map(requests_frame, allocation_map_1D.flatten(),
                                                         twilight_map_remaining_flat,
                                                         default_access_maps, custom_access_maps,
                                                         zero_out_names, nonqueue_map_file_slots_ints,
                                                         today_starting_slot, database_info_dict,
                                                         slots_needed_for_exposure_dict,
                                                         all_dates_dict, current_day, slot_size,
                                                         n_slots_in_night, n_nights_in_semester,
                                                         n_slots_in_semester)

    Aframe, Aset, schedulable_requests = hf.define_slot_index_frame(requests_frame, slots_needed_for_exposure_dict,
                                        available_indices_for_request)

    model = cf.GorubiModel(Aset, Aframe, schedulable_requests, requests_frame, database_info_dict, \
                            n_nights_in_semester, slots_needed_for_exposure_dict,
                            solve_time_limit=solve_time_limit, gurobi_output=gurobi_output, run_optimal_allocation=optional_allo,
                            semester_grid=semester_grid, quarters_grid=quarters_grid,
                            max_quarters = 10, max_unique_nights = 10) #update with flag later

    model.constraint_build_theta()
    model.constraint_one_request_per_slot()
    model.constraint_reserve_multislot_exposures()
    model.constraint_max_visits_per_night()
    model.constraint_enforce_internight_cadence()

    if optional_allo:
        model.constraint_set_max_quarters_allocated()
        model.constraint_set_max_onsky_allocated()
        model.constraint_relate_allocation_and_onsky()
        model.constraint_all_portions_of_night_represented(min_represented=1)
        model.constraint_forbidden_quarter_patterns(allow_singles=False) #update with flag value
        # model.constraint_cannot_observe_if_not_allocated(available_slots_in_each_night, n_slots_in_night)
        model.constraint_cannot_observe_if_not_allocated(twilight_map_remaining_2D)
        if os.path.exists(blackout_file):
            constraint_enforce_restricted_nights(blackout_file, limit=0)
        if os.path.exists(whiteout_file):
            constraint_enforce_restricted_nights(whiteout_file, limit=1)
        if include_aesthetic:
            model.constraint_max_consecutive_onsky(max_consecutive=6) #update with flag value
            model.constraint_minimum_consecutive_offsky(min_consecutive=10) #update with flag value
            model.constraint_maximize_baseline(max_base=5) #update with flag value

    complete_constraints_build = time.time()
    print("Total Time to build constraints: ", np.round(complete_constraints_build-start_the_clock,3))
    model.set_objective_minimize_theta()
    model.solve_model()
    complete_round1_model = time.time()
    print("Total Time to finish solver: ", np.round(complete_round1_model-start_the_clock,3))

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Retrieve data from solution
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    file = open(output_directory + "runReport.txt", "w")
    file.close()

    if optional_allo:
        print("Reading results of optimal allocation map.")
        allocation_schedule_1d = []
        for v in model.Anq.values():
            if np.round(v.X,0) == 1:
                allocation_schedule_1d.append(1)
            else:
                allocation_schedule_1d.append(0)
        allocation_schedule = np.reshape(allocation_schedule_1d, (n_nights_in_semester, n_quarters_in_night))
        holder = np.zeros(np.shape(allocation_schedule))
        allocation_map_1D, allocation_map_2D, weathered_map = mf.build_allocation_map(allocation_schedule, holder, available_slots_in_each_night, n_slots_in_night)
        mf.convert_allocation_array_to_binary(allocation_schedule, all_dates_array, output_directory + "optimal_allocation_binary_schedule.txt")

    print("Building human readable schedule.")
    combined_semester_schedule_available = hf.write_available_human_readable(
                                all_dates_dict, current_day, semester_length,
                                n_nights_in_semester, n_slots_in_night, twilight_map_remaining_2D,
                                allocation_map_2D, weathered_map, nonqueue_map_file)
    np.savetxt(output_directory + 'raw_combined_semester_schedule_available.txt',
        combined_semester_schedule_available, delimiter=',', fmt="%s")

    combined_semester_schedule_stars = hf.write_stars_schedule_human_readable( \
            combined_semester_schedule_available, model.Yrds, list(requests_frame['Starname']),
            semester_length, n_slots_in_night, n_nights_in_semester,
            all_dates_dict, slots_needed_for_exposure_dict, current_day)
    np.savetxt(output_directory + 'raw_combined_semester_schedule_Round1.txt',
        combined_semester_schedule_stars, delimiter=',', fmt="%s")

    round = 'Round 1'
    rf.build_fullness_report(combined_semester_schedule_stars, allocation_map_2D, requests_frame,
                                slot_size, round, output_directory)

    print("Writing Report.")
    filename = open(output_directory + "runReport.txt", "a")
    theta_n_var = []
    counter = 0
    for v in model.theta.values():
        varname = v.VarName
        varval = v.X
        counter += varval
    print("Sum of Theta: " + str(counter))
    filename.write("Sum of Theta: " + str(counter) + "\n")

    plot_results = False
    if plot_results:
        print("Writing cadence plot files.")
        turn_on_off_frame = pd.read_csv(turn_off_on_file)
        all_starmaps = {}
        for i in range(len(requests_frame)):
            if database_info_dict != {}:
                starmap = rf.build_observed_map_past( \
                    database_info_dict[requests_frame['Starname'][i]], semester_template_file)
            else:
                starmap = rf.build_observed_map_past([[],[],[],[]], semester_template_file)

            starmap_updated = rf.build_observed_map_future(combined_semester_schedule_stars,
                                requests_frame['Starname'][i], starmap,
                                slots_needed_for_exposure_dict,
                                np.array(allocation_all).flatten(),
                                np.array(weather_diff_remaining).flatten(),
                                all_dates_dict[current_day])

            all_starmaps[requests_frame['Starname'][i]] = starmap_updated
            future_unique_days_forecasted = 0
            for k, item in enumerate(combined_semester_schedule_stars):
                if requests_frame['Starname'][i] in combined_semester_schedule_stars[k]:
                    future_unique_days_forecasted += 1

            try:
                past_unique_dates_for_star = database_info_dict[requests_frame['Starname'][i]][1]
            except:
                past_unique_dates_for_star = []
            rf.write_cadence_plot_file(requests_frame['Starname'][i], starmap_updated,
                                        turn_on_off_frame, requests_frame,
                                        future_unique_days_forecasted,
                                        past_unique_dates_for_star,
                                        current_day, output_directory)

    complete_round1_plots = time.time()
    print("Total Time to complete Round 1: " + \
        str(np.round(complete_round1_plots-start_the_clock,3)))
    filename.write("Total Time to complete Round 1: " + \
        str(np.round(complete_round1_plots-start_the_clock,3)) + "\n")
    filename.write("\n")
    filename.write("\n")
    filename.close()
    print("Round 1 complete.")


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Initiate Round 2 Scheduling
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------

    if run_round_two:
        print("Beginning Round 2 Scheduling.")
        model.constraint_fix_previous_objective()
        model.set_objective_maximize_slots_used()
        model.solve_model()


        complete_round2_model = time.time()
        print("Total Time to finish round 2 solver: ",
            np.round(complete_round2_model-start_the_clock,3))

        combined_semester_schedule_stars = hf.write_stars_schedule_human_readable(
                combined_semester_schedule_available, Yrs, requests_frame['Starname'],
                semester_length, n_slots_in_night, n_nights_in_semester,
                all_dates_dict, slots_needed_for_exposure_dict, current_day)
        np.savetxt(output_directory + 'raw_combined_semester_schedule_Round2.txt',
            combined_semester_schedule_stars, delimiter=',', fmt="%s")

        round = 'Round 2'
        rf.build_fullness_report(combined_semester_schedule_stars, allocation_map_2D,
                                    requests_frame, slot_size, round, output_directory)

        scheduleR1 = np.loadtxt(output_directory + 'raw_combined_semester_schedule_Round1.txt',
            delimiter=',', dtype=str)
        scheduleR2 = np.loadtxt(output_directory + 'raw_combined_semester_schedule_Round2.txt',
            delimiter=',', dtype=str)
        R2_requests = rf.get_gap_filler_targets(scheduleR1, scheduleR2, all_dates_dict[current_day])
        np.savetxt(output_directory + 'Round2_Requests.txt', R2_requests, delimiter=',', fmt="%s")

        filename = open(output_directory + "runReport.txt", "a")
        theta_n_var = []
        counter = 0
        for v in theta.values():
            varname = v.VarName
            varval = v.X
            counter += varval
        filename.write("Sum of Theta: " + str(counter) + "\n")

        complete_round2_plots = time.time()
        filename.write("Total Time to complete Round 2: " + \
            str(np.round(complete_round2_plots-start_the_clock,3)) + "\n")
        filename.write("\n")
        filename.write("\n")
        filename.close()
        print("Round 2 complete.")
    else:
        print("Not running Round 2. Duplicating Raw Schedule as dummy file.")
        np.savetxt(output_directory + 'raw_combined_semester_schedule_Round2.txt', \
            combined_semester_schedule_stars, delimiter=',', fmt="%s")
        R2_requests = []
        np.savetxt(output_directory + 'Round2_Requests.txt', R2_requests, delimiter=',', fmt="%s")

    print("The optimal semester schedule is found, clear skies!")
