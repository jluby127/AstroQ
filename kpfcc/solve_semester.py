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

import kpfcc.helper_functions as hf
import kpfcc.twilight_functions as tw
import kpfcc.reporting_functions as rf
import kpfcc.processing_functions as pf
import kpfcc.mapping_functions as mf
import kpfcc.constraint_functions as cf
import kpfcc.admin_functions as af

# import line_profiler
# @profile
def run_kpfcc(manager):
    """
    Run the core kpfcc algorithm

    Args:
        manager (obj): a data_admin object that has already run the function "run_admin()"

    Returns:
        None
    """
    start_the_clock = time.time()

    print("Reading inputs and prepping files.")
    default_access_maps = mf.construct_access_dict(manager.accessibilities_file, manager.requests_frame)
    custom_access_maps = mf.construct_custom_map_dict(manager.special_map_file)
    zero_out_names = mf.construct_zero_out_arr(manager.zero_out_file)
    twilight_map_remaining_flat, twilight_map_remaining_2D, available_slots_in_each_night \
                                 = tw.construct_twilight_map(manager)
    nonqueue_map_file_slots_ints = mf.construct_nonqueue_arr(manager.nonqueue_map_file, manager.today_starting_slot)

    # When running a normal schedule, include the observatory's allocation map
    if manager.run_optimal_allocation == False:
        weather_diff_remaining, allocation_map_1D, allocation_map_2D, weathered_map = \
                                    mf.prepare_allocation_map(manager, available_slots_in_each_night)
    # When running the optimal allocation, all dates are possible except for those specifically blacked out
    else:
        # Weather arrays are zeros since no weather losses are modeled
        weather_diff_remaining = np.zeros(manager.n_nights_in_semester, dtype='int')
        weathered_map = np.zeros((manager.n_nights_in_semester, manager.n_slots_in_night), dtype='int')

        # allocation maps are ones because all nights are possible to be allocated
        allocation_map_1D = np.ones(manager.n_slots_in_semester, dtype='int')
        allocation_map_2D = np.ones((manager.n_nights_in_semester, manager.n_slots_in_night), dtype='int')
    manager.weather_diff_remaining = weather_diff_remaining

    available_indices_for_request = mf.produce_ultimate_map(manager,
                                                            allocation_map_1D.flatten(),
                                                            twilight_map_remaining_flat,
                                                            default_access_maps,
                                                            custom_access_maps,
                                                            zero_out_names,
                                                            nonqueue_map_file_slots_ints)

    Aframe, Aset, schedulable_requests = hf.define_slot_index_frame(manager.requests_frame,
                                        manager.slots_needed_for_exposure_dict, available_indices_for_request)

    model = cf.GorubiModel(manager, Aset, Aframe, schedulable_requests)

    model.constraint_build_theta()
    model.constraint_one_request_per_slot()
    model.constraint_reserve_multislot_exposures()
    model.constraint_max_visits_per_night()
    model.constraint_enforce_internight_cadence()

    if manager.run_optimal_allocation:
        model.constraint_set_max_quarters_allocated()
        model.constraint_set_max_onsky_allocated()
        model.constraint_relate_allocation_and_onsky()
        model.constraint_all_portions_of_night_represented(min_represented=1)
        model.constraint_forbidden_quarter_patterns(manager.allow_single_quarters)
        model.constraint_cannot_observe_if_not_allocated(twilight_map_remaining_2D)
        if os.path.exists(manager.blackout_file):
            constraint_enforce_restricted_nights(manager.blackout_file, limit=0)
        if os.path.exists(manager.whiteout_file):
            constraint_enforce_restricted_nights(manager.whiteout_file, limit=1)
        if manager.run_with_aesthetic:
            model.constraint_max_consecutive_onsky(manager.max_consecutive)
            model.constraint_minimum_consecutive_offsky(manager.min_consecutive)
            model.constraint_maximize_baseline(manager.max_baseline)

    complete_constraints_build = time.time()
    print("Total Time to build constraints: ", np.round(complete_constraints_build-start_the_clock,3))
    model.set_objective_minimize_theta()
    model.solve_model()
    complete_round1_model = time.time()
    print("Total Time to finish solver: ", np.round(complete_round1_model-start_the_clock,3))

    if manager.run_optimal_allocation:
        print("Reading results of optimal allocation map.")
        allocation_schedule_1d = []
        for v in model.Anq.values():
            if np.round(v.X,0) == 1:
                allocation_schedule_1d.append(1)
            else:
                allocation_schedule_1d.append(0)
        allocation_schedule = np.reshape(allocation_schedule_1d, (manager.n_nights_in_semester, manager.n_quarters_in_night))
        holder = np.zeros(np.shape(allocation_schedule))
        allocation_map_1D, allocation_map_2D, weathered_map = mf.build_allocation_map(allocation_schedule, holder, available_slots_in_each_night, manager.n_slots_in_night)
        mf.convert_allocation_array_to_binary(allocation_schedule, manager.all_dates_array, manager.output_directory + "optimal_allocation_binary_schedule.txt")

    print("Building human readable schedule.")
    combined_semester_schedule_available = hf.write_available_human_readable(
                                manager.all_dates_dict, manager.current_day, manager.semester_length,
                                manager.n_nights_in_semester, manager.n_slots_in_night, twilight_map_remaining_2D,
                                allocation_map_2D, weathered_map, manager.nonqueue_map_file)
    np.savetxt(manager.output_directory + 'raw_combined_semester_schedule_available.txt',
        combined_semester_schedule_available, delimiter=',', fmt="%s")

    combined_semester_schedule_stars = hf.write_stars_schedule_human_readable( \
            combined_semester_schedule_available, model.Yrds, list(manager.requests_frame['Starname']),
            manager.semester_length, manager.n_slots_in_night, manager.n_nights_in_semester,
            manager.all_dates_dict, manager.slots_needed_for_exposure_dict, manager.current_day)
    np.savetxt(manager.output_directory + 'raw_combined_semester_schedule_Round1.txt',
        combined_semester_schedule_stars, delimiter=',', fmt="%s")

    round = 'Round 1'
    rf.build_fullness_report(combined_semester_schedule_stars, allocation_map_2D, manager.requests_frame,
                                manager.slot_size, round, manager.output_directory)

    file = open(manager.output_directory + "runReport.txt", "w")
    file.close()
    rf.write_out_results(manager, model, round, start_the_clock)

    if manager.run_round_two:
        print("Beginning Round 2 Scheduling.")
        model.constraint_fix_previous_objective()
        model.set_objective_maximize_slots_used()
        model.solve_model()

        complete_round2_model = time.time()
        print("Total Time to finish round 2 solver: ",
            np.round(complete_round2_model-start_the_clock,3))

        combined_semester_schedule_stars = hf.write_stars_schedule_human_readable(
                combined_semester_schedule_available, model.Yrds, manager.requests_frame['Starname'],
                manager.semester_length, manager.n_slots_in_night, manager.n_nights_in_semester,
                manager.all_dates_dict, manager.slots_needed_for_exposure_dict, manager.current_day)
        np.savetxt(manager.output_directory + 'raw_combined_semester_schedule_Round2.txt',
            combined_semester_schedule_stars, delimiter=',', fmt="%s")

        round = 'Round 2'
        rf.build_fullness_report(combined_semester_schedule_stars, allocation_map_2D,
                                    manager.requests_frame, manager.slot_size, round, manager.output_directory)

        scheduleR1 = np.loadtxt(output_directory + 'raw_combined_semester_schedule_Round1.txt',
            delimiter=',', dtype=str)
        scheduleR2 = np.loadtxt(output_directory + 'raw_combined_semester_schedule_Round2.txt',
            delimiter=',', dtype=str)
        R2_requests = rf.get_gap_filler_targets(scheduleR1, scheduleR2, manager.all_dates_dict[manager.current_day])
        np.savetxt(manager.output_directory + 'Round2_Requests.txt', R2_requests, delimiter=',', fmt="%s")

        rf.write_out_results(manager, model, round, start_the_clock)

    else:
        print("Not running Round 2. Duplicating Raw Schedule as dummy file.")
        np.savetxt(manager.output_directory + 'raw_combined_semester_schedule_Round2.txt', \
            combined_semester_schedule_stars, delimiter=',', fmt="%s")
        R2_requests = []
        np.savetxt(manager.output_directory + 'Round2_Requests.txt', R2_requests, delimiter=',', fmt="%s")

    manager.combined_semester_schedule_stars = combined_semester_schedule_stars
    print("The optimal semester schedule is found, clear skies!")
