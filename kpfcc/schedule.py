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

import kpfcc.io as io
import kpfcc.access as ac
import kpfcc.management as mn
import kpfcc.maps as mp
import kpfcc.constraints as cf

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
    round_info = 'Round1'

    print("Reading inputs and prepping files.")
    twilight_map_remaining_2D, available_slots_in_each_night = ac.construct_twilight_map(manager)

    # When running a normal schedule, include the observatory's allocation map
    if manager.run_optimal_allocation == False:
        weather_diff_remaining, allocation_map_1D, allocation_map_2D, weathered_map = \
                                    mp.prepare_allocation_map(manager, available_slots_in_each_night)
    # When running the optimal allocation, all dates are possible except for those specifically blacked out
    else:
        # Weather arrays are zeros since no weather losses are modeled
        weather_diff_remaining = np.zeros(manager.n_nights_in_semester, dtype='int')
        weathered_map = np.zeros((manager.n_nights_in_semester, manager.n_slots_in_night), dtype='int')

        # allocation maps are ones because all nights are possible to be allocated
        allocation_map_1D = np.ones(manager.n_slots_in_semester, dtype='int')
        allocation_map_2D = np.ones((manager.n_nights_in_semester, manager.n_slots_in_night), dtype='int')
    manager.weather_diff_remaining = weather_diff_remaining

    print("Computing intersection of all maps. ")
    available_indices_for_request = mp.produce_ultimate_map(manager, allocation_map_1D.flatten(),
                                                            twilight_map_remaining_2D.flatten())

    print("Building Gorubi model.")
    Aframe, Aset, schedulable_requests = cf.define_slot_index_frame(manager, available_indices_for_request)
    model = cf.GorubiModel(manager, Aset, Aframe, schedulable_requests)

    model.constraint_build_theta()
    # model.constraint_build_theta_time_normalized()
    # model.constraint_build_theta_program_normalized()

    model.constraint_one_request_per_slot()
    model.constraint_reserve_multislot_exposures()
    model.constraint_max_visits_per_night()
    model.constraint_enforce_internight_cadence()
    if manager.run_optimal_allocation:
        model.constraint_set_max_quarters_allocated()
        model.constraint_set_max_onsky_allocated()
        model.constraint_relate_allocation_and_onsky()
        model.constraint_all_portions_of_night_represented()
        model.constraint_forbidden_quarter_patterns()
        model.constraint_cannot_observe_if_not_allocated(twilight_map_remaining_2D)
        if os.path.exists(manager.blackout_file):
            constraint_enforce_restricted_nights(limit=0)
        if os.path.exists(manager.whiteout_file):
            constraint_enforce_restricted_nights(limit=1)
        if manager.include_aesthetic:
            model.constraint_max_consecutive_onsky()
            model.constraint_minimum_consecutive_offsky()
            model.constraint_maximize_baseline()

    print("Total Time to build constraints: ", np.round(time.time()-start_the_clock,3))
    model.set_objective_minimize_theta()
    # model.set_objective_minimize_theta_time_norm()
    # model.set_objective_minimize_theta_prog_norm()

    model.solve_model()
    print("Total Time to finish solver: ", np.round(time.time()-start_the_clock,3))

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
        allocation_map_1D, allocation_map_2D, weathered_map = mp.build_allocation_map(allocation_schedule, holder, available_slots_in_each_night, manager.n_slots_in_night)
        mp.convert_allocation_array_to_binary(allocation_schedule, manager.all_dates_array, manager.output_directory + "optimal_allocation_binary_schedule.txt")
        manager.allocation_all = allocation_map_1D

    print("Building human readable schedule.")
    combined_semester_schedule_available = io.write_available_human_readable(manager,
                                        twilight_map_remaining_2D, allocation_map_2D, weathered_map)
    combined_semester_schedule_stars = io.write_stars_schedule_human_readable(
                                                    combined_semester_schedule_available, model.Yrds,
                                                    manager, round_info)
    io.build_fullness_report(combined_semester_schedule_stars, allocation_map_2D, manager, round_info)
    io.write_out_results(manager, model, round_info, start_the_clock)

    if manager.run_round_two:
        print("Beginning Round 2 Scheduling.")
        round_info = 'Round2'
        model.constraint_fix_previous_objective()
        model.set_objective_maximize_slots_used()
        model.solve_model()

        print("Total Time to finish Round 2 solver: ", np.round(time.time()-start_the_clock,3))

        combined_semester_schedule_stars = io.write_stars_schedule_human_readable(
                combined_semester_schedule_available, model.Yrds, manager, round_info)
        io.build_fullness_report(combined_semester_schedule_stars, allocation_map_2D, manager, round_info)
        mn.get_gap_filler_targets(manager)
        io.write_out_results(manager, model, round_info, start_the_clock)

    else:
        print("Not running Round 2. Duplicating Raw Schedule as dummy file.")
        combined_semester_schedule_stars = io.write_stars_schedule_human_readable(
                combined_semester_schedule_available, model.Yrds, manager, "Round2")
        np.savetxt(manager.output_directory + 'Round2_Requests.txt', [], delimiter=',', fmt="%s")

    manager.combined_semester_schedule_stars = combined_semester_schedule_stars

    print("Total Time to complete autoscheduler: ", np.round(time.time()-start_the_clock,3))
    print("The optimal semester schedule is found, clear skies!")
