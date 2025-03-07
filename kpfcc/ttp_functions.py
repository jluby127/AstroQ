"""
Module for running the TTP.
Designed to be only run as a function call from the generateScript.py script.

See https://github.com/lukehandley/ttp/tree/main for more info about the TTP

Example usage:
    import ttp_functions as ttp
"""
import sys
import os

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta
import gurobipy as gp
from gurobipy import GRB

import ttp.formatting as formatting
import ttp.telescope as telescope
import ttp.plotting as plotting
import ttp.model as model

# path2modules = os.path.dirname(os.path.abspath(__file__))
# sys.path.append(path2modules)
import kpfcc.backup_star_functions as bsf
import kpfcc.helper_functions as hf
import kpfcc.processing_functions as pf

named_colors = ['blue', 'red', 'green', 'gold', 'maroon', 'gray', 'orange', 'magenta', 'purple']

def run_ttp(manager):
    """
    Produce the TTP solution given the results of the autoscheduler

    Args:
        schedule_file (str) - the path and filename of the outputs of the autoscheduler
        nightly_start_stop_times_file (str) - the path and file name to the start/stop times
        request_sheet (str) - the path and filename to the requests csv
        current_day (str) - today's date in format YYYY-MM-DD
        output_dir (str) - the path to save the results

    Returns:
        None
    """

    print("Preparing for the TTP.")
    observers_path = manager.semester_directory + 'reports/observer/' + str(manager.current_day) + '/'
    check1 = os.path.isdir(observers_path)
    if not check1:
        os.makedirs(observers_path)

    observatory = telescope.Keck1()
    the_schedule = np.loadtxt(manager.schedule_file, delimiter=',', dtype=str)
    nightly_start_stop_times = pd.read_csv(manager.nightly_start_stop_times_file)

    # Determine time bounds of the night
    day_in_semester = manager.all_dates_dict[manager.current_day]
    idx = nightly_start_stop_times[nightly_start_stop_times['Date'] == str(manager.current_day)].index[0]
    night_start_time = nightly_start_stop_times['Start'][idx]
    night_stop_time = nightly_start_stop_times['Stop'][idx]
    observation_start_time = Time(str(mangaer.current_day) + "T" + str(night_start_time),
        format='isot')
    observation_stop_time = Time(str(manager.current_day) + "T" + str(night_stop_time),
        format='isot')
    total_time = np.round((observation_stop_time.jd-observation_start_time.jd)*24,3)
    print("Time in Night for Observations: " + str(total_time) + " hours.")

    # Prepare relevant files
    round_two_requests = np.loadtxt(output_path + 'Round2_Requests.txt', dtype=str)
    send_to_ttp = pf.prepare_for_ttp(manager.request_frame, the_schedule[day_in_semester],
                                        round_two_requests)
    filename = output_path + 'Selected_' + str(manager.current_day) + ".txt"
    send_to_ttp.to_csv(filename, index=False)
    target_list = formatting.theTTP(filename)

    print("Initializing the optimizer.")
    solution = model.TTPModel(observation_start_time, observation_stop_time, target_list,
                                observatory, observers_path)

    print("Processing the solution.")
    plotting.writeStarList(solution.plotly, observation_start_time, manager.current_day,
                        outputdir=observers_path)
    plotting.plot_path_2D(solution,outputdir=observers_path)
    plotting.nightPlan(solution.plotly, manager.current_day, outputdir=observers_path)
    obs_and_times = pd.read_csv(observers_path + 'ObserveOrder_' + str(manager.current_day) + ".txt")
    pf.write_starlist(manager.request_frame, solution.plotly, observation_start_time, solution.extras,
                        round_two_requests, str(manager.current_day), observers_path)
    print("The optimal path through the sky for the selected stars is found. Clear skies!")

def produce_bright_backups(backup_file, backup_observability_file, observatory,
                            observation_start_time, observation_stop_time,
                            current_day, output_dir):
    """
    Produce the TTP solution given the results of the autoscheduler

    Args:
        backup_file (str) - the path and filename request sheet equivalent for bright backups
        backup_observability_file (str) - the path and filename to the spreadsheet of how many
                                         hours each night each backup star is observable for
        observation_start_time (str) - the start time of the night in format HH:MM
        observation_stop_time (str) - the stop time of the night in format HH:MM
        current_day (str) - today's date in format YYYY-MM-DD
        output_dir (str) - the path to save the results

    Returns:
        None
    """
    print("Generating bright star backup weather script.")
    backups_path = output_dir + 'outputs/' + current_day + '/Backups/'
    check = os.path.isdir(backups_path)
    if not check:
        os.makedirs(backups_path)

    backup_starlist = pd.read_csv(backup_file)
    backup_observability_frame = pd.read_csv(backup_observability_file)

    backups = bsf.get_stars_for_tonight(backup_starlist, backup_observability_frame,
            current_day[5:], minimum_up_time=4.0)
    backups.drop(columns='index', inplace=True)

    filename_backups = backups_path + '/Backups_' + str(current_day) + ".txt"
    backups_full = bsf.get_times_worth(backups, total_time+1)
    backups_full.to_csv(filename_backups, index=False)
    targlist_backups = formatting.theTTP(filename_backups)

    solution_b = model.TTPModel(observation_start_time, observation_stop_time,
                                targlist_backups, observatory, backups_path,
                                runtime=300, optgap=0.05)

    plotting.writeStarList(solution_b.plotly, observation_start_time, current_day,
                                outputdir = backups_path)
    plotting.plot_path_2D(solution_b, outputdir = backups_path)
    plotting.nightPlan(solution_b.plotly, current_day, outputdir = backups_path)
    obs_and_times_b = pd.read_csv(backups_path + 'ObserveOrder_' + current_day + ".txt")
    pf.write_starlist(backup_starlist, solution_b.plotly, observation_start_time,
                        solution_b.extras, [], str(current_day), backups_path)
    print("Bright backups script created.")
