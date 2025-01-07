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

sys.path.append(os.environ["TTP_PATH"])
import formatting
import telescope
import plotting
import model

path2modules = os.path.dirname(os.path.abspath(__file__))
sys.path.append(path2modules)
import backup_star_functions as bsf
import helper_functions as hf
import processing_functions as pf

named_colors = ['blue', 'red', 'green', 'gold', 'maroon', 'gray', 'orange', 'magenta', 'purple']

def build_directories(current_day, output_dir):
    """
    Make necessary directories

    Args:
        current_day (str) - today's date in format YYYY-MM-DD
        output_dir (str) - the path to save the results

    Returns:
        None
    """
    # Build required directories if necessary
    savepath = output_dir + 'reports/observer/' + str(current_day) + '/'
    check1 = os.path.isdir(savepath)
    if not check1:
        os.makedirs(savepath)
    output_path = output_dir + 'outputs/' + str(current_day) + "/"
    check2 = os.path.isdir(output_path)
    if not check2:
        os.makedirs(output_path)

    return savepath, output_path

def run_ttp(schedule_file, nightly_start_stop_times_file, request_sheet, current_day, output_dir):
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
    savepath, output_path = build_directories(current_day, output_dir)
    observatory = telescope.Keck1()
    semester_start_date, semester_end_date, semester_length, semester_year, semester_letter = \
            hf.get_semester_info(current_day)
    all_dates_dict = hf.build_date_dictionary(semester_start_date, semester_length)
    the_schedule = np.loadtxt(schedule_file, delimiter=',', dtype=str)
    nightly_start_stop_times = pd.read_csv(nightly_start_stop_times_file)

    # Determine time bounds of the night
    day_in_semester = all_dates_dict[current_day]
    idx = nightly_start_stop_times[nightly_start_stop_times['Date'] == str(current_day)].index[0]
    night_start_time = nightly_start_stop_times['Start'][idx]
    night_stop_time = nightly_start_stop_times['Stop'][idx]
    observation_start_time = Time(str(current_day) + "T" + str(night_start_time),
        format='isot')
    observation_stop_time = Time(str(current_day) + "T" + str(night_stop_time),
        format='isot')
    total_time = np.round((observation_stop_time.jd-observation_start_time.jd)*24,3)
    print("Time in Night for Observations: " + str(total_time) + " hours.")

    # Prepare relevant files
    round_two_requests = np.loadtxt(output_path + 'Round2_Requests.txt', dtype=str)
    request_frame = pd.read_csv(request_sheet)
    send_to_ttp = pf.prepare_for_ttp(request_frame, the_schedule[day_in_semester],
                                        round_two_requests)
    filename = output_path + 'Selected_' + str(current_day) + ".txt"
    send_to_ttp.to_csv(filename, index=False)
    target_list = formatting.theTTP(filename)

    print("Initializing the optimizer.")
    solution = model.TTPModel(observation_start_time, observation_stop_time, target_list,
                                observatory, savepath)

    print("Processing the solution.")
    plotting.writeStarList(solution.plotly, observation_start_time, current_day,
                        outputdir=savepath)
    plotting.plot_path_2D(solution,outputdir=savepath)
    plotting.nightPlan(solution.plotly, current_day, outputdir=savepath)
    obs_and_times = pd.read_csv(savepath + 'ObserveOrder_' + str(current_day) + ".txt")
    pf.write_starlist(request_frame, solution.plotly, observation_start_time, solution.extras,
                        round_two_requests, str(current_day), savepath)
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
    savepath = output_dir + 'outputs/' + current_day + '/Backups/'
    check = os.path.isdir(output_dir + 'outputs/' + current_day + '/Backups/')
    if not check:
        os.makedirs(savepath)

    backup_starlist = pd.read_csv(backup_file)
    backup_observability_frame = pd.read_csv(backup_observability_file)

    backups = bsf.get_stars_for_tonight(backup_starlist, backup_observability_frame,
            current_day[5:], minimum_up_time=4.0)
    backups.drop(columns='index', inplace=True)

    filename_backups = savepath + '/Backups_' + str(current_day) + ".txt"
    backups_full = bsf.get_times_worth(backups, total_time+1)
    backups_full.to_csv(filename_backups, index=False)
    targlist_backups = formatting.theTTP(filename_backups)

    solution_b = model.TTPModel(observation_start_time, observation_stop_time,
                                targlist_backups, observatory, savepath,
                                runtime=300, optgap=0.05)

    plotting.writeStarList(solution_b.plotly, observation_start_time, current_day,
                                outputdir = savepath)
    plotting.plot_path_2D(solution_b, outputdir = savepath)
    plotting.nightPlan(solution_b.plotly, current_day, outputdir = savepath)
    obs_and_times_b = pd.read_csv(savepath + 'Backups/ObserveOrder_' + current_day + ".txt")
    pf.write_starlist(backup_starlist, solution_b.plotly, observation_start_time,
                        solution_b.extras, [], str(current_day),
                        savepath + "Backups/")
    print("Bright backups script created.")
