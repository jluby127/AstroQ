"""
Module for running the TTP.
Designed to be only run as a function call from the generateScript.py script.

See https://github.com/lukehandley/ttp/tree/main for more info about the TTP

Example usage:
    import ttp_functions as ttp
"""
import sys
import os
import pickle

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta
import gurobipy as gp
from gurobipy import GRB

sys.path.append('/Users/jack/Documents/github/ttp/')
import ttp.formatting as formatting
import ttp.telescope as telescope
import ttp.plotting as plotting
import ttp.model as model

import astroq.weather as wt
import astroq.io as io
import astroq.maps as mp

named_colors = ['blue', 'red', 'green', 'gold', 'maroon', 'gray', 'orange', 'magenta', 'purple']

def run_ttp(manager, include_specmatch=False):
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
    observers_path = manager.semester_directory + 'reports/observer/' + str(manager.current_day) + '/'
    check1 = os.path.isdir(observers_path)
    if not check1:
        os.makedirs(observers_path)

    observatory = telescope.Keck1()
    the_schedule = np.loadtxt(manager.upstream_path + "/outputs/" + str(manager.current_day) + "/raw_combined_semester_schedule_Round2.txt", delimiter=',', dtype=str)
    nightly_start_stop_times = pd.read_csv(manager.nightly_start_stop_times_file)

    day_in_semester = manager.all_dates_dict[manager.current_day]
    idx = nightly_start_stop_times[nightly_start_stop_times['Date'] == str(manager.current_day)].index[0]
    night_start_time = nightly_start_stop_times['Start'][idx]
    night_stop_time = nightly_start_stop_times['Stop'][idx]
    observation_start_time = Time(str(manager.current_day) + "T" + str(night_start_time),
        format='isot')
    observation_stop_time = Time(str(manager.current_day) + "T" + str(night_stop_time),
        format='isot')
    total_time = np.round((observation_stop_time.jd-observation_start_time.jd)*24,3)
    print("Time in Night for Observations: " + str(total_time) + " hours.")

    # Prepare relevant files
    round_two_requests = np.loadtxt(manager.upstream_path + "/outputs/" + str(manager.current_day) + '/Round2_Requests.txt', dtype=str, delimiter="\n")
    send_to_ttp = prepare_for_ttp(manager.requests_frame, the_schedule[day_in_semester],
                                        round_two_requests)
    if include_specmatch:
        all_specmatch = pd.read_csv(manager.DATADIR + '/SPOCS_all_stars.csv')
        all_rows_for_script_lines = pd.concat([manager.requests_frame, all_specmatch], ignore_index=True)
        all_rows_for_script_lines.reset_index(inplace=True, drop=True)

        specmatch_fillers = fill_gaps_specmatch(manager.DATADIR + '/SPOCS_all_stars.csv', manager.past_database_file, ra1=13, ra2=17)
        send_to_ttp = pd.concat([send_to_ttp, specmatch_fillers], ignore_index=True)
        send_to_ttp.reset_index(inplace=True, drop=True)

    filename = manager.upstream_path + "/outputs/" + str(manager.current_day) + '/Selected_' + str(manager.current_day) + ".txt"
    send_to_ttp.to_csv(filename, index=False)
    target_list = formatting.theTTP(filename)

    solution = model.TTPModel(observation_start_time, observation_stop_time, target_list,
                                observatory, observers_path, runtime=300, optgap=0.01, useHighEl=False)

    gurobi_model_backup = solution.gurobi_model  # backup the attribute, probably don't need this
    del solution.gurobi_model                   # remove attribute so pickle works
    save_data = [solution]
    with open(manager.reports_directory + '/observer/' + manager.current_day + '/ttp_data.pkl', 'wb') as f:
        pickle.dump(save_data, f)

    plotting.writeStarList(solution.plotly, observation_start_time, manager.current_day,
                        outputdir=observers_path)
    plotting.plot_path_2D(solution,outputdir=observers_path)
    plotting.nightPlan(solution.plotly, manager.current_day, outputdir=observers_path)
    obs_and_times = pd.read_csv(observers_path + 'ObserveOrder_' + str(manager.current_day) + ".txt")

    if include_specmatch:
        io.write_starlist(all_rows_for_script_lines, solution.plotly, observation_start_time, solution.extras,
                        round_two_requests, str(manager.current_day), observers_path)
    else:
        io.write_starlist(manager.requests_frame, solution.plotly, observation_start_time, solution.extras,
                        round_two_requests, str(manager.current_day), observers_path)
    print("The optimal path through the sky for the selected stars is found. Clear skies!")

def fill_gaps_specmatch(filename, past_history_file, ra1=13, ra2=17):
    library = pd.read_csv(filename)
    pasthistory = pd.read_csv(past_history_file)
    mask = ~library["starname"].isin(pasthistory["star_id"])

    filtered_df = library[mask]
    ra1 = filtered_df['ra'] > ra1*15
    ra2 = filtered_df['ra'] < ra2*15
    short = filtered_df['exptime'] < 1800

    filtered_df = filtered_df[ra1&ra2&short]
    filtered_df.sort_values(by='ra', inplace=True)
    filtered_df = filtered_df[-50:]
    filtered_df.reset_index(inplace=True, drop=True)

    ttp_frame = pd.DataFrame({
        "Starname":filtered_df['starname'],
        "RA":filtered_df['ra'],
        "Dec":filtered_df['dec'],
        "Exposure Time":filtered_df['exptime'],
        "Exposures Per Visit":[1]*len(filtered_df),
        "Visits In Night":[1]*len(filtered_df),
        "Intra_Night_Cadence":[0]*len(filtered_df),
        "Priority":[1]*len(filtered_df)
    })

    return ttp_frame

def produce_bright_backups(manager, nstars_max=100):
    """
    Produce the TTP solution given the results of the autoscheduler

    Args:
        manager (obj) - an instance of the manager object

    Returns:
        None
    """
    backups_path = manager.reports_directory + '/observer/' + manager.current_day + "/backups/"
    check = os.path.isdir(backups_path)
    if not check:
        os.makedirs(backups_path)

    nightly_start_stop_times = pd.read_csv(manager.nightly_start_stop_times_file)
    # Determine time bounds of the night
    day_in_semester = manager.all_dates_dict[manager.current_day]
    idx = nightly_start_stop_times[nightly_start_stop_times['Date'] == str(manager.current_day)].index[0]
    night_start_time = nightly_start_stop_times['Start'][idx]
    observation_start_time = Time(str(manager.current_day) + "T" + str(night_start_time),
        format='isot')
    night_stop_time = nightly_start_stop_times['Stop'][idx]
    observation_stop_time = Time(str(manager.current_day) + "T" + str(night_stop_time),
        format='isot')
    diff_minutes = int(abs((observation_stop_time - observation_start_time).to('min').value))
    print("Minutes on sky: ", diff_minutes)

    backup_starlist = pd.read_csv(manager.backup_file)
    manager.requests_frame = backup_starlist
    available_indices = mp.produce_ultimate_map(manager, manager.requests_frame, running_backup_stars=True)
    slots_available_tonight_for_star = {k: len(v[0]) for k, v in available_indices.items()}
    stars_with_sufficient_availability_tonight = [k for k, v in slots_available_tonight_for_star.items() if v > int(0.25*int(diff_minutes/5))]

    manager.requests_frame = backup_starlist
    isTonight = backup_starlist['starname'].isin(stars_with_sufficient_availability_tonight)
    hasDR3name = backup_starlist['gaia_id'].str.startswith('Gaia DR2')
    pool_tonight = manager.requests_frame[isTonight&hasDR3name]
    pool_tonight = pool_tonight.sample(frac=1).reset_index(drop=True)
    pool_tonight = pool_tonight[:nstars_max]

    ready_for_ttp = prepare_for_ttp(pool_tonight, list(pool_tonight['starname']), [])
    ready_for_ttp.to_csv(backups_path + "selected_stars.csv", index=False)
    target_list = formatting.theTTP(backups_path + "selected_stars.csv")

    observatory = telescope.Keck1()
    solution_b = model.TTPModel(observation_start_time, observation_stop_time,
                                target_list, observatory, backups_path,
                                runtime=120, optgap=0.05)

    plotting.writeStarList(solution_b.plotly, observation_start_time, manager.current_day,
                                outputdir = backups_path)
    plotting.plot_path_2D(solution_b, outputdir = backups_path)
    plotting.nightPlan(solution_b.plotly, manager.current_day, outputdir = backups_path)
    obs_and_times_b = pd.read_csv(backups_path + 'ObserveOrder_' + manager.current_day + ".txt")
    io.write_starlist(pool_tonight, solution_b.plotly, observation_start_time,
                        solution_b.extras, [], manager.current_day, backups_path, "backups")
    print("Bright backups script created.")

def prepare_for_ttp(request_frame, night_plan, round_two_targets):
    """
    Prepare tonight's scheduled stars for their run through the TTP (separate software package)

    Args:
        request_sheet (dataframe): the csv of PI requests
        night_plan (array): the n'th row of combined_semester_schedule array
        round_two_targets (array): a 1D list of the stars that were added in the bonus round

    Returns:
        to_ttp (dataframe): the data on the stars to be observed tonight, formatted in the way that
                            the TTP software requires as an input
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
