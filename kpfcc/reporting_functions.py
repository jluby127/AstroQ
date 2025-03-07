"""
Module for building reports after the solution is found. Designed to be only run as a function
call from the generateScript.py script.

Example usage:
    import reporting_functions as rf
"""
import os
import math
import time
from astropy.time import Time
from astropy.time import TimeDelta

import numpy as np
import pandas as pd

def build_fullness_report(combined_semester_schedule, allocation_map_2D, manager, round_info):
    """
    Determine how full the schedule is: slots available, slots scheduled, and slots required

    Args:
        combined_semester_schedule (array): a 2D array of dimensions n_nights_in_semester by
                                            n_slots_in_night where elements denote how the slot is
                                            used: target, twilight, weather, not scheduled.
        allocation_map_2D (array): a 2D array where rows represent a night and columns represent
                                   the quarter within that night. Values are 1 if that
                                   night/quarter is allocated and 0 if not.
        manager (obj): a data_admin object

    Returns:
        None
    """
    file = open(manager.output_directory + "runReport.txt", "a")
    file.write("Stats for " + str(round_info) + "\n")
    file.write("------------------------------------------------------" + "\n")
    listnames = list(manager.requests_frame['Starname'])
    unavailable = 0
    unused = 0
    used = 0
    for b, item1 in enumerate(combined_semester_schedule):
        for c, item2 in enumerate(combined_semester_schedule[b]):
            if "X" in combined_semester_schedule[b][c] or "*" in combined_semester_schedule[b][c] \
                    or "W" in combined_semester_schedule[b][c]:
                unavailable += 1
            if combined_semester_schedule[b][c] == "":
                unused += 1
            if combined_semester_schedule[b][c] in listnames or "RM___" in combined_semester_schedule[b][c]:
                used += 1
    available = unused + used
    allocated = np.sum(allocation_map_2D.flatten())
    file.write("N slots in semester:" + str(np.prod(combined_semester_schedule.shape)) + "\n")
    file.write("N available slots:" + str(allocated) + "\n")
    file.write("N slots scheduled: " + str(used) + "\n")
    file.write("N slots left empty: " + str(allocated-used) + "\n")

    total_slots_requested = 0
    for i in range(len(manager.requests_frame)):
        total_slots_requested += manager.requests_frame['# of Nights Per Semester'][i]* \
            math.ceil(manager.requests_frame['Nominal Exposure Time [s]'][i]/(manager.slot_size*60.))
    file.write("N slots requested (total): " + str(total_slots_requested) + "\n")
    percentage = np.round((used*100)/allocated,3)
    file.write("Percent full: " + str(percentage) + "%." + "\n")
    file.close()

def report_allocation_stats(allocation_map_2D, all_dates_dict, current_day, outputdir):
    """
    Return the number of Q1, Q2, Q3, and Q4 nights, as well as single, half, three quarter
    and full nights. Further produce a plot showing the allocation.

    Args:
        allocation_map_2D (array): a 2D array where rows represent a night and columns represent
                                   the quarter within that night. Values are 1 if that
                                   night/quarter is allocated and 0 if not.
        all_dates_dict (dict): a dictionary where keys are calendar dates in format (YYYY-MM-DD)
                               and values of the day number in the semester
        current_day (str): today's date in format 'YYYY-MM-DD'
        outputdir (str): the path to save the output file

    Returns:
        None
    """
    all_dates_array = list(all_dates_dict.keys())
    current_day_number = all_dates_dict[current_day]
    ff = open(outputdir + "allocation_stats.txt", "w")

    q1s = 0
    q2s = 0
    q3s = 0
    q4s = 0
    count0 = 0
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0
    for j, item in enumerate(allocation_map_2D):
        holder = [0, 0, 0, 0]
        if allocation_map_2D[j][0] == 1.:
            q1s += 1
            date_info = all_dates_array[current_day_number + j] + " - q0"
            ff.write(date_info + "\n")
            holder[0] = 1
        if allocation_map_2D[j][1] == 1.:
            q2s += 1
            date_info = all_dates_array[current_day_number + j] + " - q1"
            ff.write(date_info + "\n")
            holder[1] = 1
        if allocation_map_2D[j][2] == 1.:
            q3s += 1
            date_info = all_dates_array[current_day_number + j] + " - q2"
            ff.write(date_info + "\n")
            holder[2] = 1
        if allocation_map_2D[j][3] == 1.:
            q4s += 1
            date_info = all_dates_array[current_day_number + j] + " - q3"
            ff.write(date_info + "\n")
            holder[3] = 1
        allocated_quarters = np.sum(allocation_map_2D[j])
        if allocated_quarters == 0:
            count0 += 1
        if allocated_quarters == 1:
            count1 += 1
        if allocated_quarters == 2:
            count2 += 1
        if allocated_quarters == 3:
            count3 += 1
        if allocated_quarters == 4:
            count4 += 1

    ff.write("\n")
    ff.write("There are " + str(q1s) + " first quarters." + "\n")
    ff.write("There are " + str(q2s) + " second quarters." + "\n")
    ff.write("There are " + str(q3s) + " third quarters." + "\n")
    ff.write("There are " + str(q4s) + " fourth quarters." + "\n")
    ff.write("\n")
    ff.write("There are " + str(count0) + " no quarter nights." + "\n")
    ff.write("There are " + str(count1) + " 1 quarter nights."+ "\n")
    ff.write("There are " + str(count2) + " 2 quarter nights."+ "\n")
    ff.write("There are " + str(count3) + " 3 quarter nights."+ "\n")
    ff.write("There are " + str(count4) + " 4 quarter nights."+ "\n")
    ff.write("\n")
    total_quarters = count1 + 2*count2 + 3*count3 + 4*count4
    total_nights = count1 + count2 + count3 + count4
    ff.write("Total quarters allocated: " + str(total_quarters) + "\n")
    ff.write("Total unique nights allocated: " + str(total_nights) + "\n")
    ff.close()

# def build_observed_map_past(unique_hst_dates_observed, quarter_observed, n_obs_on_date,
#                             starmap_template_filename):
#     """
#     Construct stage one (past) of the starmap which is then used to build the cadence plot.
#
#     Args:
#         unique_hst_dates_observed (array): the dates a given target was observed in the past
#         quarter_observed (array): corresponding quarter of night when target was observed
#         n_obs_on_date (array): corresponding number of obs on of nights when target was observed
#         starmap_template_filename (str): the path and filename of the starmap template
#
#     Returns:
#         starmap (dataframe): an updated starmap which includes the past history of the target
#     """
#     starmap = pd.read_csv(starmap_template_filename)
#     observed = [False]*len(starmap)
#     n_observed = [0]*len(starmap)
#     for i, item in enumerate(unique_hst_dates_observed):
#         ind = list(starmap['Date']).index(unique_hst_dates_observed[i])
#         observed[ind + int(quarter_observed[i]-0.5)] = True
#         n_observed[ind + int(quarter_observed[i]-0.5)] = n_obs_on_date[i]
#     starmap['Observed'] = observed
#     starmap['N_obs'] = n_observed
#     return starmap

def build_observed_map_past(past_info, starmap_template_filename):
    """
    Construct stage one (past) of the starmap which is then used to build the cadence plot.

    Args:
        past_info (array): 2D array following schema of the values for the database_info_dict
        starmap_template_filename (str): the path and filename of the starmap template

    Returns:
        starmap (dataframe): an updated starmap which includes the past history of the target
    """
    starmap = pd.read_csv(starmap_template_filename)
    observed = [False]*len(starmap)
    n_observed = [0]*len(starmap)
    for i, item in enumerate(past_info[1]):
        ind = list(starmap['Date']).index(past_info[1][i])
        observed[ind + int(past_info[2][i]-0.5)] = True
        n_observed[ind + int(past_info[2][i]-0.5)] = past_info[3][i]
    starmap['Observed'] = observed
    starmap['N_obs'] = n_observed
    return starmap

def build_observed_map_future(combined_semester_schedule, starname, starmap,
                            slots_needed_for_exposure, allocation_map_1D, weather_diff_1D,
                            days_into_semester):
    """
    Construct stage two (future) of the starmap which is then used to build the cadence plot.

    Args:
        combined_semester_schedule (array): a 2D array of dimensions n_nights_in_semester by
                                            n_slots_in_night where elements denote how the slot is
                                            used: target, twilight, weather, not scheduled.
        starname (str): the request in question
        starmap (dataframe): an updated starmap which includes the past history of the target
        slots_needed_for_exposure (int): number of slots needed to complete the exposure
        allocation_map_1D (array): an array of if each quarter is allocated or not
        weather_diff_1D (array): an array of if each quarter is simulated as weathered or not
        days_into_semester (int): number of days elapsed from start of semester

    Returns:
        starmap (dataframe): an updated starmap which includes the past history of the target
    """
    starmap['Allocated'] = [False]*len(starmap)
    starmap['Weathered'] = [False]*len(starmap)

    for i, item in enumerate(combined_semester_schedule):
        if starname in combined_semester_schedule[i]:
            ind = list(combined_semester_schedule[i]).index(starname)
            quarter = determine_quarter(ind, np.shape(combined_semester_schedule)[1])
            starmap['Observed'][i*4 + quarter] = True
            starmap['N_obs'][i*4 + quarter] = list(combined_semester_schedule[i]).count(
                    starname)/slots_needed_for_exposure[starname]

    for a, item in enumerate(allocation_map_1D):
        # mulitply by 4 because each element of allocation_map_1D represents one quarter
        # starmap['Allocated'][days_into_semester*4 + a] = bool(allocation_map_1D[a])
        starmap['Allocated'][a] = bool(allocation_map_1D[a])
    for w, item in enumerate(weather_diff_1D):
        starmap['Weathered'][days_into_semester*4+ w] = bool(weather_diff_1D[w])

    return starmap

def determine_quarter(value, n_slots_in_night):
    """
    Given a slot number, determine what quarter it is to be observed in.

    Args:
        value (int): the slot number to be observed
        n_slots_in_night (int): the number of slots in the night

    Returns:
        quarter (int): the corresponding quarter (1, 2, 3, 4) of the night
    """
    if value <= int(n_slots_in_night*(1/4.)):
        quarter = 0
    elif value <= int(n_slots_in_night*(2/4.)) and value > int(n_slots_in_night*(1/4.)):
        quarter = 1
    elif value <= int(n_slots_in_night*(3/4.)) and value > int(n_slots_in_night*(2/4.)):
        quarter = 2
    elif value <= n_slots_in_night and value > int(n_slots_in_night*(3/4.)):
        quarter = 3
    elif value <= int(n_slots_in_night*(5/4.)) and value > n_slots_in_night:
        quarter = 3
    else:
        quarter = 0
        print("Houston, we've had a problem. No valid quarter: ", value, n_slots_in_night)
    return quarter

def write_cadence_plot_file(starname, starmap, turn_frame, requests_frame, future_obs,
                            unique_hst_dates_observed, current_day, outputdir):
    """
    Write the file from which later we can produce the cadence plot

    Args:
        starname (int): the string of the star's name
        starmap (dataframe): contains information on observation history and forecast
        turn_frame (dataframe): contains information on the first and last day a target
                                 is observable in each quarter of the night. Pre-computed.
        requests_frame (dataframe): the information on the request strategies
        unique_hst_dates_observed (array): list of dates previously observed
        current_day (str): today's date, format YYYY-MM-DD
        outputdir (str): the path to the save directory
    Returns:
        None
    """
    request_id = requests_frame.index[requests_frame['Starname']==str(starname)][0]
    program_code = requests_frame.loc[request_id,'Program_Code']
    save_path = outputdir + "/cadences/" + str(program_code) + "/"
    os.makedirs(save_path, exist_ok = True)

    n_obs_desired = requests_frame.loc[request_id,'# of Nights Per Semester']
    n_obs_taken = len(unique_hst_dates_observed)
    n_obs_scheduled = future_obs #np.sum(starmap['N_obs']) #n_obs_desired - n_obs_taken #
    cadence = requests_frame.loc[request_id,'Minimum Inter-Night Cadence']
    turn_index = turn_frame.index[turn_frame['Starname']==str(starname)][0]
    turns = [[turn_frame['Q1_on_date'][turn_index], turn_frame['Q1_off_date'][turn_index]],
             [turn_frame['Q2_on_date'][turn_index], turn_frame['Q2_off_date'][turn_index]],
             [turn_frame['Q3_on_date'][turn_index], turn_frame['Q3_off_date'][turn_index]],
             [turn_frame['Q4_on_date'][turn_index], turn_frame['Q4_off_date'][turn_index]]]

    commentsfile = open(save_path + starname + "_Cadence_Interactive.csv", 'w')
    commentsfile.write('#starname:' + str(starname) + '\n')
    commentsfile.write('#programcode:' + str(program_code) + '\n')
    commentsfile.write('#Nobs_scheduled:' + str(n_obs_scheduled) + '\n')
    commentsfile.write('#Nobs_desired:' + str(n_obs_desired) + '\n')
    commentsfile.write('#Nobs_taken:' + str(n_obs_taken) + '\n')
    commentsfile.write('#cadence:' + str(cadence) + '\n')
    commentsfile.write('#q1_start:' + str(turns[0][1]) + '\n')
    commentsfile.write('#q1_end:' + str(turns[0][0]) + '\n')
    commentsfile.write('#q2_start:' + str(turns[1][1]) + '\n')
    commentsfile.write('#q2_end:' + str(turns[1][0]) + '\n')
    commentsfile.write('#q3_start:' + str(turns[2][1]) + '\n')
    commentsfile.write('#q3_end:' + str(turns[2][0]) + '\n')
    commentsfile.write('#q4_start:' + str(turns[3][1]) + '\n')
    commentsfile.write('#q4_end:' + str(turns[3][0]) + '\n')
    commentsfile.write('#current_day:' + str(current_day) + '\n')
    starmap = starmap[starmap.Allocated == True]
    starmap.to_csv(commentsfile, index=False)
    commentsfile.close()

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

def get_unique_nights(star_past_obs, twilight_frame):
    """
    Parse the Jump database for the previous nights where a given star was observed

    Args:
        star_past_obs (dataframe): a subset of the Jump database dataframe, filtered to only
                                    include observations for a specific target name
        twilight_frame (dataframe): the dataframe of morning and evening twilight times for each
                                    night of the semester

    Returns:
        star_past_obs (dataframe): an updated version of the original variable which now includes
                                   columns for UTC date and hst date
        unique_hst_dates_observed (array): a 1D array length # of unique dates observed where each
                                   element is the HST date when the star was observed at least once
        quarter_observed (array): a 1D array length # of unique dates observed where each element
                                 the quarter of the night when the first of the observations
                                 for that night was taken.
    """
    unique_hst_dates_observed = []
    unique_utc_dates_observed = []
    quarter_observed = []
    all_hst_dates = []
    for i in range(len(star_past_obs)):
        timestamp = star_past_obs['utctime'][i][:-6]
        timestamp = Time(timestamp, format = 'iso')
        timestamp.format = 'jd'
        utc_date = star_past_obs['utctime'][i][:10]
        unique_utc_dates_observed.append(utc_date)

        # Note this is arbitrary 24hr subtraction to get from UT date to Hawaii date
        convert_to_hst_civil = TimeDelta(60*24*60,format='sec').jd
        hsttimestamp = timestamp - convert_to_hst_civil
        hsttimestamp.format = 'iso'
        hstdate = str(hsttimestamp)[:10]
        all_hst_dates.append(hstdate)
        hsttimestamp.format = 'jd'
        if hstdate not in unique_hst_dates_observed:
            unique_hst_dates_observed.append(hstdate)
            quarter_observed.append(get_quarter_observed(timestamp, utc_date, twilight_frame))
    star_past_obs['utc_date'] = unique_utc_dates_observed
    star_past_obs['hstDate'] = all_hst_dates

    return star_past_obs, unique_hst_dates_observed, quarter_observed

def get_quarter_observed(timestamp, utc_date, twilight_frame):
    """
    Determine which quarter of the night an observation was taken in

    Args:
        timestamp (float): a BJD value indicating the time of the observation
        utc_date (str): the calendar date of the observation in UTC, format 'YYYY-MM-DD'
        twilight_frame (dataframe): the dataframe of morning and evening twilight times

    Returns:
        quarter (int): the quarter corresponding to the timestamp
    """
    date_index = twilight_frame.index.get_loc(twilight_frame[twilight_frame['time_utc'] == utc_date].index[0])
    start_jd = twilight_frame['12_evening'][date_index]
    end_jd = twilight_frame['12_morning'][date_index]
    length_of_quarter = (end_jd - start_jd)/4.
    rel_timestamp = float(str(timestamp - start_jd))

    quarter = -10
    if rel_timestamp < length_of_quarter:
        quarter = 0.5
    elif rel_timestamp < 2*length_of_quarter and rel_timestamp > length_of_quarter:
        quarter = 1.5
    elif rel_timestamp < 3*length_of_quarter and rel_timestamp > 2*length_of_quarter:
        quarter = 2.5
    elif rel_timestamp < 4*length_of_quarter and rel_timestamp > 3*length_of_quarter:
        quarter = 3.5
    elif rel_timestamp < 5*length_of_quarter and rel_timestamp > 4*length_of_quarter:
        # allow a little bit of leeway, even if this doesn't really make sense
        quarter = 3.5
    else:
        quarter = 0.5
        print("Houston, we've had a problem: target observed in an invalid quarter.")

    return quarter

def get_nobs_on_night(star_past_obs, unique_hst_dates_observed):
    """
    Determine how many exposures were taken of a target on a given night

    Args:
        star_past_obs (dataframe): the updated version of the original variable which now includes
                                   columns for UTC date and hst date
        unique_hst_dates_observed (array): a 1D array of HST dates where the star was observed

    Returns:
        n_obs_on_date (array): a 1D array of length equal to length unique_hst_dates_observed where
                             each element is an integer of the number of exposures taken that night
    """
    n_obs_on_date = []
    for i in range(len(unique_hst_dates_observed)):
        datemask = star_past_obs['hstDate'] == unique_hst_dates_observed[i]
        n_obs_on_date.append(np.sum(datemask))
    return n_obs_on_date

def write_out_results(manager, model, round, start_the_clock):

    print("Writing Report.")
    filename = open(manager.output_directory + "runReport.txt", "a")
    theta_n_var = []
    counter = 0
    for v in model.theta.values():
        varname = v.VarName
        varval = v.X
        counter += varval
    print("Sum of Theta: " + str(counter))
    filename.write("Sum of Theta: " + str(counter) + "\n")

    print("Total Time to complete " + round + ": " + str(np.round(time.time()-start_the_clock,3)))
    filename.write("Total Time to complete " + round +  ": " + str(np.round(time.time()-start_the_clock,3)) + "\n")
    filename.write("\n")
    filename.write("\n")
    filename.close()
