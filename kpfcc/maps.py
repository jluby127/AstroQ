"""
Module for working with the various maps

This module organizes and builds the maps used by the semester solver to block out large portions
of time for one of many reasons. Designed to be only run as a function call from
the generateScript.py script.

Example usage:
    import mapping_functions as mf
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as pt

from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astropy.units as u
import astroplan as apl

import kpfcc.access as ac
import kpfcc.weather as wh

# Define list of observatories which are currently supported.
# To add an obseratory/telescope, add the Astroplan resolvable name to the list in generate_night_plan
# Then add the same name to the appropriate element of the locations dictionary.
# If necessary, add a location to the locations dictionary, if so, add the location to each of the
# pre_sunrise and post_sunrise dictionaries. Ensure times are 14 hours apart, at least one hour
# before the earliest sunset and one hour after the latest sunrise of the year.
locations = {'Hawaii':['Keck Observatory'], 'Arizona':['Kitt Peak National Observatory']}
pre_sunset = {'Hawaii':'07:30', 'Arizona':'05:00'}
post_sunrise = {'Hawaii':'17:30', 'Arizona':'14:00'}
# each telescope's pointing limits are defined by list with elements in following order:
# maximum altitude
# minimum azimuth of nasmyth deck (when no nasmyth deck, use 0)
# maximum azimuth of nasmyth deck (when no nasmyth deck, use 360)
# minimum altitude of nasmyth deck (when no nasmyth deck, use same as below)
# minimum alitude of non-nasmyth deck (when no nasmyth deck, use same as above)
# preferred minimum altitude
# maximum northern declination to enforce preferred minimum altitude
# maximum southern declination to enforce preferred minimum altitude
pointing_limits = {'Keck Observatory': [85.0, 5.3, 146.2, 33.3, 18.0, 30.0, 75.0, -35.0]}

def produce_ultimate_map(manager, allocation_map_1D, twilight_map_remaining_flat):
    """
    Combine all maps for a target to produce the final map

    Args:
        requests_frame (dataframe): the pandas dataframe containing request information

    Returns:
        available_slots_for_request (dictionary): keys are the starnames and values are the 1D array
                                                  of the final map
        available_indices_for_request (dictionary): keys are the starnames and values are a 1D array
                                                  the indices where available_slots_for_request is 1.

    """
    print("Build unique star available slot indices.")
    default_access_maps = ac.construct_access_dict(manager)
    custom_access_maps = construct_custom_map_dict(manager.special_map_file)
    zero_out_names = construct_zero_out_arr(manager.zero_out_file)

    available_slots_for_request = {}
    available_indices_for_request = {}
    for i,row in manager.requests_frame.iterrows():
        name = row['Starname']
        # add an if statement here: if star is a single shot and if running opt allo, skip
        if manager.run_optimal_allocation and row['# of Nights Per Semester'] == 1:
            print("Removing star: " + str(name) + " from model. No single shots in optimal allocation.")
        else:
            accessibility_r = default_access_maps[name]
            access = accessibility_r[manager.today_starting_slot:]

            if name in list(custom_access_maps.keys()):
                custom_map = custom_access_maps[name][manager.today_starting_slot:]
            else:
                custom_map = np.array([1]*manager.n_slots_in_semester)

            zero_out_map = np.array([1]*manager.n_slots_in_semester)
            if name in zero_out_names:
                zero_out_map[:manager.n_slots_in_night] = np.array([0]*manager.n_slots_in_night)

            respect_past_cadence = np.ones(manager.n_slots_in_semester, dtype=np.int64)
            if manager.database_info_dict != {}:
                date_last_observed = manager.database_info_dict[name][0]
                if date_last_observed != '0000-00-00':
                    date_last_observed_number = manager.all_dates_dict[date_last_observed]
                    today_number = manager.all_dates_dict[manager.current_day]
                    diff = today_number - date_last_observed_number
                    if diff < int(row['Minimum Inter-Night Cadence']):
                        block_upcoming_days = int(row['Minimum Inter-Night Cadence']) - diff
                        respect_past_cadence[:block_upcoming_days*manager.n_slots_in_night] = 0

            # Determine which nights a multi-visit request is allowed to be attempted to be scheduled.
            # This equation is a political decision and can be modified.
            # It states that for each visit, after the intra-night cadence time has elapsed,
            # we require a 90 minute window within which to allow for scheduling the next visit.
            # We then assume the next visit is scheduled at the very end of this 90 minute window,
            # which then restarts the clock for any additional visits.
            minimum_time_required = ((int(row['Desired Visits per Night']) - 1)* \
                (int(row['Minimum Intra-Night Cadence']) + 1.5))*3600 #convert hours to seconds
            minimum_slots_required = manager.slots_needed_for_exposure_dict[name]
            no_multi_visit_observations = []
            for d in range(manager.n_nights_in_semester):
                start = d*manager.n_slots_in_night
                end = start + manager.n_slots_in_night
                possible_open_slots = np.sum(allocation_map_1D[start:end] & \
                                            twilight_map_remaining_flat[start:end] & access[start:end])
                if possible_open_slots < minimum_slots_required:
                    no_multi_visit_observations.append([0]*manager.n_slots_in_night)
                else:
                    no_multi_visit_observations.append([1]*manager.n_slots_in_night)
            no_multi_visit_observations = np.array(no_multi_visit_observations)

            nonqueue_map_file_slots_ints = construct_nonqueue_arr(manager)

            # Construct the penultimate intersection of maps for the given request.
            penultimate_map = allocation_map_1D & twilight_map_remaining_flat & \
                nonqueue_map_file_slots_ints & access & custom_map & zero_out_map & \
                respect_past_cadence

            # find when target goes from available to unavailable, for any reason is not available a
            fit_within_night = np.array([1]*manager.n_slots_in_semester)
            slots_needed = manager.slots_needed_for_exposure_dict[name]
            if slots_needed > 1:
                for s in range(manager.n_slots_in_semester - 1):
                    if penultimate_map[s] == 1 and penultimate_map[s+1] == 0:
                        # The -1 below is because target can be started if just fits before unavailable
                        for e in range(slots_needed - 1):
                            fit_within_night[s - e] = 0

            # Construct the ultimate intersection of maps for the given request.
            # Define the slot indices that are available to the request for scheduling.
            available_slots_for_request[name] = penultimate_map & fit_within_night

            # reshape into n_nights_in_semester by n_slots_in_night
            available_slots_for_request[name] = np.reshape(available_slots_for_request[name], \
                                                        (manager.n_nights_in_semester, manager.n_slots_in_night))
            nightly_available_slots = []
            for d in range(len(available_slots_for_request[name])):
                 nightly_available_slots.append(list(np.where( \
                                                        available_slots_for_request[name][d] == 1)[0]))
            available_indices_for_request[name] = nightly_available_slots
    return available_indices_for_request

def construct_custom_map_dict(special_map_file):
    """
    Read in the customize acccessibility maps for unique targets, if exists.

    Args:
        special_map_file (str): path and filename of the precomputed custom maps

    Returns:
        custom_access_maps (dictionary): keys are the starnames and values are the 1D custom maps
    """
    if os.path.exists(special_map_file):
        custom_access_maps = read_accessibilty_map_dict(special_map_file)
    else:
        custom_access_maps = {}
    return custom_access_maps

def construct_zero_out_arr(zero_out_file):
    """
    Read in the list of targets to "zero out", i.e. not allowed to be scheduled only tonight.

    Args:
        zero_out_file (str): path and filename of the stars to be zeroed out

    Returns:
        zero_out_names (array): list of stars to be zeroed out
    """
    if os.path.exists(zero_out_file):
        zero_out = pd.read_csv(zero_out_file)
        zero_out_names = list(zero_out['Target'])
    else:
        zero_out_names = []
    return zero_out_names

def construct_nonqueue_arr(manager):
    """
    Build the 1D non-queue map

    Args:
        manager (obj): a data_admin object

    Returns:
        nonqueue_map_file_slots_ints (array): the 1D array of 1's and 0's indicating nonqueue map
    """
    #print("Incorporating non-queue observations.")
    if os.path.exists(manager.nonqueue_map_file):
        #print("Accommodating time-sensative non-queue observations.")
        nonqueue_map_file_slots_strs = np.loadtxt(manager.nonqueue_map_file, delimiter=',', dtype=str)
        nonqueue_map_file_slots_ints = []
        for i, item in enumerate(nonqueue_map_file_slots_strs):
            holder = []
            for j, item2 in enumerate(nonqueue_map_file_slots_strs[i]):
                if nonqueue_map_file_slots_strs[i][j] == '':
                    holder.append(1)
                else:
                    holder.append(0)
            nonqueue_map_file_slots_ints.append(holder)
        nonqueue_map_file_slots_ints = np.array(nonqueue_map_file_slots_ints).flatten()
        nonqueue_map_file_slots_ints = nonqueue_map_file_slots_ints[manager.today_starting_slot:]
    else:
        nonqueue_map_file_slots_ints = np.array(n_slots_in_semester)
        print("No non-queue observations are scheduled.")
    return nonqueue_map_file_slots_ints

def prepare_allocation_map(manager):
    """
    When not in optimal allocation mode, prepare and construct the allocation map, as well as
    perform the weather loss modeling.

    Args:
        allocation_file (str): the path and filename to the binary map of allocated quarters/nights
        current_day (str): today's date, format YYYY-MM-DD
        semester_length (int): the number of days in a full semester
        DATADIR (str): the path to the pre-computed data files
        all_dates_array (array): the list of calendar dates in the semester
        run_weather_loss (boolean): If True, run the weather loss model. If False, no weather loss.
        output_directory (str): the path to which to save
        available_slots_in_night (array): a 1D array of length n_nights_in_semester where each
                              element is an integer representing the number of slots within that
                              night that during night time
        today_starting_night (int): the day number of today's date in the semester
        n_slots_in_night (int): the number of slots in a single night

    Returns:
        weather_diff_remaining (array): 1D array for remainder of semester where elements are 1 if
                                        the night was allocated, but modelled as a weather loss
        allocation_map_1D (array): a 1D array of length equal to n_slots_in_semester,
                                1's are allocated and 0's are non-allocated slots
        allocation_map_2D (array): a 2D array of shape n_nights_in_semester by n_slots_in_night
                                    with same information as allocation_map_1D
        weathered_map (array): a 2D array of shape n_nights_in_semester by
                                    n_slots_in_night where 1's are nights that were
                                    allocated but weathered out (for plotting purposes)
    """
    # Convert allocation info from human to computer-readable
    allocation_raw = np.loadtxt(manager.allocation_file, dtype=str)
    allocation_remaining = []
    allocation_all = []
    for a in range(manager.semester_length):
        convert = list(map(int, allocation_raw[a][2:]))
        allocation_all.append(convert)
        if a >= manager.all_dates_dict[manager.current_day]:
            allocation_remaining.append(convert)
    manager.allocation_all = allocation_all
    allocation_schedule = np.reshape(allocation_all, (manager.n_nights_in_semester, manager.n_quarters_in_night))
    manager.allocation_map_2D_NQ = allocation_schedule
    manager.allocation_remaining = allocation_remaining

    # Sample out future allocated nights to simulate weather loss based on empirical weather data.
    print("Sampling out weather losses")
    loss_stats_remaining = wh.get_loss_stats(manager)
    allocation_remaining_post_weather_loss, weather_diff_remaining, weather_diff_remaining_1D, \
        days_lost = wh.simulate_weather_losses(manager.allocation_remaining, loss_stats_remaining, \
        covariance=0.14, dont_lose_nights=manager.run_weather_loss, plot=True, outputdir=manager.output_directory)
    allocation_map_1D, allocation_map_2D, weathered_map = \
        build_allocation_map(manager, allocation_remaining_post_weather_loss, weather_diff_remaining)

    manager.weather_diff_remaining = weather_diff_remaining
    wh.write_out_weather_stats(manager, days_lost, manager.allocation_remaining)
    return weather_diff_remaining, allocation_map_1D, allocation_map_2D, weathered_map

def build_allocation_map(manager, allocation_schedule, weather_diff):
    """
    Create the 1D allocation map where allocated slots are designated with a 1
    and non-allocated slots designated with a 0.

    Args:
        allocation_schedule (array): a 2D array of shape n_nights_in_semester by n_quarters_in_night
                                    elements are 1 if that night/quarter is allocated to the queue
                                    and 0 otherwise
        weather_diff (array): same as allocation_schedule but 1 if that night/quarter has been
                              declared to be weathered out (for plotting/tracking purposes)
        available_slots_in_night (array): a 1D array of length n_nights_in_semester where each
                              element is an integer representing the number of slots within that
                              night that during night time
        n_slots_in_night (int): the number of slots in a single night

    Returns:
        allocation_map_1D (array): a 1D array of length equal to n_slots_in_semester,
                                1's are allocated and 0's are non-allocated slots
        allocation_map_2D (array): a 2D array of shape n_nights_in_semester by n_slots_in_night
                                    with same information as allocation_map_1D
        allocation_map_weather_diff_2D (array): a 2D array of shape n_nights_in_semester by
                                            n_slots_in_night where 1's are nights that were
                                            allocated but weathered out (for plotting purposes)
    """
    allocation_map_1D = []
    allocation_map_weathered = []
    manager.available_slots_in_each_night_short = manager.available_slots_in_each_night[manager.today_starting_night:]

    for n in range(len(manager.available_slots_in_each_night_short)):
        allo_night_map = ac.single_night_allocated_slots(allocation_schedule[n],
                                                manager.available_slots_in_each_night_short[n], manager.n_slots_in_night)
        allocation_map_1D.append(allo_night_map)
        weather_night_map = ac.single_night_allocated_slots(weather_diff[n],
                                                manager.available_slots_in_each_night_short[n], manager.n_slots_in_night)
        allocation_map_weathered.append(weather_night_map)

    allocation_map_2D = np.reshape(allocation_map_1D,
                                                (len(manager.available_slots_in_each_night_short), manager.n_slots_in_night))
    allocation_map_weather_diff_2D = np.reshape(allocation_map_weathered,
                                                (len(manager.available_slots_in_each_night_short), manager.n_slots_in_night))
    allocation_map_1D = np.array(allocation_map_1D).flatten()

    return allocation_map_1D, allocation_map_2D, allocation_map_weather_diff_2D

def convert_allocation_array_to_binary(manager):
    """
    Used in optimal allocation, write out the results of the calendar to human readable format.

    Args:
        allocation_map_2D (array): a 2D array of length n_nights_in_semester by n_quarters_in_night,
                                [Q1, Q2, Q3, Q4] where binary elements indicate allocation.
        all_dates_array (dict): the calendar dates of the semester
        filename (str): the path and filename to save the outputted file
    Returns:
        None
    """
    file = open(manager.output_directory + "optimal_allocation_binary_schedule.txt", 'w')
    for i in range(len(manager.all_dates_array)):
        line = manager.all_dates_array[i] + " : " + str(" ".join(map(str, manager.allocation_map_2D_NQ[i-manager.today_starting_night])))
        file.write(str(line) + "\n")
    file.close()

def format_keck_allocation_info(allocation_file):
    """

    Read in allocation file, as downloaded from this site, filter by KPF and the relevant semester:
    https://www2.keck.hawaii.edu/observing/keckSchedule/queryForm.php

    Thanks to H. Isaacson for portion of code which I adapted further.
    Here we are parsing the Time column for which quarters of the night are allocated to us.
    Note: this assumes time is awarded in 0.25 night increments.
    If non-0.25 increment of time exists, assume the whole night is allocated and then
         the best way to deal with this is to later define a new map,
         set the relevant slots equal to 0, then apply map to all targets.

    Args:
        allocation_file (str): the path and filename to the downloaded csv
    Returns:
        allocation (dataframe): a dataframe containing each night's award and start/stop times
    """
    allocation = pd.read_csv(allocation_file)
    # Indicates that all awards are part of the KPF-CC queue
    allocation['Queue'] = ['True']*len(allocation)

    allocation['Time'] = allocation['Time'].str.replace(r'\s+(\d+%)', r'\1', regex=True)
    # Define a regular expression pattern to extract start time, stop time, and percentage
    pattern = r'(\d{2}:\d{2}) - (\d{2}:\d{2}) \((\d{2,3})%\)'
    # Extract the start time, stop time, and percentage using the pattern
    allocation[['Start', 'Stop', 'Percentage']] = allocation['Time'].str.extract(pattern)

    for i, item in enumerate(allocation['Percentage']):
        time_string = allocation['Start'].iloc[i]

        if item == '100':
            allocation['Start'].iloc[i] = 0
            allocation['Stop'].iloc[i]  = 1

        elif item == '75':
            if (time_string.startswith("07")) | (time_string.startswith("08")):
                allocation['Start'].iloc[i] = 0.25
                allocation['Stop'].iloc[i]  = 1
            elif (time_string.startswith("04")) | (time_string.startswith("05")):
                allocation['Start'].iloc[i] = 0
                allocation['Stop'].iloc[i]  = 0.75
            else:
                print("We have a problem, error code 1.")
                print(allocation['Date'].iloc[i],
                            allocation['Start'].iloc[i], allocation['Stop'].iloc[i])

        elif item == '50':
            if (time_string.startswith("04")) | (time_string.startswith("05")):
                allocation['Start'].iloc[i] = 0
                allocation['Stop'].iloc[i]  = 0.5
            elif (time_string.startswith("07")) | (time_string.startswith("08")):
                allocation['Start'].iloc[i] = 0.25
                allocation['Stop'].iloc[i]  = 0.75
            elif (time_string.startswith("10")) | (time_string.startswith("11")):
                allocation['Start'].iloc[i] = 0.5
                allocation['Stop'].iloc[i]  = 1
            else:
                print("We have a problem, error code 2.")
                print(allocation['Date'].iloc[i],
                            allocation['Start'].iloc[i], allocation['Stop'].iloc[i])

        elif item == '25':
            if (time_string.startswith("04")) | (time_string.startswith("05")):
                allocation['Start'].iloc[i] = 0
                allocation['Stop'].iloc[i]  = 0.25
            elif (time_string.startswith("06")) | (time_string.startswith("07")) | \
                            (time_string.startswith("08")):
                allocation['Start'].iloc[i] = 0.25
                allocation['Stop'].iloc[i]  = 0.5
            elif (time_string.startswith("09")) | (time_string.startswith("10")):
                allocation['Start'].iloc[i] = 0.5
                allocation['Stop'].iloc[i]  = 0.75
            elif (time_string.startswith("11")) | (time_string.startswith("12")) | \
                            (time_string.startswith("13")):
                allocation['Start'].iloc[i] = 0.75
                allocation['Stop'].iloc[i]  = 1
            else:
                print("We have a problem, error code 3.")
                print(allocation['Date'].iloc[i],
                                allocation['Start'].iloc[i], allocation['Stop'].iloc[i])
        else:
            print("Non-25% of night increment. Implementing whole night as precaution.")
            print("Date: ", allocation['Date'].iloc[i])
            print("True allocation amount: ", item)
            allocation['Start'].iloc[i] = 0
            allocation['Stop'].iloc[i]  = 1

    return allocation

def quarter_translator(start, stop):
    """
    Map the start/stop fractions to binary map.
    I can't think of a better/automated way to do this other than brute force, sorry.

    Args:
        start (float): the fraction of night where the night begins
        stop (float): the fraction of night where the night ends
    Returns:
        night_map (str): a string representation of the quarters of the night that are allocated
    """

    # Map start/stop to allocation
    if start == 0. and stop == 0.25:
        night_map = "1 0 0 0"
    elif start == 0.25 and stop == 0.5:
        night_map = "0 1 0 0"
    elif start == 0.5 and stop == 0.75:
        night_map = "0 0 1 0"
    elif start == 0.75 and stop == 1.:
        night_map = "0 0 0 1"

    elif start == 0. and stop == 0.5:
        night_map = "1 1 0 0"
    elif start == 0.25 and stop == 0.75:
        night_map = "0 1 1 0"
    elif start == 0.5 and stop == 1.:
        night_map = "0 0 1 1"

    elif start == 0. and stop == 0.75:
        night_map = "1 1 1 0"
    elif start == 0.25 and stop == 1.:
        night_map = "0 1 1 1"

    elif start == 0. and stop == 1.:
        night_map = "1 1 1 1"

    else:
        night_map = "0 0 0 0"
        print("We have a problem, error code 4.")

    return night_map

def quarters_observable(observability_matrix):
    """
    For a single night, using the full night's slots accessibility map,
    determine which of the quarters the given target is observable.

    Args:
        observability_matrix (array): a 1D array of length equal to n_slots_in_night representing a
                                      single night where 1 indicates target is accessible in that
                                      slot, 0 otherwise.
    Returns:
        observable_quarters (array): a 1D array of length n_quarters_in_night where each element
                                     represents if the target is at all observerable in that quarter
                                     (1st to 4th running left to right).
    """
    observable_quarters = []
    quarter = int(len(observability_matrix) / 4)

    q1 = observability_matrix[0:quarter]
    q1_up = np.sum(q1)
    if q1_up > 4:
        observable_quarters.append(1)
    else:
        observable_quarters.append(0)

    q2 = observability_matrix[quarter:2*quarter]
    q2_up = np.sum(q2)
    if q2_up > 4:
        observable_quarters.append(1)
    else:
        observable_quarters.append(0)

    q3 = observability_matrix[2*quarter:3*quarter]
    q3_up = np.sum(q3)
    if q3_up > 4:
        observable_quarters.append(1)
    else:
        observable_quarters.append(0)

    q4 = observability_matrix[3*quarter:]
    q4_up = np.sum(q4)
    if q4_up > 4:
        observable_quarters.append(1)
    else:
        observable_quarters.append(0)

    return observable_quarters

def compute_on_off_for_quarter(quarter_map, quarter):
    """
    Determine the first and last day a target is accessible in a given quarter of the night.
    Primarily for easy lookup plotting purposes later.

    Args:
        observability_matrix (array): a 1D array of length equal to n_slots_in_night representing a
                                      single night where 1 indicates target is accessible in that
                                      slot, 0 otherwise.
    Returns:
        observable_quarters (array): a 1D array of length n_quarters_in_night where each element
                                     represents if the target is at all observerable in that quarter
                                     (1st to 4th running left to right).
    """
    quarter_map = np.array(quarter_map)
    quarter_map_transpose = quarter_map.T
    quarter_map_transpose_tonight = list(quarter_map_transpose[quarter])
    try:
        on = quarter_map_transpose_tonight.index(1)
        quarter_map_transpose_tonight.reverse()
        off = len(quarter_map_transpose_tonight) - quarter_map_transpose_tonight.index(1) - 1
    except:
        on = 0
        off = 0
    return [on, off]

def round_time_to_slot(input_time, slot_size):
    """
    From a given timestamp, determine the timestamp of the nearest slot's start.

    Args:
        input_time (int): timestamp in question, format HH:MM
        slot_size (int): the size of the slots in minutes

    Returns:
        output_time (array): the time, in HH:MM format, of the start of the nearest slot
    """
    hour = int(input_time[:2])
    minute = int(input_time[3:])

    if minute%slot_size == 0:
        return input_time
    else:
        rounded = round(minute,-1)
        if rounded == 0:
            rounded = "00"
        if rounded == 60:
            rounded = "00"
            hour += 1
        if hour < 10:
            hour = "0" + str(hour)
        if hour == 24:
            hour = "00"
        output_time = str(hour) + ":" + str(rounded)
        return output_time

def convert_utc_to_hst(timestring):
    """
    Convert between the UTC time and the HST time. Note Hawaii does not change clocks with
    daylight savings so this is always a straight-forward transformation.

    Args:
        timestring (str): timestamp in question, format HH:MM
        slot_size (int): the size of the slots in minutes

    Returns:
        output_time (array): the time, in HH:MM format, of the start of the nearest slot
    """
    offset = 10 # Timezone difference between UTC and HST.
    hour = int(timestring[:2]) - offset
    if hour < 0:
        hour = 24 + hour
    if hour < 10:
        hour = '0' + str(hour)
    else:
        hour = str(hour)
    return hour + timestring[2:]
