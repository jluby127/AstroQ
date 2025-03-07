"""
Module defining functions to perform specific tasks in a helping role. Designed to be only run as
a function call from the generateScript.py script.

Example usage:
    import helper_functions as hf
"""
import os
import math

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta

def get_semester_info(current_day):
    """
    Given today's date, return information about the semester we are currently in.

    Args:
        current_day (string): today's date in format "YYYY-MM-DD"

    Returns:
        semester_start_date (string): first day of the semester (civil HST) in format "YYYY-MM-DD"
        semester_end_date (string): last day of the semester (civil HST) in format "YYYY-MM-DD"
        semester_length (int): number of days in the semester
        semester_year (int): the four digit year
        semester_letter (string): the "A" or "B" semester designation
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
            print("Time to update the leap year array!!! See line 59 in helper_functions.py!!!")
            semester_length = 0
        elif int(semester_year) in np.arange(2024, year_limit, 4):
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
        semester_start_date (string): first day of the semester (civil HST) in format "YYYY-MM-DD"
        semester_length (int): number of days in the semester

    Returns:
        all_dates (dict): a dictionary connecting calendar date to day in the semester.
    """
    all_dates = {}
    date_formal = Time(semester_start_date, format='iso',scale='utc')
    date = str(date_formal)[:10]
    all_dates[date] = 0
    for i in range(1, semester_length):
        date_formal += TimeDelta(1,format='jd')
        date = str(date_formal)[:10]
        all_dates[date] = i
    return all_dates

def current_day_tracker(current_day, all_dates):
    """
    Computes how many days in the semester have already gone by and updates n_nights_in_semester

    Args:
        current_day (string): today's date in format "YYYY-MM-DD"
        all_dates (dict): a dictionary connecting calendar date to day in the semester.

    Returns:
        days_remaining (int): the number of days remaining in the semester
    """
    remove_days = all_dates[current_day] - 1
    days_remaining = len(all_dates) - 1 - remove_days
    return days_remaining

def build_slots_required_dictionary(requests_frame, slot_size, always_round_up_flag=False):
    """
    Computes the slots needed for a given exposure for all requests.
    When always_round_up_flag is false, we can round up or down.
    Example: with 5 minute slot sizes, we want a 6 minute exposure to only require one slot
            (round down) as opposed to 2 slots (round up)

    Args:
        requests_frame (dataframe): the pandas dataframe containing the request information
        slot_size (int): the slot size in minutes
        always_round_up_flag (boolean): if true, slots needed is always larger than exposure_time

    Returns:
        slots_needed_for_exposure_dict (dictionary): keys are the request names, values are the
                        number of slots required for one exposure of that request
    """
    print("Determining slots needed for exposures.")
    slot_size = slot_size * 60 # converting to seconds

    # schedule multi-shots and multi-visits as if a single, long exposure.
    # When n_shots and n_visits are both 1, this reduces down to just the stated exposure time.
    slots_needed_for_exposure_dict = {}
    for n,row in requests_frame.iterrows():
        name = row['Starname']
        exposure_time = row['Nominal Exposure Time [s]']*row['# of Exposures per Visit'] + \
            45*(row['# Visits per Night'] - 1)
        if always_round_up_flag:
            slots_needed_for_exposure = math.ceil(exposure_time/slot_size)
        else:
            if exposure_time > slot_size:
                slots_needed_for_exposure = int(round(exposure_time/slot_size))
            else:
                slots_needed_for_exposure = 1
        slots_needed_for_exposure_dict[name] = slots_needed_for_exposure
    return slots_needed_for_exposure_dict

def slots_required_for_exposure(exposure_time, slot_size, always_round_up_flag=True):
    slot_size = slot_size * 60 # converting to seconds
    if always_round_up_flag:
        slots_needed_for_exposure = math.ceil(exposure_time/slot_size)
    else:
        if exposure_time > slot_size:
            slots_needed_for_exposure = int(round(exposure_time/slot_size))
        else:
            slots_needed_for_exposure = 1
    return slots_needed_for_exposure

def find_indices(arr, start, end):
    """
    Determine the indices in one array that are between two bounds of another array. Used to find
    the indices within arr which contain all the available indices of a given date.

    Args:
        arr (array): a 1D array of all the indices that are available to a given target. The
                    elements of this array are the indices that are available.
        start (int): the first bound
        end (boolean): the last bound

    Returns:
        first_index (int): the starting index number of arr
        last_index (int): the finishing index number of arr
    """
    # Find the first index where the value is greater than "start"
    left = 0
    right = len(arr) - 1
    first_index = -1

    while left <= right:
        mid = (left + right) // 2
        if arr[mid] > start:
            first_index = mid
            right = mid - 1  # Continue searching in the left half
        else:
            left = mid + 1  # Continue searching in the right half

    # Find the last index where the value is less than "end"
    left = 0
    right = len(arr) - 1
    last_index = -1

    while left <= right:
        mid = (left + right) // 2
        if arr[mid] < end:
            last_index = mid
            left = mid + 1  # Continue searching in the right half
        else:
            right = mid - 1  # Continue searching in the left half

    return first_index, last_index

def enforce_dates(filename, all_dates_dict):
    """
    Process the dates that need to be enforced, from the standard format in file
    Args:
        filename (str): the path and filename to open
        all_dates_dict (dict): a dictionary where keys are calendar dates in format (YYYY-MM-DD)
                                       and values of the day number in the semester
    Returns:
        enforced_dates (array): a 2D array where each element contains a 2 element list where
                            element 0 is the date number and element 1 is the quarter number
    """
    enforced_dates = []
    selections = pd.read_csv(filename)
    for s in range(len(selections)):
        night = all_dates_dict[selections['Date'][s]]
        pair = [night, selections['Quarter'][s]]
        enforced_dates.append(pair)
    return enforced_dates

def write_stars_schedule_human_readable(combined_semester_schedule, Yrds, manager, round_info):
    """
    Turns the non-square matrix of the solution into a square matrix and starts the human readable
    solution by filling in the slots where a star's exposre is started.

    Args:
        combined_semester_schedule (array): the human readable solution
        Yns (array): the Gurobi solution with keys of (starname, slot_number) and values 1 or 0.
        manager (obj): a data_admin object
        round_info (str): the solution round, ie "Round1"

    Returns:
        combined_semester_schedule (array): the updated human readable solution
    """

    end_past = manager.all_dates_dict[manager.current_day]*manager.n_slots_in_night
    all_star_schedules = {}
    for name in list(manager.requests_frame['Starname']):
        star_schedule = []
        # buffer the past with zeros
        for p in range(end_past):
            star_schedule.append(0)
        for d in range(manager.n_nights_in_semester):
            for s in range(manager.n_slots_in_night):
                try:
                    value = int(np.round(Yrds[name, d, s].x))
                except KeyError:
                    value = 0.0
                except:
                    print("Error: helper_functions.py line 224: ", name, d, s)
                star_schedule.append(value)
        all_star_schedules[name] = star_schedule

    combined_semester_schedule = combined_semester_schedule.flatten()
    for s in range(manager.n_slots_in_semester):
        slotallocated = ''
        for name in list(manager.requests_frame['Starname']):
            if all_star_schedules[name][s] == 1:
                slotallocated += str(name)
        combined_semester_schedule[s] += str(slotallocated)
    combined_semester_schedule = np.reshape(combined_semester_schedule,
            (manager.semester_length, manager.n_slots_in_night))

    # The semester solver puts a 1 only in the slot that starts the exposure for a target.
    # Therefore, many slots are empty because they are part of a multi-slot visit.
    # Here fill in the multi-slot exposures appropriately for ease of human reading and accounting.
    for n in range(manager.n_nights_in_semester-1-manager.all_dates_dict[manager.current_day], -1, -1):
        for s in range(manager.n_slots_in_night-1, -1, -1):
            if combined_semester_schedule[n+manager.all_dates_dict[manager.current_day]][s] in list(manager.requests_frame['Starname']):
                target_name = combined_semester_schedule[n+manager.all_dates_dict[manager.current_day]][s]
                slots_needed_for_exposure = manager.slots_needed_for_exposure_dict[target_name]
                if slots_needed_for_exposure > 1:
                    for e in range(1, slots_needed_for_exposure):
                        combined_semester_schedule[n+manager.all_dates_dict[manager.current_day]][s+e] += \
                                target_name
    for m in range(len(combined_semester_schedule)):
        # convert the holder string to meaningful string
        if combined_semester_schedule[m][1] == 'supercalifragilisticexpialidocious':
            for l in range(len(combined_semester_schedule[m])):
                combined_semester_schedule[m][l] = 'Past'

    np.savetxt(manager.output_directory + 'raw_combined_semester_schedule_' + round_info + '.txt',
        combined_semester_schedule, delimiter=',', fmt="%s")
    return combined_semester_schedule

def write_available_human_readable(manager, twilight_map, allocation_map_2D, weathered_map):
    """
    Fill in the human readable solution with the non-observation information: non-allocated slots,
    weather loss slots, non-queue slots, twilight slots.

    Args:
        manager (obj): a data_admin object
        twilight_map (array): the 1D array of length n_slots_in_semester where 1's represent slots
                            not in night time and 0's represent slots that are during day/twilight
        allocation_map_2D (array): a 2D array where rows represent a night and columns represent
                                   the quarter within that night. Values are 1 if that
                                   night/quarter is allocated and 0 if not.
        weathered_map (array): a 1D array of length s slots in semester where elements are 1 if
                                that slot has been modeled as lost to weather and 0 if not

    Returns:
        combined_semester_schedule (array): a 2D array of dimensions n_nights_in_semester by
                                            n_slots_in_night where elements denote how the slot is
                                            used: target, twilight, weather, not scheduled.
    """
    if os.path.exists(manager.nonqueue_map_file):
        nonqueuemap_slots_strs = np.loadtxt(manager.nonqueue_map_file, delimiter=',', dtype=str)

    # The past does not matter to us here, so specify the days/slots that are to be ignored.
    end_past = manager.all_dates_dict[manager.current_day]*manager.n_slots_in_night
    combined_semester_schedule = ['']*manager.semester_length*manager.n_slots_in_night
    combined_semester_schedule[0] = 'longwordhereformakingspace'
    for c in range(end_past):
        # for some reason when adding strings within an array, the max length of new string is the
        # length of the longest string in the whole array. So choosing an arbitrary long word
        # as a placeholder. Later I post-process this out.
        combined_semester_schedule[c] += 'supercalifragilisticexpialidocious'
    combined_semester_schedule = np.reshape(combined_semester_schedule,
            (manager.semester_length, manager.n_slots_in_night))

    for n in range(manager.semester_length - manager.all_dates_dict[manager.current_day]):
        for s in range(manager.n_slots_in_night):
            slotallocated = ''
            # remember that twilight map is "inverted": the 1's are time where it is night and the
            # 0's are time where it is day/twilight.
            if twilight_map[n][s] == 0:
                slotallocated += '*'
            if allocation_map_2D[n][s] == 0:
                slotallocated += 'X'
            if weathered_map[n][s] == 1:# and slotallocated == '':
                slotallocated += 'W'
            if os.path.exists(manager.nonqueue_map_file):
                slotallocated += str(nonqueuemap_slots_strs[n + manager.all_dates_dict[manager.current_day], ][s])
            combined_semester_schedule[n + manager.all_dates_dict[manager.current_day], ][s] += str(slotallocated)

    np.savetxt(manager.output_directory + 'raw_combined_semester_schedule_available.txt',
        combined_semester_schedule, delimiter=',', fmt="%s")
    return combined_semester_schedule

def define_slot_index_frame(manager, available_indices_for_request):
    """
    Using the dictionary of indices where each request is available, define a dataframe for which
    we will use to cut/filter/merge r,d,s tuples

    Args:
        requests_frame (dataframe): the pandas dataframe containing the request information
        slots_needed_for_exposure_dict (dictionary): keys are the request names, values are the
                        number of slots required for one exposure of that request
        available_indices_for_request (dictionary): keys are the starnames and values are a 1D array
                                                  the indices where available_slots_for_request is 1.
    Returns:
        Aframe (dataframe): a pandas dataframe containing one row for each valid r,d,s tuple
    """
    # Define the tuples of request and available slot for each request.
    # This becomes the grid over which the Gurobi variables are defined.
    # Now, slots that were never possible for scheduling are not included in the model.
    Aset = []
    Aframe_keys = []
    for n,row in manager.requests_frame.iterrows():
        name = row['Starname']
        n_visits = int(row['# Visits per Night'])
        intra = int(row['Minimum Intra-Night Cadence'])
        inter = int(row['Minimum Inter-Night Cadence'])
        slots_needed = manager.slots_needed_for_exposure_dict[name]
        for d in range(len(available_indices_for_request[name])):
            for s in available_indices_for_request[name][d]:
                Aset.append((name, d, s))
                Aframe_keys.append([name, d, s, slots_needed, n_visits, intra, inter])

    Aframe = pd.DataFrame(Aframe_keys, columns =['r', 'd', 's', 'e', 'v', 'tra', 'ter'])
    schedulable_requests = list(Aframe['r'].unique())
    for name in list(manager.requests_frame['Starname']):
        if name not in schedulable_requests:
            print("WARNING: Target " + name + " has no valid day/slot pairs and therefore is effectively removed from the model.")
    # duplicate columns for easy indexing later
    Aframe['rr'] = Aframe['r']
    Aframe['dd'] = Aframe['d']
    Aframe['ss'] = Aframe['s']

    return Aframe, Aset, schedulable_requests
