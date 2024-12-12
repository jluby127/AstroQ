"""
Module defining functions to perform specific tasks in a helping role. Designed to be only run as
a function call from the generateScript.py script.

Example usage:
    import helper_functions as hf
"""
import os
import math

import numpy as np
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

def slots_required_for_exposure(exposure_time, slot_size, always_round_up_flag=False):
    """
    Computes the slots needed for a given exposure.When always_round_up_flag is false,
    we can round up or down.
    Example: with 5 minute slot sizes, we want a 6 minute exposure to only require one slot
            (round down) as opposed to 2 slots (round up)

    Args:
        exposure_time (int): the exposure time in seconds
        slot_size (int): the slot size in minutes
        always_round_up_flag (boolean): if true, slots needed is always larger than exposure_time

    Returns:
        slots_needed_for_exposure (int): the number of slots required for this exposure
    """
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

def write_stars_schedule_human_readable(combined_semester_schedule, Yns, starnames, semester_length, n_slots_in_night,
                                        n_nights_in_semester, all_dates_dict, slots_needed_dict,
                                        current_day):
    """
    Turns the non-square matrix of the solution into a square matrix and starts the human readable
    solution by filling in the slots where a star's exposre is started.

    Args:
        Yns (array): the Gurobi solution with keys of (starname, slot_number) and values 1 or 0.
        starnames (array): a 1D list of the names of stars that make up one of the keys to Yns
        semester_length (int): the number of days in the full semester
        n_slots_in_night (int): the number of slots in a night
        n_nights_in_semester (int): the number of nights remaining in the semester
        all_dates_dict (dict): a dictionary where keys are calendar dates in format (YYYY-MM-DD)
                               and values of the day number in the semester
        current_day (str): today's date in format YYYY-MM-DD

    Returns:
        first_index (int): the starting index number of arr
        last_index (int): the finishing index number of arr
    """
    end_past = all_dates_dict[current_day]*n_slots_in_night
    n_slots_in_semester = semester_length*n_slots_in_night
    all_star_schedules = {}
    for name in starnames:
        star_schedule = []
        # buffer the past with zeros
        for p in range(end_past):
            star_schedule.append(0)
        for s in range(n_slots_in_semester):
            try:
                value = np.round(Yns[name, s].x)
            except:
                value = 0.0
            star_schedule.append(value)
        all_star_schedules[name] = star_schedule

    combined_semester_schedule = combined_semester_schedule.flatten()
    for s in range(end_past, n_slots_in_semester - end_past):
        slotallocated = ''
        for name in starnames:
            if all_star_schedules[name][s] == 1:
                slotallocated += str(name)
        combined_semester_schedule[s] += str(slotallocated)
    combined_semester_schedule = np.reshape(combined_semester_schedule,
            (semester_length, n_slots_in_night))

    # The semester solver puts a 1 only in the slot that starts the exposure for a target.
    # Therefore, many slots are empty because they are part of a multi-slot visit.
    # Here fill in the multi-slot exposures appropriately for ease of human reading and accounting.
    for n in range(n_nights_in_semester-1-all_dates_dict[current_day], -1, -1):
        for s in range(n_slots_in_night-1, -1, -1):
            if combined_semester_schedule[n+all_dates_dict[current_day]][s] in starnames:
                target_name = combined_semester_schedule[n+all_dates_dict[current_day]][s]
                slots_needed_for_exposure = slots_needed_dict[target_name]
                if slots_needed_for_exposure > 1:
                    for e in range(1, slots_needed_for_exposure):
                        combined_semester_schedule[n+all_dates_dict[current_day]][s+e] += \
                                target_name
    for m in range(len(combined_semester_schedule)):
        # convert the holder string to meaningful string
        if combined_semester_schedule[m][0] == 'supercalifragilisticexpialidocious':
            for l in range(len(combined_semester_schedule[m])):
                combined_semester_schedule[m][l] = 'Past'
    return combined_semester_schedule

def write_available_human_readable(all_dates_dict, current_day, semester_length,
                                n_nights_in_semester, n_slots_in_night,
                                twilight_map, allocation_map_2D, weathered_map, nonqueue_map):
    """
    Fill in the human readable solution with the non-observation information: non-allocated slots,
    weather loss slots, non-queue slots, twilight slots.

    Args:
        combined_semester_schedule (array): a 2D array of dimensions n_nights_in_semester by
                                            n_slots_in_night where elements denote how the slot is
                                            used: target, twilight, weather, not scheduled.
        all_dates_dict (dict): a dictionary where keys are calendar dates and values are the day
                                number within the semester
        current_day (string): today's date in format "YYYY-MM-DD"
        n_nights_in_semester (int): the number of calendar days remaining in the semester
        n_slots_in_night (int): the number of slots within a single night


        twilight_map (array): the 1D array of length n_slots_in_semester where 1's represent slots
                            not in night time and 0's represent slots that are during day/twilight
        allocation_map_2D (array): a 2D array where rows represent a night and columns represent
                                   the quarter within that night. Values are 1 if that
                                   night/quarter is allocated and 0 if not.
        weathered_map (array): a 1D array of length s slots in semester where elements are 1 if
                                that slot has been modeled as lost to weather and 0 if not
        nonqueue_map (string): the path/name of the file which has an array of dimensions n nights
                               in semester by s slots in night where elements denote if the slot is
                               reserved for a non-queue RM observation

    Returns:
        combined_semester_schedule (array): a 2D array of dimensions n_nights_in_semester by
                                            n_slots_in_night where elements denote how the slot is
                                            used: target, twilight, weather, not scheduled.
    """
    if os.path.exists(nonqueue_map):
        nonqueuemap_slots_strs = np.loadtxt(nonqueue_map, delimiter=',', dtype=str)

    # The past does not matter to us here, so specify the days/slots that are to be ignored.
    end_past = all_dates_dict[current_day]*n_slots_in_night
    combined_semester_schedule = ['']*semester_length*n_slots_in_night
    for c in range(end_past):
        # for some reason when adding strings within an array, the max length of new string is the
        # length of the longest string in the whole array. So choosing an arbitrary long word
        # as a placeholder. Later I post-process this out.
        combined_semester_schedule[c] += 'supercalifragilisticexpialidocious'
    combined_semester_schedule = np.reshape(combined_semester_schedule,
            (semester_length, n_slots_in_night))

    for n in range(semester_length - all_dates_dict[current_day]):
        for s in range(n_slots_in_night):
            slotallocated = ''
            # remember that twilight map is "inverted": the 1's are time where it is night and the
            # 0's are time where it is day/twilight.
            if twilight_map[n][s] == 0:
                slotallocated += '*'
            if allocation_map_2D[n][s] == 0:
                slotallocated += 'X'
            if weathered_map[n][s] == 1:# and slotallocated == '':
                slotallocated += 'W'
            if os.path.exists(nonqueue_map):
                slotallocated += str(nonqueuemap_slots_strs[n + all_dates_dict[current_day], ][s])
            combined_semester_schedule[n + all_dates_dict[current_day], ][s] += str(slotallocated)

    return combined_semester_schedule
