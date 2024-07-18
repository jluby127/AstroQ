import numpy as np
import matplotlib.pyplot as pt
import pandas as pd
import sys
import math
import time
import pickle
from collections import defaultdict
from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astroplan as apl
import astropy.units as u
dirpath = '/Users/jack/Documents/Github/optimalAllocation/'
sys.path.append(dirpath)
import twilightFunctions as tw


def buildHumanReadableSchedule(Yns, twilightMap, all_targets_frame, nNightsInSemester, nSlotsInNight, AvailableSlotsInGivenNight, nSlotsInQuarter, all_dates_dict, current_day, allocation_map_NS, weathered_map, slotsNeededDict, nonqueueMap_str):
    """
    Write the solved schedule to a CSV file.

    Args:
        Yns (gorubi variable): the solution of the model of shape t targets by s slots in semester. Filled with 1's and 0's if target t is scheduled to slot s.
        twilightMap (array): the 1D array of length s slots in semester where 1's represent slots not in night time and 0's represent slots that are during day/twilight
        all_targets_frame (dataframe): the requests sheet csv with information the targets and their observational strategies
        nNightsInSemester (int): the number of calendar days remaining in the semester
        nSlotsInNight (int): the number of slots within a single night
        AvailableSlotsInGivenNight (array): a 1D array of length n nights in semester where each element is an integer representing how many slots are in night time for that night
        nSlotsInQuarter (int): the number of slots in a quarter. Equal to nSlotsInNight/4
        all_dates_dict (dict): a dictionary where keys are calendar dates and values are the day number within the semester (indexed from starting date to end of semester)
        current_day (string): today's date in format "YYYY-MM-DD"
        allocation_map_NS (array): a 2D array where rows represent a night and columns represent the quarter within that night. Values are 1 if that night/quarter is allocated and 0 if not.
        weathered_map (array): a 1D array of length s slots in semester where elements are 1 if that slot has been modeled as lost to weather and 0 if not
        slotsNeededDict (dict): a dictionary where keys are target names and values are the number of slots needed to complete one exposure
        nonqueueMap_str (string): the path/name of the file which has an array of dimensions n nights in semester by s slots in night where elements denote if the slot is reserved for a non-queue RM observation

    Returns:
        combined_semester_schedule (array): a 2D array of dimensions n nights in semester by s slots in night where elements denote how the slot is used: target, twilight, weather, not scheduled.
    """

    if nonqueueMap_str != 'nofilename.csv':
        nonqueuemap_slots_strs = np.loadtxt(nonqueueMap_str, delimiter=',', dtype=str)

    # Retrieve the solution schedule from the semester solver
    # Track the fullness of the schedule
    semester_schedule = []
    fullslots = 0
    for v in Yns.values():
        if np.round(v.X,0) == 1:
            name = v.VarName[19:][:-1].split(',')[0] # From trial and error, the variable has a name that is always in front of the actual target name, use care to only get the actual target name.
            semester_schedule.append(name)
            fullslots += 1
        else:
            semester_schedule.append("")
    semester_schedule = np.reshape(semester_schedule, (len(all_targets_frame), nNightsInSemester, nSlotsInNight))

    # Initialize empty array of strings so that we are always adding to what already exists and can easily see if slots are over-committed (they shouldn't be)
    combined_semester_schedule = np.empty((len(all_dates_dict),nSlotsInNight), dtype=object)

    # The past does not matter to us here, so specify the days/slots that are to be ignored.
    # In future I might change the size of the array so that it only is from the current day to the end of the semester in size.
    # Then could remove this loop.
    for c in range(all_dates_dict[current_day]):
        for d in range(len(combined_semester_schedule[c])):
            combined_semester_schedule[c][d] = 'Past'

    # Fill the slots with the appropriate committment.
    for n in range(nNightsInSemester):
        for s in range(nSlotsInNight):
            slotallocated = ''
            for t in range(len(all_targets_frame)):
                slotallocated += semester_schedule[t][n][s]
            if nonqueueMap_str != 'nofilename.csv':
                slotallocated += str(nonqueuemap_slots_strs[n + all_dates_dict[current_day]][s])
            if twilightMap[n][s] == 0: # remember that twilight map is "inverted" aka the 1's are time where it is night and the 0's are time where it is day/twilight.
                slotallocated += '*'
            if weathered_map[n][s] == 1 and slotallocated == '':
                slotallocated += 'W'
            if allocation_map_NS[n][s] == 0: #and slotallocated == '':
                slotallocated += 'X'
            combined_semester_schedule[n+all_dates_dict[current_day]][s] = str(slotallocated)

    # Recall that the semester solver puts a 1 only in the slot that starts the exposure for a target.
    # Therefore, many slots are empty because they are part of a multi-slot visit. We need to fill those in apprpriately for ease of human reading and later accounting purposes.
    listnames = list(all_targets_frame['Starname'])
    for n in range(nNightsInSemester-1, -1, -1):
        for s in range(nSlotsInNight-1, -1, -1):
            if combined_semester_schedule[n+all_dates_dict[current_day]][s] in listnames:
                target_name = combined_semester_schedule[n+all_dates_dict[current_day]][s]
                slotsneededperExposure = slotsNeededDict[target_name]
                if slotsneededperExposure > 1:
                    for e in range(1, slotsneededperExposure):
                        combined_semester_schedule[n+all_dates_dict[current_day]][s+e] += target_name

    return combined_semester_schedule


def getGapFillerTargets(scheduleR1, scheduleR2, dayNumber):
    """
    Using the results of the two rounds of scheduling, determine what is different, ie which targets were added in Round 2

    Args:
        scheduleR1 (array): the Round 1 solution for the current day
        scheduleR2 (array): the Round 2 solution for the current day
        dayNumber (int): the day in the semester to investigate (ie. current day calendar date converted to days from start of semester)

    Returns:
        gapFillers (array): a 1D list of the target names for those that were gap fillers
    """
    new = scheduleR2[dayNumber]
    old = scheduleR1[dayNumber]
    gapFillers = [x for x in new if x not in old]
    return gapFillers


def getSemesterInfo(current_day):
    """
    Given today's date, return information about the semester we are currently in.

    Args:
        current_day (string): today's date in format "YYYY-MM-DD"

    Returns:
        semester_start_date (string): date that this semester begins on (civil HST) in format "YYYY-MM-DD"
        semester_end_date (string): date that this semester ends on (civil HST) in format "YYYY-MM-DD"
        semesterLength (int): number of days in the semester
        semesterYear (int): the four digit year
        semesterLetter (string): the "A" or "B" semester designation
    """

    # "A" semester runs from Feb 01 through July 31
    if current_day[5:7] in ['02', '03', '04', '05', '06', '07']:
        semesterLetter = 'A'
    # "B" semester runs from Aug 1 through Jan 01
    elif current_day[5:7] in ['08', '09', '10', '11', '12', '01']:
        semesterLetter = 'B'
    else:
        print("invalid date")
        return None
    semesterYear = current_day[:4]

    if semesterLetter == 'A':
        semester_start_date = semesterYear + '-02-01'
        semester_end_date = semesterYear + '-07-31'
        # check if this is a leap year
        year_limit = 2124
        # Note from Jack Lubin in the year 2024: The year 2124 is arbitrary. In this year, you'll have to update this line for another 100 years.
        # I could have extended this thousands of years in the future,
        # but thought it would be more fun if one day this line breaks every 100 years and someone has to manually update it.
        if int(semesterYear) > year_limit:
            print("Time to update the leap year array!")
        if int(semesterYear) in np.arange(2024, year_limit, 4):
            semesterLength = 182
        else:
            semesterLength = 181
    elif semesterLetter == 'B':
        semester_start_date = semesterYear + '-08-01'
        semester_end_date = str(int(semesterYear) + 1) + '-01-31'
        semesterLength = 184
    return semester_start_date, semester_end_date, semesterLength, semesterYear, semesterLetter


def buildDayDateDictionary(semester_start_date, semesterLength):
    """
    Builds a dictionary where keys are the calendar dates within the semester and values are the corresponding day numbers of the semester.

    Args:
        semester_start_date (string): date that this semester begins on (civil HST) in format "YYYY-MM-DD"
        semesterLength (int): number of days in the semester

    Returns:
        all_dates (dict): a dictionary connecting calendar date to day in the semester.
    """

    all_dates = {}
    date_formal = Time(semester_start_date, format='iso',scale='utc')
    date = str(date_formal)[:10]
    all_dates[date] = 0
    for i in range(1, semesterLength):
        date_formal += TimeDelta(1,format='jd')
        date = str(date_formal)[:10]
        all_dates[date] = i
    return all_dates


def currentDayTracker(current_day, all_dates):
    """
    Computes how many days in the semester have already gone by and updates the nNightsInSemester parameter

    Args:
        current_day (string): today's date in format "YYYY-MM-DD"
        all_dates (dict): a dictionary connecting calendar date to day in the semester.

    Returns:
        daysRemaining (int): the number of days remaining in the semester
    """
    remove_days = all_dates[current_day] - 1
    daysRemaining = len(all_dates) - 1 - remove_days
    return daysRemaining


def roundSlots(exptime, STEP):
    """
    Computes the slots needed for a given exposure but without always rounding up, rather rounding to nearest number of slots needed
    For example: under 5 minute slot sizes, we want a 6 minute exposure to only require one slot (round down) as opposed to 2 slots (round up)
    Note: do not apply this function to any exposure times that are less than 1 slot size.

    Args:
        exptime (int): the exposure time in seconds
        STEP (int): the slot size in seconds

    Returns:
        val (int): the number of slots required for this exposure
    """
    val = int(round(exptime/STEP))
    return val
