import numpy as np
import pandas as pd
import sys
import math
import time
# import pickle
import json
import os
import sys
from collections import defaultdict
from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astroplan as apl
import astropy.units as u

def buildAllocationMap(allocation_schedule, weatherDiff, AvailableSlotsInGivenNight, nSlotsInNight):
    """
    Create the 1D allocation map where allocated/available slots are designated with a 1 and non-allocated slots designated with a 0

    Args:
        allocation_schedule (array): a 2D array of shape nNightsInSemester by nQuartersInNight where each row is a night and the four elements in that row are 1 or 0 if that quarter is allocated to the queue
        weatherDiff (array): a 2D array of shape nNightsInSemester by nQuartersInNight where each row is a night and the four elements in that row are 1 or 0 if that quarter has been declared to be weathered out (for plotting/tracking purposes)
        AvailableSlotsInGivenNight (array): a 1D array of length nNightsInSemester where each element has an integer representing the number of slots within that night that are between twilight times
        nSlotsInNight (int): the number of slots in a single night

    Returns:
        allocation_map (array): a 1D array of length equal to nSlotsInSemester, 1's are allocated and 0's are non-allocated slots
        allocation_map_NS (array): a 2D array of shape nNightsInSemester by nSlotsInNight with same information as allocation_map
        allocation_map_weathered_NS (array): a 2D array of shape nNightsInSemester by nSlotsInNight where 1's are nights that were allocated but weathered out (for plotting purposes)
    """

    allocation_map = []
    allocation_map_weathered = []
    for n in range(len(AvailableSlotsInGivenNight)):
        allo_night_map = buildAllocationMapSingleNight(allocation_schedule[n], AvailableSlotsInGivenNight[n], nSlotsInNight)
        allocation_map.append(allo_night_map)
        weather_night_map = buildAllocationMapSingleNight(weatherDiff[n], AvailableSlotsInGivenNight[n], nSlotsInNight)
        allocation_map_weathered.append(weather_night_map)
    #The NS stands for Night Slot, so this version is made to be 2D whereas the allocation_map itself is a 1D list
    allocation_map_NS = np.reshape(allocation_map, (len(AvailableSlotsInGivenNight), nSlotsInNight))
    allocation_map_weathered_NS = np.reshape(allocation_map_weathered, (len(AvailableSlotsInGivenNight), nSlotsInNight))
    allocation_map = np.array(allocation_map).flatten()
    return allocation_map, allocation_map_NS, allocation_map_weathered_NS


def buildAllocationMapSingleNight(An, AvailableSlotsInTheNight, nSlotsInNight):
    """
    Determine the allocation schedule within a single night

    Args:
        An (array): a 1D array of length nQuartersInNight, 1's indicate that quarter is allocated, 0's otherwise. In order of first to last quarter in the night
        AvailableSlotsInTheNight (array): the number of slots within that night that are that are between twilight times
        nSlotsInNight (int): the number of slots in a single night

    Returns:
        allomap (array): a 1D array of length nSlotsInNight where 1's are the allocated slots that on night n and 0's are not allocated
    """

    extra = AvailableSlotsInTheNight%4
    AvailableSlotsInTheQuarter = int(AvailableSlotsInTheNight/4)
    edge = int((nSlotsInNight - AvailableSlotsInTheNight)/2)

    allomap = [0]*nSlotsInNight
    for i in range(len(An)):
        if An[i] == 1:
            start = edge + i*int(AvailableSlotsInTheQuarter)
            stop = start + int(AvailableSlotsInTheQuarter)
            if i == 3: # ensure allocation goes to up to twilight time
                stop = nSlotsInNight - edge
            for j in range(start, stop):
                allomap[j] = 1
        else:
            start = edge + i*AvailableSlotsInTheQuarter
            stop = start + AvailableSlotsInTheQuarter
            if i == 3: # prevent the last slot from being scheduled (one too many)
                stop -= 1
            for j in range(start, stop):
                allomap[j] = 0
    return allomap

def buildTwilightMap(AvailableSlotsInGivenNight, nSlotsInNight, invert=False):
    """
    Create the 1D twilight times map

    Args:
        AvailableSlotsInGivenNight (array): a 1D array of length nNightsInSemester where each element has an integer representing the number of slots within that night that are between twilight times
        nSlotsInNight (int): the number of slots in a single night
        invert (boolean): if true, flip the 1's and 0's designation. Useful because by default 0's mean we are not in day/twilight and 1's mean we are in night but for adding purposes it is useful to invert this and then np.sum
    Returns:
        nightly_twilight_map (array): a 1D array of length nSlotsInSemester where 1's indicate the slots during night time and 0's during the day time/twilight
    """

    nightly_twilight_map = []
    for i in range(len(AvailableSlotsInGivenNight)):
        if invert:
            quarterslots = [1]*nSlotsInNight
        else:
            quarterslots = [0]*nSlotsInNight
        edge = int((nSlotsInNight - AvailableSlotsInGivenNight[i])/2)
        for j in range(nSlotsInNight):
            if j >= edge and j < nSlotsInNight - edge:
                if invert:
                    quarterslots[j] = 0
                else:
                    quarterslots[j] = 1
        nightly_twilight_map.append(quarterslots)
    return nightly_twilight_map


# def writeAccessibilityMapsDict(accessDict, filename):
#     """
#     Add entries to the accessbility map dictionary
#
#     Args:
#         accessDict (dict): the original accessibility map dictionary. Keys are target names, values are 1D arrays of lenght nSlotsInSemester which indicate if a target is accessible (1) or inaccessible (0) due to telescope pointing limits, seasonal rise/set, and moon-safe distance
#         filename (str): the name of the file to re-save the updated dicationary. Best to match the original name and overwrite.
#     Returns:
#         None
#     """
#
#     # Serialize the dictionary and write it to a file
#     with open(filename, 'wb') as file:
#         pickle.dump(accessDict, file)
#     print("All star accessibility maps writing to pickle.")

# def readAccessibilityMapsDict(filename):
#     """
#     Get the Accessibilty maps for each target from pre-written file
#
#     Args:
#         filename (str): filename where the saved dictionary is stored.
#     Returns:
#         loaded_dict (dict): the python dictionary version of the saved pickle file. A dictionary where keys are target names, values are 1D arrays of lenght nSlotsInSemester which indicate if a target is accessible (1) or inaccessible (0) due to telescope pointing limits, seasonal rise/set, and moon-safe distance
#     """
#
#     # Read the serialized data from the file and deserialize it
#     with open(filename, 'rb') as file:
#         loaded_dict = pickle.load(file)
#     return loaded_dict

def writeAccessibilityMapsDict(accessDict, filename):
    """
    Add entries to the accessbility map dictionary

    Args:
        accessDict (dict): the original accessibility map dictionary. Keys are target names, values are 1D arrays of lenght nSlotsInSemester which indicate if a target is accessible (1) or inaccessible (0) due to telescope pointing limits, seasonal rise/set, and moon-safe distance
        filename (str): the name of the file to re-save the updated dicationary. Best to match the original name and overwrite.
    Returns:
        None
    """
    starnames = list(accessDict.keys())
    stringAccessDict = {}
    for s in range(len(starnames)):
        flat_access_map = np.array(accessDict[starnames[s]]).flatten()
        # because you can't json serialize an array/list or integers
        # must create a string "list". Later decode it back to a real array
        stringmap = '['
        for e in range(len(flat_access_map)):
            stringmap += str(flat_access_map[e]) + ','
        stringmap += ']'
        stringAccessDict[starnames[s]] = stringmap

    with open(filename, 'w') as convert_file:
         convert_file.write(json.dumps(stringAccessDict))
    print("All star accessibility maps written to txt file.")

def readAccessibilityMapsDict(filename):
    """
    Get the Accessibilty maps for each target from pre-written file. Need to reformat the data to array because the array is technically saved as a string.

    Args:
        filename (str): filename where the saved dictionary is stored.
    Returns:
        newDict (dict): the python dictionary version of the saved txt file. A dictionary where keys are target names, values are 1D arrays of lenght nSlotsInSemester which indicate if a target is accessible (1) or inaccessible (0) due to telescope pointing limits, seasonal rise/set, and moon-safe distance
    """
    with open(filename) as f:
        data = f.read()
    # reconstructing the data as a dictionary
    js = json.loads(data)

    newDict = {}
    starnames = list(js.keys())
    for i in range(len(starnames)):
        reformat = [int(x) for x in js[starnames[i]][1:-2].split(',')]
        newDict[starnames[i]] = reformat
    return newDict


def singleTargetAccessible(starname, ra, dec, startdate, nNightsInSemester, STEP, turnonoff=False):
    """
    Compute a target's accessibility map

    Args:
        starname (str): the name of the star (simbad resolvable)
        ra (str): the right ascension in hour angle degrees (ie. 14.33)
        dec (str): the declination in degrees (ie. +75.2)
        startdate (str): the calendar date to begin the calculation from (format is Astropy isot)
        nNightsInSemester (str): the number of calendar nights remaining in the current semester
    Returns:
        slotwindows_accessible (array): a 1D array of length nSlotsInSemester where 1 indicates the target is accessible in that slot and 0 otherwise. Will become the value within a keyed dictionary.
    """

    coords = apy.coordinates.SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
    target = apl.FixedTarget(name=starname, coord=coords)

    date_formal = Time(startdate,format='iso',scale='utc')
    date = str(date_formal)[:10]
    slotwindows_accessible = []
    quartermap = []
    for d in range(nNightsInSemester):
        oneNightAccess = isObservable(date, target, STEP)
        slotwindows_accessible.append(oneNightAccess)
        date_formal += TimeDelta(1,format='jd')
        date = str(date_formal)[:10]
        if turnonoff:
            quartermap.append(quarterObservable(oneNightAccess))

    if turnonoff:
        turns = []
        for q in range(4):
            onoffvals = quarterObservabilityCalculator(quartermap, q)
            turns.append(onoffvals)
        return slotwindows_accessible, turns
    else:
        return slotwindows_accessible


def isObservable(date, target, STEP):
    """
    Compute a target's accessibility map, checking telescope pointing limits and moon-safe distance at the beginning of every slot in a specified night.
    Later we may choose to add additional checks, like airmass below a certain limit.

    Args:
        date (str): the calendar date to compute accessibilty in format 'YYYY-MM-DD'
        target (str): name of the star
        STEP (str): the size of the slots in minutes
    Returns:
        observability_matrix (array): a 1D array of length nSlotsInNight where 1 indicates the target is accessible in that slot and 0 otherwise.
    """

    min_az = 5.3 #naysmith deck direction limits
    max_az = 146.2
    min_alt = 33.3 # Naysmith deck height
    else_min_alt = 30. #non-Naysmith deck height
    max_alt = 85.
    # This is ~20 min before earliest sunset of the year in Hawaii
    # And ~20 min after the latest sunrise of the year in Hawaii
    # This is UTC time.
    # For real runs, ensure that this is 03:30 and 17:30.
    start = date + "T03:30:00"
    startime = Time(start)
    end = date + "T17:30:00"
    endtime = Time(end)
    step = TimeDelta(STEP*60.,format='sec')
    ttemp = np.arange(startime.jd, endtime.jd, step.jd)
    t = Time(ttemp,format='jd')

    keck = apl.Observer.at_site('W. M. Keck Observatory')
    AZ = keck.altaz(t, target, grid_times_targets=True)

    keckapy = apy.coordinates.EarthLocation.of_site('Keck Observatory')
    moon = apy.coordinates.get_moon(Time(t[int(len(t)/2)],format='jd'), keckapy) # test moon/star distance at middle of night

    if moon_safe(moon, (target.ra.rad, target.dec.rad)):
        observability_matrix = []
        for i in range(len(AZ)):
            alt=AZ[i].alt.deg
            az=AZ[i].az.deg
            deck = np.where((az >= min_az) & (az <= max_az))
            deck_height = np.where((alt <= max_alt) & (alt >= min_alt))
            first = np.intersect1d(deck,deck_height)

            not_deck_1 = np.where((az < min_az))
            not_deck_2 = np.where((az > max_az))
            not_deck_height = np.where((alt <= max_alt) & (alt >= else_min_alt))
            second = np.intersect1d(not_deck_1,not_deck_height)
            third = np.intersect1d(not_deck_2,not_deck_height)

            good = np.concatenate((first,second,third))

            observability_matrix = np.zeros(len(t),dtype=int)
            observability_matrix[good] = 1
    else:
        observability_matrix = np.zeros(len(t),dtype=int)

    return observability_matrix


def moon_safe(moon,target_tuple):
    """
    Check that a coordinate is not too close to the moon. Returns true if target is sufficiently far from the moon to allow for observations.

    Args:
        moon (str): a "moon object" from Astropy
        target_tuple (tuple): the RA and Dec of a target star in hourangle and degrees format
    Returns:
        Unnamed boolean
    """

    ang_dist = apy.coordinates.angular_separation(moon.ra.rad,moon.dec.rad,target_tuple[0],target_tuple[1])
    if ang_dist*180/(np.pi) >= 30:
        return True
    else:
        return False


def getStats(accessMap):
    """
    A helper function to return the stats about a single target's accessiblity. These are useful for checking feasibilty and for generating the cadence plot of an individual target

    Args:
        accessMap (array): a 2D array of shape nNightsInSemester by nSlotsInNight where 1's indicate the target is accessible
    Returns:
        days_observable (int): the number of calendar days in the semester where the target is observable
        riseday (int): the number of days from the first day in the semester where the target is first available
        setday (int): the number of days from the first day in the semester where the target is no longer available
    """

    sumAlongDays = np.sum(accessMap, axis=1)
    timeUP = 30 # minutes
    gridpoints = int(timeUP/5)
    observableMask = sumAlongDays > gridpoints
    days_observable = np.sum(observableMask)

    riseday = -1
    i = 0
    while riseday < 0 and i < len(sumAlongDays):
        if sumAlongDays[i] > gridpoints:
            riseday = i
        i += 1
    if i == len(sumAlongDays):
        riseday = i

    setday = -1
    j = len(sumAlongDays)-1
    while setday < 0 and j > 0:
        if sumAlongDays[j] > gridpoints:
            setday = j
        j -= 1
    if j == len(sumAlongDays):
        setday = len(sumAlongDays)

    return days_observable, riseday, setday



def reformatKeckAllocationData(allocation_file):
    """
    Before a new semester, download the allocation schedule for KPF from https://www2.keck.hawaii.edu/observing/keckSchedule/queryForm.php

    Args:
        allocation_file (str): the path and filename to the downloaded csv from website above
    Returns:
        allocation (dataframe): a dataframe containing each night's award and start/stop times
    """

# Read in allocation file, as downloaded from this site filtering by KPF and the relevant semester:
# https://www2.keck.hawaii.edu/observing/keckSchedule/queryForm.php

    # Thanks to H. Isaacson for portion of code which I adapted further.
    # Here we are parsing the Time column for which quarters of the night are allocated to us.
    # Note: this assumes time is awarded in 0.25 night increments.
    # If non-0.25 increment of time exists, assume the whole night is allocated and then
    #      the best way to deal with this is to later define a new map,
    #      set the relevant slots equal to 0, then apply map to all targets.
    # -----------------------------------------------------------------------------------------

    allocation = pd.read_csv(allocation_file)

    # Manually go through this file and add a column called 'Queue' where a 'y'
    # indicates this award is part of the KPF-CC queue and a 'n' indicates it is not.
    # Or automatically add 'y' to all awards.
    allocation['Queue'] = ['y']*len(allocation)

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

        elif (item == '75'):
            if (time_string.startswith("07")) | (time_string.startswith("08")):
                allocation['Start'].iloc[i] = 0.25
                allocation['Stop'].iloc[i]  = 1
            elif (time_string.startswith("04")) | (time_string.startswith("05")):
                allocation['Start'].iloc[i] = 0
                allocation['Stop'].iloc[i]  = 0.75
            else:
                print("We have a problem, error code 1.")
                print(allocation['Date'].iloc[i], allocation['Start'].iloc[i], allocation['Stop'].iloc[i])

        elif (item == '50'):
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
                print(allocation['Date'].iloc[i], allocation['Start'].iloc[i], allocation['Stop'].iloc[i])

        elif (item == '25'):
            if (time_string.startswith("04")) | (time_string.startswith("05")):
                allocation['Start'].iloc[i] = 0
                allocation['Stop'].iloc[i]  = 0.25
            elif (time_string.startswith("06")) | (time_string.startswith("07")) | (time_string.startswith("08")):
                allocation['Start'].iloc[i] = 0.25
                allocation['Stop'].iloc[i]  = 0.5
            elif (time_string.startswith("09")) | (time_string.startswith("10")):
                allocation['Start'].iloc[i] = 0.5
                allocation['Stop'].iloc[i]  = 0.75
            elif (time_string.startswith("11")) | (time_string.startswith("12")) | (time_string.startswith("13")):
                allocation['Start'].iloc[i] = 0.75
                allocation['Stop'].iloc[i]  = 1
            else:
                print("We have a problem, error code 3.")
                print(allocation['Date'].iloc[i], allocation['Start'].iloc[i], allocation['Stop'].iloc[i])
        else:
            print("Special case: non-25% of night increment. Implementing whole night as precaution.")
            print("Date: ", allocation['Date'].iloc[i])
            print("True allocation amount: ", item)
            allocation['Start'].iloc[i] = 0
            allocation['Stop'].iloc[i]  = 1

    return allocation

def translator(start, stop):
    """
    Map the start/stop fractions to binary map. I can't think of a better/automated way to do this other than brute force, sorry.

    Args:
        start (float): the fraction of night where the night begins
        stop (float): the fraction of night where the night ends
    Returns:
        nightMap (str): a string representation of the quarters of the night that are allocated to the queue
    """

    # Map start/stop to allocation
    if start == 0. and stop == 0.25:
        nightMap = "1 0 0 0"
    elif start == 0.25 and stop == 0.5:
        nightMap = "0 1 0 0"
    elif start == 0.5 and stop == 0.75:
        nightMap = "0 0 1 0"
    elif start == 0.75 and stop == 1.:
        nightMap = "0 0 0 1"

    elif start == 0. and stop == 0.5:
        nightMap = "1 1 0 0"
    elif start == 0.25 and stop == 0.75:
        nightMap = "0 1 1 0"
    elif start == 0.5 and stop == 1.:
        nightMap = "0 0 1 1"

    elif start == 0. and stop == 0.75:
        nightMap = "1 1 1 0"
    elif start == 0.25 and stop == 1.:
        nightMap = "0 1 1 1"

    elif start == 0. and stop == 1.:
        nightMap = "1 1 1 1"

    else:
        nightMap = "0 0 0 0"
        print("We have a problem, error code 4.")

    return nightMap

def quarterObservable(observability_matrix):
    """
    For a single night, determine which of the quarters the given target is observable

    Args:
        observability_matrix (array): a 1D array of length equal to nSlotsInNight representing a single night where 1 indicates target is accessible in that slot, 0 otherwise.
    Returns:
        observable_quarters (array): a 1D array of length 4 where each element represents if the target is at all observerable in that quarter (1st to 4th running left to right).
    """
    # Note: this is only a rough calculation. Especially because the division of the 4 quarters is the same for all nights of the year, regardless of true twilight times

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


def quarterObservabilityCalculator(quartermap, quarter):
    quartermap = np.array(quartermap)
    quartermapT = quartermap.T
    quartermapT0 = list(quartermapT[quarter])
    try:
        on0 = quartermapT0.index(1)
        quartermapT0.reverse()
        off0 = len(quartermapT0) - quartermapT0.index(1) - 1
    except:
        on0 = 0
        off0 = 0
    return [on0, off0]

def roundToSlot(time, slotsize):
    hour = int(time[:2])
    minute = int(time[3:])

    if minute%slotsize == 0:
        return time
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
        newtime = str(hour) + ":" + str(rounded)
        return newtime

def convertUTC2HST(timestring):
    # assume timestring is a string of format "HH:MM"
    offset = 10 # Timezone difference between UTC and HST.
    # Offset is always the same over the course of the year because UTC doesn't do daylight savings and neither does HST.
    hour = int(timestring[:2]) - offset
    if hour < 0:
        hour = 24 + hour

    if hour < 10:
        hour = '0' + str(hour)
    else:
        hour = str(hour)
    return hour + timestring[2:]
