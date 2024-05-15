import numpy as np
import pandas as pd
import sys
import math
import time
import pickle
import os
import sys
from collections import defaultdict
from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astroplan as apl
import astropy.units as u
sys.path.append("/Users/jack/Desktop/")
sys.path.append("/Users/jack/Documents/Github/optimalAllocation/")
# import helperFunctions as hf
# import twilightFunctions as tw

def buildAllocationMap(allocation_schedule, weatherDiff, AvailableSlotsInGivenNight, twilightMaps):
    allocation_map = []
    allocation_map_weathered = []
    for n in range(len(twilightMaps)):
        allo_night_map = buildAllocationMapSingleNight(allocation_schedule[n], AvailableSlotsInGivenNight[n], twilightMaps[n])
        allocation_map.append(allo_night_map)
        weather_night_map = buildAllocationMapSingleNight(weatherDiff[n], AvailableSlotsInGivenNight[n], twilightMaps[n])
        allocation_map_weathered.append(weather_night_map)
    #The NS stands for Night Slot, so this version is made to be 2D whereas the allocation_map itself is a 1D list
    allocation_map_NS = np.reshape(allocation_map, (len(twilightMaps), len(twilightMaps[0])))
    allocation_map_weathered_NS = np.reshape(allocation_map_weathered, (len(twilightMaps), len(twilightMaps[0])))
    return allocation_map, allocation_map_NS, allocation_map_weathered_NS


def buildAllocationMapSingleNight(An, AvailableSlotsInTheNight, twilightMapNight):
    # An (list) = The allocation plan for a single night, 4 elements long and filled with 1's and 0's indicating which quarters are allocated or not.
    extra = AvailableSlotsInTheNight%4
    AvailableSlotsInTheQuarter = int(AvailableSlotsInTheNight/4)
    edge = int((len(twilightMapNight) - AvailableSlotsInTheNight)/2)

    allomap = [0]*len(twilightMapNight) #.copy()
    for i in range(len(An)):
        if An[i] == 1:
            start = edge + i*int(AvailableSlotsInTheQuarter)
            stop = start + int(AvailableSlotsInTheQuarter)
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


def buildTwilightMap(AvailableSlotsInTheNight, nSlotsInNight, invert=False, reorder=False):
    nightly_twilight_map = []
    for i in range(len(AvailableSlotsInTheNight)):
        if invert:
            quarterslots = [1]*nSlotsInNight
        else:
            quarterslots = [0]*nSlotsInNight
        edge = int((nSlotsInNight - AvailableSlotsInTheNight[i])/2)
        for j in range(nSlotsInNight):
            if j > edge and j <= nSlotsInNight - edge:
                if invert:
                    quarterslots[j] = 0
                else:
                    quarterslots[j] = 1
        nightly_twilight_map.append(quarterslots)
    return nightly_twilight_map


def writeAccessibilityMapsDict(accessDict, filename):
    # Serialize the dictionary and write it to a file
    with open(filename, 'wb') as file:
        pickle.dump(accessDict, file)
    print("All star accessibility maps writing to pickle.")

def readAccessibilityMapsDict(filename):
    # Read the serialized data from the file and deserialize it
    with open(filename, 'rb') as file:
        loaded_dict = pickle.load(file)
    return loaded_dict


def singleTargetAccessible(starname, ra, dec, startdate, nNightsInSemester):
    coords = apy.coordinates.SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
    target = apl.FixedTarget(name=starname, coord=coords)

    date_formal = Time(startdate,format='iso',scale='utc')
    date = str(date_formal)[:10]
    slotwindows_accessible = []
    for d in range(nNightsInSemester):
        slotwindows_accessible.append(isObservable(date, target))
        date_formal += TimeDelta(1,format='jd')
        date = str(date_formal)[:10]
    return slotwindows_accessible


def isObservable(date, target):
    min_az = 5.3 #naysmith deck direction limits
    max_az = 146.2
    min_alt = 33.3 # Naysmith deck height
    else_min_alt = 25. #non-Naysmith deck height
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

    return observability_matrix


def getStats(accessMap):

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
