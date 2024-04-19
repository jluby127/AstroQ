import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors

from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astroplan as apl
import astropy.units as u

def isObservable(date, target, STEP):
    # This function determines the observability of a given target on a given calendar date
    # Returns a 1D list of length equal to number of slots in a day where each element is
    # either 1 if the target is observable in that slot or 0 if it is not observable.

    min_az = 5.3 #naysmith deck direction limits
    max_az = 146.2
    min_alt = 33.3 # Naysmith deck height
    else_min_alt = 25. #non-Naysmith deck height
    max_alt = 85.

    # ~20 min before earliest sunset of the year in Hawaii
    # And ~20 min after the latest sunrise of the year in Hawaii

    # change back to 03:30 and 17:30 for real data! A 14 hour span.
    # change to 06:30 and 14:30 for simple 8 hour nights, no twilight
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



def defineTargetObservabilityForSemester(all_targets_frame, starname, STEP):

    request_id = all_targets_frame.index[all_targets_frame['Starname']==starname][0]
    ra = all_targets_frame.loc[request_id,'RA']
    dec = all_targets_frame.loc[request_id,'Dec']
    coords = apy.coordinates.SkyCoord(ra * u.deg, dec * u.deg, frame='icrs') #hourangle
    target = apl.FixedTarget(name=starname, coord=coords)

    date_formal = Time('2024-02-01',format='iso',scale='utc')
    date = str(date_formal)[:10]

    slotwindows_accessible = []
    for d in range(183):
        slotwindows_accessible.append(isObservable(date, target, STEP))
        date_formal += TimeDelta(1,format='jd')
        date = str(date_formal)[:10]

    return slotwindows_accessible



def constructAccessibilityDictionary(all_targets_frame, STEP):
    # This function loops through all the targets and constructs the accessibility map for each day
    # then saves this information in a dictionary with key value the star name
    accessDict = {}
    for x in range(len(all_targets_frame)):
        star = all_targets_frame['Starname'][x]
        #print(str(x) + " - star: ", star)
        ra = all_targets_frame['RA'][x]
        dec = all_targets_frame['Dec'][x]
        accessMap = defineTargetObservabilityForSemester(all_targets_frame, star, STEP)
        accessDict[star] = accessMap
    return accessDict
