import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import astroplan as apl
import astropy.units as u
import warnings
warnings.filterwarnings('ignore')
import time


def singleTargetAccessible(starname, ra, dec, date, startTime, stopTime, stepsize):

    coords = apy.coordinates.SkyCoord(ra * 15.0 * u.deg, dec * u.deg, frame='icrs')
    target = apl.FixedTarget(name=starname, coord=coords)

    min_az = 5.3 #naysmith deck direction limits
    max_az = 146.2
    min_alt = 33.3 # Naysmith deck height
    else_min_alt = 30. #non-Naysmith deck height
    max_alt = 85.

    start = date + "T" + startTime
    startime = Time(start)
    end = date + "T" + stopTime
    endtime = Time(end)
    step = TimeDelta(stepsize*60.,format='sec')
    ttemp = np.arange(startime.jd, endtime.jd, step.jd)
    t = Time(ttemp,format='jd')

    keck = apl.Observer.at_site('W. M. Keck Observatory')
    AZ = keck.altaz(t, target, grid_times_targets=True)

    keckapy = apy.coordinates.EarthLocation.of_site('Keck Observatory')
    moon = apy.coordinates.get_moon(Time(t[int(len(t)/2)],format='jd'), keckapy) # test moon/star distance at middle of night

    observability_matrix = np.zeros(len(t),dtype=int)
    if moon_safe(moon, (target.ra.rad, target.dec.rad)):
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
            observability_matrix[good] = 1
    nHoursObservable = np.round((np.sum(observability_matrix)*minReq)/60,1)
    return nHoursObservable

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

def computeYearlyAccessibility(starframe, calendar, starttimes, endtimes, stepsize=10):
    starframe.reset_index(inplace=True)
    backUpObservability = pd.DataFrame({'Starname':starframe['Starname']})
    for n in range(len(calendar)):
        single_night = []
        for s in range(len(starframe)):
            upTime = singleTargetAccessible(starframe['Starname'][s], starframe['RA'][s], starframe['Dec'][s], calendar[n], starttimes[n], endtimes[n], stepsize=stepsize)
            single_night.append(upTime)
        backUpObservability[calendar[n][5:]] = single_night
    return backUpObservability

def getDate(backUpObservability, current_date, minUp):
    upMask = backUpObservability[current_date] > minUp
    starlist = list(backUpObservability['Starname'][upMask])
    return starlist

def getStar(backUpObservability, starname):

    ind = backUpObservability.index[backUpObservability['Starname'] == starname].tolist()[0]
    upTimes = []
    colnames = list(backUpObservability.columns)
    for i in range(1, len(colnames)):
        upTimes.append(backUpObservability[colnames[i]][ind])
    return upTimes

def getStarsForTonight(backup_starlist, backUpObservability, current_date, minimumUpTime):
    starlist = getDate(backUpObservability, current_date, minimumUpTime)

    name = []
    ra = []
    dec = []
    exptime = []
    nshots = []
    totaltime = 0
    for i in range(len(starlist)):
        ind = backup_starlist.index[backup_starlist['Starname'] == starlist[i]].tolist()[0]
        name.append(backup_starlist['Starname'][ind])
        ra.append(backup_starlist['RA'][ind]*15.0)
        dec.append(backup_starlist['Dec'][ind])
        exptime.append(backup_starlist['ExpTime'][ind])
        if backup_starlist['ExpTime'][ind] <= 150.0 and backup_starlist['ExpTime'][ind] > 72.0:
            nshot = 2
        elif backup_starlist['ExpTime'][ind] <= 72.0 and backup_starlist['ExpTime'][ind] > 45.0:
            nshot = 3
        elif backup_starlist['ExpTime'][ind] <= 45.0:
            nshot = 5
        else:
            nshot = 1
        nshots.append(nshot)
        totaltime += backup_starlist['ExpTime'][ind]*nshot + 45*(nshot-1)

    stars4tonight = pd.DataFrame({'Starname':name, "RA":ra, "Dec":dec, 'Exposure Time':exptime,
                             'Exposures Per Visit':nshots, 'Visits In Night':[1]*len(starlist),
                             'Intra_Night_Cadence':[0]*len(starlist), 'Priority':[1]*len(starlist)})

    stars4tonight.sort_values(by='Exposure Time', ascending=False, inplace=True)
    stars4tonight.reset_index(inplace=True)

    print("There are " + str(len(stars4tonight)) + " available backup stars for tonight.")
    print("Amounting to total possible time added (not accounting for slew) of " + str(np.round(totaltime/3600,1)) + " hours.")
    return stars4tonight

def getTimesWorth(starlist, nHours):

    save_stars = list(starlist['Starname'])
    selected_stars = []

    timeUsed = 0.0
    while timeUsed < nHours and len(save_stars) != 0:
        newstar = np.random.choice(save_stars)
        ind1 = save_stars.index(newstar)
        save_stars.pop(ind1)

        ind2 = starlist.index[starlist['Starname'] == newstar].tolist()[0]
        star_row = starlist.loc[starlist['Starname'] == newstar]
        selected_stars.append(star_row)
        timeUsed += (starlist['Exposure Time'][ind2]*starlist['Exposures Per Visit'][ind2] + 45*(starlist['Exposures Per Visit'][ind2]-1))/3600

    selected_stars_frame = pd.concat(selected_stars)
    selected_stars_frame.reset_index(inplace=True, drop=True)
    return selected_stars_frame
