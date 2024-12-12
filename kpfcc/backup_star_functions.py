"""
Module for generating backup bright star list for a night. Designed to be only run as a function
call from the generateScript.py script.

Example usage:
    import backup_star_functions as bsf
"""
from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import astroplan as apl
import astropy.units as u

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

import mapping_functions as mf

def compute_hours_observable_tonight(date, target, slot_size):
    """
    Determine how many hours a star is up on a given date

    Args:
        date (str): the calendar date to compute accessibilty in format 'YYYY-MM-DD'
        target (str): an astroplan FixedTarget object
        slot_size (int): the size of the slots in minutes
    Returns:
        n_hours_observable (int): the number of hours the star is observable
    """
    observability_matrix = mf.is_observable(date, target, slot_size)
    n_hours_observable = np.round((np.sum(observability_matrix)*slot_size)/60,1)
    return n_hours_observable

def compute_yearly_accessibility(star_frame, all_dates_array, slot_size):
    """
    Pre-compute how many hours a star is up each night of the year

    Args:
        star_frame (dateframe): the equivalent of request frame but for the backup stars
        all_dates_array (array): a list of dates in format YYYY-MM-DD
        slot_size (int): the size of the slots in minutes

    Returns:
        backup_observability_frame (dataframe): contains information on each star's number of hours
                                                accessible on each night
    """
    backup_observability_frame = pd.DataFrame({'Starname':star_frame['Starname']})
    for n in range(len(calendar)):
        single_night = []
        for s in range(len(star_frame)):
            coords = apy.coordinates.SkyCoord(star_frame['RA'][s] * u.deg,
                                              star_frame['Dec'][s] * u.deg, frame='icrs')
            target = apl.FixedTarget(name=star_frame['Starname'][s], coord=coords)
            time_up = compute_hours_observable_tonight(all_dates_array[n], target, slot_size)
            single_night.append(time_up)
        backup_observability_frame[calendar[n][5:]] = single_night
    return backup_observability_frame

def get_backup_stars_for_tonight(backup_observability_frame, current_date, minimum_up_time):
    """
    Pre-compute how many hours a star is up each night of the year

    Args:
        backup_observability_frame (dataframe): contains information on each star's number of hours
                                                accessible on each night
        current_date (str): today's date in format YYYY-MM-DD
        minimum_up_time (int): minimum number of hours accessible to be considered observable

    Returns:
        starlist (array): contains the names of the stars that are observable tonight
    """
    up_mask = backup_observability_frame[current_date] > minimum_up_time
    starlist = list(backup_observability_frame['Starname'][up_mask])
    return starlist

def get_stars_for_tonight(backup_frame, backup_observability_frame, current_date,
                          minimum_up_time):
    """
    Get a list of stars for backup tonight

    Args:
        backup_frame (dataframe): contains the info on all the stars from the backup pool
        backup_observability_frame (dataframe): contains information on each star's number of hours
                                                accessible on each night
        current_date (str): today's date in format YYYY-MM-DD
        minimum_up_time (int): minimum number of hours accessible to be considered observable

    Returns:
        starlist (array): contains the names of the stars that are observable tonight
    """
    starlist = get_backup_stars_for_tonight(backup_observability_frame, current_date,
                                            minimum_up_time)
    name = []
    ra = []
    dec = []
    exposure_time = []
    n_shots = []
    total_time = 0
    for i in enumerate(starlist):
        ind = backup_frame.index[backup_frame['Starname'] == starlist[i]].tolist()[0]
        name.append(backup_frame['Starname'][ind])
        ra.append(backup_frame['RA'][ind]*15.0)
        dec.append(backup_frame['Dec'][ind])
        exposure_time.append(backup_frame['exposure_time'][ind])
        if backup_frame['exptime'][ind] <= 150.0 and backup_frame['exptime'][ind] > 72.0:
            n_shot = 2
        elif backup_frame['exptime'][ind] <= 72.0 and backup_frame['exptime'][ind] > 45.0:
            n_shot = 3
        elif backup_frame['exptime'][ind] <= 45.0:
            n_shot = 5
        else:
            n_shot = 1
        n_shots.append(n_shot)
        total_time += backup_frame['exposure_time'][ind]*n_shot + 45*(n_shot-1)

    stars_for_tonight = pd.DataFrame({'Starname':name, "RA":ra, "Dec":dec,
                             'Exposure Time':exposure_time, 'Exposures Per Visit':n_shots,
                             'Visits In Night':[1]*len(starlist),
                             'Intra_Night_Cadence':[0]*len(starlist), 'Priority':[1]*len(starlist)})

    stars_for_tonight.sort_values(by='Exposure Time', ascending=False, inplace=True)
    stars_for_tonight.reset_index(inplace=True)

    print("There are " + str(len(stars_for_tonight)) + " available backup stars for tonight.")
    print("Amounting to total possible time added (not accounting for slew) of " + \
        str(np.round(total_time/3600,1)) + " hours.")
    return stars_for_tonight

def get_times_worth(backup_list_tonight, n_hours_needed):
    """
    Get a list of stars for backup tonight that fills a given amount of hours worth of time

    Args:
        backup_list_tonight (array): contains the stars from the wider backup pool selected as
                                     observable for tonight
        n_hours_needed (int): the number of hours of open shutter time to fill

    Returns:
        selected_stars_frame (array): contains the subset of stars from the wider backup pool
                                      selected as  observable for tonight
    """
    save_stars = list(starlist['Starname'])
    selected_stars = []
    time_used = 0.0
    while time_used < n_hours_needed and len(save_stars) != 0:
        chosen_star = np.random.choice(save_stars)
        ind1 = save_stars.index(chosen_star)
        save_stars.pop(ind1)

        ind2 = starlist.index[starlist['Starname'] == chosen_star].tolist()[0]
        star_row = starlist.loc[starlist['Starname'] == chosen_star]
        selected_stars.append(star_row)
        time_used += (starlist['Exposure Time'][ind2]*starlist['Exposures Per Visit'][ind2] + \
                        45*(starlist['Exposures Per Visit'][ind2]-1))/3600
    selected_stars_frame = pd.concat(selected_stars)
    selected_stars_frame.reset_index(inplace=True, drop=True)

    return selected_stars_frame
