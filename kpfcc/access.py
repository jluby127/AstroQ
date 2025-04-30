"""
Module for computing acccessibility of requests, including from telescope pointing and twilight
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as pt
from pathlib import Path

from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astropy.units as u
import astroplan as apl

# Define list of observatories which are currently supported.
# To add an obseratory/telescope, add the Astroplan resolvable name to the list in generate_night_plan
# Then add the same name to the appropriate element of the locations dictionary.
# If necessary, add a location to the locations dictionary, if so, add the location to each of the
# pre_sunrise and post_sunrise dictionaries. Ensure times are 14 hours apart, at least one hour
# before the earliest sunset and one hour after the latest sunrise of the year.
locations = {"Keck Observatory":"Hawaii", "Kitt Peak National Observatory":"Arizona"}
pre_sunset = {'Hawaii':'03:30', 'Arizona':'05:00'}
post_sunrise = {'Hawaii':'17:30', 'Arizona':'14:00'}

def construct_access_dict(manager):
    """
    Define the dictionary of accessibility maps

    Args:
        manager (obj): a data_admin object
    Returns:
        default_access_maps (dictionary): keys are the starnames and values are the 1D access maps
    """
    print("Reading pre-computed accessibility maps.")
    rewrite_flag = False
    if os.path.exists(manager.accessibilities_file) == False:
        default_access_maps = {}
        rewrite_flag = True 
    else:
        default_access_maps = read_accessibilty_map_dict(manager.accessibilities_file)
    for n,row in manager.requests_frame.iterrows():
        name = row['starname']
        # check that this target has a pre-computed accessibility map,
        # if not, make one and add it to the file
        try:
            try_read = default_access_maps[name]
        except:
            print(name + " not found in precomputed accessibilty maps. Running now.")
            # Note: the -1 is to account for python indexing
            new_written_access_map = build_single_target_accessibility(manager, name, row['ra'], row['dec'])
            default_access_maps[name] = np.array(new_written_access_map).flatten()
            rewrite_flag = True
    if rewrite_flag:
        # overwrite with the updated file
        write_accessibilty_map_dict(default_access_maps, manager.accessibilities_file)

    return default_access_maps

def single_night_allocated_slots(allocated_quarters_tonight, available_slots_in_night,
                                n_slots_in_night):
    """
    Determine the slots of a single night that are allocated.

    Args:
        allocated_quarters_tonight (array): a 1D array of length n_quarters_in_night,
                                [Q1, Q2, Q3, Q4] where binary elements indicate allocation.
        available_slots_in_night (array): a 1D array of length n_nights_in_semester where each
                              element is an integer representing the number of slots within that
                              night that during night time
        n_slots_in_night (int): the number of slots in a single night

    Returns:
        allocated_slots_tonight (array): a 1D array of length n_slots_in_night
                        where 1's are the allocated slots that on night n and 0's are not allocated
    """

    available_slots_in_tonights_quarter = int(available_slots_in_night/4)
    edge = int((n_slots_in_night - available_slots_in_night)/2)

    allocated_slots_tonight = [0]*n_slots_in_night
    for i, item in enumerate(allocated_quarters_tonight):
        if allocated_quarters_tonight[i] == 1:
            start = edge + i*int(available_slots_in_tonights_quarter)
            stop = start + int(available_slots_in_tonights_quarter)
            if i == 3: # ensure allocation goes to up to twilight time at the end of the night.
                stop = n_slots_in_night - edge
            for j in range(start, stop):
                allocated_slots_tonight[j] = 1
        else:
            start = edge + i*available_slots_in_tonights_quarter
            stop = start + available_slots_in_tonights_quarter
            if i == 3: # prevent the last slot from being scheduled (one too many)
                stop -= 1
            for j in range(start, stop):
                allocated_slots_tonight[j] = 0
    return allocated_slots_tonight

def write_accessibilty_map_dict(access_dict, filename):
    """
    Write entries to the accessbility map dictionary

    Args:
        access_dict (dict): the pre-generated accessibility map dictionary. Keys are target names,
                            values are 1D arrays of length n_slots_in_semester indicating if a
                            target is accessible (1) or inaccessible (0) due to default telescope
                            pointing limits, seasonal rise/set, and moon-safe distance.
        filename (str): the name of the file to re-save the updated dicationary.
                        Best to match the original name and overwrite.
    Returns:
        None
    """
    starnames = list(access_dict.keys())
    string_access_dict = {}
    for s, item in enumerate(starnames):
        flat_access_map = np.array(access_dict[starnames[s]]).flatten()
        # because you can't json serialize an array/list or integers
        # must create a string "list". Later decode it back to a real array
        stringmap = '['
        for e, itme in enumerate(flat_access_map):
            stringmap += str(flat_access_map[e]) + ','
        stringmap += ']'
        string_access_dict[starnames[s]] = stringmap
    with open(filename, 'w') as convert_file:
        convert_file.write(json.dumps(string_access_dict))
    print("All star accessibility maps written to txt file: ", filename)

def read_accessibilty_map_dict(filename):
    """
    Read the accessibilty maps for each target from pre-written file. Must reformat the data to
    array because the array is technically saved as a string.

    Args:
        filename (str): filename where the saved dictionary is stored.
    Returns:
        access_dict (dict): A dictionary where keys are target names, values are 1D arrays of
                            length n_slots_in_semester indicating if a target is accessible (1)
                            or inaccessible (0) due to default telescope pointing limits, seasonal
                            rise/set, and moon-safe distance.
    """
    with open(filename) as f:
        data = f.read()
    # reconstruct the data as a dictionary
    js = json.loads(data)
    access_dict = {}
    starnames = list(js.keys())
    for i, item in enumerate(starnames):
        reformat = [int(x) for x in js[starnames[i]][1:-2].split(',')]
        access_dict[starnames[i]] = reformat
    return access_dict

def build_single_target_accessibility(manager, starname, ra, dec, min_moon_sep = 30, compute_turn_on_off=False):
    """
    Compute a target's accessibility map for the entire semester

    Args:
        starname (str): the name of the star (simbad resolvable)
        ra (str): the right ascension in hour angle degrees (ie. 14.33)
        dec (str): the declination in degrees (ie. +75.2)
        start_date (str): the calendar date to begin the calculation from (format is Astropy isot)
        n_nights_in_semester (str): the number of calendar nights remaining in the current semester
        slot_size (int): the size of the slots in minutes
        observatory (str): the Astropy and Astroplan resolvable name for an observatory
        compute_turn_on_off (boolean): a flag to determine whether or not to compute

    Returns:
        target_accessibility (array): a 1D array of length n_slots_in_semester where 1 indicates
                                     the target is accessible in that slot and 0 otherwise.
    """
    coords = apy.coordinates.SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
    target = apl.FixedTarget(name=starname, coord=coords)
    keck = apl.Observer.at_site(manager.observatory)

    date_formal = Time(manager.current_day,format='iso',scale='utc')
    date = str(date_formal)[:10]
    target_accessibility = []
    quarter_map = []

    # Set up time grid for one night, first night of the semester
    start = date + "T" + manager.daily_starting_time
    daily_start = Time(start, location=keck.location)
    daily_end = daily_start + TimeDelta(1.0, format='jd') # full day from start of first night
    tmp_slot_size = TimeDelta(5.0*u.min) 
    t = Time(np.arange(daily_start.jd, daily_end.jd, tmp_slot_size.jd), format='jd',location=keck.location)

     # Compute base alt/az pattern
    base_altaz = keck.altaz(t, target)

    df = pd.DataFrame(dict(alt=base_altaz.alt.deg, az=base_altaz.az.deg, lst=t.sidereal_time('mean').value))
    df = df.sort_values(by='lst')
    df['is_observable'] = 1 # default to observable
    idx = df.query('5.3 < az < 146.2 and alt < 33.3').index
    df.loc[idx, 'is_observable'] = 0

    # Enforce maximum elvation 
    df.loc[df.alt > 85.0, 'is_observable'] = 0
    # Enforce minumum elevation
    if target.dec.deg < 75 and target.dec.deg > -35:
        df.loc[df.alt < 30.0, 'is_observable'] = 0 # exclude targets below 30 degrees
    else:
        df.loc[df.alt < 18.0, 'is_observable'] = 0 # exclude targets below 18 degrees
    # slot midpoint of first day
    slotmidpoint0 = daily_start + (np.arange(manager.n_slots_in_night) + 0.5) *  manager.slot_size * u.min
    days = np.arange(manager.n_nights_in_semester) * u.day
    slotmidpoint = (slotmidpoint0[:,np.newaxis] + days[np.newaxis,:])

    from scipy.interpolate import interp1d
    f_is_observable = interp1d(
        df.lst, df.is_observable, kind='nearest', fill_value='extrapolate',
        assume_sorted=True
    )
    is_observable = f_is_observable(slotmidpoint.sidereal_time('mean').value)

    # Compute moon accessibility
    moon = apy.coordinates.get_moon(slotmidpoint[0,:] , keck.location)
    ang_dist = apy.coordinates.angular_separation(moon.ra.rad, moon.dec.rad, target.ra.rad, target.dec.rad)
    is_moonsafe = ang_dist*180/(np.pi) >= min_moon_sep

    # compute intersection of observable and moonsafe, should be seperate functions!
    is_accessible = is_observable.astype(bool) & is_moonsafe.astype(bool) # note that this is an element-wise AND
    return is_accessible.flatten()


def get_accessibility_stats(access_map, time_up=30, slot_size=5):
    """
    Compute stats useful for checking feasibilty and for generating the cadence plot of a target

    Args:
        access_map (array): a 2D array of shape n_nights_in_semester by n_slots_in_night where 1's
                            indicate the target is accessible
        time_up (int): the minimum number of minutes a target must be accessible in the night to be
                       considered observable in that night
        slot_size (int): the size of the slots in minutes

    Returns:
        days_observable (int): the number of days in the semester where the target achieves a
                               minimum level of observablility
        rise_day (int): the number of days from semester start where the target is first available
        set_day (int): the number of days from semester start where the target is last available
    """

    sum_along_days = np.sum(access_map, axis=1)
    gridpoints = int(time_up/slot_size)
    observable_mask = sum_along_days > gridpoints
    days_observable = np.sum(observable_mask)

    rise_day = -1
    i = 0
    while rise_day < 0 and i < len(sum_along_days):
        if sum_along_days[i] > gridpoints:
            rise_day = i
        i += 1
    if i == len(sum_along_days):
        rise_day = i

    set_day = -1
    j = len(sum_along_days)-1
    while set_day < 0 and j > 0:
        if sum_along_days[j] > gridpoints:
            set_day = j
        j -= 1
    if j == len(sum_along_days):
        set_day = len(sum_along_days)

    return days_observable, rise_day, set_day


def construct_twilight_map(manager):
    """
    Compute the number of slots per night available, based on strictly twilight times

    Args:
        current_day (str): today's date, format YYYY-MM-DD
        twilight_frame (dataframe): the pandas dataframe containing twilight information
        slot_size (int): the size of a single slot, in minutes
        all_dates_dict (dictionary): keys are calendar days, values are day of the semester
        n_slots_in_night (int): the number of slots in single night, regardless of twilight
        n_nights_in_semester (int): the number of nigths remaining in the semester

    Returns:
        twilight_map_remaining_flat (array): a 1D array where elements are 1 or 0 based if outside or inside of twilight time
        twilight_map_remaining_2D (array): a 2D version of the same information
    """
    print("Determine available slots in each night.")
    # available_slots_in_each_night is a 1D matrix of length nights n
    # This is not a gorubi variable, but a regular python variable
    # Each element will hold an integer which represents the number of slots are available in each
    # quarter of a given night, after accounting for non-observable times due to day/twilight.
    available_slots_in_each_night = []
    for date in list(manager.all_dates_array):
        slots_tonight = determine_twilight_edge(date, manager.twilight_frame, manager.slot_size)
        available_slots_in_each_night.append(slots_tonight)
    twilight_map_all = np.array(build_twilight_map(available_slots_in_each_night,
                                manager.n_slots_in_night, invert=False))
    twilight_map_remaining = twilight_map_all[manager.all_dates_dict[manager.current_day]:]
    twilight_map_remaining_1D = twilight_map_remaining.copy().flatten()
    twilight_map_remaining_2D = np.reshape(twilight_map_remaining_1D, (manager.n_nights_in_semester, manager.n_slots_in_night))
    return twilight_map_remaining_2D, available_slots_in_each_night

def generate_twilight_times(all_dates_array):
    """generate_twilight_times

    Precompute a dataframe of the morning/evening twilight times for each day in the semester

    Args:
        all_dates_array (list): the calendar dates of the semester.
                                Format: YYYY-MM-DD.
    Returns:
        twilight_frame (dataframe): the precomputed twilight times
    """
    keck = apl.Observer.at_site('W. M. Keck Observatory')

    twilight_frame = pd.DataFrame({'time_utc':all_dates_array})
    eighteen_deg_evening = []
    twelve_deg_evening = []
    six_deg_evening = []
    eighteen_deg_morning = []
    twelve_deg_morning = []
    six_deg_morning = []
    twilight_frame['timestamp'] = pd.to_datetime(twilight_frame['time_utc'])
    twilight_frame = twilight_frame.set_index('timestamp')
    # for day in twilight_frame.index.strftime('%Y-%m-%d').tolist():
    for day in twilight_frame.index.strftime(date_format='%Y-%m-%d').tolist():
        as_day = Time(day,format='iso',scale='utc')
        eighteen_deg_evening.append(keck.twilight_evening_astronomical(as_day,which='next'))
        twelve_deg_evening.append(keck.twilight_evening_nautical(as_day,which='next'))
        six_deg_evening.append(keck.twilight_evening_civil(as_day,which='next'))
        eighteen_deg_morning.append(keck.twilight_morning_astronomical(as_day,which='next'))
        twelve_deg_morning.append(keck.twilight_morning_nautical(as_day,which='next'))
        six_deg_morning.append(keck.twilight_morning_civil(as_day,which='next'))

    twilight_frame['18_evening'] = eighteen_deg_evening
    twilight_frame['12_evening'] = twelve_deg_evening
    twilight_frame['6_evening'] = six_deg_evening
    twilight_frame['18_morning'] = eighteen_deg_morning
    twilight_frame['12_morning'] = twelve_deg_morning
    twilight_frame['6_morning'] = six_deg_morning

    return twilight_frame

def get_available_slots_per_night(all_dates_array, twilight_frame, slot_size=5, start_time='03:30',
                                  n_hours_in_night=14):
    """get_available_slots_per_night

    Generates a 1D list of length equal to the number of nights in the semester where each element
    is the integer number of slots that are available to the scheduler in that night

    Args:
        all_dates_array (list): the calendar dates of the semester.
                            Format: YYYY-MM-DD.
        twilight_frame (dataframe): the precomputed twilight times
    Returns:
        None
    """
    available_slots_per_night = []
    for date in range(all_dates_array):
        edge = determine_twilight_edge(date, twilight_frame, slot_size, start_time,
                                       n_hours_in_night)
        available_slots_per_night.append(edge)
    return available_slots_per_night

def determine_twilight_edge(date, twilight_frame, slot_size, start_time='03:30',
                            n_hours_in_night=14):
    """determine_twilight_edge

    Determine how many slots within each night actually take place during night time.

    Args:
        - date (string): the calendar date of the semester. Format: YYYY-MM-DD.
        - twilight_frame (dataframe): the dataframe containing pre-computed
        - slot_size (int): the size of the a slot in minutes
        - start_time (str): the time representing the first slot of the night, format HH:MM
        - n_hours_in_night (int): number of hours in a night
    Returns:
        available_slots_in_night (int): the number of slots in night that are night time
    """
    ind = twilight_frame.index[twilight_frame['time_utc'] == date].tolist()
    twilight_evening_time_isot = twilight_frame.loc[ind,'12_evening'].values[0]
    twilight_evening_time_bjd = Time(twilight_evening_time_isot,format='jd')
    start_timestamp = date + "T" + start_time
    start_timestamp_time = Time(start_timestamp)
    n_slots = int((n_hours_in_night*60)/slot_size)

    bin_edges = []
    bin_edges.append(str(start_timestamp_time))
    for i in range(n_slots):
        start_timestamp_time += TimeDelta(slot_size/(24*60),format='jd')
        if twilight_evening_time_bjd < str(start_timestamp_time):
            edge = i
            break

    available_slots_in_night = n_slots - 2*edge
    return available_slots_in_night

def build_twilight_map(available_slots_in_night, n_slots_in_night, invert=False):
    """
    Create the 1D twilight times map for all slots in the semester

    Args:
        available_slots_in_night (array): a 1D array of length n_nights_in_semester where each
                                        element has an integer representing the number of slots
                                        within that night that are between twilight times
        n_slots_in_night (int): the number of slots in a single night
        invert (boolean): if true, flip the 1's and 0's designation. Useful because by default 0's
                          mean we are not in day/twilight and 1's mean we are in night time but for
                          adding purposes it is useful to invert this schema and then np.sum
    Returns:
        nightly_twilight_map (array): 1D array of length n_slots_in_semester where 1's indicate the
                                    slots during night time and 0's during the day time/twilight
    """
    nightly_twilight_map = []
    for i, item in enumerate(available_slots_in_night):
        if invert:
            quarterslots = [1]*n_slots_in_night
        else:
            quarterslots = [0]*n_slots_in_night
        edge = int((n_slots_in_night - available_slots_in_night[i])/2)
        for j in range(n_slots_in_night):
            if j >= edge and j < n_slots_in_night - edge:
                if invert:
                    quarterslots[j] = 0
                else:
                    quarterslots[j] = 1
        nightly_twilight_map.append(quarterslots)
    return nightly_twilight_map
