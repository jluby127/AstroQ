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
import logging
logs = logging.getLogger(__name__)

from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astropy.units as u
import astroplan as apl
from matplotlib.pylab import *
import gzip
import pickle
def save_dict_compressed(dictionary, filename):
    with gzip.open(filename, 'wb') as f:
        pickle.dump(dictionary, f)
            # Load compressed
def load_dict_compressed(filename):
    with gzip.open(filename, 'rb') as f:
        return pickle.load(f)
from scipy.interpolate import interp1d

import kpfcc.access as ac
import kpfcc.weather as wh

def produce_ultimate_map(manager, rs, running_backup_stars=False, mod=False):
    """
    Combine all maps for a target to produce the final map

    Args:
        requests_frame (dataframe): the pandas dataframe containing request information
        running_backup_stars (bool): if true, then do not run the extra map of stepping back in time to account for the starting slot fitting into the night
    Returns:
        available_indices_for_request (dictionary): keys are the starnames and values are a 1D array
                                                  the indices where available_slots_for_request is 1.

    """
    # Prepatory work
    # rs = manager.requests_frame
    date_formal = Time(manager.current_day,format='iso',scale='utc')
    date = str(date_formal)[:10]

    # Define size of grid
    ntargets = len(rs)
    nnights = manager.n_nights_in_semester
    nslots = manager.n_slots_in_night

    # Determine observability
    coords = apy.coordinates.SkyCoord(rs.ra * u.deg, rs.dec * u.deg, frame='icrs')
    targets = apl.FixedTarget(name=rs.starname, coord=coords)
    keck = apl.Observer.at_site(manager.observatory)

    # Set up time grid for one night, first night of the semester
    start = date + "T" + manager.daily_starting_time
    daily_start = Time(start, location=keck.location)
    daily_end = daily_start + TimeDelta(1.0, format='jd') # full day from start of first night
    tmp_slot_size = TimeDelta(5.0*u.min)
    t = Time(np.arange(daily_start.jd, daily_end.jd, tmp_slot_size.jd), format='jd',location=keck.location)
    t = t[np.argsort(t.sidereal_time('mean'))] # sort by lst

    # Compute base alt/az pattern, shape = (ntargets, nslots)
    coord0 = keck.altaz(t, targets, grid_times_targets=True)
    alt0 = coord0.alt.deg
    az0 = coord0.az.deg

    # 2D mask (n targets, n slots))
    is_altaz0 = np.ones_like(alt0, dtype=bool)
    is_altaz0 &= ~((5.3 < az0 ) & (az0 < 146.2) & (alt0 < 33.3)) # remove nasymth deck
    # remove min elevation for mid declination stars
    ismiddec = ((-30 < targets.dec.deg) & (targets.dec.deg < 75))
    fail = ismiddec[:,np.newaxis] & (alt0 < 30) # broadcast declination array
    is_altaz0 &= ~fail
    # all stars must be between 18 and 85 deg
    fail = (alt0 < 18) | (alt0 > 85)
    is_altaz0 &= ~fail
    # computing slot midpoint for all nights in semester 2D array (slots, nights)
    slotmidpoint0 = daily_start + (np.arange(nslots) + 0.5) *  manager.slot_size * u.min
    # days = np.arange(manager.n_nights_in_semester) * u.day
    days = np.arange(nnights) * u.day
    slotmidpoint = (slotmidpoint0[np.newaxis,:] + days[:,np.newaxis])
    # 3D mask
    is_altaz = np.empty((ntargets, nnights, nslots),dtype=bool)

    # Pre-compute the sidereal times for interpolation
    x = t.sidereal_time('mean').value
    x_new = slotmidpoint.sidereal_time('mean').value
    idx = np.searchsorted(x, x_new, side='left')
    idx = np.clip(idx, 0, len(x)-1) # Handle edge cases
    is_altaz = is_altaz0[:,idx]

    is_future = np.ones((ntargets, nnights, nslots),dtype=bool)
    # inight_current = manager.all_dates_dict[manager.current_day]
    # is_future[:,:inight_current,:] = False

    # Compute moon accessibility
    is_moon = np.ones_like(is_altaz, dtype=bool)
    moon = apy.coordinates.get_moon(slotmidpoint[:,0] , keck.location)
    # Reshaping uses broadcasting to achieve a (ntarget, night) array
    ang_dist = apy.coordinates.angular_separation(
        targets.ra.reshape(-1,1), targets.dec.reshape(-1,1),
        moon.ra.reshape(1,-1), moon.dec.reshape(1,-1),
    ) # (ntargets)
    is_moon = is_moon & (ang_dist.to(u.deg) > 30*u.deg)[:, :, np.newaxis]

    # TODO add in logic for custom maps

    # TODO add in logic to remove stars that are not observable, currently code is a no-op
    # Set to False if internight cadence is violated

    is_inter = np.ones((ntargets, nnights, nslots),dtype=bool)
    for itarget in range(ntargets):
        name = rs.iloc[itarget]['starname']
        if name in manager.database_info_dict:
            if manager.database_info_dict[name][0] != '0000-00-00' and rs.iloc[itarget]['tau_inter'] > 1: # default value if no history, and only valid for cadences beyond every night
                inight_start = manager.all_dates_dict[manager.current_day] - manager.today_starting_night
                inight_stop = min(inight_start + rs.iloc[itarget]['tau_inter'],nnights)
                is_inter[itarget,inight_start:inight_stop,:] = False

    # True if obseravtion occurs at night
    is_night = manager.twilight_map_remaining_2D.astype(bool) # shape = (nnights, nslots)
    is_night = np.ones_like(is_altaz, dtype=bool) & is_night[np.newaxis,:,:]

    is_alloc = manager.allocation_map_2D.astype(bool) # shape = (nnights, nslots)
    # is_alloc = np.repeat(is_alloc[np.newaxis, :, :], ntargets, axis=0)
    is_alloc = np.ones_like(is_altaz, dtype=bool) & is_alloc[np.newaxis,:,:] # shape = (ntargets, nnights, nslots)

    if mod:
        # turn these off when computing semesterly access for fishbowl plot
        is_inter = np.ones_like(is_altaz, dtype=bool)
        is_future = np.ones_like(is_altaz, dtype=bool)
        is_alloc = np.ones_like(is_altaz, dtype=bool)

    is_observable_now = np.logical_and.reduce([
        is_altaz,
        is_moon,
        is_night,
        is_inter,
        is_future,
        is_alloc
    ])

    is_altaz_sums = np.sum(is_altaz, axis=2)
    is_moon_sums = np.sum(is_moon, axis=2)
    is_night_sums = np.sum(is_night, axis=2)
    is_inter_sums = np.sum(is_inter, axis=2)
    is_future_sums = np.sum(is_future, axis=2)
    is_alloc_sums = np.sum(is_alloc, axis=2)
    is_intersection = np.sum(is_observable_now, axis=2)
    summed_array_names = ['is_altaz', 'is_moon', 'is_inter', 'is_future', 'is_alloc', 'is_night', 'intersection']
    summed_arrays = [is_altaz_sums, is_moon_sums, is_inter_sums, is_future_sums, is_alloc_sums, is_night_sums, is_intersection]
    summing_data = {
        str(summed_array_names[i]): arr[:, 0]
        for i, arr in enumerate(summed_arrays)
    }
    sumframe = pd.DataFrame(summing_data, index=list(rs.starname))
    sumframe['Starname'] = list(rs.starname)
    sumframe.index.name = 'starname'
    sumframe.to_csv(manager.output_directory + "/tonight_map_values.csv", index=False)

    # the target does not violate any of the observability limits in that specific slot, but
    # it does not mean it can be started at the slot. retroactively grow mask to accomodate multishot exposures.

    # Is observable now,
    is_observable = is_observable_now.copy()
    if running_backup_stars == False:
        for itarget in range(ntargets):
            e_val = manager.slots_needed_for_exposure_dict[rs.iloc[itarget]['starname']]
            if e_val == 1:
                continue

            for shift in range(1, e_val):
                # shifts the is_observable_now array to the left by shift
                # for is_observable to be true, it must be true for all shifts
                is_observable[itarget, :, :-shift] &= is_observable_now[itarget, :, shift:]

    # specify indeces of 3D observability array
    itarget, inight, islot = np.mgrid[:ntargets,:nnights,:nslots]

    # define flat table to access maps
    df = pd.DataFrame(
        {'itarget':itarget.flatten(),
         'inight':inight.flatten(),
         'islot':islot.flatten()}
    )
    df['is_observable'] = is_observable.flatten()
    available_indices_for_request = {}
    for itarget in range(ntargets):
        temp = []
        for inight in range(nnights):
            temp.append(list(islot[itarget,inight,is_observable[itarget,inight,:]]))

        available_indices_for_request[rs.iloc[itarget]['starname']] = temp
    return available_indices_for_request

def compute_twilight_map(manager):
    date_formal = Time(manager.current_day,format='iso',scale='utc')
    date = str(date_formal)[:10]
    date_number = manager.all_dates_dict[manager.current_day]

    keck = apl.Observer.at_site(manager.observatory)
    # Set up time grid for one night, first night of the semester
    daily_start = Time(date + "T" + manager.daily_starting_time, location=keck.location)
    # slotmidpoint0 = daily_start + (np.arange(manager.semester_length) + 0.5) *  manager.slot_size * u.min
    slotmidpoint0 = daily_start + (np.arange(manager.n_slots_in_night) + 0.5) *  manager.slot_size * u.min
    # days = np.arange(manager.semester_length) * u.day
    days = np.arange(manager.n_nights_in_semester) * u.day

    slotmidpoint = (slotmidpoint0[np.newaxis,:] + days[:,np.newaxis])

    evening = Time(manager.twilight_frame['12_evening'][date_number:], format='jd')[:, None]
    morning = Time(manager.twilight_frame['12_morning'][date_number:], format='jd')[:, None]
    is_night = (evening < slotmidpoint) & (slotmidpoint < morning)
    # is_night = np.ones_like(is_altaz, dtype=bool) & is_night[np.newaxis,:,:]
    twilight_2D = is_night.astype(int)
    return twilight_2D

def mod_produce_ultimate_map(manager, starname):
    """
    Compute maps quickly for plotting purposes only

    Args:
        manager
        starname

    Returns:
        the maps

    """
    # Prepatory work
    rs = manager.requests_frame
    ind = rs.index[rs['starname'] == starname].tolist()

    date_formal = Time(manager.current_day,format='iso',scale='utc')
    date = str(date_formal)[:10]

    # Define size of grid
    ntargets = 1
    nnights = manager.n_nights_in_semester
    nslots = manager.n_slots_in_night

    # Determine observability
    coords = apy.coordinates.SkyCoord(rs.ra[ind] * u.deg, rs.dec[ind] * u.deg, frame='icrs')
    targets = apl.FixedTarget(name=rs.starname[ind], coord=coords)
    keck = apl.Observer.at_site(manager.observatory)

    # Set up time grid for one night, first night of the semester
    start = date + "T" + manager.daily_starting_time
    daily_start = Time(start, location=keck.location)
    daily_end = daily_start + TimeDelta(1.0, format='jd') # full day from start of first night
    tmp_slot_size = TimeDelta(5.0*u.min)
    t = Time(np.arange(daily_start.jd, daily_end.jd, tmp_slot_size.jd), format='jd',location=keck.location)
    t = t[np.argsort(t.sidereal_time('mean'))] # sort by lst

    # Compute base alt/az pattern, shape = (ntargets, nslots)
    coord0 = keck.altaz(t, targets, grid_times_targets=True)
    alt0 = coord0.alt.deg
    az0 = coord0.az.deg

    # 2D mask (n targets, n slots))
    is_altaz0 = np.ones_like(alt0, dtype=bool)
    is_altaz0 &= ~((5.3 < az0 ) & (az0 < 146.2) & (alt0 < 33.3)) # remove nasymth deck
    # remove min elevation for mid declination stars
    ismiddec = ((-30 < targets.dec.deg) & (targets.dec.deg < 75))
    fail = ismiddec[:,np.newaxis] & (alt0 < 30) # broadcast declination array
    is_altaz0 &= ~fail
    # all stars must be between 18 and 85 deg
    fail = (alt0 < 18) | (alt0 > 85)
    is_altaz0 &= ~fail
    # computing slot midpoint for all nights in semester 2D array (slots, nights)
    slotmidpoint0 = daily_start + (np.arange(nslots) + 0.5) *  manager.slot_size * u.min
    # days = np.arange(manager.n_nights_in_semester) * u.day
    days = np.arange(nnights) * u.day
    slotmidpoint = (slotmidpoint0[np.newaxis,:] + days[:,np.newaxis])
    # 3D mask
    is_altaz = np.empty((ntargets, nnights, nslots),dtype=bool)

    # Pre-compute the sidereal times for interpolation
    x = t.sidereal_time('mean').value
    x_new = slotmidpoint.sidereal_time('mean').value
    idx = np.searchsorted(x, x_new, side='left')
    idx = np.clip(idx, 0, len(x)-1) # Handle edge cases
    is_altaz = is_altaz0[:,idx]

    is_future = np.ones((ntargets, nnights, nslots),dtype=bool)
    inight_current = manager.all_dates_dict[manager.current_day]
    is_future[:,:inight_current,:] = False

    # Compute moon accessibility
    is_moon = np.ones_like(is_altaz, dtype=bool)
    moon = apy.coordinates.get_moon(slotmidpoint[:,0] , keck.location)
    # Reshaping uses broadcasting to achieve a (ntarget, night) array
    ang_dist = apy.coordinates.angular_separation(
        targets.ra.reshape(-1,1), targets.dec.reshape(-1,1),
        moon.ra.reshape(1,-1), moon.dec.reshape(1,-1),
    ) # (ntargets)
    is_moon = is_moon & (ang_dist.to(u.deg) > 30*u.deg)[:, :, np.newaxis]

    # TODO add in logic for custom maps

    # TODO add in logic to remove stars that are not observable, currently code is a no-op
    # Set to False if internight cadence is violated

    is_inter = np.ones((ntargets, nnights, nslots),dtype=bool)
    for itarget in range(ntargets):
        name = rs.iloc[itarget]['starname']
        if name in manager.database_info_dict:
            if manager.database_info_dict[name][0] != '0000-00-00': # default value if no history
                inight_start = manager.all_dates_dict[manager.database_info_dict[name][0]]
                inight_stop = min(inight_start + rs.iloc[itarget]['tau_inter'],nnights)
                is_inter[itarget,inight_start:inight_stop,:] = False

    # True if obseravtion occurs at night
    is_night = manager.twilight_map_remaining_2D.astype(bool) # shape = (nnights, nslots)
    is_night = np.ones_like(is_altaz, dtype=bool) & is_night[np.newaxis,:,:]

    is_alloc = manager.allocation_map_2D.astype(bool) # shape = (nnights, nslots)
    is_alloc = np.ones_like(is_altaz, dtype=bool) & is_alloc[np.newaxis,:,:] # shape = (ntargets, nnights, nslots)

    is_observable_now = np.logical_and.reduce([
        is_altaz[0],
        is_moon[0],
        is_night[0],
        is_inter[0],
        is_future[0],
        is_alloc[0]
    ])

    return is_altaz, is_moon, is_night, is_inter, is_future, is_alloc

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
        nonqueue_map_file_slots_ints = np.array([1]*manager.n_slots_in_semester)
    return nonqueue_map_file_slots_ints

def prepare_allocation_map(manager):
    """
    When not in optimal allocation mode, prepare and construct the allocation map, as well as
    perform the weather loss modeling.

    Args:
        manager (obj) - a data_admin object

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
    allo = pull_allocation(manager.semester_start_date, manager.semester_length, 'KPF') # NOTE: later change instrument to 'KPF-CC'
    allocation_2D_all = build_allocation_map(allo)

    # Sample out future allocated nights to simulate weather loss based on empirical weather data.
    logs.info("Sampling out weather losses")
    loss_stats_this_semester = wh.get_loss_stats(manager)

    allocation_2D_post_weather, weather_diff, days_lost = wh.simulate_weather_losses(allocation_2D_all, loss_stats_this_semester, \
        covariance=0.14, run_weather_loss=manager.run_weather_loss, plot=True, outputdir=manager.output_directory)
    wh.write_out_weather_stats(manager, days_lost)

    return allocation_2D_all, allocation_2D_post_weather, weather_diff

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
    manager.available_slots_in_each_night_short = manager.available_slots_in_each_night#[manager.today_starting_night:]

    for n in range(len(manager.available_slots_in_each_night_short)):
        allo_night_map = ac.single_night_allocated_slots(manager.twilight_map_remaining_2D[n], allocation_schedule[n],
                                                manager.available_slots_in_each_night_short[n], manager.n_slots_in_night)
        allocation_map_1D.append(allo_night_map)
        weather_night_map = ac.single_night_allocated_slots(manager.twilight_map_remaining_2D[n], weather_diff[n],
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
                logs.error("We have a problem, error code 1.")
                logs.error(allocation['Date'].iloc[i],
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
                logs.error("We have a problem, error code 2.")
                logs.error(allocation['Date'].iloc[i],
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
                logs.error("We have a problem, error code 3.")
                logs.error(allocation['Date'].iloc[i],
                                allocation['Start'].iloc[i], allocation['Stop'].iloc[i])
        else:
            logs.error("Non-25% of night increment. Implementing whole night as precaution.")
            logs.error("Date: ", allocation['Date'].iloc[i])
            logs.error("True allocation amount: ", item)
            allocation['Start'].iloc[i] = 0
            allocation['Stop'].iloc[i]  = 1

    return allocation

def convert_allocation_info_to_binary(manager, allocation):
    # Generate the binary map for allocations this semester
    # -----------------------------------------------------------------------------------------
    logs.info("Generating binary map of allocated nights/quarters.")
    allocationMap = []
    allocationMap_ints = []
    uniqueDays = 0
    for j in range(len(manager.all_dates_array)):
        datemask = allocation['Date'] == manager.all_dates_array[j]
        oneNight = allocation[datemask]
        if np.sum(datemask) == 0:
            # for night that is not allocated
            map1 = "0 0 0 0"
            map2 = [0, 0, 0, 0]
        elif np.sum(datemask) == 1:
            # for night where only one program is allocated (regardless of length of allocation)
            oneNight.reset_index(inplace=True)
            map1 = quarter_translator(oneNight['Start'][0], oneNight['Stop'][0])
            map2 = [int(map1[0]), int(map1[2]), int(map1[4]), int(map1[6])]
            uniqueDays += 1
        elif np.sum(datemask) >= 1:
            # for night where multiple programs are allocated (regardless of their lengths)
            oneNight.reset_index(inplace=True)
            last = len(oneNight)
            map1 = quarter_translator(oneNight['Start'][0], oneNight['Stop'][last-1])
            map2 = [int(map1[0]), int(map1[2]), int(map1[4]), int(map1[6])]
            uniqueDays += 1
        else:
            logs.error("We have a problem, error code 5.")
            map1 = "0 0 0 0"
            map2 = [0, 0, 0, 0]
        allocationMap.append(map1)
        allocationMap_ints.append(map2)
    logs.info("Total number of quarters allocated: ", np.sum(allocationMap_ints))
    logs.info("Total unique nights allocated: ", uniqueDays)

    # Write the binary allocation map to file
    filename = manager.upstream_path + "inputs/allocation_schedule.txt"
    file = open(filename, 'w')
    for a in range(len(allocationMap)):
        line = manager.all_dates_array[a] + " : " + str(allocationMap[a])
        file.write(str(line) + "\n")
    file.close()

    # Produce and write the start and stop times of each night to file.
    # For the TTP.
    # -----------------------------------------------------------------------------------------
    logs.info("Generate the nightly start/stop times for observing.")
    listdates = list(allocation['Date'])
    processed_dates = []
    starts = []
    stops = []
    for t in range(len(allocation)):
        date = allocation['Date'][t]
        if date in processed_dates:
            continue
        start = allocation['Time'][t][:5]
        datecount = listdates.count(date)
        if datecount > 1:
            stop = allocation['Time'][t+datecount-1][8:13]
        else:
            stop = allocation['Time'][t][8:13]
        processed_dates.append(date)
        starts.append(start)
        stops.append(stop)
    allocation_frame = pd.DataFrame({'Date':processed_dates, 'Start':starts, 'Stop':stops})
    allocation_frame.to_csv(manager.upstream_path + 'inputs/nightly_start_stop_times.csv', index=False)

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
        logs.error("We have a problem, error code 4.")

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
