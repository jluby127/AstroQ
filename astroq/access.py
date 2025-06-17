"""
Module for computing accessibility and observability of targets, including telescope pointing,
twilight, moon constraints, and allocation maps.

This module combines functionality for:
- Basic accessibility calculations
- Comprehensive observability maps
- Allocation and scheduling constraints
"""
import json
import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as pt
from pathlib import Path
from scipy.interpolate import interp1d

from astropy.time import Time, TimeDelta
import astropy as apy
import astropy.units as u
import astroplan as apl
import astroq.weather as wh

# Set up logging
logs = logging.getLogger(__name__)

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

# Core mapping functions from maps.py
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

def compute_twilight_map(manager):
    """Compute the twilight map for the semester"""
    date_formal = Time(manager.current_day,format='iso',scale='utc')
    date = str(date_formal)[:10]
    date_number = manager.all_dates_dict[manager.current_day]

    keck = apl.Observer.at_site(manager.observatory)
    daily_start = Time(date + "T" + manager.daily_starting_time, location=keck.location)
    slotmidpoint0 = daily_start + (np.arange(manager.n_slots_in_night) + 0.5) *  manager.slot_size * u.min
    days = np.arange(manager.n_nights_in_semester) * u.day

    slotmidpoint = (slotmidpoint0[np.newaxis,:] + days[:,np.newaxis])

    evening = Time(manager.twilight_frame['12_evening'][date_number:], format='jd')[:, None]
    morning = Time(manager.twilight_frame['12_morning'][date_number:], format='jd')[:, None]
    is_night = (evening < slotmidpoint) & (slotmidpoint < morning)
    twilight_2D = is_night.astype(int)
    return twilight_2D

# Allocation and scheduling functions
def construct_custom_map_dict(special_map_file):
    """
    Construct a dictionary of custom maps from a file.
    
    This function reads a file containing custom observation maps for specific targets.
    The file should be formatted with one target per line, where each line contains:
    - The target name
    - A comma-separated list of integers representing the custom map
    
    Args:
        special_map_file (str): Path to the file containing custom maps
        
    Returns:
        dict: A dictionary mapping target names to their custom observation maps
    """
    custom_map_dict = {}
    if os.path.exists(special_map_file):
        with open(special_map_file) as f:
            for line in f:
                if line[0] == '#':
                    continue
                starname = line.split(',')[0]
                custom_map = line.split(',')[1:]
                custom_map = [int(x) for x in custom_map]
                custom_map_dict[starname] = custom_map
    return custom_map_dict

def construct_zero_out_arr(zero_out_file):
    """
    Construct an array of slots to zero out.
    
    This function reads a file containing dates that should be marked as unavailable
    for observation. The file should contain one date per line, with optional comments
    starting with '#'.
    
    Args:
        zero_out_file (str): Path to the file containing dates to zero out
        
    Returns:
        list: A list of dates that should be marked as unavailable
    """
    zero_out_arr = []
    if os.path.exists(zero_out_file):
        with open(zero_out_file) as f:
            for line in f:
                if line[0] == '#':
                    continue
                zero_out_arr.append(line.split(',')[0])
    return zero_out_arr

def construct_nonqueue_arr(manager):
    """
    Construct an array of non-queue slots.
    
    This function creates a 2D array indicating which time slots are not available
    for queue observations. It reads this information from a CSV file specified in
    the manager's nonqueue_file attribute.
    
    Args:
        manager (obj): A data_admin object containing the current state of the scheduler
        
    Returns:
        numpy.ndarray: A 2D array of shape (n_nights, n_slots) where 1 indicates
            a non-queue slot and 0 indicates a queue slot
    """
    nonqueue_arr = np.zeros((manager.n_nights_in_semester, manager.n_slots_in_night), dtype=int)
    if os.path.exists(manager.nonqueue_file):
        nonqueue_frame = pd.read_csv(manager.nonqueue_file)
        for n,row in nonqueue_frame.iterrows():
            date = row['date']
            if date in manager.all_dates_dict:
                inight = manager.all_dates_dict[date] - manager.today_starting_night
                if inight >= 0 and inight < manager.n_nights_in_semester:
                    nonqueue_arr[inight,:] = 1
    return nonqueue_arr

def prepare_allocation_map(manager):
    """
    When not in optimal allocation mode, prepare and construct the allocation map, as well as
    perform the weather loss modeling.

    Args:
        allocation_file (str): the path and filename to the binary map of allocated quarters/nights
        current_day (str): today's date, format YYYY-MM-DD
        semester_length (int): the number of days in a full semester
        DATADIR (str): the path to the pre-computed data files
        all_dates_array (array): the list of calendar dates in the semester
        run_weather_loss (boolean): If True, run the weather loss model. If False, no weather loss.
        output_directory (str): the path to which to save
        available_slots_in_night (array): a 1D array of length n_nights_in_semester where each
                              element is an integer representing the number of slots within that
                              night that during night time
        today_starting_night (int): the day number of today's date in the semester
        n_slots_in_night (int): the number of slots in a single night

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
    # Convert allocation info from human to computer-readable
    allocation_raw = np.loadtxt(manager.allocation_file, dtype=str)
    allocation_remaining = []
    allocation_all = []
    for a in range(manager.semester_length):
        convert = list(map(int, allocation_raw[a][2:]))
        allocation_all.append(convert)
        if a >= manager.all_dates_dict[manager.current_day]:
            allocation_remaining.append(convert)
    manager.allocation_all = allocation_all
    allocation_schedule_all = np.reshape(allocation_all, (manager.semester_length, manager.n_quarters_in_night))
    manager.allocation_map_2D_NQ = allocation_schedule_all
    manager.allocation_remaining = allocation_remaining

    # Sample out future allocated nights to simulate weather loss based on empirical weather data.
    logs.info("Sampling out weather losses")
    loss_stats_remaining = wh.get_loss_stats(manager)
    allocation_remaining_post_weather_loss, weather_diff_remaining, weather_diff_remaining_1D, \
        days_lost = wh.simulate_weather_losses(manager.allocation_remaining, loss_stats_remaining, \
        covariance=0.14, run_weather_loss=manager.run_weather_loss, plot=True, outputdir=manager.output_directory)
    allocation_map_1D, allocation_map_2D, weathered_map = \
        build_allocation_map(manager, allocation_remaining_post_weather_loss, weather_diff_remaining)

    manager.weather_diff_remaining = weather_diff_remaining
    wh.write_out_weather_stats(manager, days_lost, manager.allocation_remaining)
    return weather_diff_remaining, allocation_map_1D, allocation_map_2D, weathered_map

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
        allo_night_map = single_night_allocated_slots(manager.twilight_map_remaining_2D[n], allocation_schedule[n],
                                                manager.available_slots_in_each_night_short[n], manager.n_slots_in_night)
        allocation_map_1D.append(allo_night_map)
        weather_night_map = single_night_allocated_slots(manager.twilight_map_remaining_2D[n], weather_diff[n],
                                                manager.available_slots_in_each_night_short[n], manager.n_slots_in_night)
        allocation_map_weathered.append(weather_night_map)

    allocation_map_2D = np.reshape(allocation_map_1D,
                                                (len(manager.available_slots_in_each_night_short), manager.n_slots_in_night))
    allocation_map_weather_diff_2D = np.reshape(allocation_map_weathered,
                                                (len(manager.available_slots_in_each_night_short), manager.n_slots_in_night))
    allocation_map_1D = np.array(allocation_map_1D).flatten()

    return allocation_map_1D, allocation_map_2D, allocation_map_weather_diff_2D



# Utility functions
def convert_allocation_array_to_binary(manager):
    """
    Convert allocation array to binary format.
    
    This function converts the 2D allocation map into a 1D binary array where:
    - 1 indicates an allocated slot
    - 0 indicates an unallocated slot
    
    Args:
        manager (obj): A data_admin object containing the current state of the scheduler
        
    Returns:
        numpy.ndarray: A 1D binary array representing allocated slots
    """
    # Convert the allocation map to binary format
    allocation_binary = np.zeros(manager.n_slots_in_semester, dtype=int)
    for inight in range(manager.n_nights_in_semester):
        for islot in range(manager.n_slots_in_night):
            if manager.allocation_map_2D[inight,islot] == 1:
                allocation_binary[inight*manager.n_slots_in_night + islot] = 1
    return allocation_binary

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


def single_night_allocated_slots(twilight_tonight, allocated_quarters_tonight, available_slots_in_night,
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

    edge_start = np.argmax(twilight_tonight)
    edge_stop = np.argmax(twilight_tonight[::-1])

    allocated_slots_tonight = [0]*n_slots_in_night
    for i, item in enumerate(allocated_quarters_tonight):
        if allocated_quarters_tonight[i] == 1:
            start = edge_start + i*int(available_slots_in_tonights_quarter)
            stop = start + int(available_slots_in_tonights_quarter)
            if i == 3: # ensure allocation goes to up to twilight time at the end of the night.
                stop = n_slots_in_night - edge_stop
            for j in range(start, stop):
                allocated_slots_tonight[j] = 1
        else:
            start = edge_start + i*available_slots_in_tonights_quarter
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
    twilight_map_all = np.array(build_twilight_map(available_slots_in_each_night,manager.n_slots_in_night, invert=False))
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
        edge = determine_twilight_edge(date, twilight_frame, slot_size, start_time, n_hours_in_night)
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
