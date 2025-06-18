"""
Module for computing accessibility and observability of targets, including telescope pointing,
twilight, moon constraints, and allocation maps.

This module combines functionality for:
- Basic accessibility calculations
- Comprehensive observability maps
- Allocation and scheduling constraints
"""
import logging
import numpy as np
import pandas as pd
from pathlib import Path
import requests

from astropy.time import Time, TimeDelta
import astropy as apy
import astropy.units as u
import astroplan as apl
import astroq.weather as wh

# Set up logging
logs = logging.getLogger(__name__)

import astroq.management as mn

# Define list of observatories which are currently supported.
# To add an obseratory/telescope, add the Astroplan resolvable name to the list in generate_night_plan
# Then add the same name to the appropriate element of the locations dictionary.
# If necessary, add a location to the locations dictionary, if so, add the location to each of the
# pre_sunrise and post_sunrise dictionaries. Ensure times are 14 hours apart, at least one hour
# before the earliest sunset and one hour after the latest sunrise of the year.
locations = {"Keck Observatory":"Hawaii", "Kitt Peak National Observatory":"Arizona"}
pre_sunset = {'Hawaii':'03:30', 'Arizona':'05:00'}
post_sunrise = {'Hawaii':'17:30', 'Arizona':'14:00'}

# Core mapping functions from maps.py
def produce_ultimate_map(manager, rf, running_backup_stars=False):
    """
    Combine all maps for a target to produce the final map

    Args:
        rf (dataframe): request frame
        running_backup_stars (bool): if true, then do not run the extra map of stepping back in time to account for the starting slot fitting into the night

    Returns:
        available_indices_for_request (dictionary): keys are the starnames and values are a 1D array
                                                  the indices where available_slots_for_request is 1.

    """
    # Prepatory work
    date_formal = Time(manager.semester_start_date,format='iso',scale='utc')
    date = str(date_formal)[:10]
    # Define size of grid
    ntargets = len(rf)
    nnights = manager.semester_length # total nights in semester
    nslots = manager.n_slots_in_night # slots in the night


    # Determine observability
    coords = apy.coordinates.SkyCoord(rf.ra * u.deg, rf.dec * u.deg, frame='icrs')
    targets = apl.FixedTarget(name=rf.starname, coord=coords)
    keck = apl.Observer.at_site(manager.observatory)

    # Set up time grid for one night, first night of the semester
    start = date + "T" + manager.daily_starting_time #  
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
        name = rf.iloc[itarget]['starname']
        if name in manager.database_info_dict:
            if manager.database_info_dict[name][0] != '0000-00-00' and rf.iloc[itarget]['tau_inter'] > 1: # default value if no history, and only valid for cadences beyond every night
                inight_start = manager.all_dates_dict[manager.current_day] - manager.today_starting_night
                inight_stop = min(inight_start + rf.iloc[itarget]['tau_inter'],nnights)
                is_inter[itarget,inight_start:inight_stop,:] = False

    # True if obseravtion occurs at night
    tf = generate_twilight_times(manager.all_dates_array)
    evening = Time(tf['12_evening'], format='jd')[:, None]
    morning = Time(tf['12_morning'], format='jd')[:, None]
    is_night = (evening < slotmidpoint) & (slotmidpoint < morning)
    is_night = is_night.astype(bool)
    is_night = np.ones_like(is_altaz, dtype=bool) & is_night[np.newaxis,:,:]

    is_alloc = manager.allocation_map_2D.astype(bool) # shape = (nnights, nslots) # TODO compute on the fly
    # is_alloc = np.repeat(is_alloc[np.newaxis, :, :], ntargets, axis=0)
    is_alloc = np.ones_like(is_altaz, dtype=bool) & is_alloc[np.newaxis,:,:] # shape = (ntargets, nnights, nslots)

    is_observable_now = np.logical_and.reduce([
        is_altaz,
        is_moon,
        is_night,
        is_inter,
        is_future,
        is_alloc
    ])

    # the target does not violate any of the observability limits in that specific slot, but
    # it does not mean it can be started at the slot. retroactively grow mask to accomodate multishot exposures.
    # Is observable now,
    is_observable = is_observable_now.copy()
    if running_backup_stars == False:
        for itarget in range(ntargets):
            e_val = manager.slots_needed_for_exposure_dict[rf.iloc[itarget]['starname']]
            if e_val == 1:
                continue

            for shift in range(1, e_val):
                # shifts the is_observable_now array to the left by shift
                # for is_observable to be true, it must be true for all shifts
                is_observable[itarget, :, :-shift] &= is_observable_now[itarget, :, shift:]

    access = {
        'is_altaz': is_altaz,
        'is_moon': is_moon,
        'is_night': is_night,
        'is_inter': is_inter,
        'is_future': is_future,
        'is_alloc': is_alloc,
        'is_observable_now': is_observable_now,
        'is_observable': is_observable
    }
    access = np.rec.fromarrays(list(access.values()), names=list(access.keys()))
    return access

def extract_available_indices_from_record(access, manager):
    """
    Extract available indices dictionary from the record array returned by produce_ultimate_map
    
    Args:
        record_array: Record array from produce_ultimate_map containing observability masks
        manager: Manager object containing requests_frame and other data
        
    Returns:
        dict: Dictionary where keys are target names and values are lists of available slots per night
    """
    available_indices_for_request = {}
    ntargets, nnights, nslots = access.shape
    
    # specify indeces of 3D observability array
    itarget, inight, islot = np.mgrid[:ntargets,:nnights,:nslots]

    # define flat table to access maps
    df = pd.DataFrame(
        {'itarget':itarget.flatten(),
         'inight':inight.flatten(),
         'islot':islot.flatten()}
    )
    df['is_observable'] = access.is_observable.flatten()
    available_indices_for_request = {}
    for itarget in range(ntargets):
        temp = []
        for inight in range(nnights):
            temp.append(list(islot[itarget,inight,access.is_observable[itarget,inight,:]]))

        available_indices_for_request[manager.requests_frame.iloc[itarget]['starname']] = temp
    return available_indices_for_request



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

# No-op, we will specify the custom maps in the request_set object.
def construct_custom_map_dict(special_map_file):
    pass

# This will be superceeded by keck text file -> is_alloc function. 
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
        time_string = allocation.loc[i, 'Start']

        if item == '100':
            allocation.loc[i, 'Start'] = 0
            allocation.loc[i, 'Stop']  = 1

        elif item == '75':
            if (time_string.startswith("07")) | (time_string.startswith("08")):
                allocation.loc[i, 'Start'] = 0.25
                allocation.loc[i, 'Stop']  = 1
            elif (time_string.startswith("04")) | (time_string.startswith("05")):
                allocation.loc[i, 'Start'] = 0
                allocation.loc[i, 'Stop']  = 0.75
            else:
                logs.error("We have a problem, error code 1.")
                logs.error(allocation.loc[i, 'Date'],
                            allocation.loc[i, 'Start'], allocation.loc[i, 'Stop'])

        elif item == '50':
            if (time_string.startswith("04")) | (time_string.startswith("05")):
                allocation.loc[i, 'Start'] = 0
                allocation.loc[i, 'Stop']  = 0.5
            elif (time_string.startswith("07")) | (time_string.startswith("08")):
                allocation.loc[i, 'Start'] = 0.25
                allocation.loc[i, 'Stop']  = 0.75
            elif (time_string.startswith("10")) | (time_string.startswith("11")):
                allocation.loc[i, 'Start'] = 0.5
                allocation.loc[i, 'Stop']  = 1
            else:
                logs.error("We have a problem, error code 2.")
                logs.error(allocation.loc[i, 'Date'],
                            allocation.loc[i, 'Start'], allocation.loc[i, 'Stop'])

        elif item == '25':
            if (time_string.startswith("04")) | (time_string.startswith("05")):
                allocation.loc[i, 'Start'] = 0
                allocation.loc[i, 'Stop']  = 0.25
            elif (time_string.startswith("06")) | (time_string.startswith("07")) | \
                            (time_string.startswith("08")):
                allocation.loc[i, 'Start'] = 0.25
                allocation.loc[i, 'Stop']  = 0.5
            elif (time_string.startswith("09")) | (time_string.startswith("10")):
                allocation.loc[i, 'Start'] = 0.5
                allocation.loc[i, 'Stop']  = 0.75
            elif (time_string.startswith("11")) | (time_string.startswith("12")) | \
                            (time_string.startswith("13")):
                allocation.loc[i, 'Start'] = 0.75
                allocation.loc[i, 'Stop']  = 1
            else:
                logs.error("We have a problem, error code 3.")
                logs.error(allocation.loc[i, 'Date'],
                                allocation.loc[i, 'Start'], allocation.loc[i, 'Stop'])
        else:
            logs.error("Non-25% of night increment. Implementing whole night as precaution.")
            logs.error("Date: ", allocation.loc[i, 'Date'])
            logs.error("True allocation amount: ", item)
            allocation.loc[i, 'Start'] = 0
            allocation.loc[i, 'Stop']  = 1

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

def generate_twilight_times(all_dates_array):
    """generate_twilight_times

    Precompute a dataframe of the morning/evening twilight times for each day in the semester

    Args:
        all_dates_array (list): the calendar dates of the semester.
                                Format: YYYY-MM-DD.
    Returns:
        dataframe: 12 deg twilight times
    """

    df = {}
    keck = apl.Observer.at_site('W. M. Keck Observatory')
    times = Time(all_dates_array)
    df['time_utc'] = all_dates_array
    df['12_evening'] = keck.twilight_evening_nautical(times,which='next').jd
    df['12_morning'] = keck.twilight_morning_nautical(times,which='next').jd
    df = pd.DataFrame(df)
    return df

def pull_allocation(date, numdays, instrument):
    base_url = "https://www3.keck.hawaii.edu/api/getSchedule/?cmd=getSchedule"
    params = {
        'date': date,
        'numdays': numdays,
        'instrument':instrument
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        data = response.json()
    else:
        print(f"Error: {response.status_code}, {response.text}")
        return None

    dates = []
    starts = []
    ends = []
    fractions = []
    for i in range(len(data)):
        dates.append(data[i]['Date'])
        starts.append(data[i]['StartTime'])
        ends.append(data[i]['EndTime'])
        fractions.append(data[i]['FractionOfNight'])
    allo = pd.DataFrame({"Date":dates, "Start":starts, "End":ends, "Fraction":fractions})

    all_dates_dict, all_dates_array = mn.build_date_dictionary(date, numdays)
    return allo, all_dates_dict

def make_time_grid(start="03:30", end="17:30", step_min=5, base_date = "2000-01-01"):
    # Convert start and end to minutes since midnight
    h_start, m_start = map(int, start.split(":"))
    h_end, m_end = map(int, end.split(":"))

    t_start_min = h_start * 60 + m_start
    t_end_min = h_end * 60 + m_end

    # Make array of minute offsets
    minute_offsets = np.arange(t_start_min, t_end_min, step_min)
    time_grid = Time(f"{base_date}T00:00:00") + TimeDelta(minute_offsets * 60, format='sec')
    return time_grid

def closest_time_index(input_time, time_grid, base_date = "2000-01-01"):
    # Convert input_time (e.g., "12:43") to a Time object with dummy date
    input_time_obj = Time(f"{base_date}T{input_time}")

    # Find the closest time in the grid
    deltas = np.abs(input_time_obj - time_grid)
    return np.argmin(deltas)

def build_allocation_map(all_dates_dict, allo):
    grid = make_time_grid()
    allo_map = []
    allo_by_date = dict(zip(allo['Date'], allo.index))
    for date_str, day_index in sorted(all_dates_dict.items(), key=lambda x: x[1]):
        one_night = [0] * len(grid)
        if date_str in allo_by_date:
            row = allo.loc[allo_by_date[date_str]]
            start_index = closest_time_index(row['Start'], grid)
            end_index = closest_time_index(row['End'], grid)
            one_night[start_index:end_index + 1] = [1] * (end_index - start_index + 1)
        allo_map.append(one_night)
    allo_map = np.array(allo_map)
    return allo_map

def define_indices_for_requests(manager):
    """
    Using the dictionary of indices where each request is available, define a dataframe for which
    we will use to cut/filter/merge r,d,s tuples
    """
    # Get the record array from produce_ultimate_map
    record_array = produce_ultimate_map(manager, manager.requests_frame)
    
    # Convert to the expected dictionary format
    available_indices_for_request = extract_available_indices_from_record(record_array, manager)

    observability_keys = []
    strategy_keys = []
    for n,row in manager.requests_frame.iterrows():
        name = row['starname']
        if name in list(available_indices_for_request.keys()):
            max_n_visits = int(row['n_intra_max'])
            min_n_visits = int(row['n_intra_min'])
            intra = int(row['tau_intra'])
            nnights = int(row['n_inter_max'])
            inter = int(row['tau_inter'])
            slots_needed = manager.slots_needed_for_exposure_dict[name]
            strategy_keys.append([name, slots_needed, min_n_visits, max_n_visits, intra, nnights, inter])
            for d in range(len(available_indices_for_request[name])):
                for s in available_indices_for_request[name][d]:
                    observability_keys.append((name, d, s))
    strategy = pd.DataFrame(strategy_keys, columns =['id', 't_visit', 'n_intra_min', 'n_intra_max',
                                                     'tau_intra', 'n_inter_max', 'tau_inter'])
    observability = pd.DataFrame(observability_keys, columns =['id', 'd', 's'])
    return strategy, observability


