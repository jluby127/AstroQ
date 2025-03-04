"""
Module for working with the various maps

This module organizes and builds the maps used by the semester solver to block out large portions
of time for one of many reasons. Designed to be only run as a function call from
the generateScript.py script.

Example usage:
    import mapping_functions as mf
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as pt

from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astropy.units as u
import astroplan as apl

def produce_ultimate_map(requests_frame, allocation_map_1D, twilight_map_remaining_flat, 
                         default_access_maps, custom_access_maps, zero_out_names,
                         nonqueue_map_file_slots_ints, database_info_dict,
                         slots_needed_for_exposure_dict, all_dates_dict, current_day, slot_size):
    """
    Combine all maps for a target to produce the final map

    Args:
        requests_frame (dataframe): the pandas dataframe containing request information

    Returns:
        available_slots_for_request (dictionary): keys are the starnames and values are the 1D array
                                                  of the final map
        available_indices_for_request (dictionary): keys are the starnames and values are a 1D array
                                                  the indices where available_slots_for_request is 1.

    """
    print("Build unique star available slot indices.")
    available_slots_for_request = {}
    available_indices_for_request = {}
    for name in requests_frame['Starname']:
        accessibility_r = default_access_maps[name]
        access = accessibility_r[today_starting_slot:]

        if name in list(custom_access_maps.keys()):
            custom_map = custom_access_maps[name][today_starting_slot:]
        else:
            custom_map = np.array([1]*n_slots_in_semester)

        zero_out_map = np.array([1]*n_slots_in_semester)
        if name in zero_out_names:
            zero_out_map[:n_slots_in_night] = np.array([0]*n_slots_in_night)

        respect_past_cadence = np.ones(n_slots_in_semester, dtype=np.int64)
        if database_info_dict != {}:
            date_last_observed = database_info_dict[name][0]
            if date_last_observed != '0000-00-00':
                date_last_observed_number = all_dates_dict[date_last_observed]
                today_number = all_dates_dict[current_day]
                diff = today_number - date_last_observed_number
                if diff < int(row['Minimum Inter-Night Cadence']):
                    block_upcoming_days = int(row['Minimum Inter-Night Cadence']) - diff
                    respect_past_cadence[:block_upcoming_days*n_slots_in_night] = 0

        # Determine which nights a multi-visit request is allowed to be attempted to be scheduled.
        # This equation is a political decision and can be modified.
        # It states that for each visit, after the intra-night cadence time has elapsed,
        # we require a 90 minute window within which to allow for scheduling the next visit.
        # We then assume the next visit is scheduled at the very end of this 90 minute window,
        # which then restarts the clock for any additional visits.
        minimum_time_required = ((int(row['# Visits per Night']) - 1)* \
            (int(row['Minimum Intra-Night Cadence']) + 1.5))*3600 #convert hours to seconds
        minimum_slots_required = hf.slots_required_for_exposure(minimum_time_required, slot_size)
        no_multi_visit_observations = []
        for d in range(n_nights_in_semester):
            start = d*n_slots_in_night
            end = start + n_slots_in_night
            possible_open_slots = np.sum(allocation_map_1D[start:end] & \
                                        twilight_map_remaining_flat[start:end] & access[start:end])
            if possible_open_slots < minimum_slots_required:
                no_multi_visit_observations.append([0]*n_slots_in_night)
            else:
                no_multi_visit_observations.append([1]*n_slots_in_night)
        no_multi_visit_observations = np.array(no_multi_visit_observations)

        # Construct the penultimate intersection of maps for the given request.
        penultimate_map = allocation_map_1D & twilight_map_remaining_flat & \
            nonqueue_map_file_slots_ints & access & custom_map & zero_out_map & \
            respect_past_cadence

        # find when target goes from available to unavailable, for any reason is not available a
        fit_within_night = np.array([1]*n_slots_in_semester)
        slots_needed = slots_needed_for_exposure_dict[name]
        if slots_needed > 1:
            for s in range(n_slots_in_semester - 1):
                if penultimate_map[s] == 1 and penultimate_map[s+1] == 0:
                    # The -1 below is because target can be started if just fits before unavailable
                    for e in range(slots_needed - 1):
                        fit_within_night[s - e] = 0

        # Construct the ultimate intersection of maps for the given request.
        # Define the slot indices that are available to the request for scheduling.
        available_slots_for_request[name] = penultimate_map & fit_within_night

        # reshape into n_nights_in_semester by n_slots_in_night
        available_slots_for_request[name] = np.reshape(available_slots_for_request[name], \
                                                    (n_nights_in_semester, n_slots_in_night))
        nightly_available_slots = []
        for d in range(len(available_slots_for_request[name])):
             nightly_available_slots.append(list(np.where( \
                                                    available_slots_for_request[name][d] == 1)[0]))
        available_indices_for_request[name] = nightly_available_slots

        return available_indices_for_request

def construct_access_dict(accessibilities_file, requests_frame):
    """
    Define the dictionary of accessibility maps

    Args:
        accessibilities_file (str): path and filename of the precomputed access maps
        requests_frame (dataframe): the pandas dataframe containing request information

    Returns:
        default_access_maps (dictionary): keys are the starnames and values are the 1D access maps
    """
    print("Reading pre-comupted accessibility maps.")
    rewrite_flag = False
    default_access_maps = read_accessibilty_map_dict(accessibilities_file)
    for n,row in requests_frame.iterrows():
        name = row['Starname']
        # check that this target has a pre-computed accessibility map,
        # if not, make one and add it to the file
        try:
            try_read = default_access_maps[name]
        except:
            print(name + " not found in precomputed accessibilty maps. Running now.")
            # Note: the -1 is to account for python indexing
            new_written_access_map = build_single_target_accessibility(name, row['RA'], row['Dec'],
                                               semester_start_date, semester_length-1, slot_size)
            default_access_maps[name] = np.array(new_written_access_map).flatten()
            rewrite_flag = True
    if rewrite_flag:
        # overwrite with the updated file
        mf.write_accessibilty_map_dict(default_access_maps, accessibilities_file)

    return default_access_maps

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

def construct_nonqueue_arr(nonqueue_map_file):
    """
    Build the 1D non-queue map

    Args:
        nonqueue_map_file (str): path and filename of the non-queue map

    Returns:
        nonqueue_map_file_slots_ints (array): the 1D array of 1's and 0's indicating nonqueue map
    """
    print("Incorporating non-queue observations.")
    if os.path.exists(nonqueue_map_file):
        print("Constraint: accommodate time-sensative non-queue observations.")
        nonqueue_map_file_slots_strs = np.loadtxt(nonqueue_map_file, delimiter=',', dtype=str)
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
        nonqueue_map_file_slots_ints = nonqueue_map_file_slots_ints[today_starting_slot:]
    else:
        nonqueue_map_file_slots_ints = np.array(n_slots_in_semester)
        print("No non-queue observations are scheduled.")
    return nonqueue_map_file_slots_ints

def prepare_allocation_map(allocation_file, current_day, semester_length, DATADIR, \
            all_dates_array, output_directory):
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
    print("Preparing allocation map.")
    # Convert allocation info from human to computer-readable
    allocation_raw = np.loadtxt(allocation_file, dtype=str)
    allocation_remaining = []
    allocation_all = []
    for a in range(semester_length):
        convert = list(map(int, allocation_raw[a][2:]))
        allocation_all.append(convert)
        if a >= all_dates_dict[current_day]:
            allocation_remaining.append(convert)

    # Sample out future allocated nights to simulate weather loss based on empirical weather data.
    print("Sampling out weather losses")
    fn = os.path.join(DATADIR,"maunakea_weather_loss_data.csv")
    historical_weather_data = pd.read_csv(fn)
    loss_stats_remaining = []
    for i, item in enumerate(all_dates_array):
        ind = historical_weather_data.index[historical_weather_data['Date'] == \
            all_dates_array[i][5:]].tolist()[0]
        loss_stats_remaining.append(historical_weather_data['% Total Loss'][ind])

    allocation_remaining_post_weather_loss, weather_diff_remaining, weather_diff_remaining_1D, \
        days_lost = mf.simulate_weather_losses(allocation_remaining, loss_stats_remaining, \
        covariance=0.14, dont_lose_nights=run_weather_loss, plot=True, outputdir=output_directory)
    allocation_map_1D, allocation_map_2D, weathered_map = \
        mf.build_allocation_map(allocation_remaining_post_weather_loss, weather_diff_remaining,
        available_slots_in_each_night[today_starting_night:], n_slots_in_night)

    mf.write_out_weather_stats(all_dates_dict, current_day, days_lost, allocation_remaining, \
                                output_directory)
    return weather_diff_remaining, allocation_map_1D, allocation_map_2D, weathered_map

def build_allocation_map(allocation_schedule, weather_diff, available_slots_in_night,
                        n_slots_in_night):
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
    for n, item in enumerate(available_slots_in_night):
        allo_night_map = single_night_allocated_slots(allocation_schedule[n],
                                                available_slots_in_night[n], n_slots_in_night)
        allocation_map_1D.append(allo_night_map)
        weather_night_map = single_night_allocated_slots(weather_diff[n],
                                                available_slots_in_night[n], n_slots_in_night)
        allocation_map_weathered.append(weather_night_map)

    allocation_map_2D = np.reshape(allocation_map_1D,
                                                (len(available_slots_in_night), n_slots_in_night))
    allocation_map_weather_diff_2D = np.reshape(allocation_map_weathered,
                                                (len(available_slots_in_night), n_slots_in_night))
    allocation_map_1D = np.array(allocation_map_1D).flatten()

    return allocation_map_1D, allocation_map_2D, allocation_map_weather_diff_2D

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

def simulate_weather_losses(allocation_remaining, loss_stats, covariance=0.14, \
                            dont_lose_nights=False, plot=False, outputdir=None):
    """
    Simulate nights totally lost to weather usine historical data

    Args:
        allocation_remaining (array): a 1D array of length n_nights_in_semester where 1's represent
                                      allocated night and 0's represent non-allocated night
        loss_stats (array): 1D array of length n_nights_in_semester where elements are the
                            percent of the time that night is totally lost to weather
        covariance (float): the added percent that tomorrow will be lost if today is lost
        dont_lose_nights (boolean): a flag that turns off weather simulation entirely
        plot (boolean): a flag to visualize the allocation and weather loss
        outputdir (str): path to save the plot

    Returns:
        allocation_remaining_post_losses (array): 1's represent the night still is allocated
        weather_diff_remaining_2D (array): 1's represent the night still weathered
        weather_diff_remaining_1D (array): 1's represent the night still weathered
    """
    previous_day_was_lost = False
    allocation_remaining_post_losses = allocation_remaining.copy()
    counter = 0
    if dont_lose_nights == False:
        days_lost = []
        # start at 1 because we never want tonight to be simulated as total loss
        for i in range(1, len(allocation_remaining_post_losses)):
            value_to_beat = loss_stats[i]
            if previous_day_was_lost:
                value_to_beat += covariance
            roll_the_dice = np.random.uniform(0.0,1.0)

            if roll_the_dice < value_to_beat:
                # the night is simulated a total loss
                allocation_remaining_post_losses[i] = [0,0,0,0]
                previous_day_was_lost = True
                counter += 1
                days_lost.append(1)
                if plot:
                    pt.axvline(i, color='r')
            else:
                previous_day_was_lost = False
                days_lost.append(0)
                if plot:
                    pt.axvline(i, color='k')
    else:
        print('Pretending weather is always good!')
        days_lost = [0]*(len(allocation_remaining_post_losses)-1)

    weather_diff_remaining_2D = np.array(allocation_remaining) - \
                                np.array(allocation_remaining_post_losses)
    weather_diff_remaining_1D = weather_diff_remaining_2D.flatten()
    print("Total nights simulated as weathered out: " + str(counter) + " of " + \
                str(len(allocation_remaining_post_losses)) + " nights remaining.")
    if plot:
        size=15
        pt.xlabel("Days in Semester from Current Day", fontsize=size)
        pt.tick_params(axis="both", labelsize=size)
        pt.savefig(outputdir + "weather_loss_visualization.png", dpi=200,
                                    bbox_inches='tight', facecolor='w')

    return allocation_remaining_post_losses, weather_diff_remaining_2D, weather_diff_remaining_1D, days_lost

def write_out_weather_stats(all_dates_dict, current_day, days_lost, allocation_remaining, output_directory):
    """
    Write out data on the results of the weather simulation. For comparing acrossing MCMC simulations

    Args:
        all_dates_dict (dict): keys are the dates of the semester, values are the day number of semester
        current_day (star): today's date in format YYYY-MM-DD
        days_lost (array): a binary 1D array designating if the night was lost. Note that the length
                            of this array is n_nights_in_semester-1 which is smaller than len(all_dates_dict)
        allocation_remaining (array): a binary 2D array of dates by quarters where 1 indicates
                                quarter allocated and 0 not. Note that the length
                                of this array is n_nights_in_semester-1 which is smaller
                                than len(all_dates_dict)
        output_directory (str): path to save the plot

    Returns:
        None
    """
    all_dates_array = list(all_dates_dict.keys())
    sim_results = []
    allocation_statii = []
    results = []
    filename_weather = output_directory + 'Weather_Simulation_Results.csv'
    for a in range(len(all_dates_array)):
        if a < all_dates_dict[current_day]:
            sim_result = 'Past'
            allocation_status = "Past"
            result = "Past"
        elif a == all_dates_dict[current_day]:
            sim_result = '???'
            allocation_status = "True"
            result = "???"
        else:
            # Extra -1 because the days_lost array does not include today
            if days_lost[a - all_dates_dict[current_day] - 1] == 1:
                sim_result = 'Poor'
            else:
                sim_result = 'Clear'
            if np.sum(allocation_remaining[a - all_dates_dict[current_day]]) >= 1:
                allocation_status = 'True'
            else:
                allocation_status = 'False'

            if allocation_status == 'True' and sim_result == 'Poor':
                result = 'Lost'
            elif allocation_status == 'True' and sim_result == 'Clear':
                result = 'Go'
            elif allocation_status == 'False' and sim_result == 'Clear':
                result = 'Miss'
            else:
                result = 'Eh'
        sim_results.append(sim_result)
        allocation_statii.append(allocation_status)
        results.append(result)
    weather_frame = pd.DataFrame({'Date':all_dates_array, 'Allocated':allocation_statii,
                                    'SimWeather':sim_results, 'Designation':results})
    weather_frame.to_csv(filename_weather, index=False)

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

def build_single_target_accessibility (starname, ra, dec, start_date, n_nights_in_semester,
                                      slot_size, compute_turn_on_off=False):
    """
    Compute a target's accessibility map for the entire semester

    Args:
        starname (str): the name of the star (simbad resolvable)
        ra (str): the right ascension in hour angle degrees (ie. 14.33)
        dec (str): the declination in degrees (ie. +75.2)
        start_date (str): the calendar date to begin the calculation from (format is Astropy isot)
        n_nights_in_semester (str): the number of calendar nights remaining in the current semester
        slot_size (int): the size of the slots in minutes
        compute_turn_on_off (boolean): a flag to determine whether or not to compute

    Returns:
        target_accessibility (array): a 1D array of length n_slots_in_semester where 1 indicates
                                     the target is accessible in that slot and 0 otherwise.
    """
    coords = apy.coordinates.SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
    target = apl.FixedTarget(name=starname, coord=coords)

    date_formal = Time(start_date,format='iso',scale='utc')
    date = str(date_formal)[:10]
    target_accessibility = []
    quarter_map = []
    for d in range(n_nights_in_semester):
        tonights_access = is_observable(date, target, slot_size)
        target_accessibility.append(tonights_access)
        date_formal += TimeDelta(1,format='jd')
        date = str(date_formal)[:10]
        if compute_turn_on_off:
            quarter_map.append(quarters_observable(tonights_access))

    if compute_turn_on_off:
        # compute the first and last calendar date that the target is at all accessibile.
        # do so for each quarter of the night independently.
        turns = []
        for q in range(4):
            on_off = compute_on_off_for_quarter(quarter_map, q)
            turns.append(on_off)
        return target_accessibility, turns
    else:
        return target_accessibility

def is_observable(date, target, slot_size):
    """
    Compute a target's accessibility map on a given date, taking into account the
    telescope pointing limits and moon-safe distance at the beginning of every slot.

    Args:
        date (str): the calendar date to compute accessibilty in format 'YYYY-MM-DD'
        target (str): an astroplan FixedTarget object
        slot_size (int): the size of the slots in minutes
    Returns:
        observability_matrix (array): a 1D array of length n_slots_in_night where 1 indicates the
                                      target is accessible in that slot and 0 otherwise.
    """
    # Can't observe too close to zenith
    max_alt = 85.
    # Naysmith deck azimuth direction limits
    min_az = 5.3
    max_az = 146.2
    # Naysmith deck
    min_alt = 33.3
    #non-Naysmith deck minimum elevation
    else_min_alt = 18.0
    # Prefer to observe at least 30 degree altitude if they are not too far north/south
    else_min_alt_alt = 30.0
    max_north = 75.0
    max_south = -35.0

    # This is ~20 min before earliest sunset of the year in Hawaii
    # And ~20 min after the latest sunrise of the year in Hawaii
    # Both are UTC time zone.
    start = date + "T03:30:00"
    daily_start = Time(start)
    end = date + "T17:30:00"
    daily_end = Time(end)
    slot_size = TimeDelta(slot_size*60.,format='sec')
    t = Time(np.arange(daily_start.jd, daily_end.jd, slot_size.jd),format='jd')

    keck = apl.Observer.at_site('W. M. Keck Observatory')
    altaz_coordinates = keck.altaz(t, target, grid_times_targets=True)

    keckapy = apy.coordinates.EarthLocation.of_site('Keck Observatory')
    moon = apy.coordinates.get_moon(Time(t[int(len(t)/2)],format='jd'), keckapy)

    if moon_safe(moon, (target.ra.rad, target.dec.rad)):
        observability_matrix = []
        for i, item in enumerate(altaz_coordinates):
            alt=altaz_coordinates[i].alt.deg
            az=altaz_coordinates[i].az.deg
            deck = np.where((az >= min_az) & (az <= max_az))
            deck_height = np.where((alt <= max_alt) & (alt >= min_alt))
            first = np.intersect1d(deck,deck_height)

            not_deck_1 = np.where((az < min_az))
            not_deck_2 = np.where((az > max_az))

            # for targets sufficiently north or south in declination, allow access map to compute
            # any time they are above telescope pointing limits as OK. For more equitorial targets,
            # require that they be above the preferred minimum elevation.
            if target.dec.deg > max_north or target.dec.deg < max_south:
                not_deck_height = np.where((alt <= max_alt) & (alt >= else_min_alt))
            else:
                not_deck_height = np.where((alt <= max_alt) & (alt >= else_min_alt_alt))

            second = np.intersect1d(not_deck_1,not_deck_height)
            third = np.intersect1d(not_deck_2,not_deck_height)

            good = np.concatenate((first,second,third))

            observability_matrix = np.zeros(len(t),dtype=int)
            observability_matrix[good] = 1
    else:
        observability_matrix = np.zeros(len(t),dtype=int)

    return observability_matrix

def moon_safe(moon, target_tuple):
    """
    Check that a coordinate is not too close to the moon.
    Returns true if target is sufficiently far from the moon to allow for observations.

    Args:
        moon (str): a "moon object" from Astropy
        target_tuple (tuple): the RA and Dec of a target star in hourangle and degrees format
    Returns:
        Unnamed boolean
    """
    ang_dist = apy.coordinates.angular_separation(moon.ra.rad, moon.dec.rad,
                                                    target_tuple[0], target_tuple[1])
    if ang_dist*180/(np.pi) >= 30:
        return True
    else:
        return False

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
                print("We have a problem, error code 1.")
                print(allocation['Date'].iloc[i],
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
                print("We have a problem, error code 2.")
                print(allocation['Date'].iloc[i],
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
                print("We have a problem, error code 3.")
                print(allocation['Date'].iloc[i],
                                allocation['Start'].iloc[i], allocation['Stop'].iloc[i])
        else:
            print("Non-25% of night increment. Implementing whole night as precaution.")
            print("Date: ", allocation['Date'].iloc[i])
            print("True allocation amount: ", item)
            allocation['Start'].iloc[i] = 0
            allocation['Stop'].iloc[i]  = 1

    return allocation

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
        print("We have a problem, error code 4.")

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
