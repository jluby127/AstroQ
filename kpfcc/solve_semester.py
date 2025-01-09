"""
Module for running the semester solver.

This module organizes, builds, and solves the gurobi model. It also writes and organizes the
outputs of the schedule into both machine human readable forms. It is designed to produce info
that is in correct formatting for running the TTP. Further designed to be only run as a function
call from the generateScript.py script.

Example usage:
    import solve_semester as ssm
"""
import sys
import time
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB

DIR_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(DIR_PATH)
# KPF-CC specific files
import helper_functions as hf
import twilight_functions as tw
import reporting_functions as rf
import processing_functions as pf
import mapping_functions as mf
import set_functions as sf

# import line_profiler
#
# @profile
def run_kpfcc(current_day,
               requests_file,
               allocation_file,
               accessibilities_file,
               twilight_file,
               output_directory,
               slot_size,
               run_round_two,
               past_observations_file,
               semester_template_file,
               turn_off_on_file,
               nonqueue_map_file,
               special_map_file,
               zero_out_file,
               run_weather_loss,
               gurobi_output = True,
               plot_results = True,
               solve_time_limit = 300):

    """run_kpfcc
    Args:
        - current_day (str) = the calendar date of the night to produce a script.
                              Sets the "first" day of the semester from which to compute the
                            semester schedule solution from this day forward. Format: YYYY-MM-DD.
        - requests_file (str) = the path and file name to the CSV with all the PI requests.
                                Confirm that column names are correct.
        - allocation_file (str) = the path and file name to the binary map of allocated nights.
        - accessibilities_file (str) = the path and file name to the pickle file containing a
                                       dictionary of target names and associated pre-computed 1D
                                       accessibility maps of length equal to n_slots_in_semester.
        - twilight_file (str) = the path and file name to the CSV with precomputed twilight times.
        - output_directory (str) = the path where all outputs of this function should be saved.
                                    It is recommended that the path be outside the git repo.
        - slot_size (int) = the time, in minutes, for a single slot.
        - run_round_two (boolean) = when True, run the bonus round.
                                    When False, do not run the bonus round.
        - past_observations_file (str) = the path and file name of the CSV containing information
                                         on all previous observations in the semester. If file
                                         does not exist, then we are ignoring prior observations.
        - semester_template_file (str) = the path and file name of the CSV containing the visual
                                         template of the semester. For plotting purposes only.
        - turn_off_on_file (str) = the path and file name of the CSV containing the pre-computed
                                   first and last day of accessiblity for each target.
                                   For plotting purposes only.
        - nonqueue_map_file (str) = the path and file name of the CSV containining a grid of
                                    n_nights_in_semester by n_slots_in_night elements where slots
                                    reserved for non-queue observations are filled with target name.
        - special_map_file (str) = the path and file name of the CSV containining a grid of
                                    n_nights_in_semester by n_slots_in_night elements which contains
                                    information on the custom set of slots a request can be
                                    scheduled into for various reasons of the PI
        - zero_out_file (str) = the path and file name of list of stars that cannot be scheduled
                                tonight for any reason. Often this is empty.
        - run_weather_loss (boolean) = if False, then no nights are lost to weather.
        - gurobi_output (boolean) = a flag to turn off or on the feature of Gurobi printing
                                    to the terminal as it solves the model.
        - plot_results (boolean) = a flag to turn off or on the plotting outputs.
        - solve_time_limit (int) = the maximum time, in seconds, to allow Gurobi to solve the model.
    Returns:
        None
    """
    start_the_clock = time.time()

    # I suggest your output directory be something so that it doesn't autosave
    # to the same directory as the run files and crowds up the GitHub repo.
    # Note to self: move this to the generateScript.py file.
    check = os.path.isdir(output_directory)
    if not check:
        os.makedirs(output_directory)

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set up logistics parameters
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------

    # Only today is important to the semester solver. Any additional dates are only for the TTP.
    current_day = current_day[0]

    # Get semester parameters and define important quantities
    semester_start_date, semester_end_date, semester_length, semester_year, semester_letter = \
        hf.get_semester_info(current_day)
    all_dates_dict = hf.build_date_dictionary(semester_start_date, semester_length)
    all_dates_array = list(all_dates_dict.keys())
    n_nights_in_semester = hf.current_day_tracker(current_day, all_dates_dict)
    print("Total semester length: ", semester_length)
    print("There are " + str(n_nights_in_semester) + " calendar nights remaining in the semester.")

    n_quarters_in_night = 4
    n_hours_in_night = 14
    n_slots_in_quarter = int(((n_hours_in_night*60)/n_quarters_in_night)/slot_size)
    n_slots_in_night = n_slots_in_quarter*n_quarters_in_night
    n_slots_in_semester = n_slots_in_night*n_nights_in_semester

    # Define the slot and night represents the today's date
    today_starting_slot = all_dates_dict[current_day]*n_slots_in_night
    today_starting_night =  all_dates_dict[current_day]
    print("There are " + str(n_slots_in_semester) + " slots remaining in the semester.")


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Read in files and prep targets
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    print("Reading inputs and prepping files.")
    requests_frame = pd.read_csv(requests_file)
    twilight_frame = pd.read_csv(twilight_file, parse_dates=True)

    print("Compiling past observation history.")
    database_info_dict = {}
    # pf.get_kpf_past_database(past_observations_file)
    if os.path.exists(past_observations_file):
        print("Pulled database of past observations this semester.")
        database = pd.read_csv(past_observations_file)
        for i in range(len(requests_frame['Starname'])):
            starmask = database['star_id'] == requests_frame['Starname'][i]
            star_past_obs = database[starmask]
            star_past_obs.sort_values(by='utctime', inplace=True)
            star_past_obs.reset_index(inplace=True)
            total_past_observations = int(len(star_past_obs)/
                                    (requests_frame['# Visits per Night'][i]*
                                    requests_frame['# of Exposures per Visit'][i]))
            star_past_obs, unique_hst_dates_observed, quarter_observed = \
                rf.get_unique_nights(star_past_obs, twilight_frame)
            n_obs_on_date = rf.get_nobs_on_night(star_past_obs, unique_hst_dates_observed)

            if len(unique_hst_dates_observed) > 0:
                most_recent_observation_date = unique_hst_dates_observed[-1]
            else:
                # If request has not been observed, set a dummy value.
                most_recent_observation_date = '0000-00-00'
            # Within the database_info_dict dictionary, each target's data is always in order:
            # element 0 = the calendar date of the most recent observation (HST date)
            # element 1 = a list of calendar dates with at least one observation
            # element 2 = a list of the quarters where the observation took place on the
            #             corresponding the nights in element 1.
            #             If multiple visits in one night, then this is quarter of the first visit.
            # element 3 = a list of the # of observations on each past night,
            #             corresponding the nights in element 1.
            database_info_dict[requests_frame['Starname'][i]] = \
                [most_recent_observation_date, unique_hst_dates_observed, quarter_observed, \
                n_obs_on_date]
    else:
        print("No past observation history to parse.")

    print("Determining slots needed for exposures.")
    # schedule multi-shots and multi-visits as if a single, long exposure.
    # When n_shots and n_visits are both 1, this reduces down to just the stated exposure time.
    slots_needed_for_exposure_dict = {}
    for n,row in requests_frame.iterrows():
        name = row['Starname']
        singlevisit = row['Nominal Exposure Time [s]']*row['# of Exposures per Visit'] + \
            45*(row['# Visits per Night'] - 1)
        slots_needed_for_exposure_dict[name] = \
            hf.slots_required_for_exposure(singlevisit*row['# Visits per Night'], slot_size)

    print("Determine available slots in each night.")
    # available_slots_in_each_night is a 1D matrix of length nights n
    # This is not a gorubi variable, but a regular python variable
    # Each element will hold an integer which represents the number of slots are available in each
    # quarter of a given night, after accounting for non-observable times due to day/twilight.
    available_slots_in_each_night = []
    for date in all_dates_array:
        slots_tonight = tw.determine_twilight_edge(date, twilight_frame, slot_size)
        available_slots_in_each_night.append(slots_tonight)
    twilight_map_all = np.array(mf.build_twilight_map(available_slots_in_each_night,
                                n_slots_in_night, invert=False))
    twilight_map_remaining = twilight_map_all[all_dates_dict[current_day]:]
    twilight_map_remaining_flat = twilight_map_remaining.copy().flatten()
    twilight_map_remaining_2D = np.reshape(twilight_map_remaining,
                                    (n_nights_in_semester, n_slots_in_night))

    print("Reading pre-comupted accessibility maps.")
    rewrite_flag = False
    default_access_maps = mf.read_accessibilty_map_dict(accessibilities_file)
    for n,row in requests_frame.iterrows():
        name = row['Starname']
        # check that this target has a pre-computed accessibility map,
        # if not, make one and add it to the file
        try:
            try_read = default_access_maps[name]
        except:
            print(name + " not found in precomputed accessibilty maps. Running now.")
            # Note: the -1 is to account for python indexing
            new_written_access_map = mf.build_single_target_accessibility(name, row['RA'],
                                               row['Dec'], semester_start_date, semester_length-1,
                                                slot_size)
            default_access_maps[name] = np.array(new_written_access_map).flatten()
            rewrite_flag = True
    if rewrite_flag:
        # overwrite with the updated file
        mf.write_accessibilty_map_dict(default_access_maps, accessibilities_file)

    # Read in the customize acccessibility maps for unique targets, if exists.
    if os.path.exists(special_map_file):
        custom_access_maps = mf.read_accessibilty_map_dict(special_map_file)
    else:
        custom_access_maps = {}

    # Read in the list of targets to "zero out", i.e. not allowed to be scheduled only tonight.
    if os.path.exists(zero_out_file):
        zero_out = pd.read_csv(zero_out_file)
        zero_out_names = list(zero_out['Target'])
    else:
        zero_out_names = []

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
    historical_weather_data = pd.read_csv(DIR_PATH[:-5] + "data/maunakea_weather_loss_data.csv")
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

    print("Incorporating non-queue observations.")
    # Exclude slots that must be assigned to time-sensative observations
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

    # Define the tuples of request and available slot for each request.
    # This becomes the grid over which the Gurobi variables are defined.
    # Now, slots that were never possible for scheduling are not included in the model.
    Aset = []
    Aframe_keys = []
    for n,row in requests_frame.iterrows():
        name = row['Starname']
        n_visits = int(row['# Visits per Night'])
        intra = int(row['Minimum Intra-Night Cadence'])
        inter = int(row['Minimum Inter-Night Cadence'])
        slots_needed = slots_needed_for_exposure_dict[name]
        for d in range(len(available_indices_for_request[name])):
            for s in available_indices_for_request[name][d]:
                Aset.append((name, d, s))
                Aframe_keys.append([name, d, s, slots_needed, n_visits, intra, inter])

    Aframe = pd.DataFrame(Aframe_keys, columns =['r', 'd', 's', 'e', 'v', 'tra', 'ter'])
    schedulable_requests = list(Aframe['r'].unique())
    for name in list(requests_frame['Starname']):
        if name not in schedulable_requests:
            print("WARNING: Target " + name + " has no valid day/slot pairs and therefore is effectively removed from the model.")
    # duplicate columns for easy indexing later
    Aframe['rr'] = Aframe['r']
    Aframe['dd'] = Aframe['d']
    Aframe['ss'] = Aframe['s']

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set Gorubi model variables
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    print("Building variables")
    m = gp.Model('Semester_Scheduler')

    # Yrs is technically a 1D matrix indexed by tuples.
    # But in practice best think of it as a 2D square matrix of requests r and slots s, with gaps.
    # Slot s for request r will be 1 to indicate starting an exposure for that request in that slot
    Yrds = m.addVars(Aset, vtype = GRB.BINARY, name = 'Requests_Slots')

    # theta is the "shortfall" variable, continous in natural numbers.
    # theta = m.addVars(requests_frame['Starname'], name='Shortfall')
    theta = m.addVars(schedulable_requests, name='Shortfall')

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set Gorubi constraints
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------

    print("Constraint 0: Build theta variable")
    max_bonus_observations_pct = 0.5
    max_bonus_observations = 10
    aframe_slots_for_request = Aframe.groupby(['r'])[['d', 's']].agg(list)
    for name in schedulable_requests:
        idx = requests_frame.index[requests_frame['Starname']==name][0]

        if database_info_dict == {}:
            past_nights_observed = 0
        else:
            past_nights_observed = len(database_info_dict[name][1])

        # Safety valve for if the target is over-observed for any reason
        # Example: a cadence target that is also an RM target will have many more past
        # observations than is requested.
        # When we move over to parsing the Keck database for a unique request ID, instead of
        # parsing the Jump database on non-unique star name, this can be removed.
        if past_nights_observed > requests_frame['# of Nights Per Semester'][idx] + \
                    int(requests_frame['# of Nights Per Semester'][idx]*max_bonus_observations_pct):
            true_max_obs = past_nights_observed
        else:
            true_max_obs = (requests_frame['# of Nights Per Semester'][idx] - past_nights_observed)\
                   + int(requests_frame['# of Nights Per Semester'][idx]*max_bonus_observations_pct)

        m.addConstr(theta[name] >= 0, 'greater_than_zero_shortfall_' + str(name))
        # Get all (d,s) pairs for which this request is valid.
        available = list(zip(list(aframe_slots_for_request.loc[name].d), \
                                                        list(aframe_slots_for_request.loc[name].s)))
        m.addConstr(theta[name] >= ((requests_frame['# of Nights Per Semester'][idx] - \
                    past_nights_observed) - gp.quicksum(Yrds[name, d, s] for d,s in available)), \
                    'greater_than_nobs_shortfall_' + str(name))
        m.addConstr(gp.quicksum(Yrds[name, d, s] for d,s in available) <=
                    true_max_obs, 'max_unique_nights_for_request_' + str(name))

    print("Constraint 1: Enforce one request per slot.")
    # Get all requests which are valid in slot (d, s)
    Aframe_ds = pd.merge(
        Aframe.drop_duplicates(['d', 's']),
        Aframe[['r', 'd', 's']],
        suffixes=['', '2'],
        on=['d', 's'])
    ds_requests = Aframe_ds.groupby(['d','s'])[['r2']].agg(list)
    Aset_ds_no_duplicates = Aframe.copy()
    Aset_ds_no_duplicates = Aframe.drop_duplicates(subset=['d', 's'])
    # Construct the constraint
    for i, row in Aset_ds_no_duplicates.iterrows():
        requests_valid_in_slot = list(ds_requests.loc[(row.d, row.s)])[0]
        m.addConstr((gp.quicksum(Yrds[name,row.d,row.s] for name in requests_valid_in_slot) <= 1),
                        'one_request_per_slot_' + str(row.d) + "d_" + str(row.s) + "s")

    print("Constraint 2: Reserve slots for for multi-slot exposures.")
    # Get all requests that are  valid in (d,s+e) pair for a given (d,s,1..e)
    requests_valid_in_reserved_slots = pd.merge(Aframe.query('e > 1 ')['r d s e'.split()] \
        ,Aframe['r d s'.split()],on=['d'],suffixes=['','2']) \
        .query('s < s2 < s + e').groupby('r d s'.split()).agg(list)
    # If request requires only 1 slot to complete, then no constraint on reserving additional slots
    Aframe_multislots = Aframe[Aframe.e > 1]
    for i, row in Aframe_multislots.iterrows():
        # construct list of (r,d,s) indices to be constrained. These are all requests that are
        # valid in slots (d, s+1) through (d, s + e)
        # Ensuring slot (d, s) is not double filled already taken care of in Constraint 1.
        allr = list(requests_valid_in_reserved_slots.loc[row.r, row.d, row.s]['r2'])
        alls = list(requests_valid_in_reserved_slots.loc[row.r, row.d, row.s]['s2'])
        all_reserved_slots = list(zip(allr, [row.d]*len(allr), alls))
        m.addConstr((row.e*(1 - Yrds[row.r,row.d,row.s])) >= gp.quicksum(Yrds[c]
                                                        for c in all_reserved_slots),
                        'reserve_multislot_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")

    print("Constraint 3: Schedule request's maximum observations per night.")
    # Get all valid slots s for request r on day d
    aframe_slots_on_day = pd.merge(
        Aframe.drop_duplicates(['r','d',]),
        Aframe[['r','d','s']],
        suffixes=['','2'],on=['r']
    ).query('d == d2')
    slots_on_day = aframe_slots_on_day.groupby(['r','d'])[['s2']].agg(list)
    unique_request_day_pairs = Aframe.drop_duplicates(['r','d'])
    # Build the constraint
    for i, row in unique_request_day_pairs.iterrows():
        constrained_slots_tonight = np.array(slots_on_day.loc[(row.r, row.d)][0])
        m.addConstr((gp.quicksum(Yrds[row.r,row.d,ss] for ss in constrained_slots_tonight) <= row.v),
                'max_observations_per_night_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")

    print("Constraint 4: Enforce inter-night cadence.")
    aframe_intercadence = pd.merge(
        Aframe.drop_duplicates(['r','d',]),
        Aframe[['r','d','s']], #
        suffixes=['','2'],on=['r']
    ).query('d + 0 < d2 < d + ter') # When +1, it excludes the first day of the semester always
    intercadence = aframe_intercadence.groupby(['r','d'])[['d2','s2']].agg(list)
    # When inter-night cadence is 1, there will be no keys to constrain so skip
    # While the if/else statement would catch these, by shrinking the list here we do fewer
    # total steps in the loop.
    mask_inter_1 = Aframe['ter'] > 1
    # mask_inter_2 = # ensure that even if request r has at least one slot on day d,
    #                   if there are no slots for request r on days d + delta for all delta in
    #                   internight cadence, then don't include day d
    Aset_inter = Aframe[mask_inter_1]
    # We don't want duplicate slots on day d because we only need this constraint once per day
    # With duplicates, the same constraint would be applied to (r, d, s) and (r, d, s+1) which
    # is superfluous since we are summing over tonight's slots
    Aset_inter_no_duplicates = Aset_inter.copy()
    Aset_inter_no_duplicates = Aset_inter_no_duplicates.drop_duplicates(subset=['r', 'd'])
    for i, row in Aset_inter_no_duplicates.iterrows():
        constrained_slots_tonight = np.array(slots_on_day.loc[(row.r, row.d)][0])
        # Get all slots for pair (r, d) where valid
        if (row.r, row.d) in intercadence.index:
            slots_to_constrain_future = intercadence.loc[(row.r, row.d)]
            ds_pairs = list(zip(slots_to_constrain_future.d2, slots_to_constrain_future.s2))
            m.addConstr((gp.quicksum(Yrds[row.r,row.d,s2] for s2 in constrained_slots_tonight)/row.v \
                 <= (1 - (gp.quicksum(Yrds[row.r,d3,s3] for d3, s3 in ds_pairs)))), \
                'enforce_internight_cadence_' + row.r + "_" + str(row.d) + "d_" + str(row.s) + "s")
        else:
            # For request r, there are no (d,s) pairs that are within "inter cadence" days of the
            # given day d, therefore nothing to constrain. If I can find a way to filter out these
            # rows as a "mask_inter_2", then the if/else won't be needed
            continue

    complete_constraints_build = time.time()
    print("Total Time to build constraints: ",
        np.round(complete_constraints_build-start_the_clock,3))

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set Gorubi Objective and Run
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    print("Begin model solve.")
    m.setObjective(gp.quicksum(theta[name] for name in schedulable_requests), GRB.MINIMIZE)

    m.params.TimeLimit = solve_time_limit
    m.Params.OutputFlag = gurobi_output
    # Allow stop at 5% gap to prevent from spending lots of time on marginally better solution
    m.params.MIPGap = 0.05
    # More aggressive presolve gives better solution in shorter time
    m.params.Presolve = 2
    m.update()
    m.optimize()

    if m.Status == GRB.INFEASIBLE:
        print('Model remains infeasible. Searching for invalid constraints')
        search = m.computeIIS()
        print("Printing bad constraints:")
        for c in m.getConstrs():
            if c.IISConstr:
                print('%s' % c.ConstrName)
        for c in m.getGenConstrs():
            if c.IISGenConstr:
                print('%s' % c.GenConstrName)
    else:
        print("Round 1 Model Solved.")

    complete_round1_model = time.time()
    print("Total Time to finish solver: ", np.round(complete_round1_model-start_the_clock,3))
    print()
    print()
    print()

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Retrieve data from solution
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    file = open(output_directory + "runReport.txt", "w")
    file.close()

    print("Building human readable schedule.")
    combined_semester_schedule_available = hf.write_available_human_readable(
                                all_dates_dict, current_day, semester_length,
                                n_nights_in_semester, n_slots_in_night, twilight_map_remaining_2D,
                                allocation_map_2D, weathered_map, nonqueue_map_file)
    np.savetxt(output_directory + 'raw_combined_semester_schedule_available.txt',
        combined_semester_schedule_available, delimiter=',', fmt="%s")

    combined_semester_schedule_stars = hf.write_stars_schedule_human_readable( \
            combined_semester_schedule_available, Yrds, list(requests_frame['Starname']),
            semester_length, n_slots_in_night, n_nights_in_semester,
            all_dates_dict, slots_needed_for_exposure_dict, current_day)
    np.savetxt(output_directory + 'raw_combined_semester_schedule_Round1.txt',
        combined_semester_schedule_stars, delimiter=',', fmt="%s")

    round = 'Round 1'
    rf.build_fullness_report(combined_semester_schedule_stars, allocation_map_2D, requests_frame,
                                slot_size, round, output_directory)

    print("Writing Report.")
    filename = open(output_directory + "runReport.txt", "a")
    theta_n_var = []
    counter = 0
    for v in theta.values():
        varname = v.VarName
        varval = v.X
        counter += varval
    print("Sum of Theta: " + str(counter))
    filename.write("Sum of Theta: " + str(counter) + "\n")

    plot_results = False
    if plot_results:
        print("Writing cadence plot files.")
        turn_on_off_frame = pd.read_csv(turn_off_on_file)
        all_starmaps = {}
        for i in range(len(requests_frame)):
            if database_info_dict != {}:
                starmap = rf.build_observed_map_past( \
                    database_info_dict[requests_frame['Starname'][i]], semester_template_file)
            else:
                starmap = rf.build_observed_map_past([[],[],[],[]], semester_template_file)

            starmap_updated = rf.build_observed_map_future(combined_semester_schedule_stars,
                                requests_frame['Starname'][i], starmap,
                                slots_needed_for_exposure_dict,
                                np.array(allocation_all).flatten(),
                                np.array(weather_diff_remaining).flatten(),
                                all_dates_dict[current_day])

            all_starmaps[requests_frame['Starname'][i]] = starmap_updated
            future_unique_days_forecasted = 0
            for k, item in enumerate(combined_semester_schedule_stars):
                if requests_frame['Starname'][i] in combined_semester_schedule_stars[k]:
                    future_unique_days_forecasted += 1

            try:
                past_unique_dates_for_star = database_info_dict[requests_frame['Starname'][i]][1]
            except:
                past_unique_dates_for_star = []
            rf.write_cadence_plot_file(requests_frame['Starname'][i], starmap_updated,
                                        turn_on_off_frame, requests_frame,
                                        future_unique_days_forecasted,
                                        past_unique_dates_for_star,
                                        current_day, output_directory)

    complete_round1_plots = time.time()
    print("Total Time to complete Round 1: " + \
        str(np.round(complete_round1_plots-start_the_clock,3)))
    filename.write("Total Time to complete Round 1: " + \
        str(np.round(complete_round1_plots-start_the_clock,3)) + "\n")
    filename.write("\n")
    filename.write("\n")
    filename.close()
    print("Round 1 complete.")


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Initiate Round 2 Scheduling
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------

    if run_round_two:
        print("Beginning Round 2 Scheduling.")
        first_stage_objval = m.objval
        epsilon = 5
        m.params.TimeLimit = solve_time_limit
        m.Params.OutputFlag = gurobi_output
        m.params.MIPGap = 0.05
        m.addConstr(gp.quicksum(theta[name] for name in requests_frame['Starname']) <= \
                    first_stage_objval + epsilon)
        m.setObjective(gp.quicksum(slots_needed_for_exposure_dict[r]*Yrs[r,d,s]
                        for r, d, s in Aset),
                        GRB.MAXIMIZE)
        m.update()
        m.optimize()

        if m.Status == GRB.INFEASIBLE:
            print('Model remains infeasible. Searching for invalid constraints')
            search = m.computeIIS()
            print("Printing bad constraints:")
            for c in m.getConstrs():
                if c.IISConstr:
                    print('%s' % c.ConstrName)
            for c in m.getGenConstrs():
                if c.IISGenConstr:
                    print('%s' % c.GenConstrName)
        else:
            print("")
            print("")
            print("Round 2 Model Solved.")

        complete_round2_model = time.time()
        print("Total Time to finish round 2 solver: ",
            np.round(complete_round2_model-start_the_clock,3))

        combined_semester_schedule_stars = hf.write_stars_schedule_human_readable(
                combined_semester_schedule_available, Yrs, requests_frame['Starname'],
                semester_length, n_slots_in_night, n_nights_in_semester,
                all_dates_dict, slots_needed_for_exposure_dict, current_day)
        np.savetxt(output_directory + 'raw_combined_semester_schedule_Round2.txt',
            combined_semester_schedule_stars, delimiter=',', fmt="%s")

        round = 'Round 2'
        rf.build_fullness_report(combined_semester_schedule_stars, allocation_map_2D,
                                    requests_frame, slot_size, round, output_directory)

        scheduleR1 = np.loadtxt(output_directory + 'raw_combined_semester_schedule_Round1.txt',
            delimiter=',', dtype=str)
        scheduleR2 = np.loadtxt(output_directory + 'raw_combined_semester_schedule_Round2.txt',
            delimiter=',', dtype=str)
        R2_requests = rf.get_gap_filler_targets(scheduleR1, scheduleR2, all_dates_dict[current_day])
        np.savetxt(output_directory + 'Round2_Requests.txt', R2_requests, delimiter=',', fmt="%s")

        filename = open(output_directory + "runReport.txt", "a")
        theta_n_var = []
        counter = 0
        for v in theta.values():
            varname = v.VarName
            varval = v.X
            counter += varval
        filename.write("Sum of Theta: " + str(counter) + "\n")

        complete_round2_plots = time.time()
        filename.write("Total Time to complete Round 2: " + \
            str(np.round(complete_round2_plots-start_the_clock,3)) + "\n")
        filename.write("\n")
        filename.write("\n")
        filename.close()
        print("Round 2 complete.")
    else:
        print("Not running Round 2. Duplicating Raw Schedule as dummy file.")
        np.savetxt(output_directory + 'raw_combined_semester_schedule_Round2.txt', \
            combined_semester_schedule_stars, delimiter=',', fmt="%s")
        R2_requests = []
        np.savetxt(output_directory + 'Round2_Requests.txt', R2_requests, delimiter=',', fmt="%s")

    print("The optimal semester schedule is found, clear skies!")
