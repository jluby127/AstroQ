"""
Module for building reports after the solution is found. Designed to be only run as a function
call from the generateScript.py script.

Example usage:
    import reporting_functions as rf
"""
import os
import math
import time
from astropy.time import Time
from astropy.time import TimeDelta
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy import units as u

import numpy as np
import pandas as pd
import matplotlib.pyplot as pt
import re

import kpfcc.history as hs

def build_fullness_report(combined_semester_schedule, manager, round_info):
    """
    Determine how full the schedule is: slots available, slots scheduled, and slots required

    Args:
        combined_semester_schedule (array): a 2D array of dimensions n_nights_in_semester by
                                            n_slots_in_night where elements denote how the slot is
                                            used: target, twilight, weather, not scheduled.
        allocation_map_2D (array): a 2D array where rows represent a night and columns represent
                                   the quarter within that night. Values are 1 if that
                                   night/quarter is allocated and 0 if not.
        manager (obj): a data_admin object

    Returns:
        None
    """
    file_path = manager.output_directory + "runReport.txt"
    print(f"Writing to: {file_path}")
    with open(manager.output_directory + "runReport.txt", "a") as file:
        file.write("Stats for " + str(round_info) + "\n")
        file.write("------------------------------------------------------" + "\n")
        listnames = list(manager.requests_frame['starname'])
        unavailable = 0
        unused = 0
        used = 0
        for b, item1 in enumerate(combined_semester_schedule):
            for c, item2 in enumerate(combined_semester_schedule[b]):
                if "X" in combined_semester_schedule[b][c] or "*" in combined_semester_schedule[b][c] \
                        or "W" in combined_semester_schedule[b][c]:
                    unavailable += 1
                if combined_semester_schedule[b][c] == "":
                    unused += 1
                if combined_semester_schedule[b][c] in listnames or "RM___" in combined_semester_schedule[b][c]:
                    used += 1
        available = unused + used
        allocated = np.sum(manager.allocation_map_2D.flatten())
        file.write("N slots in semester:" + str(np.prod(combined_semester_schedule.shape)) + "\n")
        file.write("N available slots:" + str(allocated) + "\n")
        file.write("N slots scheduled: " + str(used) + "\n")
        file.write("N slots left empty: " + str(allocated-used) + "\n")

        total_slots_requested = 0
        for i in range(len(manager.requests_frame)):
            total_slots_requested += manager.requests_frame['n_inter_max'][i]* \
                math.ceil(manager.requests_frame['exptime'][i]/(manager.slot_size*60.))
        file.write("N slots requested (total): " + str(total_slots_requested) + "\n")
        percentage = np.round((used*100)/allocated,3)
        file.write("Percent full: " + str(percentage) + "%." + "\n")
        file.close()

def report_allocation_stats(manager):
    """
    Return the number of Q1, Q2, Q3, and Q4 nights, as well as single, half, three quarter
    and full nights. Further produce a plot showing the allocation.

    Args:
        manager (obj): a data_admin object
        allocation_map_2D (array): a 2D array where rows represent a night and columns represent
                                   the quarter within that night. Values are 1 if that
                                   night/quarter is allocated and 0 if not.

    Returns:
        None
    """
    current_day_number = manager.all_dates_dict[manager.current_day]
    ff = open(manager.output_directory + "allocation_stats.txt", "w")

    fig = pt.figure(figsize=(10,6))

    q1s = 0
    q2s = 0
    q3s = 0
    q4s = 0
    count0 = 0
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0
    for j, item in enumerate(manager.allocation_map_2D_NQ):
        holder = [0, 0, 0, 0]
        if manager.allocation_map_2D_NQ[j][0] == 1:
            q1s += 1
            # Note: I recently had to add the "-1" to each of these lines. Not sure why, but it broke when doing optimal allocation.
            date_info = manager.all_dates_array[current_day_number + j - 1] + " - q0"
            #ff.write(date_info + "\n")
            holder[0] = 1
            pt.axvline(current_day_number + j, ymin=0, ymax=0.25, color='b')
        if manager.allocation_map_2D_NQ[j][1] == 1:
            q2s += 1
            date_info = manager.all_dates_array[current_day_number + j - 1] + " - q1"
            #ff.write(date_info + "\n")
            holder[1] = 1
            pt.axvline(current_day_number + j, ymin=0.25, ymax=0.5, color='b')
        if manager.allocation_map_2D_NQ[j][2] == 1:
            q3s += 1
            date_info = manager.all_dates_array[current_day_number + j - 1] + " - q2"
            #ff.write(date_info + "\n")
            holder[2] = 1
            pt.axvline(current_day_number + j, ymin=0.5, ymax=0.75, color='b')
        if manager.allocation_map_2D_NQ[j][3] == 1:
            q4s += 1
            date_info = manager.all_dates_array[current_day_number + j - 1] + " - q3"
            #ff.write(date_info + "\n")
            holder[3] = 1
            pt.axvline(current_day_number + j, ymin=0.75, ymax=1.0, color='b')
        allocated_quarters = np.sum(manager.allocation_map_2D_NQ[j])
        if allocated_quarters == 0:
            count0 += 1
        if allocated_quarters == 1:
            count1 += 1
        if allocated_quarters == 2:
            count2 += 1
        if allocated_quarters == 3:
            count3 += 1
        if allocated_quarters == 4:
            count4 += 1

    size=20
    pt.ylim(0,1)
    pt.ylabel('Quarter of Night', fontsize=size)
    pt.xlabel('Day in Semester', fontsize=size)
    pt.tick_params(axis='both', labelsize=size)
    pt.savefig(manager.output_directory + "allocation_visualization.png", dpi=200, bbox_inches='tight', facecolor='w')

    ff.write("\n")
    ff.write("There are " + str(q1s) + " first quarters." + "\n")
    ff.write("There are " + str(q2s) + " second quarters." + "\n")
    ff.write("There are " + str(q3s) + " third quarters." + "\n")
    ff.write("There are " + str(q4s) + " fourth quarters." + "\n")
    ff.write("\n")
    ff.write("There are " + str(count0) + " no quarter nights." + "\n")
    ff.write("There are " + str(count1) + " 1 quarter nights."+ "\n")
    ff.write("There are " + str(count2) + " 2 quarter nights."+ "\n")
    ff.write("There are " + str(count3) + " 3 quarter nights."+ "\n")
    ff.write("There are " + str(count4) + " 4 quarter nights."+ "\n")
    ff.write("\n")
    total_quarters = count1 + 2*count2 + 3*count3 + 4*count4
    total_nights = count1 + count2 + count3 + count4
    ff.write("Total quarters allocated: " + str(total_quarters) + "\n")
    ff.write("Total unique nights allocated: " + str(total_nights) + "\n")
    ff.close()

def write_out_results(manager, theta, round, start_the_clock):
    """
    Write the Run Report
    Args:
        manager (obj): a data_admin object
        theta (obj): a Scheduler.theta object
        rount (str): "Round1" or "Round2"
        start_the_clock (obj): time object of when the scheduler code began

    Returns:
        None
    """
    filename = open(manager.output_directory + "runReport.txt", "a")
    theta_n_var = []
    counter = 0
    for v in theta.values():
        varname = v.VarName
        varval = v.X
        counter += varval
    print("Sum of Theta: " + str(counter))
    filename.write("Sum of Theta: " + str(counter) + "\n")
    print("Total Time to complete " + round + ": " + str(np.round(time.time()-start_the_clock,3)))
    filename.write("Total Time to complete " + round +  ": " + str(np.round(time.time()-start_the_clock,3)) + "\n")
    filename.close()

def build_observed_map_past(past_info, starmap_template_filename):
    """
    Construct stage one (past) of the starmap which is then used to build the cadence plot.

    Args:
        past_info (array): 2D array following schema of the values for the database_info_dict
        starmap_template_filename (str): the path and filename of the starmap template

    Returns:
        starmap (dataframe): an updated starmap which includes the past history of the target
    """
    starmap = pd.read_csv(starmap_template_filename)
    observed = [False]*len(starmap)
    n_observed = [0]*len(starmap)
    for i, item in enumerate(past_info[1]):
        ind = list(starmap['Date']).index(past_info[1][i])
        observed[ind + int(past_info[2][i]-0.5)] = True
        n_observed[ind + int(past_info[2][i]-0.5)] = past_info[3][i]
    starmap['Observed'] = observed
    starmap['N_obs'] = n_observed
    return starmap

def build_observed_map_future(manager, combined_semester_schedule, starname, starmap):
    """
    Construct stage two (future) of the starmap which is then used to build the cadence plot.

    Args:
        manager (obj): a data_admin object
        combined_semester_schedule (array): a 2D array of dimensions n_nights_in_semester by
                                            n_slots_in_night where elements denote how the slot is
                                            used: target, twilight, weather, not scheduled.
        starname (str): the request in question
        starmap (dataframe): an updated starmap which includes the past history of the target

    Returns:
        starmap (dataframe): an updated starmap which includes the past history of the target
    """
    starmap['Allocated'] = [False]*len(starmap)
    starmap['Weathered'] = [False]*len(starmap)

    for i, item in enumerate(combined_semester_schedule):
        if starname in combined_semester_schedule[i]:
            ind = list(combined_semester_schedule[i]).index(starname)
            quarter = determine_quarter(ind, np.shape(combined_semester_schedule)[1])
            starmap['Observed'][i*4 + quarter] = True
            starmap['N_obs'][i*4 + quarter] = list(combined_semester_schedule[i]).count(
                    starname)/manager.slots_needed_for_exposure_dict[starname]

    for a, item in enumerate(manager.allocation_all):
        # mulitply by 4 because each element of allocation_all represents one quarter
        # starmap['Allocated'][days_into_semester*4 + a] = bool(allocation_all[a])
        starmap['Allocated'][a] = bool(manager.allocation_all[a])
    for w, item in enumerate(manager.weather_diff_remaining.flatten()):
        starmap['Weathered'][manager.all_dates_dict[manager.current_day]*4+ w] = bool(manager.weather_diff_remaining.flatten()[w])

    return starmap

def determine_quarter(value, n_slots_in_night):
    """
    Given a slot number, determine what quarter it is to be observed in.

    Args:
        value (int): the slot number to be observed
        n_slots_in_night (int): the number of slots in the night

    Returns:
        quarter (int): the corresponding quarter (1, 2, 3, 4) of the night
    """
    if value <= int(n_slots_in_night*(1/4.)):
        quarter = 0
    elif value <= int(n_slots_in_night*(2/4.)) and value > int(n_slots_in_night*(1/4.)):
        quarter = 1
    elif value <= int(n_slots_in_night*(3/4.)) and value > int(n_slots_in_night*(2/4.)):
        quarter = 2
    elif value <= n_slots_in_night and value > int(n_slots_in_night*(3/4.)):
        quarter = 3
    elif value <= int(n_slots_in_night*(5/4.)) and value > n_slots_in_night:
        quarter = 3
    else:
        quarter = 0
        print("Houston, we've had a problem. No valid quarter: ", value, n_slots_in_night)
    return quarter

def write_stars_schedule_human_readable(combined_semester_schedule, Yrds, manager, round_info):
    """
    Turns the non-square matrix of the solution into a square matrix and starts the human readable
    solution by filling in the slots where a star's exposre is started.

    Args:
        combined_semester_schedule (array): the human readable solution
        Yns (array): the Gurobi solution with keys of (starname, slot_number) and values 1 or 0.
        manager (obj): a data_admin object
        round_info (str): the solution round, ie "Round1"

    Returns:
        combined_semester_schedule (array): the updated human readable solution
    """
    end_past = manager.all_dates_dict[manager.current_day]*manager.n_slots_in_night
    all_star_schedules = {}
    for name in list(manager.requests_frame['starname']):
        star_schedule = []
        # buffer the past with zeros
        for p in range(end_past):
            star_schedule.append(0)
        for d in range(manager.n_nights_in_semester):
            for s in range(manager.n_slots_in_night):
                try:
                    value = int(np.round(Yrds[name, d, s].x))
                except KeyError:
                    value = 0.0
                except:
                    print("Error: helper_functions.py line 224: ", name, d, s)
                star_schedule.append(value)
        all_star_schedules[name] = star_schedule

    combined_semester_schedule = combined_semester_schedule.flatten()
    for s in range(manager.n_slots_in_semester):
        slotallocated = ''
        for name in list(manager.requests_frame['starname']):
            if all_star_schedules[name][manager.today_starting_slot+s] == 1:
                slotallocated += str(name)
                # print(s, combined_semester_schedule[s], "update: " + str(combined_semester_schedule[manager.today_starting_slot+s]), slotallocated)
        combined_semester_schedule[manager.today_starting_slot+s] += str(slotallocated)
    combined_semester_schedule = np.reshape(combined_semester_schedule,
            (manager.semester_length, manager.n_slots_in_night))

    # The semester solver puts a 1 only in the slot that starts the exposure for a target.
    # Therefore, many slots are empty because they are part of a multi-slot visit.
    # Here fill in the multi-slot exposures appropriately for ease of human reading and accounting.
    for n in range(manager.semester_length -1-manager.all_dates_dict[manager.current_day], -1, -1):
        for s in range(manager.n_slots_in_night-1, -1, -1):
            if combined_semester_schedule[n+manager.all_dates_dict[manager.current_day]][s] in list(manager.requests_frame['starname']):
                target_name = combined_semester_schedule[n+manager.all_dates_dict[manager.current_day]][s]
                slots_needed_for_exposure = manager.slots_needed_for_exposure_dict[target_name]
                if slots_needed_for_exposure > 1:
                    for e in range(1, slots_needed_for_exposure):
                        combined_semester_schedule[n+manager.all_dates_dict[manager.current_day]][s+e] += \
                                target_name
    for m in range(len(combined_semester_schedule)):
        # convert the holder string to meaningful string
        if combined_semester_schedule[m][1] == 'supercalifragilisticexpialidocious':
            for l in range(len(combined_semester_schedule[m])):
                combined_semester_schedule[m][l] = 'Past'

    np.savetxt(manager.output_directory + 'raw_combined_semester_schedule_' + round_info + '.txt',
        combined_semester_schedule, delimiter=',', fmt="%s")
    return combined_semester_schedule

def write_available_human_readable(manager):
    """
    Fill in the human readable solution with the non-observation information: non-allocated slots,
    weather loss slots, non-queue slots, twilight slots.

    Args:
        manager (obj): a data_admin object
        twilight_map (array): the 1D array of length n_slots_in_semester where 1's represent slots
                            not in night time and 0's represent slots that are during day/twilight
        allocation_map_2D (array): a 2D array where rows represent a night and columns represent
                                   the quarter within that night. Values are 1 if that
                                   night/quarter is allocated and 0 if not.
        weathered_map (array): a 1D array of length s slots in semester where elements are 1 if
                                that slot has been modeled as lost to weather and 0 if not

    Returns:
        combined_semester_schedule (array): a 2D array of dimensions n_nights_in_semester by
                                            n_slots_in_night where elements denote how the slot is
                                            used: target, twilight, weather, not scheduled.
    """
    if os.path.exists(manager.nonqueue_map_file):
        nonqueuemap_slots_strs = np.loadtxt(manager.nonqueue_map_file, delimiter=',', dtype=str)

    # The past does not matter to us here, so specify the days/slots that are to be ignored.
    end_past = manager.all_dates_dict[manager.current_day]*manager.n_slots_in_night
    combined_semester_schedule = ['']*manager.semester_length*manager.n_slots_in_night
    combined_semester_schedule[0] = 'longwordhereformakingspace'
    for c in range(end_past):
        # for some reason when adding strings within an array, the max length of new string is the
        # length of the longest string in the whole array. So choosing an arbitrary long word
        # as a placeholder. Later I post-process this out.
        # For some reason this word had to be different than the one above...
        combined_semester_schedule[c] += 'supercalifragilisticexpialidocious'
    combined_semester_schedule = np.reshape(combined_semester_schedule,
            (manager.semester_length, manager.n_slots_in_night))

    for n in range(manager.semester_length - manager.all_dates_dict[manager.current_day]):
        for s in range(manager.n_slots_in_night):
            slotallocated = ''
            # remember that twilight map is "inverted": the 1's are time where it is night and the
            # 0's are time where it is day/twilight.
            if manager.twilight_map_remaining_2D[n][s] == 0:
                slotallocated += '*'
            if manager.allocation_map_2D[n][s] == 0:
                slotallocated += 'X'
            if manager.weathered_map[n][s] == 1:# and slotallocated == '':
                slotallocated += 'W'
            if os.path.exists(manager.nonqueue_map_file):
                slotallocated += str(nonqueuemap_slots_strs[n + manager.all_dates_dict[manager.current_day], ][s])
            combined_semester_schedule[n + manager.all_dates_dict[manager.current_day], ][s] += str(slotallocated)

    np.savetxt(manager.output_directory + 'raw_combined_semester_schedule_available.txt',
        combined_semester_schedule, delimiter=',', fmt="%s")
    return combined_semester_schedule

def serialize_schedule(Yrds, manager):
    """
    Turns the non-square matrix of the solution into a square matrix and starts the human readable
    solution by filling in the slots where a star's exposre is started.

    Args:
        combined_semester_schedule (array): the human readable solution
        Yns (array): the Gurobi solution with keys of (starname, slot_number) and values 1 or 0.
        manager (obj): a data_admin object

    Returns:
        None
    """
    df = pd.DataFrame(Yrds.keys(),columns=['r','d','s'])
    df['value'] = [Yrds[k].x for k in Yrds.keys()]
    sparse = df.query('value>0')
    sparse.drop(columns=['value'], inplace=True)
    # sparse only has keys from Yrds that have values = 1, so only scheduled slots and only the starting slot of the observation
    sparse.to_csv(manager.output_directory + "serialized_outputs_sparse.csv", index=False, na_rep="")
    day, slot = np.mgrid[:manager.semester_length,:manager.n_slots_in_night]
    dense1 = pd.DataFrame(dict(d=day.flatten(), s=slot.flatten()))
    dense1 = pd.merge(dense1, sparse, left_on=['d','s'],right_on=['d','s'],how='left')
    dense1['r'] = dense1['r'].fillna('')
    # dense1 has keys for all days and slots, where no star was scheduled to start its observation, the r column is blank
    dense1.to_csv(manager.output_directory + "serialized_outputs_dense_v1.csv", index=False, na_rep="")

    dense2 = dense1.copy()
    isNight = np.array(manager.twilight_map_remaining_2D).flatten()
    isAlloc = np.array(manager.allocation_map_2D).flatten()
    isClear = np.array(manager.weathered_map).flatten()
    # have to go backwards otherwise you're adding stars into slots and then testing if the star is in the next slot
    for slot in range(manager.n_slots_in_semester-1, -1, -1):
        name_string = ""
        if isNight[slot] == 0:
            name_string += "*"
        if isAlloc[slot] == 0:
            name_string += "X"
        if isClear[slot] == 1:
            name_string += "W"
        dense2['r'][slot] = name_string + str(dense2['r'][slot])

        if dense2['r'][slot] in list(manager.requests_frame['starname']):
            slots_needed = manager.slots_needed_for_exposure_dict[dense2['r'][slot]]
            if slots_needed > 1:
                for t in range(1, slots_needed):
                    dense2['r'][slot + t] = str(dense2['r'][slot + t]) + str(dense2['r'][slot])

    # dense2 has keys for all days and slots, manually fill in the reserved slots for each observation and fill in Past/Twilight/Weather info
    dense2.to_csv(manager.output_directory + "serialized_outputs_dense_v2.csv", index=False, na_rep="")

def write_starlist(frame, solution_frame, night_start_time, extras, filler_stars, current_day,
                    outputdir, version='nominal'):
    """
    Generate the nightly script in the correct format.

    Args:
        frame (dataframe): the csv of PI requests for just the targets that were selected
                            to be observed tonight
        solution_frame (dataframe): the solution to the TTP model.plotly
        night_start_time (astropy time object): Beginning of observing interval
        extras (array): starnames of "extra" stars (those not fit into the script)
        filler_stars (array): star names of the stars added in the bonus round
        current_day (str): today's date in format YYYY-MM-DD
        outputdir (str): the directory to save the script file

    Returns:
        None
    """
    total_exptime = 0
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    script_file = os.path.join(outputdir,'script_{}_{}.txt'.format(current_day, version))

    lines = []
    for i, item in enumerate(solution_frame['Starname']):
        filler_flag = solution_frame['Starname'][i] in filler_stars
        row = frame.loc[frame['starname'] == solution_frame['Starname'][i]]
        row.reset_index(inplace=True)
        total_exptime += float(row['exptime'][0])

        start_exposure_hst = str(TimeDelta(solution_frame['Start Exposure'][i]*60,format='sec') + \
                                                night_start_time)[11:16]
        first_available_hst = str(TimeDelta(solution_frame['First Available'][i]*60,format='sec')+ \
                                                night_start_time)[11:16]
        last_available_hst = str(TimeDelta(solution_frame['Last Available'][i]*60,format='sec') + \
                                                night_start_time)[11:16]

        lines.append(format_kpf_row(row, start_exposure_hst, first_available_hst,last_available_hst,
                                    current_day, filler_flag = filler_flag))

    lines.append('')
    lines.append('X' * 45 + 'EXTRAS' + 'X' * 45)
    lines.append('')

    for j in range(len(extras['Starname'])):
        if extras['Starname'][j] in filler_stars:
            filler_flag = True
        else:
            filler_flag = False
        row = frame.loc[frame['starname'] == extras['Starname'][j]]
        row.reset_index(inplace=True)
        lines.append(format_kpf_row(row, '56:78', extras['First Available'][j],
                    extras['Last Available'][j], current_day, filler_flag, True))

    # add buffer lines to end of file
    lines.append("")
    lines.append("")

    with open(script_file, 'w') as f:
        f.write('\n'.join(lines))
    print("Total Open Shutter Time Scheduled: " + str(np.round((total_exptime/3600),2)) + " hours")
    return lines

# Define column names in advance
columns = [
    'star_name', 'ra', 'dec', 'equinox', 'jmag', 'gmag', 'teff',
    'gaia_dr3_id', 'program_code', 'priority', 'semester',
    'obs_time', 'first_avail', 'last_avail'
]
def parse_star_line(line):
    line = line.strip()

    # Case 1: blank line
    if not line:
        return {col: '' for col in columns}

    # Case 2: line with Xs
    if 'EXTRAS' in line and re.fullmatch(r'X*EXTRASX*', line):
        row = {col: 'XX' for col in columns}
        row['gaia_dr3_id'] = 'EXTRAS'
        return row

    # Normal case: parse structured line
    tokens = line.split()

    try:
        star_name = tokens[0]
        ra = f"{tokens[1]}h{tokens[2]}m{tokens[3]}s"
        dec = f"{tokens[4]}d{tokens[5]}m{tokens[6]}s"
        equinox = tokens[7]

        jmag = re.search(r'jmag=([\d.]+)', line)
        gmag = re.search(r'gmag=([\d.]+)', line)
        teff = re.search(r'Teff=([\d.]+)', line)

        jmag = jmag.group(1) if jmag else ''
        gmag = gmag.group(1) if gmag else ''
        teff = teff.group(1) if teff else ''

        gaia_id_match = re.search(r'DR[23]_\d+', line)
        gaia_id = gaia_id_match.group(0) if gaia_id_match else ''

        post_gaia_tokens = line.split(gaia_id)[-1].strip().split() if gaia_id else []

        program_code = post_gaia_tokens[0] if len(post_gaia_tokens) > 0 else ''
        priority = post_gaia_tokens[1] if len(post_gaia_tokens) > 1 else ''
        semester = post_gaia_tokens[2] if len(post_gaia_tokens) > 2 else ''
        obs_time      = str(post_gaia_tokens[3]) if len(post_gaia_tokens) > 3 else ''
        first_avail   = str(post_gaia_tokens[4]) if len(post_gaia_tokens) > 4 else ''
        last_avail    = str(post_gaia_tokens[5]) if len(post_gaia_tokens) > 5 else ''

        return {
            'star_name': star_name,
            'ra': ra,
            'dec': dec,
            'equinox': equinox,
            'jmag': jmag,
            'gmag': gmag,
            'teff': teff,
            'gaia_dr3_id': gaia_id,
            'program_code': program_code,
            'priority': priority,
            'semester': semester,
            'obs_time': obs_time,
            'first_avail': first_avail,
            'last_avail': last_avail
        }

    except Exception as e:
        # If parsing fails, return blank row (optional)
        return {col: ' ' for col in columns}


def format_kpf_row(row, obs_time, first_available, last_available, current_day,
                    filler_flag = False, extra=False):
    """
    Format request data in the specific way needed for the script (relates to the Keck "Magiq"
    software's data ingestion requirements).

    Args:
        row (dataframe): a single row from the requests sheet dataframe
        obs_time (str): the timestamp of the night to begin the exposure according to the TTP.
                        In format HH:MM in HST timezone
        first_available (str): the timestamp of the night where the star is first accessible.
                                In format HH:MM in HST timezone.
        last_available (str): the timestamp of the night where the star is last accessible.
                                In format HH:MM in HST timezone.
        filler_flag (boolean): True of the target was added in the bonus round
        extra (boolean): is this an "extra" target

    Returns:
        line (str): the properly formatted string to be included in the script file
    """

    equinox = '2000'
    updated_ra, updated_dec = pm_correcter(row['ra'][0], row['dec'][0],
                                row['pmra'][0], row['pmdec'][0], current_day, equinox=equinox)
    if updated_dec[0] != "-":
        updated_dec = "+" + updated_dec

    cpsname = hs.cps_star_name(row['starname'][0])
    namestring = ' '*(16-len(cpsname[:16])) + cpsname[:16]

    jmagstring = ('jmag=' + str(np.round(float(row['jmag'][0]),1)) + ' '* \
        (4-len(str(np.round(row['jmag'][0],1)))))
    exposurestring = (' '*(4-len(str(int(row['exptime'][0])))) + \
        str(int(row['exptime'][0])) + '/' + \
        str(int(row['exptime'][0])) + ' '* \
        (4-len(str(int(row['exptime'][0])))))

    ofstring = ('1of' + str(int(row['n_intra_max'][0])))
    scstring = 'sc=' + 'T'

    numstring = str(int(row['n_exp'][0])) + "x"
    gmagstring = 'gmag=' + str(np.round(float(row['gmag'][0]),1)) + \
                                                ' '*(4-len(str(np.round(row['gmag'][0],1))))
    teffstr = 'Teff=' + str(int(row['teff'][0])) + \
                                    ' '*(4-len(str(int(row['teff'][0]))))

    gaiastring = str(row['gaia_id'][0]) + ' '*(25-len(str(row['gaia_id'][0])))
    programstring = row['program_code'][0]

    if filler_flag:
        # All targets added in round 2 bonus round are lower priority
        priostring = "p3"
    else:
        priostring = "p1"

    if extra == False:
        timestring2 = str(obs_time)
    else:
        # designate a nonsense time
        timestring2 = "56:78"

    line = (namestring + ' ' + updated_ra + ' ' + updated_dec + ' ' + str(equinox) + ' '
                + jmagstring + ' ' + exposurestring + ' ' + ofstring + ' ' + scstring +  ' '
                + numstring + ' '+ gmagstring + ' ' + teffstr + ' ' + gaiastring + ' CC '
                        + priostring + ' ' + programstring + ' ' + timestring2 +
                         ' ' + first_available  + ' ' + last_available )

    # if not pd.isnull(row['Observing Notes'][0]):
    #     line += (' ' + str(row['Observing Notes'][0]))

    return line

def pm_correcter(ra, dec, pmra, pmdec, current_day, equinox="2000"):
    """
    Update a star's coordinates due to proper motion.

    Args:
        ra (float): RA in degrees
        dec (float): Dec in degrees
        pmra (float): proper motion in RA (mas/yr), including cos(Dec)
        pmdec (float): proper motion in Dec (mas/yr)
        equinox (str): original epoch (e.g. '2000.0')
        current_day (str): date to which to propagate (e.g. '2025-04-30')

    Returns:
        formatted_ra (str), formatted_dec (str): updated coordinates as strings
    """
    start_time = Time(f'J{equinox}')
    current_time = Time(current_day)
    coord = SkyCoord(
        ra=ra * u.deg,
        dec=dec * u.deg,
        pm_ra_cosdec=pmra * u.mas/u.yr,
        pm_dec=pmdec * u.mas/u.yr,
        obstime=start_time
    )
    new_coord = coord.apply_space_motion(new_obstime=current_time)
    formatted_ra = new_coord.ra.to_string(unit=u.hourangle, sep=' ', pad=True, precision=1)
    formatted_dec = new_coord.dec.to_string(unit=u.deg, sep=' ', pad=True, precision=0)

    return formatted_ra, formatted_dec
