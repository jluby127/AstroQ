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

import numpy as np
import pandas as pd

def build_fullness_report(combined_semester_schedule, allocation_map_2D, manager, round_info):
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
    file = open(manager.output_directory + "runReport.txt", "a")
    file.write("Stats for " + str(round_info) + "\n")
    file.write("------------------------------------------------------" + "\n")
    listnames = list(manager.requests_frame['Starname'])
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
    allocated = np.sum(allocation_map_2D.flatten())
    file.write("N slots in semester:" + str(np.prod(combined_semester_schedule.shape)) + "\n")
    file.write("N available slots:" + str(allocated) + "\n")
    file.write("N slots scheduled: " + str(used) + "\n")
    file.write("N slots left empty: " + str(allocated-used) + "\n")

    total_slots_requested = 0
    for i in range(len(manager.requests_frame)):
        total_slots_requested += manager.requests_frame['# of Nights Per Semester'][i]* \
            math.ceil(manager.requests_frame['Nominal Exposure Time [s]'][i]/(manager.slot_size*60.))
    file.write("N slots requested (total): " + str(total_slots_requested) + "\n")
    percentage = np.round((used*100)/allocated,3)
    file.write("Percent full: " + str(percentage) + "%." + "\n")
    file.close()

def report_allocation_stats(manager, allocation_map_2D):
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
    ff = open(outputdir + "allocation_stats.txt", "w")

    q1s = 0
    q2s = 0
    q3s = 0
    q4s = 0
    count0 = 0
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0
    for j, item in enumerate(allocation_map_2D):
        holder = [0, 0, 0, 0]
        if allocation_map_2D[j][0] == 1.:
            q1s += 1
            date_info = manager.all_dates_array[manager.current_day_number + j] + " - q0"
            ff.write(date_info + "\n")
            holder[0] = 1
        if allocation_map_2D[j][1] == 1.:
            q2s += 1
            date_info = manager.all_dates_array[manager.current_day_number + j] + " - q1"
            ff.write(date_info + "\n")
            holder[1] = 1
        if allocation_map_2D[j][2] == 1.:
            q3s += 1
            date_info = manager.all_dates_array[manager.current_day_number + j] + " - q2"
            ff.write(date_info + "\n")
            holder[2] = 1
        if allocation_map_2D[j][3] == 1.:
            q4s += 1
            date_info = manager.all_dates_array[manager.current_day_number + j] + " - q3"
            ff.write(date_info + "\n")
            holder[3] = 1
        allocated_quarters = np.sum(allocation_map_2D[j])
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

def write_out_results(manager, model, round, start_the_clock):
    """
    Write the Run Report
    Args:
        manager (obj): a data_admin object
        model (obj): a GurobiModel object
        rount (str): "Round1" or "Round2"
        start_the_clock (obj): time object of when the scheduler code began

    Returns:
        None
    """
    print("Writing Report.")
    filename = open(manager.output_directory + "runReport.txt", "a")
    theta_n_var = []
    counter = 0
    for v in model.theta.values():
        varname = v.VarName
        varval = v.X
        counter += varval
    print("Sum of Theta: " + str(counter))
    filename.write("Sum of Theta: " + str(counter) + "\n")

    print("Total Time to complete " + round + ": " + str(np.round(time.time()-start_the_clock,3)))
    filename.write("Total Time to complete " + round +  ": " + str(np.round(time.time()-start_the_clock,3)) + "\n")
    filename.write("\n")
    filename.write("\n")
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
    for name in list(manager.requests_frame['Starname']):
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
        for name in list(manager.requests_frame['Starname']):
            if all_star_schedules[name][s] == 1:
                slotallocated += str(name)
        combined_semester_schedule[s] += str(slotallocated)
    combined_semester_schedule = np.reshape(combined_semester_schedule,
            (manager.semester_length, manager.n_slots_in_night))

    # The semester solver puts a 1 only in the slot that starts the exposure for a target.
    # Therefore, many slots are empty because they are part of a multi-slot visit.
    # Here fill in the multi-slot exposures appropriately for ease of human reading and accounting.
    for n in range(manager.n_nights_in_semester-1-manager.all_dates_dict[manager.current_day], -1, -1):
        for s in range(manager.n_slots_in_night-1, -1, -1):
            if combined_semester_schedule[n+manager.all_dates_dict[manager.current_day]][s] in list(manager.requests_frame['Starname']):
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

def write_available_human_readable(manager, twilight_map, allocation_map_2D, weathered_map):
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
        combined_semester_schedule[c] += 'supercalifragilisticexpialidocious'
    combined_semester_schedule = np.reshape(combined_semester_schedule,
            (manager.semester_length, manager.n_slots_in_night))

    for n in range(manager.semester_length - manager.all_dates_dict[manager.current_day]):
        for s in range(manager.n_slots_in_night):
            slotallocated = ''
            # remember that twilight map is "inverted": the 1's are time where it is night and the
            # 0's are time where it is day/twilight.
            if twilight_map[n][s] == 0:
                slotallocated += '*'
            if allocation_map_2D[n][s] == 0:
                slotallocated += 'X'
            if weathered_map[n][s] == 1:# and slotallocated == '':
                slotallocated += 'W'
            if os.path.exists(manager.nonqueue_map_file):
                slotallocated += str(nonqueuemap_slots_strs[n + manager.all_dates_dict[manager.current_day], ][s])
            combined_semester_schedule[n + manager.all_dates_dict[manager.current_day], ][s] += str(slotallocated)

    np.savetxt(manager.output_directory + 'raw_combined_semester_schedule_available.txt',
        combined_semester_schedule, delimiter=',', fmt="%s")
    return combined_semester_schedule

def write_starlist(frame, solution_frame, night_start_time, extras, filler_stars, current_day,
                    outputdir):
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
    script_file = os.path.join(outputdir,'script_{}_{}.txt'.format(current_day,'nominal'))
    print('Writing starlist to ' + script_file)

    lines = []
    for i, item in enumerate(solution_frame['Starname']):
        filler_flag = solution_frame['Starname'][i] in filler_stars
        row = frame.loc[frame['Starname'] == solution_frame['Starname'][i]]
        row.reset_index(inplace=True)
        total_exptime += float(row['Nominal Exposure Time [s]'][0])

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
        row = frame.loc[frame['Starname'] == extras['Starname'][j]]
        row.reset_index(inplace=True)
        lines.append(format_kpf_row(row, '56:78', extras['First Available'][j],
                    extras['Last Available'][j], current_day, filler_flag, True))

    # add buffer lines to end of file
    lines.append("")
    lines.append("")

    with open(script_file, 'w') as f:
        f.write('\n'.join(lines))
    print("Total Open Shutter Time Scheduled: " + str(np.round((total_exptime/3600),2)) + " hours")

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

    epochstr = '2024'
    updated_ra, updated_dec = pm_correcter(row['RA'][0], row['Dec'][0],
                                row['Proper Motion in RA [miliarcseconds/year]'][0],
                                row['Proper Motion in Dec [miliarcseconds/year]'][0],
                                epochstr, current_day, verbose=False)
    if updated_dec[0] != "-":
        updated_dec = "+" + updated_dec

    namestring = ' '*(16-len(row['Starname'][0][:16])) + row['Starname'][0][:16]

    jmagstring = ('jmag=' + str(np.round(float(row['J Magnitude'][0]),1)) + ' '* \
        (4-len(str(np.round(row['J Magnitude'][0],1)))))
    exposurestring = (' '*(4-len(str(int(row['Nominal Exposure Time [s]'][0])))) + \
        str(int(row['Nominal Exposure Time [s]'][0])) + '/' + \
        str(int(row['Maximum Exposure Time [s]'][0])) + ' '* \
        (4-len(str(int(row['Maximum Exposure Time [s]'][0])))))

    ofstring = ('1of' + str(int(row['# Visits per Night'][0])))

    if row['Simucal'][0]:
        scval = 'T'
    else:
        scval = 'F'
    scstring = 'sc=' + scval

    numstring = str(int(row['# of Exposures per Visit'][0])) + "x"
    gmagstring = 'gmag=' + str(np.round(float(row['G Magnitude'][0]),1)) + \
                                                ' '*(4-len(str(np.round(row['G Magnitude'][0],1))))
    teffstr = 'Teff=' + str(int(row['Effective Temperature [Kelvin]'][0])) + \
                                    ' '*(4-len(str(int(row['Effective Temperature [Kelvin]'][0]))))

    gaiastring = str(row['GAIA Identifier'][0]) + ' '*(25-len(str(row['GAIA Identifier'][0])))
    programstring = row['Program_Code'][0]

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

    line = (namestring + ' ' + updated_ra + ' ' + updated_dec + ' ' + str(epochstr) + ' '
                + jmagstring + ' ' + exposurestring + ' ' + ofstring + ' ' + scstring +  ' '
                + numstring + ' '+ gmagstring + ' ' + teffstr + ' ' + gaiastring + ' CC '
                        + priostring + ' ' + programstring + ' ' + timestring2 +
                         ' ' + first_available  + ' ' + last_available )

    if not pd.isnull(row['Observing Notes'][0]):
        line += (' ' + str(row['Observing Notes'][0]))

    return line

def pm_correcter(ra, dec, pmra, pmdec, epochstr, current_day, verbose=False):
    """
    Update a star's coordinates due to proper motion

    Args:
        ra (str): star's old coordinate RA in units of degrees
        dec (str): star's old coordinate Dec in units of degrees
        pmra (str): star's proper motion in the RA dimension, units of millearcseconds per year
        pmdec (str): star's proper motion in the Dec dimension, units of millearcseconds per year
        epochstr (str): a string of the year those coordinates are updated to
        verbose (boolean): True to print out to command line

    Returns:
        formatted_ra (str): the updated RA position in units of degrees
        formatted_dec (str):  the updated Dec position in units of degrees
    """
    # requires giving RA and Dec in degrees
    # example: RA = 321.5 and Dec = 15.6
    # note that degrees are not hour angles!
    # this code converts RA from degrees to hourangle at the end

    current_time = Time(current_day)  # You can adjust the date as needed
    ra_deg = Angle(ra, unit=u.deg)  # RA in degrees
    dec_deg = Angle(dec, unit=u.deg)  # Dec in degrees
    pm_ra = pmra * u.mas/u.yr
    pm_dec = pmdec * u.mas/u.yr
    epochtime = Time('J' + epochstr)
    ra_advanced_deg = (ra_deg + (pm_ra * (current_time - epochtime).to(u.yr)).to(u.deg))/15
    dec_advanced_deg = dec_deg + (pm_dec * (current_time - epochtime).to(u.yr)).to(u.deg)

    formatted_ra = ra_advanced_deg.to_string(unit=u.deg, sep=' ', pad=True, precision=1)
    formatted_dec = dec_advanced_deg.to_string(unit=u.deg, sep=' ', pad=True, precision=0)

    return formatted_ra, formatted_dec
