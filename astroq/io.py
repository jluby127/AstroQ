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

import astroq.history as hs
import astroq.access as ac

def build_fullness_report(serialized_output_sparse, manager, round_info):
    """
    Determine how full the schedule is: slots available, slots scheduled, and slots required

    Args:
        serialized_output_sparse (DataFrame): DataFrame with columns ['r', 'd', 's'] containing scheduled targets
        manager (obj): a data_admin object
        round_info (str): Round information for the report

    Returns:
        None
    """
    import astroq.access as ac
    import math
    
    file_path = manager.output_directory + "runReport.txt"
    print(f"Writing to: {file_path}")
    
    # Get access records array to determine available slots
    access = ac.produce_ultimate_map(manager, manager.requests_frame)
    
    # Calculate total available slots using is_alloc (sum of one slice = available slots)
    total_available_slots = np.sum(access['is_alloc'][0])  # Use first target's slice since all targets have same allocation
    
    # Calculate total scheduled slots from serialized output
    total_scheduled_slots = len(serialized_output_sparse)
    
    # Add reserved slots for multi-slot exposures
    total_reserved_slots = 0
    for _, row in serialized_output_sparse.iterrows():
        target_name = row['r']
        # Find the target in requests frame
        target_request = manager.requests_frame[manager.requests_frame['starname'] == target_name]
        if len(target_request) > 0:
            # Get slots needed for this exposure
            slots_needed = manager.slots_needed_for_exposure_dict.get(target_name, 1)
            # Add reserved slots (excluding the starting slot which is already counted)
            total_reserved_slots += max(0, slots_needed - 1)
    
    # Total slots used (scheduled + reserved)
    total_slots_used = total_scheduled_slots + total_reserved_slots
    
    # Calculate total slots requested
    total_slots_requested = 0
    for _, row in manager.requests_frame.iterrows():
        slots_needed_per_visit = manager.slots_needed_for_exposure_dict[row['starname']]
        total_slots_requested += int(row['n_inter_max']) * int(slots_needed_per_visit) * int(row['n_intra_max'])
    
    # Calculate utilization percentages
    utilization_available = (total_slots_used / total_available_slots * 100) if total_available_slots > 0 else 0
    utilization_requested = (total_slots_used / total_slots_requested * 100) if total_slots_requested > 0 else 0
    
    # Write report
    with open(manager.output_directory + "runReport.txt", "a") as file:
        file.write("Stats for " + str(round_info) + "\n")
        file.write("------------------------------------------------------" + "\n")
        file.write(f"Total starting slots: {total_scheduled_slots} ({total_scheduled_slots * manager.slot_size / 60:.1f} hours)\n")
        file.write(f"Total reserved slots: {total_reserved_slots} ({total_reserved_slots * manager.slot_size / 60:.1f} hours)\n")
        file.write(f"Total scheduled slots: {total_slots_used} ({total_slots_used * manager.slot_size / 60:.1f} hours)\n")
        file.write(f"Total requested slots: {total_slots_requested} ({total_slots_requested * manager.slot_size / 60:.1f} hours)\n")
        file.write(f"Utilization of requested time: {utilization_requested:.1f}%\n")
        file.write("---\n")
        file.write(f"Total available slots: {total_available_slots} ({total_available_slots * manager.slot_size / 60:.1f} hours)\n")
        file.write(f"Empty available slots: {total_available_slots - total_slots_used} ({(total_available_slots - total_slots_used) * manager.slot_size / 60:.1f} hours)\n")
        file.write(f"Utilization of available time: {utilization_available:.1f}%\n")
        file.write("---\n")

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
    sparse = df.query('value>0').copy()
    sparse.drop(columns=['value'], inplace=True)
    sparse.to_csv(manager.output_directory + "serialized_outputs_sparse.csv", index=False, na_rep="")
    
    day, slot = np.mgrid[:manager.semester_length,:manager.n_slots_in_night]
    dense1 = pd.DataFrame(dict(d=day.flatten(), s=slot.flatten()))
    dense1 = pd.merge(dense1, sparse, left_on=['d','s'],right_on=['d','s'],how='left')
    dense1['r'] = dense1['r'].fillna('')
    # dense1 has keys for all days and slots, where no star was scheduled to start its observation, the r column is blank
    dense1.to_csv(manager.output_directory + "serialized_outputs_dense_v1.csv", index=False, na_rep="")

    dense2 = dense1.copy()
    access = ac.produce_ultimate_map(manager, manager.requests_frame) #temp until we pickle the manager and read it in
    isAlloc = access['is_alloc'].flatten()
    isClear = access['is_clear'].flatten()
    # have to go backwards otherwise you're adding stars into slots and then testing if the star is in the next slot
    for slot in range(manager.n_slots_in_semester-1, -1, -1):
        name_string = ""
        if isAlloc[slot] == 0:
            name_string += "X"
        if isClear[slot] == 1:
            name_string += "W"
        dense2.loc[slot, 'r'] = name_string + str(dense2.loc[slot, 'r'])

        if dense2.loc[slot, 'r'] in list(manager.requests_frame['starname']):
            slots_needed = manager.slots_needed_for_exposure_dict[dense2.loc[slot, 'r']]
            if slots_needed > 1:
                for t in range(1, slots_needed):
                    dense2.loc[slot + t, 'r'] = str(dense2.loc[slot + t, 'r']) + str(dense2.loc[slot, 'r'])

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

    if not pd.isnull(row['Observing Notes'][0]):
        line += (' ' + str(row['Observing Notes'][0]))

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


