"""
Module for building reports after the solution is found. Designed to be only run as a function
call from the generateScript.py script.

Example usage:
    import reporting_functions as rf
"""

# Standard library imports
import math
import os
import re
import time

# Third-party imports
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time, TimeDelta
from astropy import units as u
import numpy as np
import pandas as pd

# Local imports
import astroq.history as hs
import astroq.access as ac

def serialize_schedule(Yrds, semester_planner):
    """
    Turns the non-square matrix of the solution into a square matrix and starts the human readable
    solution by filling in the slots where a star's exposre is started.

    Args:
        Yrds (array): the Gurobi solution with keys of (starname, day, slot) and values 1 or 0.
        planner (obj): a SemesterPlanner object

    Returns:
        None
    """
    df = pd.DataFrame(Yrds.keys(),columns=['r','d','s'])
    df['value'] = [Yrds[k].x for k in Yrds.keys()]
    sparse = df.query('value>0').copy()
    sparse.drop(columns=['value'], inplace=True)
    sparse.to_csv(semester_planner.output_directory + "semester_plan.csv", index=False, na_rep="")
    semester_planner.future_forecast = "semester_plan.csv"

    day, slot = np.mgrid[:semester_planner.semester_length,:semester_planner.n_slots_in_night]
    dense1 = pd.DataFrame(dict(d=day.flatten(), s=slot.flatten()))
    dense1 = pd.merge(dense1, sparse, left_on=['d','s'],right_on=['d','s'],how='left')
    dense1['r'] = dense1['r'].fillna('')
    # dense1 has keys for all days and slots, where no star was scheduled to start its observation, the r column is blank
    dense1.to_csv(semester_planner.output_directory + "serialized_outputs_dense_v1.csv", index=False, na_rep="")

    dense2 = dense1.copy()
    # Use the stored access record from the first run (no need to recompute)
    access = semester_planner.access_record
    isAlloc = access['is_alloc'].flatten()
    isClear = access['is_clear'].flatten()
    # have to go backwards otherwise you're adding stars into slots and then testing if the star is in the next slot
    for slot in range(semester_planner.n_slots_in_semester-1, -1, -1):
        name_string = ""
        if isAlloc[slot] == 0:
            name_string += "X"
        if isClear[slot] == 1:
            name_string += "W"
        dense2.loc[slot, 'r'] = name_string + str(dense2.loc[slot, 'r'])

        if dense2.loc[slot, 'r'] in list(semester_planner.requests_frame['starname']):
            slots_needed = semester_planner.slots_needed_for_exposure_dict[dense2.loc[slot, 'r']]
            if slots_needed > 1:
                for t in range(1, slots_needed):
                    dense2.loc[slot + t, 'r'] = str(dense2.loc[slot + t, 'r']) + str(dense2.loc[slot, 'r'])

    # dense2 has keys for all days and slots, manually fill in the reserved slots for each observation and fill in Past/Twilight/Weather info
    dense2.to_csv(semester_planner.output_directory + "serialized_outputs_dense_v2.csv", index=False, na_rep="")
    
    # Generate the fullness report
    build_fullness_report(semester_planner, "Round1")

def build_fullness_report(semester_planner, round_info):
    """
    Determine how full the schedule is: slots available, slots scheduled, and slots required

    Args:
        semester_planner (obj): a SemesterPlanner object
        round_info (str): information about the optimization round

    Returns:
        None
    """
    file_path = semester_planner.output_directory + "runReport.txt"
    print(f"Writing to: {file_path}")
    
    # Read the semester plan CSV file
    semester_plan_path = os.path.join(semester_planner.output_directory, "semester_plan.csv")
    if not os.path.exists(semester_plan_path):
        print(f"Warning: semester_plan.csv not found at {semester_plan_path}")
        return
    
    schedule_df = pd.read_csv(semester_plan_path)
    
    # Get access record data
    access = semester_planner.access_record
    is_alloc = access['is_alloc']
    
    # Calculate statistics
    # is_alloc is 3D array, but all 2D arrays are the same, so use the first one
    is_alloc_2d = is_alloc[0]  # Take the first 2D array
    total_slots_in_semester = is_alloc_2d.shape[0] * is_alloc_2d.shape[1]  # Total slots = n_nights * n_slots_per_night
    allocated_slots = np.sum(is_alloc_2d)  # Slots that are allocated (not X)
    
    # Calculate scheduled slots (starting slots only)
    scheduled_starting_slots = len(schedule_df)  # Number of starting slots with stars scheduled
    
    # Calculate reserved slots (slots beyond starting slots for multi-slot exposures)
    reserved_slots = 0
    for _, row in schedule_df.iterrows():
        star_name = row['r']
        if star_name in semester_planner.slots_needed_for_exposure_dict:
            slots_needed = semester_planner.slots_needed_for_exposure_dict[star_name]
            reserved_slots += slots_needed - 1  # Subtract 1 because the starting slot is already counted
    
    total_scheduled_slots = scheduled_starting_slots + reserved_slots
    empty_slots = allocated_slots - total_scheduled_slots
    
    # Calculate total slots requested using slots_needed_for_exposure_dict and n_intra_max
    total_slots_requested = 0
    for i in range(len(semester_planner.requests_frame)):
        star_name = semester_planner.requests_frame['starname'][i]
        if star_name in semester_planner.slots_needed_for_exposure_dict:
            slots_needed = semester_planner.slots_needed_for_exposure_dict[star_name]
            total_slots_requested += slots_needed * semester_planner.requests_frame['n_intra_max'][i] * semester_planner.requests_frame['n_inter_max'][i]
    
    # Calculate percentages
    percentage_of_available = np.round((total_scheduled_slots * 100) / allocated_slots, 3) if allocated_slots > 0 else 0
    percentage_of_requested = np.round((total_scheduled_slots * 100) / total_slots_requested, 3) if total_slots_requested > 0 else 0
    
    with open(semester_planner.output_directory + "runReport.txt", "w") as file:
        file.write("Stats for " + str(round_info) + "\n")
        file.write("------------------------------------------------------" + "\n")
        file.write("N slots in semester:" + str(total_slots_in_semester) + "\n")
        file.write("N available slots:" + str(allocated_slots) + "\n")
        file.write("N starting slots scheduled: " + str(scheduled_starting_slots) + "\n")
        file.write("N reserved slots: " + str(reserved_slots) + "\n")
        file.write("N total slots scheduled: " + str(total_scheduled_slots) + "\n")
        file.write("N slots left empty: " + str(empty_slots) + "\n")
        file.write("N slots requested (total): " + str(total_slots_requested) + "\n")
        file.write("Utilization (% of available slots): " + str(percentage_of_available) + "%" + "\n")
        file.write("Utilization (% of requested slots): " + str(percentage_of_requested) + "%" + "\n")
        file.close()
        
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
    # Cast starname column to strings to ensure proper matching
    frame['starname'] = frame['starname'].astype(str)
    
    # Cast extras star names to strings to ensure proper matching
    if extras is not None and len(extras) > 0:
        if hasattr(extras, 'astype'):
            # If extras is a DataFrame
            extras['Starname'] = extras['Starname'].astype(str)
        else:
            # If extras is a list, convert each star name to string
            extras['Starname'] = [str(star) for star in extras['Starname']]
    
    total_exptime = 0
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    script_file = os.path.join(outputdir,'script_{}_{}.txt'.format(current_day, version))

    lines = []
    for i, item in enumerate(solution_frame['Starname']):
        filler_flag = solution_frame['Starname'][i] in filler_stars
        row = frame.loc[frame['starname'] == solution_frame['Starname'][i]]
        row.reset_index(inplace=True)
        total_exptime += float(row['exptime'].iloc[0])

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
    # Handle missing pmra/pmdec columns with default values
    pmra = row.get('pmra', pd.Series([0.0])).iloc[0] if 'pmra' in row else 0.0
    pmdec = row.get('pmdec', pd.Series([0.0])).iloc[0] if 'pmdec' in row else 0.0
    updated_ra, updated_dec = pm_correcter(row['ra'].iloc[0], row['dec'].iloc[0],
                                pmra, pmdec, current_day, equinox=equinox)
    if updated_dec[0] != "-":
        updated_dec = "+" + updated_dec

    cpsname = hs.crossmatch_star_name(row['starname'].iloc[0])
    namestring = ' '*(16-len(cpsname[:16])) + cpsname[:16]

    # Handle missing columns with default values
    jmag_val = row.get('jmag', [15.0])[0] if 'jmag' in row else 15.0
    gmag_val = row.get('gmag', [15.0])[0] if 'gmag' in row else 15.0
    teff_val = row.get('teff', [5000])[0] if 'teff' in row else 5000
    gaia_id_val = row.get('gaia_id', ['UNKNOWN'])[0] if 'gaia_id' in row else 'UNKNOWN'
    
    # Convert to float safely, with fallback to defaults
    try:
        jmag_val = float(jmag_val) if jmag_val is not None else 15.0
    except (ValueError, TypeError):
        jmag_val = 25.0
    
    try:
        gmag_val = float(gmag_val) if gmag_val is not None else 15.0
    except (ValueError, TypeError):
        gmag_val = 25.0
    
    try:
        teff_val = float(teff_val) if teff_val is not None else 5000
    except (ValueError, TypeError):
        teff_val = 0.0
    
    jmagstring = ('jmag=' + str(np.round(float(jmag_val),1)) + ' '* \
        (4-len(str(np.round(float(jmag_val),1)))))
    exposurestring = (' '*(4-len(str(int(row['exptime'].iloc[0])))) + \
        str(int(row['exptime'].iloc[0])) + '/' + \
        str(int(row['exptime'].iloc[0])) + ' '* \
        (4-len(str(int(row['exptime'].iloc[0])))))

    ofstring = ('1of' + str(int(row['n_intra_max'].iloc[0])))
    scstring = 'sc=' + 'T'

    numstring = str(int(row['n_exp'].iloc[0])) + "x"
    gmagstring = 'gmag=' + str(np.round(float(gmag_val),1)) + \
                                                ' '*(4-len(str(np.round(float(gmag_val),1))))
    teffstr = 'Teff=' + str(int(teff_val)) + \
                                    ' '*(4-len(str(int(teff_val))))

    gaiastring = str(gaia_id_val) + ' '*(25-len(str(gaia_id_val)))
    programstring = row['program_code'].iloc[0]

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

    # Handle missing Observing Notes column
    observing_notes = row.get('Observing Notes', [''])[0] if 'Observing Notes' in row else ''
    if observing_notes and not pd.isnull(observing_notes):
        line += (' ' + str(observing_notes))

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


