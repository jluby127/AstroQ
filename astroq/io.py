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
        # Simplified - assume all slots are allocated
        allocated = manager.n_nights_in_semester * manager.n_slots_in_night
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

def serialize_schedule(Yrds, planner):
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
    sparse.to_csv(planner.output_directory + "semester_plan.csv", index=False, na_rep="")

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
    # Handle missing pmra/pmdec columns with default values
    pmra = row.get('pmra', [0.0])[0] if 'pmra' in row else 0.0
    pmdec = row.get('pmdec', [0.0])[0] if 'pmdec' in row else 0.0
    updated_ra, updated_dec = pm_correcter(row['ra'][0], row['dec'][0],
                                pmra, pmdec, current_day, equinox=equinox)
    if updated_dec[0] != "-":
        updated_dec = "+" + updated_dec

    cpsname = hs.crossmatch_star_name(row['starname'][0])
    namestring = ' '*(16-len(cpsname[:16])) + cpsname[:16]

    # Handle missing columns with default values
    jmag_val = row.get('jmag', [15.0])[0] if 'jmag' in row else 15.0
    gmag_val = row.get('gmag', [15.0])[0] if 'gmag' in row else 15.0
    teff_val = row.get('teff', [5000])[0] if 'teff' in row else 5000
    gaia_id_val = row.get('gaia_id', ['UNKNOWN'])[0] if 'gaia_id' in row else 'UNKNOWN'
    
    jmagstring = ('jmag=' + str(np.round(float(jmag_val),1)) + ' '* \
        (4-len(str(np.round(jmag_val,1)))))
    exposurestring = (' '*(4-len(str(int(row['exptime'][0])))) + \
        str(int(row['exptime'][0])) + '/' + \
        str(int(row['exptime'][0])) + ' '* \
        (4-len(str(int(row['exptime'][0])))))

    ofstring = ('1of' + str(int(row['n_intra_max'][0])))
    scstring = 'sc=' + 'T'

    numstring = str(int(row['n_exp'][0])) + "x"
    gmagstring = 'gmag=' + str(np.round(float(gmag_val),1)) + \
                                                ' '*(4-len(str(np.round(gmag_val,1))))
    teffstr = 'Teff=' + str(int(teff_val)) + \
                                    ' '*(4-len(str(int(teff_val))))

    gaiastring = str(gaia_id_val) + ' '*(25-len(str(gaia_id_val)))
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


