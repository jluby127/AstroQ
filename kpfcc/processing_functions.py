"""
Module for processin the inputs and outputs of the autoscheduler to/from various sources.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import processing_functions as pf
"""
import os

from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.time import TimeDelta
import requests
import re
from bs4 import BeautifulSoup

import numpy as np
import pandas as pd

def access_jump():
    """
    Access the Jump database. You will need correct permissions to access this secured database
    hosted by Caltech Cadence servers. Eventually we will pull from the Keck Observatory Archive
    (KOA) directly.

    Args:
        None

    Returns:
        s (object): a session request
    """

    login_url = 'https://jump.caltech.edu/user/login/'
    s = requests.session()
    s.get(login_url)
    csrftoken = s.cookies['csrftoken']
    # you'll need to add your credentials for username and password
    payload = {'action':'login', 'username':'jblubin','password':'jumpinJackspass',
               'csrfmiddlewaretoken': csrftoken}
    new_login = s.post(login_url, data = payload, headers = dict(Referer = login_url))
    return s

def get_database_explorer(name, path_for_csv, url='https://jump.caltech.edu/explorer/', links=[]):
    """
    Pull the relevant query from the Jump database.
    Eventually we will pull from the Keck Observatory Archive (KOA) directly.

    Args:
        name (str): the name of the pre-built query in the Jump database
        path_for_csv (str): the link to the Jump database

    Returns:
        None
    """
    session = access_jump()
    response = session.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    table = soup.find("tbody", attrs={"class":"list"})
    for row in table.find_all("tr"):
        tab = row.find("td", attrs={"class":"name"})
        try:
            x = tab.a.string
        except:
            pass
        else:
            if x == name:
                for l in row.find_all("a", href = re.compile('download')):
                    links.append('/'.join(url.split('/')[:3])+l.get('href'))
    if links != []:
        response = session.get(links[0])
        open(path_for_csv, 'wb').write(response.content)
    else:
        print("Can't find that query name, try again. (needs to be exact)")
    session.keep_alive = False

def get_kpf_past_database(path_to_csv):
    """
    A user-friendly front facing function to pull past observations from the Jump database.
    Note that each semester, you will have to log into Jump, find this pre-built query,
    and manually update the start and end date fields.
    Eventually we will pull from the Keck Observatory Archive (KOA) directly.

    Args:
        path_to_csv (str): the directory to save the data

    Returns:
        None
    """
    name = 'Get all KPF Observations in Semester'
    get_database_explorer(name, path_to_csv)
    print("This semester's past KPF observations pulled from Jump. Saved to csv: " + path_to_csv)

def prepare_for_ttp(request_frame, night_plan, round_two_targets):
    """
    Prepare tonight's scheduled stars for their run through the TTP (separate software package)

    Args:
        request_sheet (dataframe): the csv of PI requests
        night_plan (array): the n'th row of combined_semester_schedule array
        round_two_targets (array): a 1D list of the stars that were added in the bonus round

    Returns:
        to_ttp (dataframe): the data on the stars to be observed tonight, formatted in the way that
                            the TTP software requires as an input
    """
    ignore = ['*', 'W', '', '*X', 'X']
    selected_stars = []
    for i, item in enumerate(night_plan):
        if night_plan[i] not in ignore and night_plan[i][:4] != "RM___":
            selected_stars.append(night_plan[i])
            ignore.append(night_plan[i])

    starnames = []
    ras = []
    decs = []
    exposure_times = []
    exposures_per_visit = []
    visits_in_night = []
    cadences = []
    priorities = []
    for j, item in enumerate(selected_stars):
        idx = request_frame.index[request_frame['Starname']==str(selected_stars[j])][0]
        starnames.append(str(request_frame['Starname'][idx]))
        ras.append(request_frame['RA'][idx])
        decs.append(request_frame['Dec'][idx])
        exposure_times.append(int(request_frame['Nominal Exposure Time [s]'][idx]))
        exposures_per_visit.append(int(request_frame['# of Exposures per Visit'][idx]))
        visits_in_night.append(int(request_frame['# Visits per Night'][idx]))
        cadences.append(int(request_frame['Minimum Intra-Night Cadence'][idx]))
        # higher numbers are higher priorities, filler targets get low priority
        if str(selected_stars[j]) in round_two_targets:
            prior = 1
        else:
            prior = 10
        priorities.append(prior)
    to_ttp = pd.DataFrame({"Starname":starnames,"RA":ras,"Dec":decs,
                          "Exposure Time":exposure_times,
                          "Exposures Per Visit":exposures_per_visit,
                          "Visits In Night":visits_in_night, "Intra_Night_Cadence":cadences,
                          "Priority":priorities})
    return to_ttp

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
