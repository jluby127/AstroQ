from collections import defaultdict
import numpy as np
import astropy.units as u
import pandas as pd
import math
import astropy as apy
import astroplan as apl
from astropy.coordinates import Angle
from astropy.time import Time
from astropy.time import TimeDelta
import requests
import re
from bs4 import BeautifulSoup
import os


def login_JUMP():
    """
    Access the Jump database. You will need correct permissions to access this secured database hosted by Caltech Cadence servers.
    Eventually this will no longer be needed, instead we will pull from the Keck Observatory Archive (KOA) directly.

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
    Eventually this will no longer be needed, instead we will pull from the Keck Observatory Archive (KOA) directly.

    Args:
        name (str): the name of the pre-built query in the Jump database
        path_for_csv (str): the link to the Jump database

    Returns:
        None
    """

    # log into JUMP and go to DataBase page
    session = login_JUMP()
    response = session.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    # find the table of queries
    table = soup.find("tbody", attrs={"class":"list"})
    # find the correct row by the query name
    for row in table.find_all("tr"):
        tab = row.find("td", attrs={"class":"name"})
        try:
            x = tab.a.string
        except:
            pass
        else:
            # once it finds the match, save the download link
            if x == name:
                for l in row.find_all("a", href = re.compile('download')):
                    links.append('/'.join(url.split('/')[:3])+l.get('href'))
    # make sure it finds it before saving
    if links != []:
        response = session.get(links[0])
        # pretty sure you can just read this directly into pandas if that's what you want
        open(path_for_csv, 'wb').write(response.content)
    else:
        print(" It isn't finding that query name, try again. (needs to be exact)")
    session.keep_alive = False
    return


def getKPFAllObservations(path_to_csv):
    """
    A user-friendly front facing function to pull past observations from the Jump database.
    Note that each semester, you will have to log into Jump, find this pre-built query, and manually update the start and end date fields.
    Eventually this will no longer be needed, instead we will pull from the Keck Observatory Archive (KOA) directly.

    Args:
        path_to_csv (str): the directory to save the resulting CSV of data pulled from the Jump database

    Returns:
        None
    """

    name = 'Get all KPF Observations in Semester'
    get_database_explorer(name, path_to_csv)
    print("All Past KPF observations for this semester pulled from Jump. Saved to csv: " + path_to_csv)


def getUniqueNights(star_past_obs, twilight):
    """
    Parse the Jump database for the previous nights where a given star was observed

    Args:
        star_past_obs (dataframe): a subset of the Jump database dataframe, filtered to only include observations for a specific target name
        twilight (dataframe): the dataframe of morning and evening twilight times for each night of the semester

    Returns:
        star_past_obs (dataframe): an updated version of the original variable which now includes columns for UTC date and hst date
        unique_hstdates_observed (array): a 1D array length # of unique dates observed where each element is the HST date when the star was observed at least once
        quarterObserved (array): (array) a 1D array length # of unique dates observed where each element the quarter of the night when the first of the observations for that night was taken.
    """

    unique_hstdates_observed = []
    unique_utcdates_observed = []
    quarterObserved = []
    all_hstdates = []
    for i in range(len(star_past_obs)):
        timestamp = star_past_obs['utctime'][i][:-6]
        timestamp = Time(timestamp, format = 'iso')
        timestamp.format = 'jd'
        utcdate = star_past_obs['utctime'][i][:10]
        unique_utcdates_observed.append(utcdate)

        convertHSTcivil = TimeDelta(60*24*60,format='sec').jd #Note this is arbitrary 24hr subtraction to get from UT date to Hawaii date
        hsttimestamp = timestamp - convertHSTcivil
        hsttimestamp.format = 'iso'
        hstdate = str(hsttimestamp)[:10]
        all_hstdates.append(hstdate)
        hsttimestamp.format = 'jd'
        if hstdate not in unique_hstdates_observed:
            unique_hstdates_observed.append(hstdate)
            quarterObserved.append(getQuartersObserved(timestamp, utcdate, twilight))

    star_past_obs['utcDate'] = unique_utcdates_observed
    star_past_obs['hstDate'] = all_hstdates

    return star_past_obs, unique_hstdates_observed, quarterObserved


def getQuartersObserved(timestamp, utcdate, twilight):
    """
    Determine which quarter of the night an observation was taken in

    Args:
        timestamp (float): a BJD value indicating the time of the observation
        utcdate (str): the calendar date of the observation in UTC, format 'YYYY-MM-DD'
        twilight (dataframe): the dataframe of morning and evening twilight times for each night of the semester

    Returns:
        star_past_obs (dataframe): an updated version of the original variable which now includes columns for UTC date and hst date
        unique_hstdates_observed (array): a 1D array length # of unique dates observed where each element is the HST date when the star was observed at least once
        quarterObserved (array): (array) a 1D array length # of unique dates observed where each element the quarter of the night when the first of the observations was taken.
    """

    dateInd = twilight.index.get_loc(twilight[twilight['time_utc'] == utcdate].index[0])
    start_jd = twilight['12_evening'][dateInd]
    end_jd = twilight['12_morning'][dateInd]
    lengthOfQuarter = (end_jd - start_jd)/4.
    rel_timestamp = float(str(timestamp - start_jd))

    quarter = -10
    if rel_timestamp < lengthOfQuarter:
        quarter = 0.5
    elif rel_timestamp < 2*lengthOfQuarter and rel_timestamp > lengthOfQuarter:
        quarter = 1.5
    elif rel_timestamp < 3*lengthOfQuarter and rel_timestamp > 2*lengthOfQuarter:
        quarter = 2.5
    elif rel_timestamp < 4*lengthOfQuarter and rel_timestamp > 3*lengthOfQuarter:
        quarter = 3.5
    elif rel_timestamp < 5*lengthOfQuarter and rel_timestamp > 4*lengthOfQuarter:
        # allow a little bit of leeway, even if this doesn't really make sense
        quarter = 3.5
    else:
        quarter = 0.5
        print("Houston, we've had a problem: target observed in an invalid quarter.")
        # Note to self: I'm not sure why but a lot of targets fall into this bin. Suspcicious.
    return quarter


def getNobs_on_Night(star_past_obs, unique_hstdates_observed):
    """
    Determine how many exposures were taken of a target on a given night

    Args:
        star_past_obs (dataframe): the updated version of the original variable which now includes columns for UTC date and hst date
        unique_hstdates_observed (array): a 1D array of unique HST dates where the star was observed

    Returns:
        Nobs_on_date (array): a 1D array of length equal to lenght unique_hstdates_observed where each element is an integer representing the number of exposures taken that night
    """

    Nobs_on_date = []
    for i in range(len(unique_hstdates_observed)):
        datemask = star_past_obs['hstDate'] == unique_hstdates_observed[i]
        Nobs_on_date.append(np.sum(datemask))
    return Nobs_on_date


def prepareTTP(request_sheet, night_plan, filltargets):
    """
    Prepare tonight's scheduled stars for their run through the TTP (separate software package)

    Args:
        request_sheet (dataframe): the csv of PI requests
        night_plan (array): the n'th row of combined_semester_schedule array representing night n where target names are put into slots
        filltargets (array): a 1D list of target names of the stars that were added in the bonus round

    Returns:
        toTTP (dataframe): the data on these stars formatted in the way that the TTP software requires as an input
    """

    ignore = ['*', 'W', '', '*X', 'X']
    selected_stars = []
    for i in range(len(night_plan)):
        if night_plan[i] not in ignore and night_plan[i][:4] != "RM__":
            selected_stars.append(night_plan[i])
            ignore.append(night_plan[i])

    # build the dataframe with correct info and correct headers
    all_targets_frame = pd.read_csv(request_sheet)
    starnames = []
    RAs = []
    Decs = []
    ExpTimes = []
    ExpsPerVisits = []
    nVisits = []
    cadences = []
    priorities = []
    for j in range(len(selected_stars)):
        idx = all_targets_frame.index[all_targets_frame['Starname']==str(selected_stars[j])][0]
        starnames.append(str(all_targets_frame['Starname'][idx]))
        RAs.append(all_targets_frame['RA'][idx])
        Decs.append(all_targets_frame['Dec'][idx])
        ExpTimes.append(int(all_targets_frame['Nominal Exposure Time [s]'][idx]))
        ExpsPerVisits.append(int(all_targets_frame['# of Exposures per Visit'][idx]))
        nVisits.append(int(all_targets_frame['# Visits per Night'][idx]))
        cadences.append(int(all_targets_frame['Minimum Intra-Night Cadence'][idx]))
        # higher numbers are higher priorities
        if str(selected_stars[j]) in filltargets:
            prior = 1 # the filler targets get low priority
        else:
            prior = 10
        priorities.append(prior)
    toTTP = pd.DataFrame({"Starname":starnames,"RA":RAs,"Dec":Decs,
                          "Exposure Time":ExpTimes,"Exposures Per Visit":ExpsPerVisits,
                          "Visits In Night":nVisits,"Intra_Night_Cadence":cadences,"Priority":priorities})
    return toTTP


def write_starlist(frame, orderedList, extras, gapFillers, condition, current_day, outputdir):
    """
    Generate the nightly script in the correct format.
    Note: this function is here only for completeness but the script generation is now fully done by the TTP standalone package

    Args:
        frame (dataframe): the csv of PI requests for just the targets that were selected to be observed tonight
        orderedList (array): a 1D array of the names of the targets selected to be observed tonight, in the order they were selected to be observed by the TTP
        extras (array): starnames of extra stars
        gapFillers (array): star names of the stars added in the bonus round
        condition (str): the weather conditions for this script. Currently only "nominal" is supported.
        current_day (str): today's date in format YYYY-MM-DD
        outputdir (str): the directory to save the script file

    Returns:
        None
    """

    total_exptime = 0

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    script_file = os.path.join(outputdir,'script_{}_{}.txt'.format(current_day,condition))
    print('Writing starlist to ' + script_file)

    lines = []
    for i in range(len(orderedList)):
        if orderedList['Target'][i] in gapFillers:
            fillerFlag = True
        else:
            fillerFlag = False
        row = frame.loc[frame['Starname'] == orderedList['Target'][i]]
        row.reset_index(inplace=True)
        total_exptime += float(row['Nominal Exposure Time [s]'][0])
        lines.append(format_kpf_row(row, orderedList['StartExposure'][i], current_day, fillerFlag = fillerFlag))

    lines.append('')
    lines.append('X' * 45 + 'EXTRAS' + 'X' * 45)
    lines.append('')

    for j in range(len(extras['Starname'])):
        if extras['Starname'][j] in gapFillers:
            fillerFlag = True
        else:
            fillerFlag = False
        row = frame.loc[frame['Starname'] == extras['Starname'][j]]
        row.reset_index(inplace=True)
        lines.append(format_kpf_row(row, '56:78', current_day, fillerFlag, False, True))

    # add buffer lines to end of file
    lines.append("")
    lines.append("")
    #Formatted starlists appear as text files in the directory
    with open(script_file, 'w') as f:
        f.write('\n'.join(lines))
    print("Total Open Shutter Time Scheduled: " + str(np.round((total_exptime/3600),2)) + " hours.")


def format_kpf_row(row, obs_time, current_day, fillerFlag = False, Jtwothousand = False, extra=False):
    """
    Format request data in the specific way needed for the script (relates to the Keck "Magiq" software's data ingestion requirements)
    Note: this function is here only for completeness but the script generation is now fully done by the TTP standalone package

    Args:
        row (dataframe): a single row from the requests sheet dataframe
        obs_time (str): the timestamp of the night to begin the exposure according to the TTP. In format HH:MM in HST timezone
        fillerFlag (boolean): true of the target was added in the bonus round
        Jtwothousand (boolean): False if the coordinates epoch is not J2000
        extra (boolean): is this an extra target

    Returns:
        line (str): the properly formatted string to be included in the script file
    """

    #Just a bunch of string formatting. This prints standard starlists as ordered by the salesman optimization
    if Jtwothousand:
        # convert the decimal RA/Dec into H:M:S format
        # these are the old J2000 coordinates
        coord = SkyCoord(ra=row['RA'][0]*u.degree, dec=row['Dec'][0]*u.degree, frame='icrs')
        temp = coord.to_string('hmsdms')
        radec = temp.split(' ')
        ra_h = radec[0][0:2]
        ra_m = radec[0][3:5]
        ra_s = radec[0][6:8]
        pm = radec[1][0]
        dec_h = radec[1][1:3]
        dec_m = radec[1][4:6]
        dec_s = radec[1][7:9]
        rastring = ra_h + " " + ra_m + " " + ra_s
        decstring = pm + dec_h + " " + dec_m + " " + dec_s
    else:
        # this is the updated coordinates, see function above

        # these are stars that coords were taken from Gaia DR3 catalog and so epoch is 2016, not the usual 2000
        specialstatus = ['TIC405230669', 'TIC302453291', 'TIC254296853', 'TIC399571743',
                         'TIC408506314', 'TIC408508003', 'TIC36646052', 'TIC441076258',
                         'TOI-5695', 'TOI-6022', 'TOI-6137']
        if row['Starname'][0] in specialstatus:
            ep = '2016'
        else:
            ep = '2000'
        # rastring, decstring = pm_correcter(row['Starname'][0], row['RA'][0], row['Dec'][0], row['pmRA'][0], row['pmDec'][0], ep, verbose=False)
        rastring, decstring = pm_correcter(row['RA'][0], row['Dec'][0], row['Proper Motion in RA [miliarcseconds/year]'][0], row['Proper Motion in Dec [miliarcseconds/year]'][0], ep, current_day, verbose=False)
        if decstring[0] != "-":
            decstring = "+" + decstring

    namestring = ' '*(16-len(row['Starname'][0][:16])) + row['Starname'][0][:16]

    epochstring = '2000'

    jmagstring = ('jmag=' + str(np.round(float(row['J Magnitude'][0]),1)) + ' '*(4-len(str(np.round(row['J Magnitude'][0],1)))))

    exposurestring = (' '*(4-len(str(int(row['Nominal Exposure Time [s]'][0])))) + str(int(row['Nominal Exposure Time [s]'][0])) + '/'
                        + str(int(row['Maximum Exposure Time [s]'][0])) + ' '*(4-len(str(int(row['Maximum Exposure Time [s]'][0])))))

    ofstring = ('1of' + str(int(row['# Visits per Night'][0])))

    if row['Simucal'][0]:
        scval = 'T'
    else:
        scval = 'F'
    scstring = 'sc=' + scval

    numstring = (str(int(row['# of Exposures per Visit'][0])) + "x")

    gmagstring = ('gmag=' + str(np.round(float(row['G Magnitude'][0]),1)) + ' '*(4-len(str(np.round(row['G Magnitude'][0],1)))))

    teffstr = 'Teff=' + str(int(row['Effective Temperature [Kelvin]'][0])) + ' '*(4-len(str(int(row['Effective Temperature [Kelvin]'][0]))))

    if str(row['GAIA Identifier'][0]) != "NoGaiaName":
        # gaiastring = str(row['GAIA Identifier'][0][5:]) + ' '*(25-len(str(row['GAIA Identifier'][0][5:])))
        gaiastring = str(row['GAIA Identifier'][0]) + ' '*(25-len(str(row['GAIA Identifier'][0])))
    else:
        gaiastring = str(row['GAIA Identifier'][0]) + ' '*(25-len(str(row['GAIA Identifier'][0])))

    programstring = row['Program_Code'][0]

    if fillerFlag:
        priostring = "p3" # All targets added in round 2 bonus round are inherently lower priority
    else:
        priostring = "p1"

    if extra == False:
        timestring2 = str(obs_time)
    else:
        timestring2 = "56:78" # choose a nonsense time

    line = (namestring + ' ' + rastring + ' ' + decstring + ' ' + str(epochstring) + ' '
                + jmagstring + ' ' + exposurestring + ' ' + ofstring + ' ' + scstring +  ' '
                + numstring + ' '+ gmagstring + ' ' + teffstr + ' ' + gaiastring + ' CC '
                        + priostring + ' ' + programstring + ' ' + timestring2)

    if not pd.isnull(row['Observing Notes'][0]):
        line += (' ' + str(row['Observing Notes'][0]))

    return line


def pm_correcter(ra, dec, pmra, pmdec, epochstr, current_day, verbose=False):
    """
    Update a star's coordinates due to proper motion
    Note: this function is here only for completeness but the script generation is now fully done by the TTP standalone package

    Args:
        ra (str): the target star's old coordinate RA in units of degrees
        dec (str): the target star's old coordinate Dec in units of degrees
        pmra (str): the target star's proper motion in the RA dimension, in units of millearcseconds per year
        pmdec (str): the target star's proper motion in the Dec dimension, in units of millearcseconds per year
        epochstr (str): a string of the year those coordinates are updated to
        verbose (boolean): true to print out to command line

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

    if verbose:
        print(ra)
        print(ra_deg)
        print(ra_advanced_deg)
        print()
        print(dec)
        print(dec_deg)
        print(dec_advanced_deg)

    # Format the coordinates in the specified format
    formatted_ra = ra_advanced_deg.to_string(unit=u.deg, sep=' ', pad=True, precision=1)
    formatted_dec = dec_advanced_deg.to_string(unit=u.deg, sep=' ', pad=True, precision=0)

    return formatted_ra, formatted_dec
