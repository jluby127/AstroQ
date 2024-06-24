from collections import defaultdict
from gurobipy import *
import numpy as np
import astropy.units as u
import pandas as pd
import gurobipy as gp
import math
import astropy as apy
import astroplan as apl
from astropy.coordinates import Angle
from astropy.time import Time
from astropy.time import TimeDelta
import requests
from bs4 import BeautifulSoup


def login_JUMP():
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


def getKPFAllObservations(path_to_csv, flag=False):
    name = 'Get all KPF Observations in Semester'
    path_to_csv = '/Users/jack/Desktop/Jump_allObs2024A.csv'
    # comment this line out when playing with synthetic schedules
    get_database_explorer(name, path_to_csv)
    print("All KPF observations pulled from Jump. Saved to csv: " + path_to_csv)


def getUniqueNights(star_past_obs, twilight):

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

        convertHSTcivil = TimeDelta(60*24*60,format='sec').jd #Note this is arbitrary 24hr subtraction
        hsttimestamp = timestamp - convertHSTcivil
        hsttimestamp.format = 'iso'
        hstdate = str(hsttimestamp)[:10]
        all_hstdates.append(hstdate)
        hsttimestamp.format = 'jd'
        if hstdate not in unique_hstdates_observed:
            unique_hstdates_observed.append(hstdate)
            quarterObserved.append(getQuartersObserved(timestamp, utcdate, twilight))
            #Nobs_on_date.append(1)

    star_past_obs['utcDate'] = unique_utcdates_observed
    star_past_obs['hstDate'] = all_hstdates

    return star_past_obs, unique_hstdates_observed, quarterObserved


def getQuartersObserved(timestamp, utcdate, twilight):

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
    elif rel_timestamp < 4.5*lengthOfQuarter and rel_timestamp > 4*lengthOfQuarter:
        # allow a little bit of leeway
        quarter = 3.5
    else:
        print("Houston, we've had a problem: target observed in an invalid quarter.")
    return quarter


def getNobs_on_Night(star_past_obs, unique_hstdates_observed):
    Nobs_on_date = []
    for i in range(len(unique_hstdates_observed)):
        datemask = star_past_obs['hstDate'] == unique_hstdates_observed[i]
        Nobs_on_date.append(np.sum(datemask))
    return Nobs_on_date


def prepareTTP(request_sheet, night_plan):
    # night_plan is the n'th row of combined_semester_schedule array representing night n. Will include * and maybe W, must be processed
    # Sometimes a W and a starname are in the same slot. Must filter out.
    ignore = ['*', 'W', '']
    selected_stars = []
    for i in range(len(night_plan)):
        if night_plan[i] not in ignore:
            selected_stars.append(night_plan[i])
            ignore.append(night_plan[i])
            # try:
            #     if night_plan[i][0] != 'W' and night_plan[i][:4] != 'RM__':
            #         selected_stars.append(night_plan[i])
            #         ignore.append(night_plan[i])
            # except:
            #     selected_stars.append(night_plan[i])
            #     ignore.append(night_plan[i])

    # build the dataframe with correct info and correct headers
    all_targets_frame = pd.read_csv(request_sheet)
    RAs = []
    Decs = []
    ExpTimes = []
    ExpsPerVisits = []
    nVisits = []
    cadences = []
    for j in range(len(selected_stars)):
        #print(selected_stars[j])
        #print(all_targets_frame.index[all_targets_frame['Starname']==str(selected_stars[j])])
        #print(all_targets_frame.index[all_targets_frame['Starname']==str(selected_stars[j])][0])
        idx = all_targets_frame.index[all_targets_frame['Starname']==str(selected_stars[j])][0]
        RAs.append(all_targets_frame['RA'][idx])
        Decs.append(all_targets_frame['Dec'][idx])
        ExpTimes.append(all_targets_frame['Nominal_ExpTime'][idx])
        ExpsPerVisits.append(all_targets_frame['N_Observations_per_Visit'][idx])
        nVisits.append(all_targets_frame['N_Visits_per_Night'][idx])
        cadences.append(all_targets_frame['Intra_Night_Cadence'][idx])
    toTTP = pd.DataFrame({"Starname":selected_stars,"RA":RAs,"Dec":Decs,
                          "Exposure Time":ExpTimes,"Exposures Per Visit":ExpsPerVisits,
                          "Visits In Night":nVisits,"Intra_Night_Cadence":cadences})
    return toTTP


def write_starlist(frame,orderedList,condition,current_day,outputdir):
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    script_file = os.path.join(outputdir,'script_{}_{}.txt'.format(current_day,condition))
    print('Writing starlist to ' + script_file)

    lines = []
    for i in range(len(orderedList)):
        #row = frame.iloc[i]
        # print(i, orderedList['Target'][i])
        row = frame.loc[frame['Starname'] == orderedList['Target'][i]]
        row.reset_index(inplace=True)
        # print(row)
        # print(row['Starname'])
        # print(row['Starname'][0])
        lines.append(format_kpf_row(row, orderedList['StartExposure'][i]))

    lines.append('')
    lines.append('X' * 45 + 'EXTRAS' + 'X' * 45)
    lines.append('')

    #Formatted starlists appear as text files in the directory
    with open(script_file, 'w') as f:
        f.write('\n'.join(lines))


def format_kpf_row(row, obs_time, Jtwothousand = False, extra=False):
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
        rastring, decstring = pm_correcter(row['Starname'][0], row['RA'][0], row['Dec'][0], row['pmRA'][0], row['pmDec'][0], ep, verbose=False)
        if decstring[0] != "-":
            decstring = "+" + decstring

    namestring = ' '*(16-len(row['Starname'][0][:16])) + row['Starname'][0][:16]

    epochstring = '2000'

    jmagstring = ('jmag=' + str(np.round(float(row['Jmag'][0]),1)) + ' '*(4-len(str(np.round(row['Jmag'][0],1)))))

    exposurestring = (' '*(4-len(str(int(row['Nominal_ExpTime'][0])))) + str(int(row['Nominal_ExpTime'][0])) + '/'
                        + str(int(row['Max_ExpTime'][0])) + ' '*(4-len(str(int(row['Max_ExpTime'][0])))))

    ofstring = ('1of' + str(int(row['N_Visits_per_Night'][0])))

    if row['Simulcal'][0]:
        scval = 'T'
    else:
        scval = 'F'
    scstring = 'sc=' + scval

    numstring = (str(int(row['N_Observations_per_Visit'][0])) + "x")

    gmagstring = ('gmag=' + str(np.round(float(row['Gmag'][0]),1)) + ' '*(4-len(str(np.round(row['Gmag'][0],1)))))

    teffstr = 'Teff=' + str(int(row['Teff'][0])) + ' '*(4-len(str(int(row['Teff'][0]))))

    if str(row['UpdatedGaia'][0]) != "NoGaiaName":
        gaiastring = str(row['UpdatedGaia'][0][5:]) + ' '*(25-len(str(row['UpdatedGaia'][0][5:])))
    else:
        gaiastring = str(row['UpdatedGaia'][0]) + ' '*(25-len(str(row['UpdatedGaia'][0])))

    programstring = row['Program_Code'][0]

    priostring = "p1" # For now, no priorities.

    if extra == False:
        timestring2 = str(obs_time)#[11:-7]
    else:
        timestring2 = "56:78" # choose a nonsense time

    line = (namestring + ' ' + rastring + ' ' + decstring + ' ' + str(epochstring) + ' '
                + jmagstring + ' ' + exposurestring + ' ' + ofstring + ' ' + scstring +  ' '
                + numstring + ' '+ gmagstring + ' ' + teffstr + ' ' + gaiastring + ' CC '
                        + priostring + ' ' + programstring + ' ' + timestring2)

    if not pd.isnull(row['Comment'][0]):
        line += (' ' + str(row['Comment'][0]))

    return line


def pm_correcter(star_in, ra, dec, pmra, pmdec, epochstr, verbose=False):
    # requires giving RA and Dec in degrees
    # example: RA = 321.5 and Dec = 15.6
    # note that degrees are not hour angles!
    # this code converts RA from degrees to hourangle at the end

    current_time = Time('2024-01-01')  # You can adjust the date as needed

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
