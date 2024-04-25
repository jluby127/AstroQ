from collections import defaultdict
from gurobipy import *
import numpy as np
import astropy.units as u
import pandas as pd
import gurobipy as gp
import math
import astropy as apy
import astroplan as apl
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
