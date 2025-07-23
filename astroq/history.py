"""
Module for processing the inputs and outputs of the autoscheduler to/from various sources.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import processing_functions as pf
"""
import os

from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.time import TimeDelta

import numpy as np
import pandas as pd
try:
    JUMP_USERNAME = os.environ['KPFCC_JUMP_USERNAME']
    JUMP_PASSWORD = os.environ['KPFCC_JUMP_PASSWORD']
except:
    print("If pulling past history from Jump, first you must set your environment variables.")
    JUMP_USERNAME = "none"
    JUMP_PASSWORD = "none"

def get_unique_nights(star_past_obs, twilight_frame):
    """
    Parse the Jump database for the previous nights where a given star was observed

    Args:
        star_past_obs (dataframe): a subset of the Jump database dataframe, filtered to only
                                    include observations for a specific target name
        twilight_frame (dataframe): the dataframe of morning and evening twilight times for each
                                    night of the semester

    Returns:
        star_past_obs (dataframe): an updated version of the original variable which now includes
                                   columns for UTC date and hst date
        unique_hst_dates_observed (array): a 1D array length # of unique dates observed where each
                                   element is the HST date when the star was observed at least once
        quarter_observed (array): a 1D array length # of unique dates observed where each element
                                 the quarter of the night when the first of the observations
                                 for that night was taken.
    """
    unique_hst_dates_observed = []
    unique_utc_dates_observed = []
    quarter_observed = []
    all_hst_dates = []
    for i in range(len(star_past_obs)):
        timestamp = star_past_obs['utctime'][i][:-6]
        timestamp = Time(timestamp, format = 'iso')
        timestamp.format = 'jd'
        utc_date = star_past_obs['utctime'][i][:10]
        unique_utc_dates_observed.append(utc_date)

        # Note this is arbitrary 24hr subtraction to get from UT date to Hawaii date
        convert_to_hst_civil = TimeDelta(60*24*60,format='sec').jd
        hsttimestamp = timestamp - convert_to_hst_civil
        hsttimestamp.format = 'iso'
        hstdate = str(hsttimestamp)[:10]
        all_hst_dates.append(hstdate)
        hsttimestamp.format = 'jd'
        if hstdate not in unique_hst_dates_observed:
            unique_hst_dates_observed.append(hstdate)
            quarter_observed.append(get_quarter_observed(timestamp, utc_date, twilight_frame))
    star_past_obs['utc_date'] = unique_utc_dates_observed
    star_past_obs['hstDate'] = all_hst_dates

    return star_past_obs, unique_hst_dates_observed, quarter_observed

def get_quarter_observed(timestamp, utc_date, twilight_frame):
    """
    Determine which quarter of the night an observation was taken in

    Args:
        timestamp (float): a BJD value indicating the time of the observation
        utc_date (str): the calendar date of the observation in UTC, format 'YYYY-MM-DD'
        twilight_frame (dataframe): the dataframe of morning and evening twilight times

    Returns:
        quarter (int): the quarter corresponding to the timestamp
    """
    date_index = twilight_frame.index.get_loc(twilight_frame[twilight_frame['time_utc'] == utc_date].index[0])
    start_jd = twilight_frame['12_evening'][date_index]
    end_jd = twilight_frame['12_morning'][date_index]
    length_of_quarter = (end_jd - start_jd)/4.
    rel_timestamp = float(str(timestamp - start_jd))

    quarter = -10
    if rel_timestamp < length_of_quarter:
        quarter = 0.5
    elif rel_timestamp < 2*length_of_quarter and rel_timestamp > length_of_quarter:
        quarter = 1.5
    elif rel_timestamp < 3*length_of_quarter and rel_timestamp > 2*length_of_quarter:
        quarter = 2.5
    elif rel_timestamp < 4*length_of_quarter and rel_timestamp > 3*length_of_quarter:
        quarter = 3.5
    elif rel_timestamp < 5*length_of_quarter and rel_timestamp > 4*length_of_quarter:
        # allow a little bit of leeway, even if this doesn't really make sense
        quarter = 3.5
    else:
        quarter = 0.5
        print("Houston, we've had a problem: target observed in an invalid quarter.")

    return quarter

def get_nobs_on_night(star_past_obs, unique_hst_dates_observed):
    """
    Determine how many exposures were taken of a target on a given night

    Args:
        star_past_obs (dataframe): the updated version of the original variable which now includes
                                   columns for UTC date and hst date
        unique_hst_dates_observed (array): a 1D array of HST dates where the star was observed

    Returns:
        n_obs_on_date (array): a 1D array of length equal to length unique_hst_dates_observed where
                             each element is an integer of the number of exposures taken that night
    """
    n_obs_on_date = []
    for i in range(len(unique_hst_dates_observed)):
        datemask = star_past_obs['hstDate'] == unique_hst_dates_observed[i]
        n_obs_on_date.append(np.sum(datemask))
    return n_obs_on_date

def real_star_name(name):
    if name.isdigit():
        return 'HD ' + name
    elif name.startswith('KIC') and name[3:].isdigit():
        return 'KIC ' + name[3:]
    elif name.startswith('T00') and name[3:].isdigit() and len(name[3:]) == 3:
        return 'TOI-' + name[3:]
    elif name.startswith('T0') and name[2:].isdigit() and len(name[2:]) == 4:
        return 'TOI-' + name[2:]
    elif name.startswith('KEPLER'):
        name2 = name.replace("KEPLER", "Kepler")
        return name2
    else:
        return name

def cps_star_name(name):
    if name.startswith('HD'):
        name2 = name.replace(" ", "")
        return name2[2:]
    elif name.startswith('KIC'):
        name2 = name.replace(" ", "")
        return name2
    elif name.startswith('TYC'):
        name2 = name.replace(" ", "-")
        return name2
    elif name.startswith('Kepler'):
        name2 = name.replace("Kepler", "KEPLER")
        return name2
    else:
        return name

def build_past_history(past_observations_file, requests_frame, twilight_frame):
    """
    Construct the past history dictionary having pulled a table of observations from Jump.
    Later this function will be replaced by one which pulls history from Keck OB database.

    Args:
        past_observations_file (str): the path and filename of the csv pulled from Jump
        requests_frame (dataframe): the dataframe containing the request information
        twilight_frame (dataframe): the dataframe containing the twilight information

    Returns:
        database_info_dict (dictionary): contains keys of starnames and values are lists where
                    elements follow order described below.
    """
    database_info_dict = {}
    if os.path.exists(past_observations_file):
        print("Pulled database of past observations this semester.")
        database = pd.read_csv(past_observations_file)
        # ensure naming conventions are consistent
        database['star_id_cps'] = database['star_id']
        database['star_id'] = database['star_id_cps'].apply(real_star_name)

        for i in range(len(requests_frame['starname'])):
            starmask = database['star_id'] == requests_frame['starname'][i]
            star_past_obs = database[starmask]
            star_past_obs.sort_values(by='utctime', inplace=True)
            star_past_obs.reset_index(inplace=True)
            total_past_observations = int(len(star_past_obs)/
                                    (requests_frame['n_intra_max'][i]*
                                    requests_frame['n_exp'][i]))
            star_past_obs, unique_hst_dates_observed, quarter_observed = \
                get_unique_nights(star_past_obs, twilight_frame)
            n_obs_on_date = get_nobs_on_night(star_past_obs, unique_hst_dates_observed)

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
            database_info_dict[requests_frame['starname'][i]] = \
                [most_recent_observation_date, unique_hst_dates_observed, quarter_observed, \
                n_obs_on_date]
    else:
        print("No past observation history to parse.")
    return database_info_dict

# NOTE: we are not yet actually using the history field of the OBs, so this is not ready yet but will be used extensively.
# -----------------------------------------------------------------------------------------------------------------------------
# def is_observation_complete(ob, entry):
#     """
#     Determine if a particular entry in the history of the OB is a completed observation. For the
#     accounting and charging of time purposes. Note, this logic is a political decision and can be
#     modified at any time.
#
#     Args:
#         ob (object) - the json object for a single OB. For OB structure info, see here: XYZ.
#         entry (dict) - one element of the history array within the OB
#     Returns:
#         good (boolean) - True if passes unit tests
#     """
#     good = True
#     # If you got at least half of your desired number of exposures, the observation is good
#     if len(entry['exposure times']) < 0.5*ob['observation']['num_exposures']:
#         good = False
#     # If each of your exposures, except the last, was at least 75% complete in time, then your observation is good
#     each_exp = []
#     for i in range(len(entry['exposure times'])):
#         if entry['exposure times'][i] < 0.75*ob['observation']['Exp Time']:
#             each_exp.append(False)
#         else:
#             each_exp.append(True)
#     if np.sum(each_exp) < len(entry['exposure times']) - 1:
#         good = False
#     return good
#
# def get_past_history(ob):
#     """
#     Convert the history field of the OB into the format expected by the autoscheduler.
#
#     Args:
#         ob (object) - the json object for a single OB. For OB structure info, see here: XYZ.
#     Returns:
#         history (list) - See definition in comment below. Intended to be a value in a dictionary
#                          where the keys are the star names.
#     """
#     unique = ob['metadata']['semid'] + "_" + str(ob['target']['target_name'])
#     submitter = ob['metadata']['submitter']
#     # Within the database_info_dict dictionary, each target's data is always in order:
#     # element 0 = the calendar date of the most recent observation (HST date)
#     # element 1 = a list of calendar dates with at least one observation
#     # element 2 = a list of the quarters where the observation took place on the
#     #             corresponding the nights in element 1.
#     #             If multiple visits in one night, then this is quarter of the first visit.
#     # element 3 = a list of the # of observations on each past night,
#     #             corresponding the nights in element 1.
#     observed_exps = []
#     observed_times = []
#     observed_dates = []
#     observed_quarters = []
#     for obs in ob['metadata']['history']:
#         if is_observation_complete(ob, entry):
#             n_exp = len(obs['exposure times'])
#             observed_exps.append(n_exp)
#
#             sumtime = np.sum(obs['exposure times'])
#             observed_times.append(sumtime)
#
#             tstamp = obs['timestamp'][:21]
#             date = obs['timestamp'][:10]
#             observed_dates.append(date)
#
#             quarter = get_quarter_observed(tstamp, date, twilight_frame)
#             observed_quarters.append(quarter)
#         else:
#             print("Observation for " + unique + " on night " + str(obs['timestamp'][:10]) + " submitted by " + submitter + " is not complete.")
#     history = [observed_dates[0], observed_dates, observed_quarters, observed_exps]
#     return history
