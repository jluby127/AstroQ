"""
Module for processing the past history of observations from OB database or user inputted file. 
Used to populate the past.csv file for the optimization.
"""

# Standard library imports
import os
from collections import namedtuple

# Third-party imports
import pandas as pd
import requests
from astropy.time import Time, TimeDelta

# Global timezone offset: hours difference between local time and UT
# Positive values mean local time is ahead of UT (e.g., +8 for UTC+8)
# Negative values mean local time is behind UT (e.g., -10 for UTC-10)
TIMEZONE_OFFSET_HOURS = -10  # Adjust this value for your timezone

# Define StarHistory namedtuple at module level for proper pickling
StarHistory = namedtuple('StarHistory', [
    'name',
    'date_last_observed',
    'total_n_exposures',
    'total_n_visits',
    'total_n_unique_nights',
    'total_open_shutter_time',
    'n_obs_on_nights',
    'n_visits_on_nights',
])

def pull_OB_histories(semester):
    """
    Pull the latest database OBs down to local.

    Args:
        semester (str) - the semester from which to query OBs, format YYYYL
        histories (bool) - if True, pull the history of OBs for the semester, if False, pull the latest OBs for the semester

    Returns:
        data (json) - the OB information in json format
    """
    url = "https://www3.keck.hawaii.edu/api/kpfcc/getObservingBlockHistory"
    params = {}
    params["semester"] = semester
    try:
        data = requests.get(url, params=params, auth=(os.environ['KECK_OB_DATABASE_API_USERNAME'], os.environ['KECK_OB_DATABASE_API_PASSWORD']))
        data = data.json()
        return data
    except:
        print("ERROR")
        return

def write_OB_histories_to_csv(histories):
    """
    Prepare dataframe of past history for writing to CSV.

    Args:
        histories (json): the OB histories in json format

    Returns:
        df (pandas DataFrame): the OB histories in dataframe format
    """
    rows = []
    for entry in histories["history"]:
        # Zip together times and durations for each exposure in this record
        for start_time, duration in zip(entry["exposure_start_times"], entry["exposure_times"]):
            rows.append({
                "id": entry.get("id", ""),
                "target": entry.get("target", ""),
                "semid": entry.get("semid", ""),
                "timestamp": entry.get("timestamp", ""),
                "exposure_start_time": start_time,
                "exposure_time": duration,
                "observer": entry.get("observer", ""),
                # "junk": entry.get("junk", ""),
            })
    df = pd.DataFrame(rows)
    df.sort_values(by='timestamp', inplace=True)
    return df

def write_OB_histories_to_csv_JUMP(requests_frame, jump_file):
    """
    Convert a JUMP-style CSV to the OB histories format and prepare dataframe for writing to CSV.

    Args:
        requests_frame (pandas DataFrame): the requests frame in dataframe format
        jump_file (str): the path to the JUMP-style CSV file

    Returns:
        df (pandas DataFrame): the OB histories in dataframe format
    """
    # Check if file is empty or doesn't exist
    if not os.path.exists(jump_file) or os.path.getsize(jump_file) == 0:
        # Create empty DataFrame with expected columns
        empty_df = pd.DataFrame(columns=['id', 'target', 'semid', 'timestamp', 'exposure_start_time', 'exposure_time', 'observer'])
        empty_df.to_csv(output_file, index=False)
        return empty_df
    
    try:
        # Load the jump file
        df = pd.read_csv(jump_file)
    except pd.errors.EmptyDataError:
        # Create empty DataFrame with expected columns
        empty_df = pd.DataFrame(columns=['id', 'target', 'semid', 'timestamp', 'exposure_start_time', 'exposure_time', 'observer'])
        empty_df.to_csv(output_file, index=False)
        return empty_df
    
    # Crossmatch star names
    df['target'] = df['star_id'].apply(crossmatch_star_name)
    
    # Prepare a lookup for id and semid from the requests_frame
    req = requests_frame
    # Create lookup dictionary manually to handle duplicate starnames
    req_lookup = {}
    for _, row in req.iterrows():
        starname = row['starname']
        unique_id = row['unique_id']
        if unique_id not in req_lookup:
            req_lookup[unique_id] = {'starname': row['starname'], 'program_code': row['program_code']}

    # Map id and semid using the crossmatched target
    df['unique_id'] = df['target'].map(lambda name: req_lookup.get(name, {}).get('unique_id', ''))
    df['semid'] = df['target'].map(lambda name: req_lookup.get(name, {}).get('program_code', ''))
    
    # Filter out rows where starname doesn't have a corresponding entry in lookup table
    df = df[df['id'] != '']
    
    # Convert datetime format from space-separated to ISO format with 'T'
    df['timestamp'] = df['utctime'].str.replace(' ', 'T')
    df['exposure_start_time'] = df['utctime'].str.replace(' ', 'T')
    df['observer'] = 'none'
    
    # Keep only the required columns
    out_df = df[['id', 'target', 'semid', 'timestamp', 'exposure_start_time', 'exposure_time', 'observer']]
    return out_df

def process_star_history(filename):
    """
    Process the past.csv file to return a dict of star histories.

    Args:
        filename (str): the path to the past.csv file

    Returns:
        result (dict): a dictionary of star histories where keys are 'id', values are objects with attributes: date_last_observed, total_n_exposures, total_n_visits, total_n_unique_nights, total_open_shutter_time
    """

    # Check if file is empty or doesn't exist
    if not os.path.exists(filename) or os.path.getsize(filename) == 0:
        # Create empty DataFrame with expected columns
        df = pd.DataFrame(columns=['id', 'target', 'semid', 'timestamp', 'exposure_start_time', 'exposure_time', 'observer'])
    else:
        try:
            df = pd.read_csv(filename)
        except pd.errors.EmptyDataError:
            # Create empty DataFrame with expected columns
            df = pd.DataFrame(columns=['id', 'target', 'semid', 'timestamp', 'exposure_start_time', 'exposure_time', 'observer'])
    result = {}

    # Group by target (star)
    for unique_id, star_df in df.groupby('id'):
        starname = df[df['id'] == unique_id]['target'].values[0]
        # Group by visit (timestamp)
        visits = []
        for ts, visit_df in star_df.groupby('timestamp'):
            # Each visit is a DataFrame of exposures
            exposure_start_time = list(visit_df['exposure_start_time'])
            exposure_time = list(visit_df['exposure_time'])
            junk = list(visit_df['junk']) if 'junk' in visit_df.columns else []
            # Convert junk to bool if needed
            junk = [bool(j) for j in junk]
            # Apply junk logic: skip if >= half are True
            if junk and sum(junk) >= len(junk) / 2:
                continue
            visits.append({
                'id': visit_df['id'].iloc[0],
                'exposure_start_time': exposure_start_time,
                'exposure_time': exposure_time,
                'timestamp': ts,
            })
        if not visits:
            continue
        star_id = visits[0]['id']
        all_times = [Time(t) for v in visits for t in v['exposure_start_time']]
        # Adjust UT time to local timezone
        ut_time = max(all_times) if all_times else None
        if ut_time:
            # Subtract timezone offset (negative offset means local time is behind UT)
            local_time = ut_time - TimeDelta(abs(TIMEZONE_OFFSET_HOURS) / 24, format='jd')
            date_last_observed = local_time.isot[:10]
        else:
            date_last_observed = ''
        total_n_exposures = sum(len(v['exposure_start_time']) for v in visits)
        total_n_visits = len(visits)
        unique_nights = set(Time(t).isot[:10] for v in visits for t in v['exposure_start_time'])
        total_n_unique_nights = len(unique_nights)
        total_open_shutter_time = int(round(sum(sum(map(float, v['exposure_time'])) for v in visits)))
        # Build n_obs_on_nights and n_visits_on_nights
        n_obs_on_nights = {}
        n_visits_on_nights = {}
        for v in visits:
            visit_date = v['timestamp'][:10]
            n_obs_on_nights[visit_date] = n_obs_on_nights.get(visit_date, 0) + len(v['exposure_start_time'])
            n_visits_on_nights[visit_date] = n_visits_on_nights.get(visit_date, 0) + 1
        result[unique_id] = StarHistory(
            name=starname,
            date_last_observed=date_last_observed,
            total_n_exposures=total_n_exposures,
            total_n_visits=total_n_visits,
            total_n_unique_nights=total_n_unique_nights,
            total_open_shutter_time=total_open_shutter_time,
            n_obs_on_nights=n_obs_on_nights,
            n_visits_on_nights=n_visits_on_nights
        )
    return result

def crossmatch_star_name(name):
    """
    Convert between canonical and CPS star naming conventions.
    If given a CPS name, returns the canonical name.
    If given a canonical name, returns the CPS name.

    Args:
        name (str): the star name to convert

    Returns:
        name (str): the converted star name
    """
    if name.startswith('KEPLER'):
        return name.replace("KEPLER", "Kepler")
    if name.startswith('Kepler'):
        return name.replace("Kepler", "KEPLER")
    return name

