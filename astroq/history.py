"""
Module for processing the inputs and outputs of the autoscheduler to/from various sources.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import processing_functions as pf
"""
import os
import numpy as np
import pandas as pd
import json
import requests

from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.time import TimeDelta
from collections import defaultdict
from collections import namedtuple

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

def write_OB_histories_to_csv(histories, savepath):
    """
    Write the OB histories to a CSV file.
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
    df.to_csv(savepath, index=False)

def write_OB_histories_to_csv_JUMP(requests_frame, jump_file):
    """
    Convert a JUMP-style CSV to the OB histories format and write to CSV.
    """
    # Load the jump file
    df = pd.read_csv(jump_file)
    
    # Crossmatch star names
    df['target'] = df['star_id'].apply(crossmatch_star_name)
    
    # Prepare a lookup for id and semid from the requests_frame
    req = requests_frame
    req_lookup = req.set_index('starname')[['unique_id', 'program_code']].to_dict(orient='index')

    # Map id and semid using the crossmatched target
    df['id'] = df['target'].map(lambda name: req_lookup.get(name, {}).get('unique_id', ''))
    df['semid'] = df['target'].map(lambda name: req_lookup.get(name, {}).get('program_code', ''))
    
    # Rename columns
    df['timestamp'] = df['utctime']
    df['exposure_start_time'] = df['utctime']
    df['observer'] = 'none'
    
    # Keep only the required columns
    out_df = df[['id', 'target', 'semid', 'timestamp', 'exposure_start_time', 'exposure_time', 'observer']]
    
    # Write to the manager's past_file location
    # out_df.to_csv(manager.upstream_path + manager.past_file, index=False)
    return out_df

def process_star_history(filename):
    """
    Given a CSV file (output of write_OB_histories_to_csv), return a dict:
      keys are 'id', values are objects with attributes:
      date_last_observed, total_n_exposures, total_n_visits, total_n_unique_nights, total_open_shutter_time
    Each key is for one target (star).
    """
    df = pd.read_csv(filename)
    StarHistory = namedtuple('StarHistory', [
        'id',
        'date_last_observed',
        'total_n_exposures',
        'total_n_visits',
        'total_n_unique_nights',
        'total_open_shutter_time',
        'n_obs_on_nights',
    ])
    result = {}

    # Group by target (star)
    for starname, star_df in df.groupby('target'):
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
        date_last_observed = max(all_times).isot[:10] if all_times else ''
        total_n_exposures = sum(len(v['exposure_start_time']) for v in visits)
        total_n_visits = len(visits)
        unique_nights = set(Time(t).isot[:10] for v in visits for t in v['exposure_start_time'])
        total_n_unique_nights = len(unique_nights)
        total_open_shutter_time = int(round(sum(sum(map(float, v['exposure_time'])) for v in visits)))
        # Build n_obs_on_nights
        n_obs_on_nights = {}
        for v in visits:
            visit_date = v['timestamp'][:10]
            n_obs_on_nights[visit_date] = n_obs_on_nights.get(visit_date, 0) + len(v['exposure_start_time'])
        result[starname] = StarHistory(
            id=star_id,
            date_last_observed=date_last_observed,
            total_n_exposures=total_n_exposures,
            total_n_visits=total_n_visits,
            total_n_unique_nights=total_n_unique_nights,
            total_open_shutter_time=total_open_shutter_time,
            n_obs_on_nights=n_obs_on_nights,
        )
    return result

def crossmatch_star_name(name):
    """
    Convert between canonical and CPS star naming conventions.
    If given a CPS name, returns the canonical name.
    If given a canonical name, returns the CPS name.
    """
    # HD <-> CPS
    if name.isdigit():
        return 'HD ' + name
    if name.startswith('HD'):
        return name.replace(" ", "")[2:]
    # KIC <-> CPS
    if name.startswith('KIC') and name[3:].isdigit():
        if " " in name:
            return name.replace(" ", "")
        else:
            return 'KIC ' + name[3:]
    # TYC <-> CPS
    if name.startswith('TYC'):
        if " " in name:
            return name.replace(" ", "-")
        else:
            return name.replace("-", " ")
    # TOI <-> CPS
    if name.startswith('T00') and name[3:].isdigit() and len(name[3:]) == 3:
        return 'TOI-' + name[3:]
    if name.startswith('T0') and name[2:].isdigit() and len(name[2:]) == 4:
        return 'TOI-' + name[2:]
    if name.startswith('TOI-'):
        if len(name) == 7:
            return 'T00' + name[4:]
        elif len(name) == 8:
            return 'T0' + name[4:]
    # Kepler <-> CPS
    if name.startswith('KEPLER'):
        return name.replace("KEPLER", "Kepler")
    if name.startswith('Kepler'):
        return name.replace("Kepler", "KEPLER")
    return name

