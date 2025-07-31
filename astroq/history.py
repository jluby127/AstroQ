"""
Module for processing the inputs and outputs of the autoscheduler to/from various sources.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import processing_functions as pf
"""

# Standard library imports
import os
from collections import namedtuple

# Third-party imports
import pandas as pd
import requests


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

def write_OB_histories_to_csv_JUMP(requests_frame, jump_file, output_file):
    """
    Convert a JUMP-style CSV to the OB histories format and write to CSV.
    """
    # Load the jump file
    df = pd.read_csv(jump_file)
    
    # Crossmatch star names
    df['target'] = df['star_id'].apply(crossmatch_star_name)
    
    # Prepare a lookup for id and semid from the requests_frame
    req = requests_frame
    # Create lookup dictionary manually to handle duplicate starnames
    req_lookup = {}
    for _, row in req.iterrows():
        starname = row['starname']
        if starname not in req_lookup:
            req_lookup[starname] = {'unique_id': row['unique_id'], 'program_code': row['program_code']}

    # Map id and semid using the crossmatched target
    df['id'] = df['target'].map(lambda name: req_lookup.get(name, {}).get('unique_id', ''))
    df['semid'] = df['target'].map(lambda name: req_lookup.get(name, {}).get('program_code', ''))
    
    # Filter out rows where starname doesn't have a corresponding entry in lookup table
    df = df[df['id'] != '']
    
    # Convert datetime format from space-separated to ISO format with 'T'
    df['timestamp'] = df['utctime'].str.replace(' ', 'T')
    df['exposure_start_time'] = df['utctime'].str.replace(' ', 'T')
    df['observer'] = 'none'
    
    # Keep only the required columns
    out_df = df[['id', 'target', 'semid', 'timestamp', 'exposure_start_time', 'exposure_time', 'observer']]
    
    # Write to the manager's past_file location
    out_df.to_csv(output_file, index=False)
    return out_df

def process_star_history_dict(history_dict, output_csv):
    """
    Given a dict with a 'history' key (list of visit dicts), write a CSV with columns:
      id, date_last_observed, total_n_exposures, total_n_visits, total_n_unique_nights, total_open_shutter_time
    Each row is for one target (star).
    """
    history = history_dict['history']
    # Group all visits by target
    from collections import defaultdict
    star_visits = defaultdict(list)
    for visit in history:
        starname = visit.get('target', 'UNKNOWN')
        star_visits[starname].append(visit)

    rows = []
    for starname, visits in star_visits.items():
        # Filter out visits where >= half of 'junk' is True
        filtered_visits = []
        for v in visits:
            junk = v.get('junk', None)
            if junk is not None and len(junk) > 0:
                if sum(junk) >= len(junk) / 2:
                    continue  # skip this visit
            filtered_visits.append(v)
        if not filtered_visits:
            continue
        star_id = filtered_visits[0].get('id', starname)
        all_times = [Time(t) for v in filtered_visits for t in v.get('exposure_start_times', [])]
        date_last_observed = max(all_times).isot[:10] if all_times else ''
        total_n_exposures = sum(len(v.get('exposure_start_times', [])) for v in filtered_visits)
        total_n_visits = len(filtered_visits)
        unique_nights = set(Time(t).isot[:10] for v in filtered_visits for t in v.get('exposure_start_times', []))
        total_n_unique_nights = len(unique_nights)
        total_open_shutter_time = int(round(sum(sum(v.get('exposure_times', [])) for v in filtered_visits)))
        rows.append({
            'id': star_id,
            'date_last_observed': date_last_observed,
            'total_n_exposures': total_n_exposures,
            'total_n_visits': total_n_visits,
            'total_n_unique_nights': total_n_unique_nights,
            'total_open_shutter_time': total_open_shutter_time
        })
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    return df

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
    if name.startswith('KEPLER'):
        return name.replace("KEPLER", "Kepler")
    if name.startswith('Kepler'):
        return name.replace("Kepler", "KEPLER")
    return name

