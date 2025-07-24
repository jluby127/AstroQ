"""
Module for processing the inputs and outputs of the autoscheduler to/from various sources.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import processing_functions as pf
"""
import os
import numpy as np
import pandas as pd

from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.time import TimeDelta

def process_star_history(history_dict, output_csv):
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

