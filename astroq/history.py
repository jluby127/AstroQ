"""
Module for processing the past history of observations from OB database or user inputted file. 
Used to populate the past.csv file for the optimization.
"""

# Standard library imports
import os
from collections import namedtuple
from datetime import timedelta

# Third-party imports
import astropy.units as u
import astroplan as apl
import pandas as pd
import requests
from astropy.time import Time, TimeDelta


def _make_observer(observatory_string):
    """Observer for sunset; None if site lookup fails."""
    if not observatory_string or not str(observatory_string).strip():
        return None
    s = str(observatory_string).strip()
    for name in (s, 'mauna kea', 'Keck'):
        try:
            return apl.Observer.at_site(name)
        except Exception:
            continue
    return None


def _local_datetime_from_ut(t_ut, utc_offset_hours):
    """Apply fixed offset (local = UT + utc_offset_hours as implemented elsewhere)."""
    return (t_ut - TimeDelta(abs(utc_offset_hours) / 24, format='jd')).to_value('datetime')


def local_time_str_from_ut(exposure_start_str, utc_offset_hours):
    """Wall-clock local string from UT exposure start (matches process_star_history formatting)."""
    t_ut = Time(exposure_start_str)
    dt = _local_datetime_from_ut(t_ut, utc_offset_hours)
    return dt.strftime('%Y-%m-%d %H:%M:%S')


def _night_of_key_from_cell(night_of_val):
    """YYYY-MM-DD from CSV NightOf (date-only or YYYY-MM-DDTHH:MM...)."""
    if night_of_val is None:
        return None
    try:
        if isinstance(night_of_val, float) and pd.isna(night_of_val):
            return None
    except (TypeError, ValueError):
        pass
    s = str(night_of_val).strip()
    if not s or s.lower() == 'nan':
        return None
    return s[:10] if len(s) >= 10 else None


def night_of_date_from_sunset(t_ut, observer, utc_offset_hours):
    """
    Local civil calendar date (YYYY-MM-DD) of the observing night: the local date of the
    most recent sunset at the observatory before this exposure (in UTC time order).
    Falls back to local-noon night boundary if sunset cannot be computed.
    """
    if observer is not None:
        try:
            sunset = observer.sun_set_time(t_ut, which='previous', horizon=-0.34 * u.deg)
            local_sunset = sunset - TimeDelta(abs(utc_offset_hours) / 24, format='jd')
            return local_sunset.to_value('datetime').date().isoformat()
        except Exception:
            pass
    dt = _local_datetime_from_ut(t_ut, utc_offset_hours)
    if dt.hour < 12:
        return (dt.date() - timedelta(days=1)).isoformat()
    return dt.date().isoformat()


def _exposure_ut_column(df):
    if 'exposure_start_time_UT' in df.columns:
        return 'exposure_start_time_UT'
    if 'exposure_start_time' in df.columns:
        return 'exposure_start_time'
    return None

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
    'exposure_start_times',  # list of "YYYY-MM-DD HH:MM" strings (local time) per exposure
])

def write_OB_histories_to_csv(histories, utc_offset_hours=-10, observatory=None):
    """
    Prepare dataframe of past history for writing to CSV.

    Args:
        histories (json): the OB histories in json format
        utc_offset_hours (float): local minus UT from config [global] UTCoffset.
        observatory (str, optional): site name for astroplan (e.g. config global observatory).

    Returns:
        df (pandas DataFrame): the OB histories in dataframe format
    """
    observer = _make_observer(observatory)
    rows = []
    for entry in histories["history"]:
        for start_time, duration in zip(entry["exposure_start_times"], entry["exposure_times"]):
            t_ut = Time(start_time)
            night_of = night_of_date_from_sunset(t_ut, observer, utc_offset_hours)
            rows.append({
                "id": entry.get("id", ""),
                "target": entry.get("target", ""),
                "semid": entry.get("semid", ""),
                "timestamp": entry.get("timestamp", ""),
                "exposure_start_time_UT": start_time,
                "exposure_start_time_local": local_time_str_from_ut(start_time, utc_offset_hours),
                "exposure_time": duration,
                "NightOf": night_of,
                "observer": entry.get("observer", ""),
            })
    empty_cols = [
        'id', 'target', 'semid', 'timestamp',
        'exposure_start_time_UT', 'exposure_start_time_local', 'exposure_time', 'NightOf', 'observer',
    ]
    df = pd.DataFrame(rows)
    if len(df) == 0:
        return pd.DataFrame(columns=empty_cols)
    df.sort_values(by='timestamp', inplace=True)
    return df
    

def process_star_history(filename, utc_offset_hours=-10, observatory=None):
    """
    Process the past.csv file to return a dict of star histories.

    Args:
        filename (str): the path to the past.csv file
        utc_offset_hours (float): hours difference local minus UT (negative = local behind UT,
            e.g. -10 for HST). Pass SemesterPlanner's value from config [global] UTCoffset.
        observatory (str, optional): site name for NightOf sunset logic (matches prep / config).

    Returns:
        result (dict): a dictionary of star histories where keys are 'id', values are objects with attributes: date_last_observed, total_n_exposures, total_n_visits, total_n_unique_nights, total_open_shutter_time
    """

    empty_cols = [
        'id', 'target', 'semid', 'timestamp',
        'exposure_start_time_UT', 'exposure_start_time_local', 'exposure_time', 'NightOf', 'observer',
    ]
    if not os.path.exists(filename) or os.path.getsize(filename) == 0:
        df = pd.DataFrame(columns=empty_cols)
    else:
        try:
            df = pd.read_csv(filename)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=empty_cols)
    result = {}
    observer = _make_observer(observatory)
    has_night_of_col = 'NightOf' in df.columns

    # Group by target (star)
    for unique_id, star_df in df.groupby('id'):
        ut_col = _exposure_ut_column(star_df)
        if ut_col is None:
            continue
        starname = df[df['id'] == unique_id]['target'].values[0]
        visits = []
        for ts, visit_df in star_df.groupby('timestamp'):
            exposure_ut = []
            exposure_time = []
            night_keys = []
            junk = list(visit_df['junk']) if 'junk' in visit_df.columns else []
            junk = [bool(j) for j in junk]
            if junk and sum(junk) >= len(junk) / 2:
                continue
            for _, row in visit_df.iterrows():
                t = row[ut_col]
                exposure_ut.append(t)
                exposure_time.append(row['exposure_time'])
                if has_night_of_col and pd.notna(row.get('NightOf', None)):
                    nk = _night_of_key_from_cell(row['NightOf'])
                else:
                    nk = None
                if not nk:
                    nk = night_of_date_from_sunset(Time(t), observer, utc_offset_hours)
                night_keys.append(nk)
            visits.append({
                'id': visit_df['id'].iloc[0],
                'exposure_start_time': exposure_ut,
                'exposure_time': exposure_time,
                'timestamp': ts,
                'night_keys': night_keys,
            })
        if not visits:
            continue
        all_night_keys = [k for v in visits for k in v['night_keys'] if k]
        date_last_observed = max(all_night_keys) if all_night_keys else ''
        total_n_exposures = sum(len(v['exposure_start_time']) for v in visits)
        total_n_visits = len(visits)
        unique_nights = set(all_night_keys)
        total_n_unique_nights = len(unique_nights)
        total_open_shutter_time = int(round(sum(sum(map(float, v['exposure_time'])) for v in visits)))
        n_obs_on_nights = {}
        n_visits_on_nights = {}
        for v in visits:
            for k in v['night_keys']:
                if k:
                    n_obs_on_nights[k] = n_obs_on_nights.get(k, 0) + 1
            for k in set(x for x in v['night_keys'] if x):
                n_visits_on_nights[k] = n_visits_on_nights.get(k, 0) + 1
        exposure_start_times = []
        for v in visits:
            for t in v['exposure_start_time']:
                ut_time = Time(t)
                local_time = ut_time - TimeDelta(abs(utc_offset_hours) / 24, format='jd')
                exposure_start_times.append(local_time.to_value('datetime').strftime('%Y-%m-%d %H:%M'))
        result[unique_id] = StarHistory(
            name=starname,
            date_last_observed=date_last_observed,
            total_n_exposures=total_n_exposures,
            total_n_visits=total_n_visits,
            total_n_unique_nights=total_n_unique_nights,
            total_open_shutter_time=total_open_shutter_time,
            n_obs_on_nights=n_obs_on_nights,
            n_visits_on_nights=n_visits_on_nights,
            exposure_start_times=exposure_start_times
        )
    return result

