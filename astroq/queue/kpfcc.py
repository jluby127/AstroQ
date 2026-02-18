"""
Module for preparing all data specific to from Keck Observatory's custom made Observing Block (OB) database.
This is specific to the KPF-CC program and the observatory's infrastructure as way to power the prep kpfcc command.
New observatories should write their own module to connect to a new "prep <your observatory>" command.
"""

# Standard library imports
import json
import logging
import os

# Third-party imports
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
import astropy.units as u
import astroplan as apl
import astropy.coordinates as apy

# Local imports
from astroq.access import Access

logs = logging.getLogger(__name__)

# Column definitions: mapping from original names to new names and data types
column_definitions = {
    '_id': {'new_name': 'unique_id', 'type': 'string'},
    'metadata.semid': {'new_name': 'program_code', 'type': 'string'},
    'target.target_name': {'new_name': 'starname', 'type': 'string'},
    'target.ra': {'new_name': 'ra', 'type': 'string'},
    'target.dec': {'new_name': 'dec', 'type': 'string'},
    'observation.exposure_time': {'new_name': 'exptime', 'type': 'Int64'},
    'observation.num_exposures': {'new_name': 'n_exp', 'type': 'Int64'},
    'schedule.num_nights_per_semester': {'new_name': 'n_inter_max', 'type': 'Int64'},
    'schedule.num_internight_cadence': {'new_name': 'tau_inter', 'type': 'Int64'},
    'schedule.desired_num_visits_per_night': {'new_name': 'n_intra_max', 'type': 'Int64'},
    'schedule.minimum_num_visits_per_night': {'new_name': 'n_intra_min', 'type': 'Int64'},
    'schedule.num_intranight_cadence': {'new_name': 'tau_intra', 'type': 'Float64'},
    'schedule.minimum_elevation': {'new_name': 'minimum_elevation', 'type': 'Float64'},
    'schedule.minimum_moon_separation': {'new_name': 'minimum_moon_separation', 'type': 'Float64'},
    'schedule.weather_band_1': {'new_name': 'weather_band_1', 'type': 'boolean'},
    'schedule.weather_band_2': {'new_name': 'weather_band_2', 'type': 'boolean'},
    'schedule.weather_band_3': {'new_name': 'weather_band_3', 'type': 'boolean'},
    'target.gaia_id': {'new_name': 'gaia_id', 'type': 'string'},
    'target.t_eff': {'new_name': 'teff', 'type': 'Float64'},
    'target.j_mag': {'new_name': 'jmag', 'type': 'Float64'},
    'target.g_mag': {'new_name': 'gmag', 'type': 'Float64'},
    'target.pm_ra': {'new_name': 'pmra', 'type': 'Float64'},
    'target.pm_dec': {'new_name': 'pmdec', 'type': 'Float64'},
    'target.epoch': {'new_name': 'epoch', 'type': 'Float64'},
    'observation.exp_meter_threshold': {'new_name': 'exp_meter_threshold', 'type': 'Float64'},
    'metadata.ob_inactive': {'new_name': 'inactive', 'type': 'boolean'},
}

# Required fields for OBs to be considered valid
# All fields listed here must be present in the OB for it to pass validation
required_fields = list(column_definitions.keys())

def pull_OBs(semester):
    """
    Pull the latest info from Keck Observatory's KPF-CC database OBs down to local machine.
    Note you must set environment variables KECK_OB_DATABASE_API_USERNAME and KECK_OB_DATABASE_API_PASSWORD to your credentials.

    Args:
        semester (str) - the semester from which to query OBs, format YYYYL

    Returns:
        data (json) - the OB information in json format
    """
    url = "https://www3.keck.hawaii.edu/api/kpfcc/getAllSemesterObservingBlocks"
    params = {}
    params["semester"] = semester
    try:
        data = requests.get(url, params=params, auth=(os.environ['KECK_OB_DATABASE_API_USERNAME'], os.environ['KECK_OB_DATABASE_API_PASSWORD']))
        data = data.json()
        return data
    except:
        print("ERROR")
        return

def _validate_datetime_format(datetime_str):
    """
    Validate that datetime string follows YYYY-MM-DDTHH:MM or YYYY-MM-DDTHH:MM:SS format.
    
    Args:
        datetime_str (str): The datetime string to validate
        
    Returns:
        bool: True if format is valid, False otherwise
    """
    import re
    if not isinstance(datetime_str, str):
        return False
    # Accept YYYY-MM-DDTHH:MM or YYYY-MM-DDTHH:MM:SS (Keck API returns the latter)
    pattern = r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}(:\d{2})?$'
    return bool(re.match(pattern, datetime_str))

def format_custom_csv(OBs):
    """
    Format the custom.csv file from the OBs.

    Args:
        OBs (json): the OB information in json format

    Returns:
        custom_frame (pandas DataFrame): a DataFrame with the custom information, equivalent to the custom.csv file.
    """
    rows = []
    for ob in OBs['observing_blocks']:
        if "custom_time_constraints" in ob.get('schedule', {}):
            ctc = ob['schedule']['custom_time_constraints']
            unique_id = ob['_id']
            starname = ob['target']['target_name']

            if isinstance(ctc, list) and len(ctc) > 0:
                for constraint in ctc:
                    if isinstance(constraint, dict):
                        start = constraint.get('start_datetime', '')
                        stop = constraint.get('end_datetime', '')
                        if start and stop and _validate_datetime_format(start) and _validate_datetime_format(stop):
                            rows.append({
                                'unique_id': unique_id,
                                'starname': starname,
                                'start': start,
                                'stop': stop
                            })
            elif isinstance(ctc, dict):
                start = ctc.get('start_datetime', '')
                stop = ctc.get('end_datetime', '')
                if start and stop and _validate_datetime_format(start) and _validate_datetime_format(stop):
                    rows.append({
                        'unique_id': unique_id,
                        'starname': starname,
                        'start': start,
                        'stop': stop
                    })

    if len(rows) > 0:
        custom_frame = pd.DataFrame(rows)
    else:
        custom_frame = pd.DataFrame(columns=['unique_id', 'starname', 'start', 'stop'])
    
    return custom_frame 
        
def pull_allocation_info(start_date, numdays, instrument, conversion_ratio=12.0):
    """
    Pull the allocation information directly from the Keck Observatory's operations schedule via the API.

    Args:
        start_date (str): the start date of the allocation (day one of the semester)
        numdays (int): the number of days beyond the start_date to pull allocation information (usually ~180)
        instrument (str): the instrument to pull allocation data, here it is "KPF-CC"
        conversion_ratio (float): factor to convert nights to hours (e.g. hours per night). Default 12.0.

    Returns:
        allocation_frame (pandas DataFrame): a DataFrame with the allocation information, equivalent to the allocation.csv file.
        hours_by_program (dict): a dictionary mapping the program code to the total hours allocated to that program
        nights_by_program (dict): a dictionary mapping the program code to the total nights allocated to that program
    """
    params = {}
    params['cmd'] = "getSchedule"
    params["date"] = start_date
    params["numdays"] = numdays
    params["instrument"] = instrument 
    url = "https://www3.keck.hawaii.edu/api/schedule/getSchedule"
    try:
        data = requests.get(url, params=params)
        data_json = json.loads(data.text)
        df = pd.DataFrame(data_json)
        awarded_programs = df['ProjCode'].unique()
        df['start'] = pd.to_datetime(df['Date'] + ' ' + df['StartTime']).dt.strftime('%Y-%m-%dT%H:%M')
        df['stop']  = pd.to_datetime(df['Date'] + ' ' + df['EndTime']).dt.strftime('%Y-%m-%dT%H:%M')

        allocation_frame = df[['start', 'stop']].copy() # TODO: add observer and comment
        
        nights_by_program = df.groupby('ProjCode')['FractionOfNight'].sum().round(3).to_dict()
        hours_by_program = {k: round(v * conversion_ratio, 3) for k, v in nights_by_program.items()}
    except:
        print("ERROR: allocation information not found. Double check date and instrument. Saving an empty file.")
        allocation_frame = pd.DataFrame(columns=['start', 'stop'])
        awarded_programs = []
        hours_by_program = {}
        nights_by_program = {}
    return allocation_frame, hours_by_program, nights_by_program

def format_keck_allocation_info(allocation_file):
    """
    An alternate way to produce the allocation.csv file. Read in a Keck operations schedule file.
    
    Args:
        allocation_file (str): the path and filename to the downloaded csv
        
    Returns:
        allocation_frame (pandas DataFrame): a DataFrame with the allocation information, equivalent to the allocation.csv file.
        hours_by_program (dict): a dictionary mapping the program code to the total hours allocated to that program
        nights_by_program (dict): a dictionary mapping the program code to the total nights allocated to that program
    """
    allocation = pd.read_csv(allocation_file)
        
    # Convert start and stop times to datetime for hour calculation
    allocation['start'] = pd.to_datetime(allocation['Date'] + ' ' + allocation['StartTime']).dt.strftime('%Y-%m-%dT%H:%M')
    allocation['stop'] = pd.to_datetime(allocation['Date'] + ' ' + allocation['EndTime']).dt.strftime('%Y-%m-%dT%H:%M')

    # Calculate hours for each row
    start_times = pd.to_datetime(allocation['start'])
    stop_times = pd.to_datetime(allocation['stop'])
    hours_per_row = (stop_times - start_times).dt.total_seconds() / 3600
    
    # Add hours to the dataframe
    allocation['hours'] = hours_per_row
    
    # Create allocation frame with start/stop times for saving
    allocation_frame = allocation[['start', 'stop']].copy()
    
    # Calculate total hours per ProjCode
    hours_by_program = allocation.groupby('ProjCode')['hours'].sum().round(3).to_dict()
    nights_by_program = allocation.groupby('ProjCode')['FractionOfNight'].sum().round(3).to_dict()

    return allocation_frame, hours_by_program, nights_by_program

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

def get_request_sheet(OBs, awarded_programs, savepath):
    """
    Produce the request.csv file from the json OBs.

    Args:
        OBs (json): the OB information in json format
        awarded_programs (list): a list of the awarded programs
        savepath (str): the path and filename where to save the request sheet

    Returns:
        good_obs (pandas DataFrame): a DataFrame with the OBs that pass the checks
        bad_obs_values (pandas DataFrame): a DataFrame with the values of the bad OBs fields
        bad_obs_hasFields (pandas DataFrame): a DataFrame with the indication of fields existing or not for thebad OBs
        bad_obs_count_by_semid (pandas DataFrame): a DataFrame with the count of bad OBs by semester, for admin plotting purposes
        bad_field_histogram (pandas DataFrame): a DataFrame with the histogram of bad OBs by field, for admin plotting purposes
    """
    good_obs, bad_obs_values, bad_obs_hasFields = sort_good_bad(OBs, awarded_programs)

    # Filter bad OBs to only those in awarded programs
    if 'metadata.semid' in bad_obs_values.columns:
        mask = bad_obs_values['metadata.semid'].isin(awarded_programs)
        bad_obs_values = bad_obs_values[mask].reset_index(drop=True)
        bad_obs_hasFields = bad_obs_hasFields[mask].reset_index(drop=True)

    bad_obs_count_by_semid, bad_field_histogram = analyze_bad_obs(good_obs, bad_obs_values, bad_obs_hasFields, awarded_programs)
    good_obs.sort_values(by='program_code', inplace=True)
    good_obs.reset_index(inplace=True, drop=True)

    # good_obs['active'] = [True] * len(good_obs)
    
    # Cast starname column to strings to ensure proper matching
    if 'starname' in good_obs.columns:
        good_obs['starname'] = good_obs['starname'].astype(str)
    
    os.makedirs(os.path.dirname(savepath), exist_ok=True)
    return good_obs, bad_obs_values, bad_obs_hasFields, bad_obs_count_by_semid, bad_field_histogram

def flatten(d, parent_key='', sep='.'):
    """
    Flatten a dictionary into a single level.

    Args:
        d (dict): the nested dictionary to flatten
        parent_key (str): the parent key
        sep (str): the separator between the parent key and the child key

    Returns:
        items (dict): a dictionary with the flattened keys and values
    """
    items = {}
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.update(flatten(v, new_key, sep=sep))
        else:
            items[new_key] = v
    return items

def apply_safety_valves(value_df, presence_df):
    """
    Apply safety valve defaults to fill in missing or empty values for certain fields.
    This function modifies value_df and presence_df in place.

    Args:
        value_df (pandas DataFrame): DataFrame with OB values
        presence_df (pandas DataFrame): DataFrame indicating field presence

    Returns:
        value_df (pandas DataFrame): Modified DataFrame with safety valve defaults applied
        presence_df (pandas DataFrame): Modified DataFrame with presence updated
    """
    # Define default values for safety valves
    safety_valve_defaults = {
        'target.gaia_id': 'NoGaiaName',
        'target.t_eff': -1000.0,
        'observation.exp_meter_threshold': 50000.0, #absurdly high so that it is not used in the computation of exposure times
        'schedule.num_intranight_cadence': 0,
        'schedule.num_intranight_cadence': 0,
        'schedule.num_inter_cadence': 0,
        'schedule.n_inter_max': 0,
        'schedule.n_intra_max': 1,
        'schedule.n_inter_min': 1,
        'schedule.n_exp': 1,
        'schedule.minimum_elevation': 33,
        'schedule.minimum_moon_separation': 33,
        'schedule.weather_band_1': True,
        'schedule.weather_band_2': True,
        'schedule.weather_band_3': False,
        'metadata.ob_inactive': False,
    }
    
    # Apply safety valves using a loop
    for col_name, default_value in safety_valve_defaults.items():
        if col_name not in value_df.columns:
            # Column doesn't exist, create it with default value
            value_df[col_name] = default_value
            presence_df[col_name] = True
        else:
            # Column exists, fill NaN values with default
            value_df[col_name] = value_df[col_name].fillna(default_value)
            # Also handle empty strings for string columns
            if isinstance(default_value, str):
                value_df[col_name] = value_df[col_name].replace('', default_value)
            presence_df[col_name] = presence_df[col_name] | value_df[col_name].notna()
            
    # Special case for weather bands based on metadata.semid
    # this was only for 2025B while weather bands were being developed
    if 'metadata.semid' in value_df.columns:
        # Check for 2025B_E473 semid and set opposite weather band values
        mask_2025B_E473 = value_df['metadata.semid'] == '2025B_E473'
        if mask_2025B_E473.any():
            # Set weather bands to opposite values for 2025B_E473
            if 'schedule.weather_band_1' in value_df.columns:
                value_df.loc[mask_2025B_E473, 'schedule.weather_band_1'] = False
            if 'schedule.weather_band_2' in value_df.columns:
                value_df.loc[mask_2025B_E473, 'schedule.weather_band_2'] = False
            if 'schedule.weather_band_3' in value_df.columns:
                value_df.loc[mask_2025B_E473, 'schedule.weather_band_3'] = True

    return value_df, presence_df

def create_checks_dataframes(OBs, required_fields):
    """
    Create the dataframes to determine the good and bad OBs.

    Args:
        OBs (json): the OB information in json format
        required_fields (list): a list of the required fields that must be present

    Returns:
        value_df (pandas DataFrame): a DataFrame with the values of the OBs
        presence_df (pandas DataFrame): a DataFrame with the indication of fields existing or not for the OBs
        all_true_mask (pandas Series): a mask indicating which OBs in the list are good. 
    """
    # Store flattened rows and collect all keys
    flat_value_rows = []
    flat_presence_rows = []
    all_keys = set()
    for entry in OBs['observing_blocks']:
        flat = flatten(entry)
        flat_value_rows.append(flat)
        presence_row = {k: True for k in flat}
        flat_presence_rows.append(presence_row)
        all_keys.update(flat.keys())

    columns = sorted(all_keys)

    value_df = pd.DataFrame([
        [row.get(col, np.nan) for col in columns]
        for row in flat_value_rows
    ], columns=columns)

    presence_df = pd.DataFrame([
        [row.get(col, False) for col in columns]
        for row in flat_presence_rows
    ], columns=columns)

    # Optional: add row labels
    index_labels = [f"entry_{i+1}" for i in range(len(value_df))]
    value_df.index = index_labels
    presence_df.index = index_labels

    # Catch values that exist, but are None, "<NA>", or blank
    for col in columns:
        for idx in index_labels:
            val = value_df.at[idx, col]
            if not pd.api.types.is_scalar(val) or pd.isna(val):
                presence_df.at[idx, col] = False
            elif isinstance(val, str) and val.strip() == "":
                # Catch empty strings or whitespace-only strings
                presence_df.at[idx, col] = False

    return value_df, presence_df

def cast_columns(df):
    """
    Cast columns to their appropriate data types based on the column_definitions dictionary.

    Args:
        df (pandas DataFrame): the DataFrame to cast the columns of

    Returns:
        df (pandas DataFrame): the DataFrame with the columns cast to the appropriate data types
    """
    df = df.copy()
    for col, col_info in column_definitions.items():
        if col in df.columns:
            dtype = col_info['type']
            if dtype in ['Int64', 'Float64']:
                df[col] = pd.to_numeric(df[col], errors='coerce').astype(dtype)
            elif dtype == 'string':
                df[col] = df[col].astype('string')
            elif dtype == 'boolean':
                # Convert to boolean, handling various representations
                df[col] = df[col].astype('boolean')
            else:
                raise ValueError(f"Unsupported dtype: {dtype}. Only 'Int64', 'Float64', 'string', and 'boolean' are allowed.")
    return df

def validate_and_convert_coordinates(df):
    """
    Validate and convert RA/Dec coordinates from string format (hourangle/deg) to degrees.
    Removes rows with invalid coordinates and prints warnings for removed targets.

    Args:
        df (pandas DataFrame): DataFrame with 'ra' and 'dec' columns as strings in hourangle/deg format

    Returns:
        df (pandas DataFrame): DataFrame with valid coordinates converted to degrees, invalid rows removed
    """
    ra_list = df['ra'].astype(str).tolist()
    dec_list = df['dec'].astype(str).tolist()
    
    # Try to create SkyCoord and handle invalid coordinates
    valid_indices = []
    invalid_targets = []
    
    for i, (ra, dec) in enumerate(zip(ra_list, dec_list)):
        try:
            # Test if this coordinate pair is valid
            test_coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
            valid_indices.append(i)
        except Exception as e:
            # Get the target info for reporting (using renamed column names)
            target_id = df.iloc[i].get('unique_id', 'Unknown ID')
            target_name = df.iloc[i].get('starname', 'Unknown Name')
            invalid_targets.append({
                'id': target_id,
                'starname': target_name,
                'ra': ra,
                'dec': dec,
                'error': str(e)
            })
    
    # Print information about invalid targets
    if invalid_targets:
        print(f"Warning: {len(invalid_targets)} targets have invalid coordinates and will be removed:")
        for target in invalid_targets:
            print(f"  ID: {target['id']}, Star: {target['starname']}, RA: {target['ra']}, Dec: {target['dec']}")
            print(f"    Error: {target['error']}")
    
    # Filter to only valid coordinates
    if len(valid_indices) < len(df):
        df = df.iloc[valid_indices].reset_index(drop=True)
        ra_list = [ra_list[i] for i in valid_indices]
        dec_list = [dec_list[i] for i in valid_indices]
    
    # Now create SkyCoord with only valid coordinates
    coords = SkyCoord(ra=ra_list, dec=dec_list, unit=(u.hourangle, u.deg))
    df['ra'] = coords.ra.deg
    df['dec'] = coords.dec.deg
    
    return df

def sort_good_bad(OBs, awarded_programs):
    """
    Sort the OBs into good and bad buckets.
    
    Args:
        OBs (json): the OB information in json format
        awarded_programs (list): a list of the awarded programs

    Returns:
        trimmed_good (pandas DataFrame): a DataFrame with the good OBs
        bad_OBs_values (pandas DataFrame): a DataFrame with the values of the bad OBs fields
        bad_OBs_hasFields (pandas DataFrame): a DataFrame with the indication of fields existing or not for the bad OBs
    """

    # Create the dataframes
    OB_values, OB_hasFields = create_checks_dataframes(OBs, required_fields)
    
    # Apply safety valves
    run_safety_valves = True
    if run_safety_valves:
        OB_values, OB_hasFields = apply_safety_valves(OB_values, OB_hasFields)
    
    # Create masks considering only the required fields
    def row_is_good(row):
        # Check that all required fields are present
        # If a required field is missing from the dataframe, row.get returns False (default)
        return all(
            row.get(col, False) for col in required_fields
        )
    
    # Apply the good/bad row determination
    pass_OBs_mask = OB_hasFields.apply(row_is_good, axis=1)

    bad_OBs_values = OB_values[~pass_OBs_mask]
    bad_OBs_values.reset_index(inplace=True, drop='True')
    bad_OBs_hasFields = OB_hasFields[~pass_OBs_mask]
    bad_OBs_hasFields.reset_index(inplace=True, drop='True')

    mask = bad_OBs_values['metadata.semid'].isin(awarded_programs)
    bad_OBs_values = bad_OBs_values[mask].reset_index(drop=True)
    bad_OBs_hasFields = bad_OBs_hasFields[mask].reset_index(drop=True)

    good_OB_values = OB_values[pass_OBs_mask]
    good_OB_values.reset_index(inplace=True, drop='True')
    good_OBs = cast_columns(good_OB_values)

    good_OBs_awarded = good_OBs[good_OBs['metadata.semid'].isin(awarded_programs)]
    good_OBs_awarded.reset_index(inplace=True, drop='True')
    
    # Create column mapping from column_definitions
    new_column_names = {col: col_info['new_name'] for col, col_info in column_definitions.items()}
    trimmed_good = good_OBs_awarded[list(column_definitions.keys())].rename(columns=new_column_names)

    # Validate and convert coordinates
    trimmed_good = validate_and_convert_coordinates(trimmed_good)

    return trimmed_good, bad_OBs_values, bad_OBs_hasFields

def recompute_exposure_times(request_frame, slowdown_factor):
    """
    Recompute the exposure times for the request frame based on the band number slowdown factor.

    Args:
        request_frame (pandas DataFrame): the request.csv in dataframe format
        slowdown_factor (float): the slowdown factor to apply to the exposure times

    Returns:
        new_exptimes (list): a list of the new exposure times based on slowdown. 
    """
    # These values determined emperically using KPF data spanning a year. 
    # Do not change unless you have good reason.
    factor = 40
    slope_median = -0.362
    intercept_median = 8.889

    rate = slope_median*request_frame['gmag'] + intercept_median
    time = (request_frame['exp_meter_threshold']*factor*10**6)/(10**rate)
    time = time.clip(lower=12)
    if slowdown_factor > 1:
        newtime = (time * slowdown_factor).clip(upper=request_frame["exptime"]).round().astype("Int64")
    else:
        newtime = (time * slowdown_factor).clip(lower=request_frame["exptime"]).round().astype("Int64")
    return newtime

def analyze_bad_obs(trimmed_good, bad_OBs_values, bad_OBs_hasFields, awarded_programs, required_fields=required_fields):
    """
    Analyze the bad OBs and produce a count of bad OBs by semester and a histogram of bad OBs by field.

    Args:
        trimmed_good (pandas DataFrame): the good OBs 
        bad_OBs_values (pandas DataFrame): the values of the fields in the bad OBs 
        bad_OBs_hasFields (pandas DataFrame): the existence of the fields in the bad OBs
        awarded_programs (list): a list of the awarded programs
        required_fields (list): a list of the required fields

    Returns:
        - bad_obs_count_by_semid: dict {metadata.semid: count of bad OBs}
        - bad_field_histogram: dict {field: count of times field was missing in a bad OB}
    """
    # 1. Count bad OBs per metadata.semid
    if 'metadata.semid' in bad_OBs_values.columns:
        bad_obs_count_by_semid = bad_OBs_values['metadata.semid'].value_counts().to_dict()
    else:
        bad_obs_count_by_semid = {}
    # Ensure all awarded_programs are present as keys
    if awarded_programs is not None:
        for semid in awarded_programs:
            if semid not in bad_obs_count_by_semid:
                bad_obs_count_by_semid[semid] = 0

    # 2. Histogram of missing fields (reasons for bad OBs) - only check required fields
    bad_field_histogram = {col: 0 for col in required_fields if col in bad_OBs_hasFields.columns}
    for idx, row in bad_OBs_hasFields.iterrows():
        for col in bad_field_histogram:
            if not bool(row.get(col, False)):
                bad_field_histogram[col] += 1

    return bad_obs_count_by_semid, bad_field_histogram

def plot_bad_obs_histograms(bad_obs_count_by_semid, bad_field_histogram):
    """
    Plots histograms for bad_obs_count_by_semid and bad_field_histogram.
    X: keys, Y: values.

    Args:
        bad_obs_count_by_semid (dict): a dictionary mapping the program code to the count of bad OBs
        bad_field_histogram (dict): a dictionary mapping the field to the count of times field was missing in a bad OB

    Returns:
        None
    """
    import matplotlib.pyplot as plt

    # Plot for bad_obs_count_by_semid
    plt.figure(figsize=(10, 4))
    plt.bar(bad_obs_count_by_semid.keys(), bad_obs_count_by_semid.values())
    plt.xlabel('Program (metadata.semid)')
    plt.ylabel('Number of Bad OBs')
    plt.title('Number of Bad OBs per Program')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    # Plot for bad_field_histogram
    plt.figure(figsize=(12, 4))
    plt.bar(bad_field_histogram.keys(), bad_field_histogram.values())
    plt.xlabel('Field')
    plt.ylabel('Count as Reason for Bad OB')
    plt.title('Frequency of Each Field as Reason for Bad OB')
    plt.xticks(rotation=90, ha='right')
    plt.tight_layout()
    plt.show()

def inspect_row(df_exists, df_values, row_num, required_fields=required_fields):
    """
    Inspect and print a summary of a specific row's key existence and requirement status.

    Args:
        df_exists (pandas DataFrame): the existence of the fields in the OBs
        df_values (pandas DataFrame): the values of the fields in the OBs
        row_num (int): the row number to inspect
        required_fields (list): a list of the required fields

    Returns:
        email_body (str): the email body for the inspection
    """
    row_exists = df_exists.iloc[row_num]
    row_values = df_values.iloc[row_num]

    lines = []
    lines.append(f"_id:         {row_values['_id']}")
    lines.append(f"target_name: {row_values['target.target_name']}")
    lines.append(f"semester_id: {row_values['metadata.semid']}")
    lines.append("-" * 40)
    lines.append(f"{'Missing But Required Fields':<40}")
    lines.append("-" * 40)

    for col in required_fields:
        if col in df_exists.columns:
            exists = bool(row_exists[col])
            if not exists:
                lines.append(f"{col:<40}")

    email_body = email_template.format(
    semid=row_values['metadata.semid'],
    starname=row_values['target.target_name'],
    _id=row_values['_id'],
    badparams="\n".join(lines)
    )
    return email_body

email_template = """
Hello,

As the PI of a KPFCC program, {semid}, you are receiving this email because at least one of the OBs
you submitted to the queue is incorrect, insufficient, or both. This email is a courtesy notice that the
following OB has been rejected by the Community Cadence system and therefore is not scheduled for any
observations. When this OB has been remedied, the system will automatically begin including it in the
optimal scheduling algorithm. Each day the algorithm pulls from the OB database and performs these checks and so
each day that the OB is not in compliance, you will receive this email. To effectively shut off this email, either
1) fix the affected OB
2) delete the OB

The following OB id/target_name was rejected for the follwowing reason(s). These fields are missing or incorrectly formatted: \n
{badparams}

Note: for fields that begin with "target", you may just need to hit the bullseye button again on the webform and then resubmit.

If you have any questions about why this OB was rejected or on the process of remedy/resubmission, please reach out
to KPF-CC Project Scientist Jack Lubin (jblubin@ucla.edu)

Very best,
Jack
"""

# This portion not ready yet, but saving for future use.
#---------------------------------------------------------
# import smtplib
# from email.mime.multipart import MIMEMultipart
# from email.mime.text import MIMEText

# Email Configuration
# SMTP_SERVER = 'smtp.example.com'
# SMTP_PORT = 587
# SENDER_EMAIL = 'youremail@example.com'
# SENDER_PASSWORD = 'yourpassword'
#
# def send_email(receiver_email, subject, body):
#     msg = MIMEMultipart()
#     msg['From'] = SENDER_EMAIL
#     msg['To'] = receiver_email
#     msg['Subject'] = subject
#     msg.attach(MIMEText(body, 'plain'))
#
#     try:
#         # Connect to SMTP server
#         server = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
#         server.starttls()  # Secure the connection
#         server.login(SENDER_EMAIL, SENDER_PASSWORD)
#         text = msg.as_string()
#         server.sendmail(SENDER_EMAIL, receiver_email, text)
#         server.quit()
#         print(f"Email sent to {receiver_email}")
#     except Exception as e:
#         print(f"Failed to send email to {receiver_email}: {e}")

# def define_email(bad_obs, i):
#     badkeys = list(bad_obs.keys())
#     email_address = bad_obs[badkeys[i]][4]
#     human_unique = str(bad_obs[badkeys[i]][0])+"_"+str(bad_obs[badkeys[i]][1])+"_"+str(bad_obs[badkeys[i]][2])
#     reasons=""
#     for j in range(len(bad_obs[badkeys[i]][5])):
#         if j == len(bad_obs[badkeys[i]][5])-1:
#             add = ""
#         else:
#             add = "\n"
#         reasons += "-- " + str(bad_obs[badkeys[i]][5][j]) + add

#     email_body = email_template.format(
#         name=bad_obs[badkeys[i]][3],
#         semester=bad_obs[badkeys[i]][0],
#         program=bad_obs[badkeys[i]][1],
#         starname=bad_obs[badkeys[i]][2],
#         _id=badkeys[i],
#         h_id=human_unique,
#         badparams=reasons
#     )
#     return email_address, email_body

def filter_request_csv(request_df, weather_band_num):
    """
    Filter request.csv file to only keep rows where weather_band_X = True
    
    Args:
        request_file_path (str): Path to the request.csv file
        weather_band_num (int): Weather band number to filter by
        
    Returns:
        bool: True if filtering was successful, False otherwise
    """
    weather_band_col = f'weather_band_{weather_band_num}'
    
    if weather_band_col in request_df.columns:
        filtered_df = request_df[request_df[weather_band_col] == True]
    else:
        print(f'Warning: Column {weather_band_col} not found in request.csv. No filtering applied.')
    return filtered_df
    
def update_allocation_file(allocation_df, current_date):
    """
    Update allocation.csv file with today's 12-degree twilight times
    
    Args:
        allocation_file_path (str): Path to the allocation.csv file
        current_date (str): Current date in YYYY-MM-DD format
        
    Returns:
        bool: True if update was successful, False otherwise
    """
    date_exists = False
    date_idx = -1
    
    # Check if current date exists in allocation file
    for idx, row in allocation_df.iterrows():
        row_date = str(row['start'])[:10]  # Get YYYY-MM-DD portion
        if row_date == current_date:
            date_exists = True
            date_idx = idx
            break
    
    # Get 12-degree twilight times for current date
    observatory = 'Keck Observatory'
    keck = apl.Observer.at_site(observatory)
    day = Time(current_date, format='isot', scale='utc')
    
    evening_12 = keck.twilight_evening_nautical(day, which='next')
    morning_12 = keck.twilight_morning_nautical(day, which='next')
    
    if not date_exists:
        print(f'Adding allocation row for current_day: {current_date}')
        # Add new row at the bottom
        new_row = pd.DataFrame({
            'start': [evening_12.strftime('%Y-%m-%dT%H:%M')],
            'stop': [morning_12.strftime('%Y-%m-%dT%H:%M')]
        })
        allocation_df = pd.concat([allocation_df, new_row], ignore_index=True)
        allocation_df.loc[len(allocation_df)-1, 'comment'] = 'added as part of full-band processing'
        print(f'Added allocation: {evening_12.iso} to {morning_12.iso}')
    else:
        print(f'Updating existing allocation row for current_day: {current_date}')
        # Update existing row
        allocation_df.loc[date_idx, 'start'] = evening_12.strftime('%Y-%m-%dT%H:%M')
        allocation_df.loc[date_idx, 'stop'] = morning_12.strftime('%Y-%m-%dT%H:%M')
        allocation_df.loc[date_idx, 'comment'] = 'added as part of full-band processing'
        print(f'Updated allocation: {evening_12.iso} to {morning_12.iso}')
    
    return allocation_df

def write_starlist(frame, solution_frame, night_start_time, extras, filler_stars, current_day,
                    outputdir, version='nominal'):
    """
    Generate the nightly script in the format required by the Keck "Magiq" software. 
    Backwards compatable to pre-KPF-CC observing.

    Args:
        frame (dataframe): the request_frame of just the targets that were selected to be observed tonight
        solution_frame (dataframe): the solution attribute from the TTP model.plotly object
        night_start_time (astropy time object): Beginning of observing interval
        extras (array): starnames of "extra" stars (those not fit into the script)
        filler_stars (array): star names of the stars added in the bonus round
        current_day (str): today's date in format YYYY-MM-DD
        outputdir (str): the directory to save the script file
        version (str): a tag for thescript (e.g. nominal, slowdown, backups, etc)

    Returns:
        lines (str): the script file as a string
    """
    # Cast starname column to strings to ensure proper matching
    frame['starname'] = frame['starname'].astype(str)
    
    # Cast extras star names to strings to ensure proper matching
    if extras is not None and len(extras) > 0:
        if hasattr(extras, 'astype'):
            # If extras is a DataFrame
            extras['Starname'] = extras['Starname'].astype(str)
        else:
            # If extras is a list, convert each star name to string
            extras['Starname'] = [str(star) for star in extras['Starname']]
    
    total_exptime = 0
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    script_file = os.path.join(outputdir,'script_{}_{}.txt'.format(current_day, version))

    lines = []
    for i, item in enumerate(solution_frame['Starname']):
        filler_flag = solution_frame['Starname'][i] in filler_stars
        row = frame.loc[frame['unique_id'] == solution_frame['Starname'][i]]
        row.reset_index(inplace=True)
        total_exptime += float(row['exptime'].iloc[0])

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
        row = frame.loc[frame['unique_id'] == extras['Starname'][j]]
        row.reset_index(inplace=True)
        lines.append(format_kpf_row(row, '24:00', extras['First Available'][j],
                    extras['Last Available'][j], current_day, filler_flag, True))

    # add buffer lines to end of file
    lines.append("")
    lines.append("")

    with open(script_file, 'w') as f:
        f.write('\n'.join(lines))
    print("Total Open Shutter Time Scheduled: " + str(np.round((total_exptime/3600),2)) + " hours")
    return lines

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

    equinox = '2000'
    # Handle missing pmra/pmdec columns with default values
    pmra = row.get('pmra', pd.Series([0.0])).iloc[0] if 'pmra' in row else 0.0
    pmdec = row.get('pmdec', pd.Series([0.0])).iloc[0] if 'pmdec' in row else 0.0
    updated_ra, updated_dec = pm_correcter(row['ra'].iloc[0], row['dec'].iloc[0],
                                pmra, pmdec, current_day, equinox=equinox)
    if updated_dec[0] != "-":
        updated_dec = "+" + updated_dec

    starname_str = str(row['starname'].iloc[0])
    namestring = ' '*(16-len(starname_str[:16])) + starname_str[:16]

    # Handle missing columns with default values
    jmag_val = row.get('jmag', [15.0])[0] if 'jmag' in row else 15.0
    gmag_val = row.get('gmag', [15.0])[0] if 'gmag' in row else 15.0
    teff_val = row.get('teff', [5000])[0] if 'teff' in row else 5000
    gaia_id_val = row.get('gaia_id', ['UNKNOWN'])[0] if 'gaia_id' in row else 'UNKNOWN'
    
    # Convert to float safely, with fallback to defaults
    try:
        jmag_val = float(jmag_val) if jmag_val is not None else 15.0
    except (ValueError, TypeError):
        jmag_val = 25.0
    
    try:
        gmag_val = float(gmag_val) if gmag_val is not None else 15.0
    except (ValueError, TypeError):
        gmag_val = 25.0
    
    try:
        teff_val = float(teff_val) if teff_val is not None else 5000
    except (ValueError, TypeError):
        teff_val = 0.0
    
    jmagstring = ('jmag=' + str(np.round(float(jmag_val),1)) + ' '* \
        (4-len(str(np.round(float(jmag_val),1)))))
    exposurestring = (' '*(4-len(str(int(row['exptime'].iloc[0])))) + \
        str(int(row['exptime'].iloc[0])) + '/' + \
        str(int(row['exptime'].iloc[0])) + ' '* \
        (4-len(str(int(row['exptime'].iloc[0])))))

    ofstring = ('1of' + str(int(row['n_intra_max'].iloc[0])))
    scstring = 'sc=' + 'T'

    numstring = str(int(row['n_exp'].iloc[0])) + "x"
    gmagstring = 'gmag=' + str(np.round(float(gmag_val),1)) + \
                                                ' '*(4-len(str(np.round(float(gmag_val),1))))
    teffstr = 'Teff=' + str(int(teff_val)) + \
                                    ' '*(4-len(str(int(teff_val))))

    gaiastring = str(gaia_id_val) + ' '*(25-len(str(gaia_id_val)))
    programstring = row['program_code'].iloc[0]

    if filler_flag:
        # All targets added in round 2 bonus round are lower priority
        priostring = "p3"
    else:
        priostring = "p1"

    if extra == False:
        timestring2 = str(obs_time)
    else:
        # designate a nonsense time
        timestring2 = "24:00"

    line = (namestring + ' ' + updated_ra + ' ' + updated_dec + ' ' + str(equinox) + ' '
                + jmagstring + ' ' + exposurestring + ' ' + ofstring + ' ' + scstring +  ' '
                + numstring + ' '+ gmagstring + ' ' + teffstr + ' ' + gaiastring + ' CC '
                        + priostring + ' ' + programstring + ' ' + timestring2 +
                         ' ' + first_available  + ' ' + last_available )

    # Handle missing Observing Notes column
    observing_notes = row.get('Observing Notes', [''])[0] if 'Observing Notes' in row else ''
    if observing_notes and not pd.isnull(observing_notes):
        line += (' ' + str(observing_notes))

    return line

def pm_correcter(ra, dec, pmra, pmdec, current_day, equinox="2000"):
    """
    Update a star's coordinates due to proper motion.

    Args:
        ra (float): RA in degrees
        dec (float): Dec in degrees
        pmra (float): proper motion in RA (mas/yr), including cos(Dec)
        pmdec (float): proper motion in Dec (mas/yr)
        equinox (str): original epoch (e.g. '2000.0')
        current_day (str): date to which to propagate (e.g. '2025-04-30')

    Returns:
        formatted_ra (str), formatted_dec (str): updated coordinates as strings
    """
    start_time = Time(f'J{equinox}')
    current_time = Time(current_day)
    coord = SkyCoord(
        ra=ra * u.deg,
        dec=dec * u.deg,
        pm_ra_cosdec=pmra * u.mas/u.yr,
        pm_dec=pmdec * u.mas/u.yr,
        obstime=start_time
    )
    new_coord = coord.apply_space_motion(new_obstime=current_time)
    formatted_ra = new_coord.ra.to_string(unit=u.hourangle, sep=' ', pad=True, precision=1)
    formatted_dec = new_coord.dec.to_string(unit=u.deg, sep=' ', pad=True, precision=0)

    return formatted_ra, formatted_dec

class Access_KPFCC(Access):
    """
    Keck Observatory-specific Access class that inherits from the base Access class.
    Overrides compute_altaz() and compute_clear() methods with Keck-specific implementations.
    """
    
    def __init__(self, 
                 semester_start_date, 
                 semester_length, 
                 n_nights_in_semester,
                 today_starting_night,  
                 current_day, 
                 all_dates_dict, 
                 all_dates_array, 
                 slot_size, 
                 slots_needed_for_exposure_dict, 
                 custom_file, 
                 allocation_file, 
                 past_history, 
                 output_directory, 
                 run_weather_loss, 
                 run_band3, 
                 observatory_string, 
                 request_frame, 
                 ):
        """
        Initialize the Access_KPFCC object with explicit parameters.
        Calls the parent Access.__init__() to set up base class attributes.
            
        Args:
            semester_start_date: Start date of the semester
            semester_length: Total number of nights in the semester
            n_nights_in_semester: Number of remaining nights in the semester
            today_starting_night: Starting night number for today
            current_day: Current day identifier
            all_dates_dict: Dictionary mapping dates to day numbers
            all_dates_array: Array of date strings for the semester
            slot_size: Size of each time slot in minutes
            slots_needed_for_exposure_dict: Dictionary mapping star names to required slots
            custom_file: Path to custom times file
            allocation_file: Path to allocation file
            past_history: Past observation history
            output_directory: Directory for output files
            run_weather_loss: Whether to run weather loss simulation
            run_band3: Whether to run band 3 (used for not peforming the is_observble step for the football plot)
            observatory_string: Observatory name/location string
            request_frame: DataFrame containing request information
        """
        # Call parent class __init__ to set up base attributes
        super().__init__(
                 semester_start_date, 
                 semester_length, 
                 n_nights_in_semester,
                 today_starting_night,  
                 current_day, 
                 all_dates_dict, 
                 all_dates_array, 
                 slot_size, 
                 slots_needed_for_exposure_dict, 
                 custom_file, 
                 allocation_file, 
                 past_history, 
                 output_directory, 
                 run_weather_loss, 
                 run_band3, 
                 observatory_string, 
                 request_frame)
        
        self.weather_loss_file = os.path.join(self.DATADIR, "maunakea_weather_loss_data.csv")

        self.use_K1 = True
        # See here: https://www2.keck.hawaii.edu/inst/common/TelLimits.html
        if self.use_K1:
            # K1 Telescope limits 
            self.nays_az_low = 5.3
            self.nays_az_high = 146.2
            self.nays_alt = 33.3
        else:   
            # K2 Telescope limits 
            self.nays_az_low = 185.3
            self.nays_az_high = 332.8
            self.nays_alt = 36.8
        self.tel_min = 18
        self.tel_max = 85

    def compute_altaz(self, tel_min):
        """
        Compute boolean mask of is_altaz for targets according to a minimum elevation. 
        Specific to Keck Observatory's K1/K2 pointing limits.
        This method overrides the base class method to include Keck-specific nasmyth deck constraints.

        Args:
            tel_min (float): the minimum elevation for the telescope (ignored, uses self.tel_min instead)
        
        Returns:
            is_altaz (array): boolean mask of is_altaz for targets
        """
        # Compute base alt/az pattern, shape = (ntargets, nslots)

        altazes = self.observatory.altaz(self.timegrid, self.targets, grid_times_targets=True)
        alts = altazes.alt.deg
        azes = altazes.az.deg
        min_elevation = self.request_frame['minimum_elevation'].values  # Get PI-desired minimum elevation values
        min_elevation = np.maximum(min_elevation, self.tel_min)  # Ensure minimum elevation is at least tel_min
        
        # 2D mask (n targets, n slots))
        is_altaz0 = np.ones_like(alts, dtype=bool)
        # Remove nasmyth deck - Keck-specific constraint
        is_altaz0 &= ~((self.nays_az_low < azes ) & (azes < self.nays_az_high) & (alts < self.nays_alt))
       
        # Remove min elevation using per-row minimum_elevation values for all stars
        fail = alts < min_elevation[:, np.newaxis]  # broadcast elevation array
        is_altaz0 &= ~fail
        
        # All stars must be between tel_min and tel_max deg
        fail = (alts < self.tel_min) | (alts > self.tel_max)
        is_altaz0 &= ~fail

        # Pre-compute the sidereal times for interpolation
        x = self.timegrid.sidereal_time('mean').value
        x_new = self.slotmidpoints.sidereal_time('mean').value
        idx = np.searchsorted(x, x_new, side='left')
        idx = np.clip(idx, 0, len(x)-1)  # Handle edge cases

        # Create 3D mask by indexing the 2D pattern
        self.is_altaz = is_altaz0[:,idx]

    def compute_clear(self):
        """
        Compute boolean mask of is_clear for all targets according to the clear times.
        """
        self.is_clear = np.ones_like(self.is_altaz, dtype=bool)
        if self.run_weather_loss:
            logs.info("Running weather loss model.")
            self.get_loss_stats(self.weather_loss_file)
            self.is_clear = self.simulate_weather_losses(covariance=0.14)
            self.is_clear = np.tile(self.is_clear[np.newaxis, :, :], (self.ntargets, 1, 1))
        else:
            logs.info("Pretending weather is always clear!")
            self.is_clear = np.ones((self.ntargets, self.nnights, self.nslots), dtype=bool)