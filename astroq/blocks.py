"""
Module for processing data from Keck Observatory's custom made Observing Block (OB) database.

Example usage:
    import ob_functions as ob
"""

# Standard library imports
import json
import os

# Third-party imports
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import astroplan as apl

from kpf_etc.etc import kpf_etc_snr

exception_fields = ['_id', 'del_flag', 'metadata.comment', 'metadata.details', 'metadata.history',
                    'metadata.instruments', 'metadata.is_approved', 'metadata.last_modification',
                    'metadata.ob_feasible', 'metadata.observer_name', 'metadata.state', 'metadata.status',
                    'metadata.submitted', 'metadata.submitter', 'metadata.tags', 'observation.auto_exp_meter',
                    'observation.auto_nd_filters', 'observation.cal_n_d_1', 'observation.cal_n_d_2',
                    'observation.exp_meter_exp_time', 'observation.exp_meter_mode', 'observation.guide_here',
                    'observation.object', 'observation.take_simulcal', 'observation.exp_meter_bin',
                    'observation.trigger_ca_h_k', 'observation.trigger_green', 'observation.trigger_red',
                    'schedule.accessibility_map', 'schedule.days_observable', 'schedule.fast_read_mode_requested',
                    'schedule.num_visits_per_night',
                    'schedule.rise_semester_day', 'schedule.scheduling_mode', 'schedule.sets_semester_day',
                    'schedule.total_observations_requested', 'schedule.total_time_for_target',
                    'schedule.total_time_for_target_hours', 'target.isNew', 'target.parallax', 'target.equinox', 'target.systemic_velocity',
                    'target.tic_id', 'target.two_mass_id', 'schedule.weather_band', 'target.catalog_comment',
                    'calibration.cal_n_d_1', 'calibration.cal_n_d_2', 'calibration.cal_source', 'calibration.exp_meter_bin',
                    'calibration.exp_meter_exp_time', 'calibration.exp_meter_mode', 'calibration.exp_meter_threshold',
                    'calibration.exposure_time', 'calibration.intensity_monitor', 'calibration.num_exposures', 'calibration.object',
                    'calibration.open_science_shutter', 'calibration.open_sky_shutter', 'calibration.take_simulcal',
                    'calibration.trigger_ca_h_k', 'calibration.trigger_green', 'calibration.trigger_red',
                    'calibration.wide_flat_pos', 'observation.block_sky', 'observation.nod_e', 'observation.nod_n',
                    'schedule.isNew', 'observation.isNew', 'schedule.comment', 'target.d_ra', 'target.d_dec', 'target.undefined',
                    'target.ra_deg', 'target.dec_deg', 'observation.undefined', 'schedule.num_visits_per_night', 'schedule.undefined',
                    'schedule.custom_time_constraints', #make this an exception field because it is handled elsewhere to make the custom.csv file
                    'schedule.desired_num_visits_per_night', 'schedule.minimum_num_visits_per_night', 'history', # NOTE: this line will be removed 
                    'schedule.weather_band_1', 'schedule.weather_band_2', 'schedule.weather_band_3'# NOTE: this line will be removed 
                    # 'observation.exp_meter_threshold', 'schedule.minimum_elevation', 'schedule.minimum_moon_separation', # for now. 
]

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
}

def pull_OBs(semester):
    """
    Pull the latest database OBs down to local.

    Args:
        semester (str) - the semester from which to query OBs, format YYYYL
        histories (bool) - if True, pull the history of OBs for the semester, if False, pull the latest OBs for the semester

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

def format_custom_csv(OBs):
    """
    Format the custom csv file for the OBs.
    """
    rows = []
    for ob in OBs['observing_blocks']:
        try:
            ctc = ob['schedule']['custom_time_constraints']
            unique_id = ob['_id']
            starname = ob['target']['target_name']

            # Handle ctc as a list with multiple constraints
            if isinstance(ctc, list) and len(ctc) > 0:
                # Process each constraint in the list
                for i, constraint in enumerate(ctc):
                    if isinstance(constraint, dict):
                        start = constraint.get('start_datetime', '')
                        stop = constraint.get('end_datetime', '')
                        
                        # Only add rows that have the required data
                        if start and stop:
                            rows.append({
                                'unique_id': unique_id,
                                'starname': starname,
                                'start': start,
                                'stop': stop
                            })
            elif isinstance(ctc, dict):
                # If it's already a dict, use it directly
                start = ctc.get('start_datetime', '')
                stop = ctc.get('end_datetime', '')
                
                if start and stop:
                    rows.append({
                        'unique_id': unique_id,
                        'starname': starname,
                        'start': start,
                        'stop': stop
                    })
        except Exception as e:
            pass

    # Create DataFrame and save to CSV
    if len(rows) > 0:
        df = pd.DataFrame(rows)
    else:
        # Create empty DataFrame with proper column headers
        df = pd.DataFrame(columns=['unique_id', 'starname', 'start', 'stop'])
    
    return df 
        
def pull_allocation_info(start_date, numdays, instrument):
    params = {}
    params['cmd'] = "getSchedule"
    params["date"] = start_date
    params["numdays"] = numdays
    params["instrument"] = instrument 
    url = "https://www3.keck.hawaii.edu/api/getSchedule/"
    try:
        data = requests.get(url, params=params)
        data_json = json.loads(data.text)
        df = pd.DataFrame(data_json)
        awarded_programs = df['ProjCode'].unique()
        df['start'] = pd.to_datetime(df['Date'] + ' ' + df['StartTime']).dt.strftime('%Y-%m-%dT%H:%M')
        df['stop']  = pd.to_datetime(df['Date'] + ' ' + df['EndTime']).dt.strftime('%Y-%m-%dT%H:%M')
        allocation_frame = df[['start', 'stop']].copy() # TODO: add observer and comment
        
        # Calculate hours for each row
        start_times = pd.to_datetime(df['start'])
        stop_times = pd.to_datetime(df['stop'])
        hours_per_row = (stop_times - start_times).dt.total_seconds() / 3600
        
        # Add ProjCode and hours to the dataframe
        df['hours'] = hours_per_row
        
        # Calculate total hours per ProjCode
        hours_by_program = df.groupby('ProjCode')['hours'].sum().round(3).to_dict()
        
    except:
        print("ERROR: allocation information not found. Double check date and instrument. Saving an empty file.")
        allocation_frame = pd.DataFrame(columns=['start', 'stop'])
        awarded_programs = []
        hours_by_program = {}
    return allocation_frame, hours_by_program

def format_keck_allocation_info(allocation_file):
    """
    Read in allocation file and parse start/stop times to calculate hours by program.
    
    Args:
        allocation_file (str): the path and filename to the downloaded csv
        savepath (str): the path and filename where to save the processed allocation data
        
    Returns:
        hours_by_program (dict): dictionary mapping ProjCode to total hours
    """
    allocation = pd.read_csv(allocation_file)
    
    # Parse the Time column to extract start and stop times
    # Pattern: "10:28 - 15:07 ( 50%)"
    pattern = r'(\d{2}:\d{2}) - (\d{2}:\d{2}) \(\s*(\d{2,3})%\)'
    allocation[['Start', 'Stop', 'Percentage']] = allocation['Time'].str.extract(pattern)
    
    # Convert start and stop times to datetime for hour calculation
    allocation['start'] = pd.to_datetime(allocation['Date'] + ' ' + allocation['Start']).dt.strftime('%Y-%m-%dT%H:%M')
    allocation['stop'] = pd.to_datetime(allocation['Date'] + ' ' + allocation['Stop']).dt.strftime('%Y-%m-%dT%H:%M')
    
    # Calculate hours for each row
    start_times = pd.to_datetime(allocation['start'])
    stop_times = pd.to_datetime(allocation['stop'])
    hours_per_row = (stop_times - start_times).dt.total_seconds() / 3600
    
    # Add hours to the dataframe
    allocation['hours'] = hours_per_row
    
    # Create allocation frame with start/stop times for saving
    allocation_frame = allocation[['start', 'stop']].copy()
    
    # Ensure the directory exists before saving
    os.makedirs(os.path.dirname(savepath), exist_ok=True)
    
    # Calculate total hours per ProjCode
    hours_by_program = allocation.groupby('ProjCode')['hours'].sum().round(3).to_dict()
    
    return allocation_frame, hours_by_program

def get_request_sheet(OBs, awarded_programs, savepath):
    good_obs, bad_obs_values, bad_obs_hasFields = sort_good_bad(OBs, awarded_programs)

    # Filter bad OBs to only those in awarded programs
    if 'metadata.semid' in bad_obs_values.columns:
        mask = bad_obs_values['metadata.semid'].isin(awarded_programs)
        bad_obs_values = bad_obs_values[mask].reset_index(drop=True)
        bad_obs_hasFields = bad_obs_hasFields[mask].reset_index(drop=True)

    bad_obs_count_by_semid, bad_field_histogram = analyze_bad_obs(good_obs, bad_obs_values, bad_obs_hasFields, awarded_programs)
    good_obs.sort_values(by='program_code', inplace=True)
    good_obs.reset_index(inplace=True, drop=True)
    
    # Cast starname column to strings to ensure proper matching
    if 'starname' in good_obs.columns:
        good_obs['starname'] = good_obs['starname'].astype(str)
    
    os.makedirs(os.path.dirname(savepath), exist_ok=True)
    return good_obs, bad_obs_values, bad_obs_hasFields, bad_obs_count_by_semid, bad_field_histogram

def flatten(d, parent_key='', sep='.'):
    items = {}
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.update(flatten(v, new_key, sep=sep))
        else:
            items[new_key] = v
    return items

def create_checks_dataframes(OBs, exception_fields):
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

    # Catch values that exist, but are None or "<NA>"
    for col in columns:
        for idx in index_labels:
            val = value_df.at[idx, col]
            if not pd.api.types.is_scalar(val) or pd.isna(val):
            # if pd.isna(value_df.at[idx, col]):
            # if value_df.at[idx, col] is None or value_df.at[idx, col] == "<NA>":
                presence_df.at[idx, col] = False

    run_safety_valves = True
    if run_safety_valves:
        # Define default values for safety valves
        safety_valve_defaults = {
            'schedule.num_intranight_cadence': 0,
            'schedule.desired_num_visits_per_night': 1,
            'schedule.minimum_num_visits_per_night': 0,  # Will be overridden by special logic below
            'target.gaia_id': 'NoGaiaName',
            'observation.exp_meter_threshold': -1.0,
            'schedule.minimum_elevation': 33,
            'schedule.minimum_moon_separation': 33,
            'schedule.weather_band_1': True,
            'schedule.weather_band_2': True,
            'schedule.weather_band_3': False
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
                
        # Special case for fixing default ExpMeterThreshold
        # Using default of 1.6, this is in MegaPhotons/A which gives SNR ~150
        if 'observation.exp_meter_threshold' in value_df.columns:
            value_df['observation.exp_meter_threshold'] = 1.6

        # Special case for weather bands based on metadata.semid
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
        
        # Special case: minimum_num_visits_per_night should use desired_num_visits_per_night if available
        if 'schedule.desired_num_visits_per_night' in value_df.columns:
            if 'schedule.minimum_num_visits_per_night' in value_df.columns:
                # Fill NaN values in minimum with corresponding desired values
                value_df['schedule.minimum_num_visits_per_night'] = value_df['schedule.minimum_num_visits_per_night'].fillna(value_df['schedule.desired_num_visits_per_night'])
            else:
                # Use desired values as minimum
                value_df['schedule.minimum_num_visits_per_night'] = value_df['schedule.desired_num_visits_per_night']
            presence_df['schedule.minimum_num_visits_per_night'] = presence_df['schedule.minimum_num_visits_per_night'] | value_df['schedule.minimum_num_visits_per_night'].notna()
        
        # Special case: if schedule.num_nights_per_semester == 1, set schedule.num_internight_cadence to 0
        if 'schedule.num_nights_per_semester' in value_df.columns and 'schedule.num_internight_cadence' in value_df.columns:
            mask = value_df['schedule.num_nights_per_semester'] == 1
            value_df.loc[mask, 'schedule.num_internight_cadence'] = 0
            presence_df.loc[mask, 'schedule.num_internight_cadence'] = True


    # Create masks considering the exception fields
    def row_is_good(row):
        # Ignore the exception fields in the presence check
        return all(
            row.get(col, False) if col not in exception_fields else True
            for col in row.index
        )

    # Apply the good/bad row determination
    all_true_mask = presence_df.apply(row_is_good, axis=1)
    some_false_mask = ~all_true_mask

    # Apply masks to both dataframes
    complete_value_df = value_df[all_true_mask]
    incomplete_value_df = value_df[some_false_mask]

    complete_presence_df = presence_df[all_true_mask]
    incomplete_presence_df = presence_df[some_false_mask]

    return value_df, presence_df, all_true_mask

def cast_columns(df):
    """
    Cast columns to their appropriate data types based on the column_definitions dictionary.
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

def sort_good_bad(OBs, awarded_programs):

    OB_values, OB_hasFields, pass_OBs_mask = create_checks_dataframes(OBs, exception_fields)

    bad_OBs_values = OB_values[~pass_OBs_mask]
    bad_OBs_values.reset_index(inplace=True, drop='True')
    bad_OBs_hasFields = OB_hasFields[~pass_OBs_mask]
    bad_OBs_hasFields.reset_index(inplace=True, drop='True')

    if 'metadata.semid' in bad_OBs_values.columns:
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
    trimmed_good.columns.values[9] = 'n_intra_max'

    ra_list = trimmed_good['ra'].astype(str).tolist()
    dec_list = trimmed_good['dec'].astype(str).tolist()
    
    # Try to create SkyCoord and handle invalid coordinates
    valid_indices = []
    invalid_targets = []
    
    for i, (ra, dec) in enumerate(zip(ra_list, dec_list)):
        try:
            # Test if this coordinate pair is valid
            test_coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
            valid_indices.append(i)
        except Exception as e:
            # Get the target info for reporting
            target_id = trimmed_good.iloc[i].get('_id', 'Unknown ID')
            target_name = trimmed_good.iloc[i].get('target.target_name', 'Unknown Name')
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
    if len(valid_indices) < len(trimmed_good):
        trimmed_good = trimmed_good.iloc[valid_indices].reset_index(drop=True)
        ra_list = [ra_list[i] for i in valid_indices]
        dec_list = [dec_list[i] for i in valid_indices]
    
    # Now create SkyCoord with only valid coordinates
    coords = SkyCoord(ra=ra_list, dec=dec_list, unit=(u.hourangle, u.deg))
    trimmed_good['ra'] = coords.ra.deg
    trimmed_good['dec'] = coords.dec.deg

    return trimmed_good, bad_OBs_values, bad_OBs_hasFields

def recompute_exposure_times(request_frame, slowdown_factor):
    """
    Recompute the exposure times for the request frame based on the band number slowdown factor.
    """
    new_exptimes = []
    for i in range(len(request_frame)):
        SNR = 120*(request_frame['exp_meter_threshold'][i]**0.5)
        try:
            new_exptime = kpf_etc_snr(request_frame['teff'][i], request_frame['gmag'][i], SNR, 604)
            if np.isnan(new_exptime):
                new_exptime = request_frame['exptime'][i]*slowdown_factor
        except Exception as e:
            print(f"Error on star {request_frame['starname'][i]} in row {i}: {e}")
            new_exptime = 3600.0
        print(request_frame['starname'][i])
        print(f"New exptime: {new_exptime}")
        print(f"Request exptime: {request_frame['exptime'][i]}")
        final_time = min([new_exptime*slowdown_factor, request_frame['exptime'][i]])
        print(f"Final time: {final_time}")
        print("--------------------------------")
        new_exptimes.append(final_time)
    return new_exptimes

def analyze_bad_obs(trimmed_good, bad_OBs_values, bad_OBs_hasFields, awarded_programs,exception_fields=exception_fields):
    """
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

    # 2. Histogram of missing fields (reasons for bad OBs)
    bad_field_histogram = {col: 0 for col in bad_OBs_hasFields.columns if col not in exception_fields}
    for idx, row in bad_OBs_hasFields.iterrows():
        for col in bad_field_histogram:
            if not bool(row.get(col, True)):
                bad_field_histogram[col] += 1

    return bad_obs_count_by_semid, bad_field_histogram

def plot_bad_obs_histograms(bad_obs_count_by_semid, bad_field_histogram):
    """
    Plots histograms for bad_obs_count_by_semid and bad_field_histogram.
    X: keys, Y: values.
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

def inspect_row(df_exists, df_values, row_num, exception_fields=exception_fields):
    """
    Inspect and print a summary of a specific row's key existence and requirement status.
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

    for col in df_exists.columns:
        exists = bool(row_exists[col])
        is_required = col not in exception_fields
        if is_required and not exists:
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
