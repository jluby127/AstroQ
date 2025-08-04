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

exception_fields = ['_id', 'del_flag', 'metadata.comment', 'metadata.details', 'metadata.history',
                    'metadata.instruments', 'metadata.is_approved', 'metadata.last_modification',
                    'metadata.ob_feasible', 'metadata.observer_name', 'metadata.state', 'metadata.status',
                    'metadata.submitted', 'metadata.submitter', 'metadata.tags', 'observation.auto_exp_meter',
                    'observation.auto_nd_filters', 'observation.cal_n_d_1', 'observation.cal_n_d_2',
                    'observation.exp_meter_exp_time', 'observation.exp_meter_mode',
                    'observation.exp_meter_threshold', 'observation.guide_here',
                    'observation.object', 'observation.take_simulcal', 'observation.exp_meter_bin',
                    'observation.trigger_ca_h_k', 'observation.trigger_green', 'observation.trigger_red',
                    'schedule.accessibility_map', 'schedule.days_observable', 'schedule.fast_read_mode_requested',
                    'schedule.minimum_elevation', 'schedule.minimum_moon_separation', 'schedule.num_visits_per_night',
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
                    'schedule.desired_num_visits_per_night', 'schedule.minimum_num_visits_per_night', 'history'# NOTE: this line will be removed 
]

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
        
def pull_allocation_info(start_date, numdays, instrument, savepath):
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
    allocation_frame.to_csv(savepath, index=False)
    return hours_by_program

def get_request_sheet(OBs, awarded_programs, savepath):
    good_obs, bad_obs_values, bad_obs_hasFields = sort_good_bad(OBs, awarded_programs)
    
    # Add "HD" prefix to star names that are only numbers
    good_obs['starname'] = good_obs['starname'].apply(lambda x: f"HD{x}" if str(x).isdigit() else x)
    
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
    
    good_obs.to_csv(savepath, index=False)
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
        # Safety valve for missing required fields
        if 'schedule.num_intranight_cadence' in value_df.columns:
            value_df['schedule.num_intranight_cadence'] = value_df['schedule.num_intranight_cadence'].fillna(0)
            presence_df['schedule.num_intranight_cadence'] = presence_df['schedule.num_intranight_cadence'] | value_df['schedule.num_intranight_cadence'].notna()
        else:
            value_df['schedule.num_intranight_cadence'] = 0
            presence_df['schedule.num_intranight_cadence'] = True

        if 'schedule.desired_num_visits_per_night' not in value_df.columns:
            value_df['schedule.desired_num_visits_per_night'] = 1
            presence_df['schedule.desired_num_visits_per_night'] = True
        else:
            value_df['schedule.desired_num_visits_per_night'] = 1
            presence_df['schedule.desired_num_visits_per_night'] = True

        if 'schedule.minimum_num_visits_per_night' in value_df.columns:
            if 'schedule.desired_num_visits_per_night' in value_df.columns:
                value_df['schedule.minimum_num_visits_per_night'] = value_df['schedule.minimum_num_visits_per_night'].fillna(value_df['schedule.desired_num_visits_per_night'])
                presence_df['schedule.minimum_num_visits_per_night'] = presence_df['schedule.minimum_num_visits_per_night'] | value_df['schedule.minimum_num_visits_per_night'].notna()
            else:
                value_df['schedule.minimum_num_visits_per_night'] = value_df['schedule.minimum_num_visits_per_night'].fillna(0)
                presence_df['schedule.minimum_num_visits_per_night'] = presence_df['schedule.minimum_num_visits_per_night'] | value_df['schedule.minimum_num_visits_per_night'].notna()
        else:
            if 'schedule.desired_num_visits_per_night' in value_df.columns:
                value_df['schedule.minimum_num_visits_per_night'] = value_df['schedule.desired_num_visits_per_night']
                presence_df['schedule.minimum_num_visits_per_night'] = value_df['schedule.minimum_num_visits_per_night'].notna()
            else:
                value_df['schedule.minimum_num_visits_per_night'] = 0
                presence_df['schedule.minimum_num_visits_per_night'] = True

        # New safety valve: if schedule.num_nights_per_semester == 1, set schedule.num_internight_cadence to 0
        if 'schedule.num_nights_per_semester' in value_df.columns and 'schedule.num_internight_cadence' in value_df.columns:
            mask = value_df['schedule.num_nights_per_semester'] == 1
            value_df.loc[mask, 'schedule.num_internight_cadence'] = 0
            presence_df.loc[mask, 'schedule.num_internight_cadence'] = True

        # Safety valve: if target.teff is missing, assign 0
        if 'target.teff' not in value_df.columns:
            value_df['target.t_eff'] = 0
            presence_df['target.t_eff'] = True

        # Safety valve: if target.gaia_id is missing or empty, assign 'NoGaiaName'
        if 'target.gaia_id' not in value_df.columns:
            value_df['target.gaia_id'] = 'NoGaiaName'
            presence_df['target.gaia_id'] = True
        else:
            value_df['target.gaia_id'] = value_df['target.gaia_id'].fillna('NoGaiaName').replace('', 'NoGaiaName')
            presence_df['target.gaia_id'] = presence_df['target.gaia_id'] | value_df['target.gaia_id'].notna()

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

    type_dict = {'observation.exposure_time':'Int64',
                'observation.num_exposures':'Int64',
                'schedule.num_internight_cadence':'Int64',
                'schedule.num_intranight_cadence':'Float64',
                'schedule.num_nights_per_semester':'Int64',
                'schedule.minimum_num_visits_per_night':'Int64',
                'schedule.desired_num_visits_per_night':'Int64',
                }

    df = df.copy()
    for col, dtype in type_dict.items():
        if col in df.columns:
            if dtype in ['Int64', 'Float64']:
                df[col] = pd.to_numeric(df[col], errors='coerce').astype(dtype)
            else:
                raise ValueError(f"Unsupported dtype: {dtype}. Only 'Int64' and 'Float64' are allowed.")
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
    

    columns_to_keep = [
        '_id',
        'metadata.semid',
        'target.target_name',
        'target.ra',
        'target.dec',
        'observation.exposure_time',
        'observation.num_exposures',
        'schedule.num_nights_per_semester',
        'schedule.num_internight_cadence',
        'schedule.desired_num_visits_per_night',
        'schedule.minimum_num_visits_per_night',
        'schedule.num_intranight_cadence',
        'schedule.weather_band',
        'target.gaia_id',
        'target.t_eff',
        'target.j_mag',
        'target.g_mag',
        'target.pm_ra',
        'target.pm_dec',
        'target.epoch',
    ]

    new_column_names = {
        '_id':'unique_id',
        'metadata.semid':'program_code',
        'target.target_name':'starname',
        'target.ra':'ra',
        'target.dec':'dec',
        'observation.exposure_time':'exptime',
        'observation.num_exposures':'n_exp',
        'schedule.num_nights_per_semester':'n_inter_max',
        'schedule.num_internight_cadence':'tau_inter',
        'schedule.desired_num_visits_per_night':'n_intra_max',
        'schedule.minimum_num_visits_per_night':'n_intra_min',
        'schedule.num_intranight_cadence':'tau_intra',
        'schedule.weather_band':'weather_band',
        'target.gaia_id':'gaia_id',
        'target.t_eff':'teff',
        'target.j_mag':'jmag',
        'target.g_mag':'gmag',
        'target.pm_ra':'pmra',
        'target.pm_dec':'pmdec',
        'target.epoch':'epoch',
    }

    trimmed_good = good_OBs_awarded[columns_to_keep].rename(columns=new_column_names)
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
