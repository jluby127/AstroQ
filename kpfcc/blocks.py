"""
Module for processing data from Keck Observatory's custom made Observing Block (OB) database.

Example usage:
    import ob_functions as ob
"""
import json
import requests
import pandas as pd
import numpy as np
import os
from astropy.time import Time

def refresh_local_data(semester):
    """
    Pull the latest database OBs down to local.

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
    
def get_request_sheet(OBs, awarded_programs, savepath):

    good_obs, bad_obs_values, bad_obs_hasFields = sort_good_bad(OBs, awarded_programs)
    good_obs.sort_values(by='program_code', inplace=True)
    good_obs.reset_index(inplace=True, drop=True)
    good_obs.to_csv(savepath, index=False)
    return good_obs, bad_obs_values, bad_obs_hasFields


def is_observation_complete(ob, entry):
    """
    Determine if a particular entry in the history of the OB is a completed observation. For the
    accounting and charging of time purposes. Note, this logic is a political decision and can be
    modified at any time.

    Args:
        ob (object) - the json object for a single OB. For OB structure info, see here: XYZ.
        entry (dict) - one element of the history array within the OB
    Returns:
        good (boolean) - True if passes unit tests
    """
    good = True
    # If you got at least half of your desired number of exposures, the observation is good
    if len(entry['exposure times']) < 0.5*ob['observation']['num_exposures']:
        good = False
    # If each of your exposures, except the last, was at least 75% complete in time, then your observation is good
    each_exp = []
    for i in range(len(entry['exposure times'])):
        if entry['exposure times'][i] < 0.75*ob['observation']['Exp Time']:
            each_exp.append(False)
        else:
            each_exp.append(True)
    if np.sum(each_exp) < len(entry['exposure times']) - 1:
        good = False
    return good

def get_past_history(ob):
    """
    Convert the history field of the OB into the format expected by the autoscheduler.

    Args:
        ob (object) - the json object for a single OB. For OB structure info, see here: XYZ.
    Returns:
        history (list) - See definition in comment below. Intended to be a value in a dictionary
                         where the keys are the star names.
    """
    unique = ob['metadata']['semid'] + "_" + str(ob['target']['target_name'])
    submitter = ob['metadata']['submitter']
    # Within the database_info_dict dictionary, each target's data is always in order:
    # element 0 = the calendar date of the most recent observation (HST date)
    # element 1 = a list of calendar dates with at least one observation
    # element 2 = a list of the quarters where the observation took place on the
    #             corresponding the nights in element 1.
    #             If multiple visits in one night, then this is quarter of the first visit.
    # element 3 = a list of the # of observations on each past night,
    #             corresponding the nights in element 1.
    observed_exps = []
    observed_times = []
    observed_dates = []
    observed_quarters = []
    for obs in ob['metadata']['history']:
        if is_observation_complete(ob, entry):
            n_exp = len(obs['exposure times'])
            observed_exps.append(n_exp)

            sumtime = np.sum(obs['exposure times'])
            observed_times.append(sumtime)

            tstamp = obs['timestamp'][:21]
            date = obs['timestamp'][:10]
            observed_dates.append(date)

            quarter = get_quarter_observed(tstamp, date, twilight_frame)
            observed_quarters.append(quarter)
        else:
            print("Observation for " + unique + " on night " + str(obs['timestamp'][:10]) + " submitted by " + submitter + " is not complete.")
    history = [observed_dates[0], observed_dates, observed_quarters, observed_exps]
    return history

def flatten(d, parent_key='', sep='.'):
    items = {}
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.update(flatten(v, new_key, sep=sep))
        else:
            items[new_key] = v
    return items

exception_fields = ['_id', 'del_flag', 'metadata.comment', 'metadata.details', 'metadata.history', 
                    'metadata.instruments', 'metadata.is_approved', 'metadata.last_modification', 
                    'metadata.ob_feasible', 'metadata.observer_name', 'metadata.state', 'metadata.status', 
                    'metadata.submitted', 'metadata.submitter', 'metadata.tags', 'observation.auto_exp_meter', 
                    'observation.auto_nd_filters', 'observation.cal_n_d_1', 'observation.cal_n_d_2', 
                    'observation.exp_meter_exp_time', 'observation.exp_meter_mode', 
                    'observation.exp_meter_threshold', 'observation.guide_here', 
                    'observation.object', 'observation.take_simulcal', 
                    'observation.trigger_ca_h_k', 'observation.trigger_green', 'observation.trigger_red', 
                    'schedule.accessibility_map', 'schedule.days_observable', 'schedule.fast_read_mode_requested',
                    'schedule.minimum_elevation', 'schedule.minimum_moon_separation',  
                    'schedule.rise_semester_day', 'schedule.scheduling_mode', 'schedule.sets_semester_day', 
                    'schedule.total_observations_requested', 'schedule.total_time_for_target', 
                    'schedule.total_time_for_target_hours', 'target.isNew', 'target.parallax', 'target.equinox', 'target.systemic_velocity',
                    'target.tic_id', 'target.two_mass_id', 'schedule.weather_band', 'target.catalog_comment']

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
                'schedule.num_visits_per_night':'Int64',
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
    
    good_OB_values = OB_values[pass_OBs_mask]
    good_OB_values.reset_index(inplace=True, drop='True')
    good_OBs = cast_columns(good_OB_values)

    good_OBs_awarded = good_OBs[good_OBs['metadata.semid'].isin(awarded_programs)]
    good_OBs_awarded.reset_index(inplace=True, drop='True')

    columns_to_keep = [
        '_id',
        'metadata.semid',
        'target.target_name',
        'target.ra_deg',
        'target.dec_deg',
        'observation.exposure_time',
        'observation.num_exposures',
        'schedule.num_nights_per_semester',
        'schedule.num_internight_cadence',
        'schedule.num_visits_per_night',
        'schedule.num_visits_per_night',
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
        'target.ra_deg':'ra',
        'target.dec_deg':'dec',
        'observation.exposure_time':'exptime',
        'observation.num_exposures':'n_exp',
        'schedule.num_nights_per_semester':'n_inter_max',
        'schedule.num_internight_cadence':'tau_inter',
        'schedule.num_visits_per_night':'n_intra_max',
        'schedule.num_visits_per_night':'n_intra_min',
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
    trimmed_good
    return trimmed_good, bad_OBs_values, bad_OBs_hasFields

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

# import smtplib
# from email.mime.multipart import MIMEMultipart
# from email.mime.text import MIMEText

# Email Configuration
SMTP_SERVER = 'smtp.example.com'
SMTP_PORT = 587
SENDER_EMAIL = 'youremail@example.com'
SENDER_PASSWORD = 'yourpassword'

def send_email(receiver_email, subject, body):
    msg = MIMEMultipart()
    msg['From'] = SENDER_EMAIL
    msg['To'] = receiver_email
    msg['Subject'] = subject
    msg.attach(MIMEText(body, 'plain'))

    try:
        # Connect to SMTP server
        server = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
        server.starttls()  # Secure the connection
        server.login(SENDER_EMAIL, SENDER_PASSWORD)
        text = msg.as_string()
        server.sendmail(SENDER_EMAIL, receiver_email, text)
        server.quit()
        print(f"Email sent to {receiver_email}")
    except Exception as e:
        print(f"Failed to send email to {receiver_email}: {e}")


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
