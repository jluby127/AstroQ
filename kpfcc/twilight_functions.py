"""
Module for performing functions relating to twilight time determination.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import twilight_functions as tf
"""
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta
import astroplan as apl

def generate_twilight_times(all_dates_array):
    """generate_twilight_times

    Precompute a dataframe of the morning/evening twilight times for each day in the semester

    Args:
        all_dates_array (list): the calendar dates of the semester.
                                Format: YYYY-MM-DD.
    Returns:
        twilight_frame (dataframe): the precomputed twilight times
    """
    keck = apl.Observer.at_site('W. M. Keck Observatory')

    twilight_frame = pd.DataFrame({'time_utc':all_dates_array})
    eighteen_deg_evening = []
    twelve_deg_evening = []
    six_deg_evening = []
    eighteen_deg_morning = []
    twelve_deg_morning = []
    six_deg_morning = []
    for day in twilight_frame.index.strftime(date_format='%Y-%m-%d').tolist():
        as_day = Time(day,format='iso',scale='utc')
        eighteen_deg_evening.append(keck.twilight_evening_astronomical(as_day,which='next'))
        twelve_deg_evening.append(keck.twilight_evening_nautical(as_day,which='next'))
        six_deg_evening.append(keck.twilight_evening_civil(as_day,which='next'))
        eighteen_deg_morning.append(keck.twilight_morning_astronomical(as_day,which='next'))
        twelve_deg_morning.append(keck.twilight_morning_nautical(as_day,which='next'))
        six_deg_morning.append(keck.twilight_morning_civil(as_day,which='next'))

    twilight_frame['18_evening'] = eighteen_deg_evening
    twilight_frame['12_evening'] = twelve_deg_evening
    twilight_frame['6_evening'] = six_deg_evening
    twilight_frame['18_morning'] = eighteen_deg_morning
    twilight_frame['12_morning'] = twelve_deg_morning
    twilight_frame['6_morning'] = six_deg_morning

    return twilight_frame

def get_available_slots_per_night(all_dates_array, twilight_frame, slot_size=5, start_time='03:30',
                                  n_hours_in_night=14):
    """get_available_slots_per_night

    Generates a 1D list of length equal to the number of nights in the semester where each element
    is the integer number of slots that are available to the scheduler in that night

    Args:
        all_dates_array (list): the calendar dates of the semester.
                            Format: YYYY-MM-DD.
        twilight_frame (dataframe): the precomputed twilight times
    Returns:
        None
    """
    available_slots_per_night = []
    for date in range(all_dates_array):
        edge = determine_twilight_edge(date, twilight_frame, slot_size, start_time,
                                       n_hours_in_night)
        available_slots_per_night.append(edge)
    return available_slots_per_night

def determine_twilight_edge(date, twilight_frame, slot_size, start_time='03:30',
                            n_hours_in_night=14):
    """determine_twilight_edge

    Determine how many slots within each night actually take place during night time.

    Args:
        - date (string): the calendar date of the semester. Format: YYYY-MM-DD.
        - twilight_frame (dataframe): the dataframe containing pre-computed
        - slot_size (int): the size of the a slot in minutes
        - start_time (str): the time representing the first slot of the night, format HH:MM
        - n_hours_in_night (int): number of hours in a night
    Returns:
        available_slots_in_night (int): the number of slots in night that are night time
    """
    ind = twilight_frame.index[twilight_frame['time_utc'] == date].tolist()
    twilight_evening_time_isot = twilight_frame.loc[ind,'12_evening'].values[0]
    twilight_evening_time_bjd = Time(twilight_evening_time_isot,format='jd')
    start_timestamp = date + "T" + start_time
    start_timestamp_time = Time(start_timestamp)
    n_slots = int((n_hours_in_night*60)/slot_size)

    bin_edges = []
    bin_edges.append(str(start_timestamp_time))
    for i in range(n_slots):
        start_timestamp_time += TimeDelta(slot_size/(24*60),format='jd')
        if twilight_evening_time_bjd < str(start_timestamp_time):
            edge = i
            break

    available_slots_in_night = n_slots - 2*edge
    return available_slots_in_night
