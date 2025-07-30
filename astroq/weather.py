"""
Module defining functions to perform weather loss simulations

Example usage:
    import weather as wh
"""
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as pt
import warnings
warnings.filterwarnings('ignore')
import logging
logs = logging.getLogger(__name__)

from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import astroplan as apl
import astropy.units as u

import astroq.access as ac

def get_loss_stats(manager):
    """
    Get the loss probabilities for this semester's dates

    Args:
        manager (obj): a data_admin object

    Returns:
        loss_stats_this_semester (array): each element is the percent of total loss for nights,
    """
    historical_weather_data = pd.read_csv(os.path.join(manager.DATADIR,"maunakea_weather_loss_data.csv"))
    loss_stats_this_semester = []
    for i, item in enumerate(manager.all_dates_array):
        ind = historical_weather_data.index[historical_weather_data['Date'] == \
            manager.all_dates_array[i][5:]].tolist()[0]
        loss_stats_this_semester.append(historical_weather_data['% Total Loss'][ind])
    return loss_stats_this_semester

def simulate_weather_losses(manager, loss_stats, covariance=0.14, run_weather_loss=False):
    """
    Simulate nights totally lost to weather usine historical data

    Args:
        loss_stats (array): 1D array of manager.semester_length where elements are the
                            percent of the time that night is totally lost to weather
        covariance (float): the added percent that tomorrow will be lost if today is lost
        run_weather_loss (boolean): a flag that turns on/off weather simulation entirely

    Returns:
        is_clear (array): Trues represent clear nights, Falses represent weathered nights
    """
    previous_day_was_lost = False
    is_clear = np.ones((manager.n_nights_in_semester, int((24*60)/ manager.slot_size)),dtype=bool)
    if run_weather_loss:
        for i in range(len(loss_stats)):
            if i != today_index:
                value_to_beat = loss_stats[i]
                if previous_day_was_lost:
                    value_to_beat += covariance
                roll_the_dice = np.random.uniform(0.0,1.0)

                if roll_the_dice < value_to_beat:
                    # the night is simulated a total loss
                    is_clear[i] = np.zeros(len(allocation_post_losses[i]))
                    previous_day_was_lost = True
                else:
                    previous_day_was_lost = False
        logs.info("Total nights simulated as weathered out: " + str(np.any(is_clear, axis=1).sum()) + " of " + str(len(is_clear)) + " nights remaining.")
    else:
        logs.info('Pretending weather is always good!')

    return is_clear

def write_out_weather_stats(manager, is_clear):
    """
    Write out data on the results of the weather simulation. For comparing acrossing MCMC simulations

    Args:
        is_clear (array): 2D array, Trues represent clear nights, Falses represent weathered nights

    Returns:
        None
    """
    sim_results = []
    allocation_statii = []
    results = []
    filename_weather = manager.output_directory + 'Weather_Simulation_Results.csv'
    for a in range(len(is_clear)):
        if a < manager.all_dates_dict[manager.current_day]:
            sim_result = 'Past'
            allocation_status = "Past"
            result = "Past"
        elif a == manager.all_dates_dict[manager.current_day]:
            sim_result = '???'
            allocation_status = "True"
            result = "???"
        else:
            # Extra -1 because the days_lost array does not include today
            if days_lost[a - manager.all_dates_dict[manager.current_day] - 1] == 1:
                sim_result = 'Poor'
            else:
                sim_result = 'Clear'
            if np.sum(allocation[a - manager.all_dates_dict[manager.current_day]]) >= 1:
                allocation_status = 'True'
            else:
                allocation_status = 'False'

            if allocation_status == 'True' and sim_result == 'Poor':
                result = 'Lost'
            elif allocation_status == 'True' and sim_result == 'Clear':
                result = 'Go'
            elif allocation_status == 'False' and sim_result == 'Clear':
                result = 'Miss'
            else:
                result = 'Eh'
        sim_results.append(sim_result)
        allocation_statii.append(allocation_status)
        results.append(result)
    weather_frame = pd.DataFrame({'Date':manager.all_dates_array, 'Allocated':allocation_statii,
                                    'SimWeather':sim_results, 'Designation':results})
    weather_frame.to_csv(filename_weather, index=False)
