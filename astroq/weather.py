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

def get_loss_stats(manager):
    """
    Simulate nights totally lost to weather usine historical data

    Args:
        manager (obj): a data_admin object

    Returns:
        loss_stats_remaining (array): each element is the percent of total loss for nights,
                                      starting with tomorrow and going to end of semester
    """
    historical_weather_data = pd.read_csv(os.path.join(manager.DATADIR,"maunakea_weather_loss_data.csv"))
    loss_stats_remaining = []
    for i, item in enumerate(manager.all_dates_array):
        ind = historical_weather_data.index[historical_weather_data['Date'] == \
            manager.all_dates_array[i][5:]].tolist()[0]
        loss_stats_remaining.append(historical_weather_data['% Total Loss'][ind])
    return loss_stats_remaining


def simulate_weather_losses(allocation_remaining, loss_stats, covariance=0.14, \
                            run_weather_loss=False, plot=False, outputdir=None):
    """
    Simulate nights totally lost to weather usine historical data

    Args:
        allocation_remaining (array): a 1D array of length n_nights_in_semester where 1's represent
                                      allocated night and 0's represent non-allocated night
        loss_stats (array): 1D array of length n_nights_in_semester where elements are the
                            percent of the time that night is totally lost to weather
        covariance (float): the added percent that tomorrow will be lost if today is lost
        run_weather_loss (boolean): a flag that turns on/off weather simulation entirely
        plot (boolean): a flag to visualize the allocation and weather loss
        outputdir (str): path to save the plot

    Returns:
        allocation_remaining_post_losses (array): 1's represent the night still is allocated
        weather_diff_remaining_2D (array): 1's represent the night still weathered
        weather_diff_remaining_1D (array): 1's represent the night still weathered
    """
    previous_day_was_lost = False
    allocation_remaining_post_losses = allocation_remaining.copy()
    counter = 0
    if run_weather_loss:
        days_lost = []
        # start at 1 because we never want tonight to be simulated as total loss
        for i in range(1, len(allocation_remaining_post_losses)):
            value_to_beat = loss_stats[i]
            if previous_day_was_lost:
                value_to_beat += covariance
            roll_the_dice = np.random.uniform(0.0,1.0)

            if roll_the_dice < value_to_beat:
                # the night is simulated a total loss
                allocation_remaining_post_losses[i] = [0,0,0,0]
                previous_day_was_lost = True
                counter += 1
                days_lost.append(1)
                if plot:
                    pt.axvline(i, color='r')
            else:
                previous_day_was_lost = False
                days_lost.append(0)
                if plot:
                    pt.axvline(i, color='k')
    else:
        logs.info('Pretending weather is always good!')
        days_lost = [0]*(len(allocation_remaining_post_losses)-1)

    weather_diff_remaining_2D = np.array(allocation_remaining) - \
                                np.array(allocation_remaining_post_losses)
    weather_diff_remaining_1D = weather_diff_remaining_2D.flatten()
    logs.info("Total nights simulated as weathered out: " + str(counter) + " of " + \
                str(len(allocation_remaining_post_losses)) + " nights remaining.") # info
    if plot:
        size=15
        pt.xlabel("Days in Semester from Current Day", fontsize=size)
        pt.tick_params(axis="both", labelsize=size)
        pt.savefig(outputdir + "weather_loss_visualization.png", dpi=200,
                                    bbox_inches='tight', facecolor='w')

    return allocation_remaining_post_losses, weather_diff_remaining_2D, weather_diff_remaining_1D, days_lost

def write_out_weather_stats(manager, days_lost, allocation_remaining):
    """
    Write out data on the results of the weather simulation. For comparing acrossing MCMC simulations

    Args:
        all_dates_dict (dict): keys are the dates of the semester, values are the day number of semester
        current_day (star): today's date in format YYYY-MM-DD
        days_lost (array): a binary 1D array designating if the night was lost. Note that the length
                            of this array is n_nights_in_semester-1 which is smaller than len(all_dates_dict)
        allocation_remaining (array): a binary 2D array of dates by quarters where 1 indicates
                                quarter allocated and 0 not. Note that the length
                                of this array is n_nights_in_semester-1 which is smaller
                                than len(all_dates_dict)
        output_directory (str): path to save the plot

    Returns:
        None
    """
    sim_results = []
    allocation_statii = []
    results = []
    filename_weather = manager.output_directory + 'Weather_Simulation_Results.csv'
    for a in range(len(manager.all_dates_array)):
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
            if np.sum(allocation_remaining[a - manager.all_dates_dict[manager.current_day]]) >= 1:
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
