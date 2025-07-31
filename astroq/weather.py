"""
Module defining functions to perform weather loss simulations

Example usage:
    import weather as wh
"""

# Standard library imports
import logging
import os
import warnings

# Third-party imports
import numpy as np
import pandas as pd

# Suppress warnings
warnings.filterwarnings('ignore')

logs = logging.getLogger(__name__)

DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),'data')

def get_loss_stats(all_dates_array):
    """
    Get the loss probabilities for this semester's dates

    Args:
        all_dates_array (array): Array of date strings for the semester

    Returns:
        loss_stats_this_semester (array): each element is the percent of total loss for nights,
    """
    historical_weather_data = pd.read_csv(os.path.join(DATADIR,"maunakea_weather_loss_data.csv"))
    loss_stats_this_semester = []
    for i, item in enumerate(all_dates_array):
        ind = historical_weather_data.index[historical_weather_data['Date'] == \
            all_dates_array[i][5:]].tolist()[0]
        loss_stats_this_semester.append(historical_weather_data['% Total Loss'][ind])
    return loss_stats_this_semester

def simulate_weather_losses(semester_length, n_nights_in_semester, slot_size, loss_stats, covariance=0.14):
    """
    Simulate nights totally lost to weather using historical data

    Args:
        semester_length (int): Total number of nights in the semester
        n_nights_in_semester (int): Number of remaining nights in the semester
        slot_size (int): Size of time slots in minutes
        loss_stats (array): 1D array of semester_length where elements are the
                            percent of the time that night is totally lost to weather
        covariance (float): the added percent that tomorrow will be lost if today is lost

    Returns:
        is_clear (array): Trues represent clear nights, Falses represent weathered nights
    """
    previous_day_was_lost = False
    is_clear = np.ones((n_nights_in_semester, int((24*60)/ slot_size)),dtype=bool)
    for i in range(len(loss_stats)):
        if i < n_nights_in_semester:  # Only process remaining nights
            value_to_beat = loss_stats[i]
            if previous_day_was_lost:
                value_to_beat += covariance
            roll_the_dice = np.random.uniform(0.0,1.0)

            if roll_the_dice < value_to_beat:
                # the night is simulated a total loss
                is_clear[i] = np.zeros(is_clear.shape[1])  # Set all slots to False
                previous_day_was_lost = True
            else:
                previous_day_was_lost = False
    logs.info("Total nights simulated as weathered out: " + str(np.sum(~np.any(is_clear, axis=1))) + " of " + str(len(is_clear)) + " nights remaining.")

    return is_clear
