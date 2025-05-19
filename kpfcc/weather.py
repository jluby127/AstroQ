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

from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import astroplan as apl
import astropy.units as u

import kpfcc.access as ac

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
        print('Pretending weather is always good!')
        days_lost = [0]*(len(allocation_remaining_post_losses)-1)

    weather_diff_remaining_2D = np.array(allocation_remaining) - \
                                np.array(allocation_remaining_post_losses)
    weather_diff_remaining_1D = weather_diff_remaining_2D.flatten()
    print("Total nights simulated as weathered out: " + str(counter) + " of " + \
                str(len(allocation_remaining_post_losses)) + " nights remaining.")
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

def compute_hours_observable_tonight(date, target, slot_size):
    """
    Determine how many hours a star is up on a given date

    Args:
        date (str): the calendar date to compute accessibilty in format 'YYYY-MM-DD'
        target (str): an astroplan FixedTarget object
        slot_size (int): the size of the slots in minutes
    Returns:
        n_hours_observable (int): the number of hours the star is observable
    """
    observability_matrix = ac.is_observable(date, target, slot_size)
    n_hours_observable = np.round((np.sum(observability_matrix)*slot_size)/60,1)
    return n_hours_observable

def compute_yearly_accessibility(manager):
    """
    Pre-compute how many hours a star is up each night of the year

    Args:
        manager (obj): a data_admin object

    Returns:
        backup_observability_frame (dataframe): contains information on each star's number of hours
                                                accessible on each night
    """
    star_frame = pd.read_csv(manager.backup_file)
    backup_observability_frame = pd.DataFrame({'Starname':star_frame['Starname']})
    for n in range(len(calendar)):
        single_night = []
        for s in range(len(star_frame)):
            coords = apy.coordinates.SkyCoord(star_frame['RA'][s] * u.deg,
                                              star_frame['Dec'][s] * u.deg, frame='icrs')
            target = apl.FixedTarget(name=star_frame['Starname'][s], coord=coords)
            time_up = compute_hours_observable_tonight(manager.all_dates_array[n], target, manager.slot_size)
            single_night.append(time_up)
        backup_observability_frame[calendar[n][5:]] = single_night
    return backup_observability_frame

def get_stars_for_tonight(manager, minimum_up_time=4):
    """
    Get a list of stars for backup tonight

    Args:
        manager (obj): a data_admin object
        minimum_up_time (int): minimum number of hours accessible to be considered observable

    Returns:
        starlist (array): contains the names of the stars that are observable tonight
    """
    backup_frame = pd.read_csv(manager.backup_file)
    backup_observability_frame = pd.read_csv(manager.backup_observability_file)
    up_mask = backup_observability_frame[manager.current_date] > minimum_up_time
    starlist = list(backup_observability_frame['Starname'][up_mask])

    name = []
    ra = []
    dec = []
    exposure_time = []
    n_shots = []
    total_time = 0
    for i, item in enumerate(starlist):
        ind = backup_frame.index[backup_frame['Starname'] == starlist[i]].tolist()
        if len(ind) != 0:
            ind = ind[0]
            name.append(backup_frame['Starname'][ind])
            ra.append(backup_frame['RA'][ind]*15.0)
            dec.append(backup_frame['Dec'][ind])
            exposure_time.append(backup_frame['Nominal Exposure Time [s]'][ind])
            if backup_frame['Nominal Exposure Time [s]'][ind] <= 150.0 and backup_frame['Nominal Exposure Time [s]'][ind] > 72.0:
                n_shot = 2
            elif backup_frame['Nominal Exposure Time [s]'][ind] <= 72.0 and backup_frame['Nominal Exposure Time [s]'][ind] > 45.0:
                n_shot = 3
            elif backup_frame['Nominal Exposure Time [s]'][ind] <= 45.0:
                n_shot = 5
            else:
                n_shot = 1
            n_shots.append(n_shot)
            total_time += backup_frame['Nominal Exposure Time [s]'][ind]*n_shot + 45*(n_shot-1)

    stars_for_tonight = pd.DataFrame({'Starname':name, "RA":ra, "Dec":dec,
                             'Exposure Time':exposure_time, 'Exposures Per Visit':n_shots,
                             'Visits In Night':[1]*len(name),
                             'Intra_Night_Cadence':[0]*len(name), 'Priority':[1]*len(name)})

    stars_for_tonight.sort_values(by='Exposure Time', ascending=False, inplace=True)
    stars_for_tonight.reset_index(inplace=True)

    print("There are " + str(len(stars_for_tonight)) + " available backup stars for tonight.")
    print("Amounting to total possible time added (not accounting for slew) of " + \
        str(np.round(total_time/3600,1)) + " hours.")
    return stars_for_tonight

def get_times_worth(backup_list_tonight, n_hours_needed):
    """
    Get a list of stars for backup tonight that fills a given amount of hours worth of time

    Args:
        backup_list_tonight (array): contains the stars from the wider backup pool selected as
                                     observable for tonight
        n_hours_needed (int): the number of hours of open shutter time to fill

    Returns:
        selected_stars_frame (array): contains the subset of stars from the wider backup pool
                                      selected as  observable for tonight
    """
    save_stars = list(backup_list_tonight['Starname'])
    selected_stars = []
    time_used = 0.0
    while time_used < n_hours_needed and len(save_stars) != 0:
        chosen_star = np.random.choice(save_stars)
        ind1 = save_stars.index(chosen_star)
        save_stars.pop(ind1)

        ind2 = backup_list_tonight.index[backup_list_tonight['Starname'] == chosen_star].tolist()[0]
        star_row = backup_list_tonight.loc[backup_list_tonight['Starname'] == chosen_star]
        selected_stars.append(star_row)
        time_used += (backup_list_tonight['Exposure Time'][ind2]*backup_list_tonight['Exposures Per Visit'][ind2] + \
                        45*(backup_list_tonight['Exposures Per Visit'][ind2]-1))/3600
    selected_stars_frame = pd.concat(selected_stars)
    selected_stars_frame.reset_index(inplace=True, drop=True)

    return selected_stars_frame
