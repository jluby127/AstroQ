"""
Module for plotting.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import plotting_functions as ptf
"""
import os
import warnings
warnings.filterwarnings('ignore')
import plotly.graph_objects as go
import plotly.express as px
from astropy.time import Time

import numpy as np
import pandas as pd

import kpfcc.management as mn
import kpfcc.history as hs


class StarTracker:
    """
    An object to easily pass around compute and pass around information on the past/future scheduling of stars.
    """
    def __init__(self, manager):
        """
        Define the StarTracker object. Set up the shared variables.

        Args:
            manager (obj): a data_admin object, see admin_functions.py

        Returns:
            None
        """
        self.manager = manager
        self.slot_size = manager.slot_size#*60 # convert to seconds
        self.forecast = manager.future_forecast
        self.first_forecast_dir = manager.folder_forecasts
        self.cadence_files_dir = manager.folder_cadences

        self.programs = self.manager.requests_frame['Program_Code'].unique()
        account_types = ["admin"]
        account_types.extend(self.programs)
        styles = ["birds_eyes", "cofs", "cadences"]
        self.create_folders(manager.reports_directory, account_types, [self.manager.current_day], styles)

        try:
            self.past = pd.read_csv(self.manager.past_database_file)
        except:
            self.past = pd.DataFrame(columns=['star_id', 'utctime'])
        try:
            self.nonqueue = pd.read_csv(self.manager.nonqueue_file)
        except:
            self.nonqueue = pd.DataFrame(columns=['Starname','Start','Stop','ProgramID','Comment'])
        forecast_by_day = []
        with open(self.forecast, 'r') as file:
            for line in file:
                splits = line.split(',')
                forecast_by_day.append(splits)
        self.forecast = forecast_by_day

        print("Data successfully loaded for calendar date: ", self.manager.current_day)

    def create_folders(self, main_path, account_types, dates, styles):
        """
        Check if folders exist, if not, create them

        Args:
            main_path (str): primary root path to save to
            account_types (str): the user acount types (admin, PI, etc)
            dates (str): the dates in question
            styles (str): the names of the plots to be generated, store them separately

        Returns:
            None
        """
        print("Checking/Creating folder structure")
        if not os.path.exists(main_path + 'observer/'):
            os.makedirs(main_path + 'observer/')
        if not os.path.exists(main_path + 'observer/' + self.manager.current_day + "/"):
            os.makedirs(main_path + 'observer/' + self.manager.current_day + "/")
        for account in account_types:
            for date in dates:
                base_dir = os.path.join(main_path, account, date)
                for style in styles:
                    style_dir = os.path.join(base_dir, style)
                    if not os.path.exists(style_dir):
                        os.makedirs(style_dir)

    def get_star_stats(self, starname):
        """
        Grab the observational stategy information for a given star

        Args:
            starname (str): the name of the star in question

        Returns:
            expected_nobs_per_night (int): how many exposures we expect to take
            total_observations_requested (int): sum of observational strategie values
            exposure_time (int): exposure time of single shot
            slots_per_night (int): number of slots required to complete all exposures in a night
            program (str): the program code
        """
        index = self.manager.requests_frame.loc[self.manager.requests_frame['Starname'] == starname].index
        program = str(self.manager.requests_frame['Program_Code'][index].values[0])
        exposure_time = int(self.manager.requests_frame['Nominal Exposure Time [s]'][index].values[0])
        slots_per_night = mn.compute_slots_required_for_exposure(
                        exposure_time, self.slot_size, False)*self.manager.requests_frame['Desired Visits per Night'][index].values[0]
        slots_per_visit = mn.compute_slots_required_for_exposure(
                        exposure_time, self.slot_size, False)
        expected_nobs_per_night = int(self.manager.requests_frame['# of Exposures per Visit'][index]) * \
                        int(self.manager.requests_frame['Desired Visits per Night'][index])
        total_observations_requested = expected_nobs_per_night * \
                        int(self.manager.requests_frame['# of Nights Per Semester'][index])

        return expected_nobs_per_night, total_observations_requested, exposure_time, \
                slots_per_night, program, slots_per_visit

    def get_star_past(self, starname):
        """
        Gather the information about a star's past observation history this semester

        Args:
            starname (str): the name of the star in question

        Returns:
            observations_past (dict): a dictionary where keys are the dates of an observation
                                      and the values are number of observations taken on that night
        """
        starmask = self.past['star_id'] == starname
        star_obs_past = self.past[starmask]
        star_obs_past.sort_values(by='utctime', inplace=True)
        star_obs_past.reset_index(inplace=True)
        star_obs_past, unique_hst_dates_past, quarters_observed_past = \
                    hs.get_unique_nights(star_obs_past, self.manager.twilight_frame)
        nobs_on_date_past = hs.get_nobs_on_night(star_obs_past, unique_hst_dates_past)
        observations_past = {}
        for i in range(len(unique_hst_dates_past)):
            observations_past[unique_hst_dates_past[i]] = nobs_on_date_past[i]
        return observations_past

    def get_star_forecast(self, starname):
        """
        Gather the information about a star's future forecast this semester

        Args:
            starname (str): the name of the star in question

        Returns:
            observations_future (dict): a dictionary where keys are the dates of an observation and
                                       the values are number of observations to be taken that night
        """
        index = self.manager.requests_frame.loc[self.manager.requests_frame['Starname'] == starname].index
        requested_nobs_per_night, total_observations_requested, exposure_time, slots_per_night, \
                    program, slots_per_visit = self.get_star_stats(starname)
        observations_future = {}
        for i in range(len(self.forecast)):
            if starname in self.forecast[i]:
                # observations_future[self.manager.all_dates_array[i]] = requested_nobs_per_night
                observations_future[self.manager.all_dates_array[i]] = list(self.forecast[i]).count(starname)/slots_per_visit
        return observations_future

    def write_star_first_forecast(self, program, starname):
        """
        After the first time running a forecast, record the information for later comparison

        Args:
            program (str): the program code
            starname (str): the name of the star in question

        Returns:
            None
        """
        observations_future = self.get_star_forecast(starname)
        forecast_frame = pd.DataFrame({"date":self.manager.all_dates_array,
                                       "pct_complete":observations_future})
        forecast_frame.to_csv(self.manager.reports_directory + program + "/first_forecast/" + starname + \
                            ".csv", index=False)

    def get_star_first_forecast(self, program, name):
        """
        Retrieve the information for later comparison

        Args:
            program (str): the program code
            name (str): the name of the star in question

        Returns:
            first_forecast (dataframe): contains information on the first forecast of the target
        """
        try:
            first_forecast = pd.read_csv(self.first_forecast_dir + program + "/" + name + ".csv")
        except:
            first_forecast = None
        return first_forecast
