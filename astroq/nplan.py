"""
Night Planning Module (nplan.py)

Module for night-level observation planning and optimization.
Uses the Target & Time Planner (TTP) to optimize nightly observation sequences.

Main Functions:
- run_ttp(manager): Optimize nightly observation sequences
- produce_bright_backups(manager): Create backup target lists for poor weather
- prepare_for_ttp(...): Prepare data for the TTP system

See https://github.com/lukehandley/ttp/tree/main for more info about the TTP
"""
import sys
import os
import pickle

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta
import gurobipy as gp
from gurobipy import GRB

sys.path.append('/Users/jack/Documents/github/ttp/')
import ttp.formatting as formatting
import ttp.telescope as telescope
import ttp.plotting as plotting
import ttp.model as model

import astroq.weather as wt
import astroq.io as io
import astroq.maps as mp

named_colors = ['blue', 'red', 'green', 'gold', 'maroon', 'gray', 'orange', 'magenta', 'purple']


class NightPlanner(object):
    """
    Night Planner for optimizing observation sequences within a single night.
    
    Uses the Target & Time Planner (TTP) to create optimal observation schedules
    and backup target lists for poor weather conditions.
    """
    
    def __init__(self, manager):
        """
        Initialize the Night Planner with a data manager.
        
        Args:
            manager: Data manager object containing configuration and data
        """
        self.manager = manager
        
    def run_ttp(self):
        """
        Produce the TTP solution given the results of the autoscheduler.
        
        Optimizes the nightly observation sequence for scheduled targets.
        """
        observers_path = self.manager.semester_directory + 'reports/observer/' + str(self.manager.current_day) + '/'
        check1 = os.path.isdir(observers_path)
        if not check1:
            os.makedirs(observers_path)

        observatory = telescope.Keck1()
        the_schedule = np.loadtxt(self.manager.upstream_path + "/outputs/" + str(self.manager.current_day) + "/raw_combined_semester_schedule_Round2.txt", delimiter=',', dtype=str)
        nightly_start_stop_times = pd.read_csv(self.manager.nightly_start_stop_times_file)

        day_in_semester = self.manager.all_dates_dict[self.manager.current_day]
        idx = nightly_start_stop_times[nightly_start_stop_times['Date'] == str(self.manager.current_day)].index[0]
        night_start_time = nightly_start_stop_times['Start'][idx]
        night_stop_time = nightly_start_stop_times['Stop'][idx]
        observation_start_time = Time(str(self.manager.current_day) + "T" + str(night_start_time),
            format='isot')
        observation_stop_time = Time(str(self.manager.current_day) + "T" + str(night_stop_time),
            format='isot')
        total_time = np.round((observation_stop_time.jd-observation_start_time.jd)*24,3)
        print("Time in Night for Observations: " + str(total_time) + " hours.")

        # Prepare relevant files
        round_two_requests = np.loadtxt(self.manager.upstream_path + "/outputs/" + str(self.manager.current_day) + '/Round2_Requests.txt', dtype=str, delimiter="\n")
        send_to_ttp = self.prepare_for_ttp(self.manager.requests_frame, the_schedule[day_in_semester], round_two_requests)

        filename = self.manager.upstream_path + "/outputs/" + str(self.manager.current_day) + '/Selected_' + str(self.manager.current_day) + ".txt"
        send_to_ttp.to_csv(filename, index=False)
        target_list = formatting.theTTP(filename)

        solution = model.TTPModel(observation_start_time, observation_stop_time, target_list,
                                    observatory, observers_path, runtime=10, optgap=0.01, useHighEl=False)

        gurobi_model_backup = solution.gurobi_model  # backup the attribute, probably don't need this
        del solution.gurobi_model                   # remove attribute so pickle works
        save_data = [solution]
        with open(self.manager.reports_directory + '/observer/' + self.manager.current_day + '/ttp_data.pkl', 'wb') as f:
            pickle.dump(save_data, f)

        plotting.writeStarList(solution.plotly, observation_start_time, self.manager.current_day,
                            outputdir=observers_path)
        plotting.plot_path_2D(solution,outputdir=observers_path)
        plotting.nightPlan(solution.plotly, self.manager.current_day, outputdir=observers_path)
        obs_and_times = pd.read_csv(observers_path + 'ObserveOrder_' + str(self.manager.current_day) + ".txt")
        io.write_starlist(self.manager.requests_frame, solution.plotly, observation_start_time, solution.extras,
                            round_two_requests, str(self.manager.current_day), observers_path)
        print("The optimal path through the sky for the selected stars is found. Clear skies!")

    def produce_bright_backups(self, nstars_max=100):
        """
        Produce backup target lists for poor weather conditions.
        
        Args:
            nstars_max (int): Maximum number of backup stars to include
        """
        backups_path = self.manager.reports_directory + '/observer/' + self.manager.current_day + "/backups/"
        check = os.path.isdir(backups_path)
        if not check:
            os.makedirs(backups_path)

        nightly_start_stop_times = pd.read_csv(self.manager.nightly_start_stop_times_file)
        # Determine time bounds of the night
        day_in_semester = self.manager.all_dates_dict[self.manager.current_day]
        idx = nightly_start_stop_times[nightly_start_stop_times['Date'] == str(self.manager.current_day)].index[0]
        night_start_time = nightly_start_stop_times['Start'][idx]
        observation_start_time = Time(str(self.manager.current_day) + "T" + str(night_start_time),
            format='isot')
        night_stop_time = nightly_start_stop_times['Stop'][idx]
        observation_stop_time = Time(str(self.manager.current_day) + "T" + str(night_stop_time),
            format='isot')
        diff_minutes = int(abs((observation_stop_time - observation_start_time).to('min').value))
        print("Minutes on sky: ", diff_minutes)

        backup_starlist = pd.read_csv(self.manager.backup_file)
        self.manager.requests_frame = backup_starlist
        available_indices = mp.produce_ultimate_map(self.manager, self.manager.requests_frame, running_backup_stars=True)
        slots_available_tonight_for_star = {k: len(v[0]) for k, v in available_indices.items()}
        stars_with_sufficient_availability_tonight = [k for k, v in slots_available_tonight_for_star.items() if v > int(0.25*int(diff_minutes/5))]

        self.manager.requests_frame = backup_starlist
        isTonight = backup_starlist['starname'].isin(stars_with_sufficient_availability_tonight)
        hasDR3name = backup_starlist['gaia_id'].str.startswith('Gaia DR2')
        pool_tonight = self.manager.requests_frame[isTonight&hasDR3name]
        pool_tonight = pool_tonight.sample(frac=1).reset_index(drop=True)
        pool_tonight = pool_tonight[:nstars_max]

        ready_for_ttp = self.prepare_for_ttp(pool_tonight, list(pool_tonight['starname']), [])
        ready_for_ttp.to_csv(backups_path + "selected_stars.csv", index=False)
        target_list = formatting.theTTP(backups_path + "selected_stars.csv")

        observatory = telescope.Keck1()
        solution_b = model.TTPModel(observation_start_time, observation_stop_time,
                                    target_list, observatory, backups_path,
                                    runtime=10, optgap=0.05)

        plotting.writeStarList(solution_b.plotly, observation_start_time, self.manager.current_day,
                                    outputdir = backups_path)
        plotting.plot_path_2D(solution_b, outputdir = backups_path)
        plotting.nightPlan(solution_b.plotly, self.manager.current_day, outputdir = backups_path)
        obs_and_times_b = pd.read_csv(backups_path + 'ObserveOrder_' + self.manager.current_day + ".txt")
        io.write_starlist(pool_tonight, solution_b.plotly, observation_start_time,
                            solution_b.extras, [], self.manager.current_day, backups_path, "backups")
        print("Bright backups script created.")

    def prepare_for_ttp(self, request_frame, night_plan, round_two_targets):
        """
        Prepare tonight's scheduled stars for their run through the TTP.

        Args:
            request_frame (dataframe): the pandas dataframe of PI requests
            night_plan (array): the n'th row of combined_semester_schedule array
            round_two_targets (array): a 1D list of the stars that were added in the bonus round

        Returns:
            to_ttp (dataframe): the data on the stars to be observed tonight, formatted for TTP
        """
        ignore = ['*', 'W', '', '*X', 'X']
        selected_stars = []
        for i, item in enumerate(night_plan):
            if night_plan[i] not in ignore and night_plan[i][:4] != "RM___":
                selected_stars.append(night_plan[i])
                ignore.append(night_plan[i])

        starnames = []
        ras = []
        decs = []
        exposure_times = []
        exposures_per_visit = []
        visits_in_night = []
        cadences = []
        priorities = []
        for j, item in enumerate(selected_stars):
            idx = request_frame.index[request_frame['starname']==str(selected_stars[j])][0]
            starnames.append(str(request_frame['starname'][idx]))
            ras.append(request_frame['ra'][idx])
            decs.append(request_frame['dec'][idx])
            exposure_times.append(int(request_frame['exptime'][idx]))
            exposures_per_visit.append(int(request_frame['n_exp'][idx]))
            visits_in_night.append(int(request_frame['n_intra_max'][idx]))
            cadences.append(int(request_frame['tau_intra'][idx]))
            # higher numbers are higher priorities, filler targets get low priority
            if str(selected_stars[j]) in round_two_targets:
                prior = 1
            else:
                prior = 10
            priorities.append(prior)
        to_ttp = pd.DataFrame({"Starname":starnames,"RA":ras,"Dec":decs,
                              "Exposure Time":exposure_times,
                              "Exposures Per Visit":exposures_per_visit,
                              "Visits In Night":visits_in_night, "Intra_Night_Cadence":cadences,
                              "Priority":priorities})
        return to_ttp


# Legacy function wrappers for backwards compatibility
def run_ttp(manager):
    """
    Legacy wrapper function for backwards compatibility.
    Use NightPlanner class for new code.
    """
    planner = NightPlanner(manager)
    planner.run_ttp()


def produce_bright_backups(manager, nstars_max=100):
    """
    Legacy wrapper function for backwards compatibility.
    Use NightPlanner class for new code.
    """
    planner = NightPlanner(manager)
    planner.produce_bright_backups(nstars_max)


def prepare_for_ttp(request_frame, night_plan, round_two_targets):
    """
    Legacy wrapper function for backwards compatibility.
    Use NightPlanner class for new code.
    """
    # This function is stateless, so we can just call it directly
    # Create a temporary planner instance
    temp_planner = NightPlanner(None)
    return temp_planner.prepare_for_ttp(request_frame, night_plan, round_two_targets)
