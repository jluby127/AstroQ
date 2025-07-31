"""
Module for plotting.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import plotting_functions as ptf
"""

# Standard library imports
from collections import defaultdict
import os
import pickle

# Third-party imports
import numpy as np
import pandas as pd
import seaborn as sns
from astropy.time import Time, TimeDelta

# Local imports
import astroq.access as ac
import astroq.management as mn


labelsize = 38

class StarPlotter(object):
    """
        Define the StarPlotter class, which contains all information about a single star and its program.
        from which we can easily produce plots dynamically.
    """

    def __init__(self, name):
        self.starname = name

    def get_stats(self, requests_frame, slot_size):
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
        index = requests_frame.loc[requests_frame['starname'] == self.starname].index
        self.ra = float(requests_frame['ra'][index].values[0])
        self.dec = float(requests_frame['dec'][index].values[0])
        self.program = str(requests_frame['program_code'][index].values[0])
        self.exptime = int(requests_frame['exptime'][index].values[0])
        self.n_exp = int(requests_frame['n_exp'][index])
        self.n_intra_max = int(requests_frame['n_intra_max'][index])
        self.n_intra_min = int(requests_frame['n_intra_min'][index])
        self.tau_intra = int(requests_frame['tau_intra'][index])
        self.n_inter_max = int(requests_frame['n_inter_max'][index])
        self.tau_inter = int(requests_frame['tau_inter'][index])
        self.slots_per_night = mn.compute_slots_required_for_exposure(self.exptime, slot_size, False)*self.n_intra_max
        self.slots_per_visit = mn.compute_slots_required_for_exposure(self.exptime, slot_size, False)
        self.expected_nobs_per_night = self.n_exp * self.n_intra_max
        self.total_observations_requested = self.expected_nobs_per_night * self.n_inter_max

    def get_past(self, past):
        """
        Gather the information about a star's past observation history this semester.
        Args:
            past (DataFrame): DataFrame of past observations, must have 'target' and 'timestamp' columns.
        Returns:
            None. Sets self.observations_past as a dict: {date: n_obs}
        """
        # Filter to this star
        star_obs_past = past[past['target'] == self.starname]
        # Parse date from timestamp and group by date
        star_obs_past = star_obs_past.copy()
        star_obs_past['date'] = star_obs_past['timestamp'].str[:10]
        observations_past = star_obs_past.groupby('date').size().to_dict()
        self.observations_past = observations_past

    def get_future(self, forecaset_file, all_dates_array):
        """
        Process a DataFrame with columns ['r', 'd', 's'] to gather the future schedule for this star.

        Args:
            forecast_df (pd.DataFrame): DataFrame with columns ['r', 'd', 's']
            all_dates_array (list): List of all dates in the semester, indexed by 'd'
        """
        forecast_df = pd.read_csv(forecaset_file)
        # Only keep rows for this star
        star_rows = forecast_df[forecast_df['r'] == self.starname]
        # Count number of slots scheduled per night (d)
        observations_future = {}
        for d, group in star_rows.groupby('d'):
            # d may be int or str; ensure it's int for indexing
            d_idx = int(d)
            date = all_dates_array[d_idx]
            n_slots = len(group)
            # Optionally, normalize by slots_per_visit or n_intra_max as before
            observations_future[date] = n_slots
        self.observations_future = observations_future

    def get_map(self, semester_planner):
        """
        Build the starmap for this star using the new schedule format (future_forecast DataFrame with columns r, d, s).
        Only set starmap[d, s] = 1 if sched['r'] == self.starname.
        """
        forecast_df = pd.read_csv(semester_planner.future_forecast)
        n_nights = semester_planner.semester_length
        n_slots = int((24 * 60) / semester_planner.slot_size)
        starmap = np.zeros((n_nights, n_slots), dtype=int)
        for _, sched in forecast_df.iterrows():
            if sched['r'] == self.starname:
                d = int(sched['d'])
                s = int(sched['s'])
                starmap[d, s] = 1
        self.starmap = starmap.T
        # assert np.shape(starmap) == (semester_planner.semester_length, int((24 * 60) / semester_planner.slot_size)), np.shape(starmap)
        # assert np.shape(starmap) == (10, 4), np.shape(starmap)

def process_stars(semester_planner):

    # Create a starmap of the times when we cannot observe due to twilight and allocation constraints
    # Used in the birdseye view plot to blackout the unavailble squares

    # Create Access object with required parameters
    access_obj = ac.Access(
        semester_start_date=semester_planner.semester_start_date,
        semester_length=semester_planner.semester_length,
        slot_size=semester_planner.slot_size,
        observatory=semester_planner.observatory,
        current_day=semester_planner.current_day,
        all_dates_dict=semester_planner.all_dates_dict,
        all_dates_array=semester_planner.all_dates_array,
        n_nights_in_semester=semester_planner.n_nights_in_semester,
        custom_file=semester_planner.custom_file,
        allocation_file=semester_planner.allocation_file,
        past_history=semester_planner.past_history,
        today_starting_night=semester_planner.today_starting_night,
        slots_needed_for_exposure_dict=semester_planner.slots_needed_for_exposure_dict,
        run_weather_loss=semester_planner.run_weather_loss,
        output_directory=semester_planner.output_directory
    )
    access = access_obj.produce_ultimate_map(semester_planner.requests_frame)
    nulltime = access['is_alloc'][0]
    nulltime = 1 - nulltime
    nulltime = np.array(nulltime).T

    # Previously, there was a unique call to star names, every row of the request frame will be unique already when we switch to "id"
    starnames = semester_planner.requests_frame['starname'].unique()
    programs = semester_planner.requests_frame['program_code'].unique()

    # Make colors consistent for all stars in each program
    colors = sns.color_palette("deep", len(programs))
    rgb_strings = [f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})" for r, g, b in colors]
    program_colors_rgb_vals = dict(zip(programs, rgb_strings))

    all_stars = []
    i = 0 
    for i, row in semester_planner.requests_frame.iterrows():
        # Create a StarPlotter object for each request, fill and compute relavant information
        newstar = StarPlotter(row['starname'])
        newstar.get_map(semester_planner)
        newstar.get_stats(semester_planner.requests_frame, semester_planner.slot_size)
        print(list(semester_planner.past_history.keys()))
        if newstar.starname in list(semester_planner.past_history.keys()):
            newstar.observations_past = semester_planner.past_history[newstar.starname].n_obs_on_nights
        else:
            newstar.observations_past = {}
        newstar.get_future(semester_planner.future_forecast, semester_planner.all_dates_array)

        # Create COF arrays for each request
        combined_set = set(list(newstar.observations_past.keys()) + list(newstar.observations_future.keys()))
        newstar.dates_observe = [newstar.n_exp*newstar.n_intra_max if date in combined_set else 0 for date in semester_planner.all_dates_array]
        # don't assume that all future observations forecast for getting all desired n_intra_max
        for b in range(len(newstar.dates_observe)):
            if semester_planner.all_dates_array[b] in list(newstar.observations_future.keys()):
                index = list(newstar.observations_future.keys()).index(semester_planner.all_dates_array[b])
                newstar.dates_observe[b] *= list(newstar.observations_future.values())[index]
        newstar.cume_observe = np.cumsum(newstar.dates_observe)
        newstar.cume_observe_pct = np.round((np.cumsum(newstar.dates_observe)/newstar.total_observations_requested)*100.,3)

        # Create consistent colors across programs, and random colors for each star within programs
        newstar.program_color_rgb = program_colors_rgb_vals[newstar.program]
        newstar.star_color_rgb = rgb_strings[np.random.randint(0, len(rgb_strings)-1)]
        newstar.draw_lines = False
        newstar.maps_names = ['is_alloc', 'is_custom', 'is_altaz', 'is_moon', 'is_inter', 'is_future']
        newstar.maps = [access[name] for name in newstar.maps_names] 
        newstar.allow_mapview = True

        all_stars.append(newstar)
        i += 1

    # Now create StarPlotter objects for each program, as it were one star.
    # These will not have all the attributes, but we only need these for the admin COF plot
    # These StarPlotter objects cannot be used to create a birdseye plot, they don't have all attributes
    unique_programs = sorted(set(star.program for star in all_stars))
    programs_as_stars = {}
    for i in range(len(unique_programs)):

        prog_indices = [j for j, star in enumerate(all_stars) if star.program == unique_programs[i]]
        prog_objs = [star for j, star in enumerate(all_stars) if star.program == unique_programs[i]]

        # This is the quasi-StarPlotter object definition
        programmatic_star = StarPlotter(all_stars[prog_indices[0]].program)
        programmatic_star.program = all_stars[prog_indices[0]].program

        # Compute the COF data for all stars in the given program
        cume_observe = [all_stars[k].cume_observe for k in prog_indices]
        programmatic_star.cume_observe = np.sum([all_stars[k].cume_observe for k in prog_indices], axis=0)
        stars_stacked = np.vstack(cume_observe)
        summed_cumulative = np.sum(stars_stacked, axis=0)
        max_value = np.sum([all_stars[k].total_observations_requested for k in prog_indices])
        programmatic_star.cume_observe_pct = summed_cumulative / max_value * 100

        # Compute sum of starmaps
        super_map = np.zeros(np.shape(all_stars[prog_indices[0]].starmap))
        for m in range(len(prog_indices)):
            super_map += all_stars[prog_indices[m]].starmap
        programmatic_star.starmap = super_map
        programmatic_star.total_observations_requested = np.sum([all_stars[k].total_observations_requested for k in prog_indices])
        programmatic_star.draw_lines = False
        programmatic_star.allow_mapview = False

        # Set colors to match program color
        programmatic_star.program_color_rgb = all_stars[prog_indices[0]].program_color_rgb
        programmatic_star.star_color_rgb = all_stars[prog_indices[0]].program_color_rgb

        # Create list of "stars" objects which are really the programmatic overview
        programs_as_stars[all_stars[prog_indices[0]].program] = programmatic_star

    # Group stars into lists by program indexed by a dictionary
    program_dict = defaultdict(list)
    for obj in all_stars:
        program_dict[obj.program].append(obj)

    return program_dict, programs_as_stars, nulltime

def write_star_objects(savepath, data):
    with open(savepath + "star_objects.pkl", "wb") as f:
        pickle.dump(data, f)
    print("Plotting objects saved to " + savepath + "star_objects.pkl")

def read_star_objects(savepath):
    with open(savepath, "rb") as f:
        loaded_obj_list = pickle.load(f)
    return loaded_obj_list

def build_semester_webapp_pkl(semester_planner):
    stars_in_program, programs_as_stars, nulltime = process_stars(semester_planner)
    save_path = semester_planner.output_directory
    os.makedirs(save_path, exist_ok = True)
    write_star_objects(save_path, [stars_in_program, programs_as_stars, nulltime])