"""
Module for plotting.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import plotting_functions as ptf
"""

import os
import pandas as pd
import numpy as np
import json
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from collections import defaultdict
import pickle

import kpfcc.management as mn
import kpfcc.history as hs
import kpfcc.maps as mp

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

    def get_past(self, past, twilight_frame):
        """
        Gather the information about a star's past observation history this semester

        Args:
            starname (str): the name of the star in question

        Returns:
            observations_past (dict): a dictionary where keys are the dates of an observation
                                      and the values are number of observations taken on that night
        """
        starmask = past['star_id'] == self.starname
        star_obs_past = past[starmask]
        star_obs_past.sort_values(by='utctime', inplace=True)
        star_obs_past.reset_index(inplace=True)
        star_obs_past, unique_hst_dates_past, quarters_observed_past = \
                    hs.get_unique_nights(star_obs_past, twilight_frame)
        nobs_on_date_past = hs.get_nobs_on_night(star_obs_past, unique_hst_dates_past)
        observations_past = {}
        for i in range(len(unique_hst_dates_past)):
            observations_past[unique_hst_dates_past[i]] = nobs_on_date_past[i]
        self.observations_past = observations_past

    def get_future(self, forecast, all_dates_array):
        """
        Gather the information about a star's future forecast this semester

        Args:
            starname (str): the name of the star in question

        Returns:
            observations_future (dict): a dictionary where keys are the dates of an observation and
                                       the values are number of observations to be taken that night
        """
        observations_future = {}
        for i in range(len(forecast)):
            if self.starname in forecast[i]:
                observations_future[all_dates_array[i]] = (list(forecast[i]).count(self.starname)/self.slots_per_visit)/self.n_intra_max#*self.n_exp*self.n_intra_max
        self.observations_future = observations_future

    def get_map(self, forecast, semester_length):
        """
        Gather the information about a star's future forecast this semester

        Args:
            starname (str): the name of the star in question

        Returns:
            None
        """
        starmap = (forecast == self.starname).astype(int).to_numpy()
        starmap =  np.array(pad_rows_top(starmap, semester_length)).T
        self.starmap = starmap

def process_stars(manager):

    # Create a starmap of the times when we cannot observe due to twilight and allocation constraints
    # Used in the birdseye view plot to blackout the unavailble squares
    nulltime = manager.twilight_map_remaining_2D&manager.allocation_map_2D
    nulltime = 1 - nulltime
    nulltime = pad_rows_top(nulltime, manager.semester_length)
    nulltime = np.array(nulltime).T

    try:
        past = pd.read_csv(self.manager.past_database_file)
    except:
        past = pd.DataFrame(columns=['star_id', 'utctime'])

    starnames = manager.requests_frame['starname'].unique()
    programs = manager.requests_frame['program_code'].unique()

    # Make colors consistent for all stars in each program
    colors = sns.color_palette("deep", len(programs))
    rgb_strings = [f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})" for r, g, b in colors]
    program_colors_rgb_vals = dict(zip(programs, rgb_strings))

    # prepare the serialized output
    forecast = pd.read_csv(manager.future_forecast, header=None).astype(str)
    forecastT = np.array(forecast)

    all_stars = []
    for name in starnames:
        # Create a StarPlotter object for each request, fill and compute relavant information
        newstar = StarPlotter(name)
        newstar.get_stats(manager.requests_frame, manager.slot_size)
        newstar.get_past(past, manager.twilight_frame)
        newstar.get_future(forecastT, manager.all_dates_array)

        # Create birdseye map for each request
        newstar.get_map(forecast, manager.semester_length)

        # Create COF arrays for each request
        combined_set = set(list(newstar.observations_past.keys()) + list(newstar.observations_future.keys()))
        newstar.dates_observe = [newstar.n_exp*newstar.n_intra_max if date in combined_set else 0 for date in manager.all_dates_array]
        # don't assume that all future observations forecast for getting all desired n_intra_max
        for b in range(len(newstar.dates_observe)):
            if manager.all_dates_array[b] in list(newstar.observations_future.keys()):
                index = list(newstar.observations_future.keys()).index(manager.all_dates_array[b])
                newstar.dates_observe[b] *= list(newstar.observations_future.values())[index]
        newstar.cume_observe = np.cumsum(newstar.dates_observe)
        newstar.cume_observe_pct = np.round((np.cumsum(newstar.dates_observe)/newstar.total_observations_requested)*100.,3)

        # Create consistent colors across programs, and random colors for each star within programs
        newstar.program_color_rgb = program_colors_rgb_vals[newstar.program]
        newstar.star_color_rgb = rgb_strings[np.random.randint(0, len(rgb_strings)-1)]
        newstar.draw_lines = False

        all_stars.append(newstar)

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

def pad_rows_top(arr, target_rows):
    current_rows, cols = arr.shape
    missing_rows = target_rows - current_rows
    if missing_rows < 0:
        raise ValueError("Array already has more rows than target.")
    padding = np.ones((missing_rows, cols), dtype=arr.dtype)
    padded_arr = np.vstack((padding, arr))
    return padded_arr

def generate_birds_eye(manager, availablity, all_stars, filename=''):

    fig = go.Figure()
    fig.update_layout(width=1200, height=800)
    fig.add_trace(go.Heatmap(
        z=availablity,
        colorscale=[[0, 'rgba(0,0,0,0)'], [1, f"rgb(200, 200, 200)"]],
        zmin=0, zmax=1,
        opacity=1.0,
        showscale=False,
        name="Not On Sky",
        showlegend=False,
    ))

    for i in range(len(all_stars)):

        fig.add_trace(go.Heatmap(
            z=all_stars[i].starmap,
            colorscale=[[0, 'rgba(0,0,0,0)'], [1, all_stars[i].program_color_rgb]],
            zmin=0, zmax=1,
            opacity=1.0,
            showscale=False,
            name=all_stars[i].starname,
            hovertemplate='<b>' + str(all_stars[i].starname) +
                '</b><br><b>Date: %{x}</b><br><b>Slot: %{y}</b><br>Forecasted N_Obs: ' + \
                str(all_stars[i].total_observations_requested) + '<extra></extra>',
            showlegend=False,
        ))

        if all_stars[i].draw_lines:
            # Add connecting line for points with value 1
            points = np.argwhere(all_stars[i].starmap == 1)
            sorted_indices = np.argsort(points[:, 1])  # sort by x (column index)
            x_coords = points[sorted_indices, 1]
            y_coords = points[sorted_indices, 0]
            fig.add_trace(go.Scatter(
                x=x_coords,
                y=y_coords,
                mode='lines+markers',
                line=dict(color=all_stars[i].program_color_rgb, width=2),
                marker=dict(size=6, color=all_stars[i].program_color_rgb),
                name='Connected Points'
            ))

    # Add vertical grid lines every slot (x)
    for x in np.arange(0.5, all_stars[i].starmap.shape[1], 1):
        fig.add_shape(
            type="line",
            x0=x, x1=x,
            y0=0, y1=all_stars[i].starmap.shape[0] - 1,
            line=dict(color="lightgray", width=1),
            layer="below"
        )

    # # Add horizontal grid lines every slot (y)
    # for y in np.arange(0.5, all_stars[i].starmap.shape[0], 1):
    #     fig.add_shape(
    #         type="line",
    #         x0=0, x1=all_stars[i].starmap.shape[1] - 1,
    #         y0=y, y1=y,
    #         line=dict(color="lightgray", width=1),
    #         layer="below"
    #     )

    fig.update_layout(
        title="Birds Eye View Semester Schedule",
        yaxis_title="Slot in Night",
        xaxis_title="Night in Semester",
        xaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=np.arange(0, manager.semester_length, 21),
            tickmode='array',
            showgrid=False,  # Turn off native grid
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=np.arange(0, manager.n_slots_in_night, 12),
            tickmode='array',
            showgrid=False,  # Turn off native grid
        ),
        template="plotly_white",
        showlegend=True,
        legend=dict(
            font=dict(size=labelsize-10)
        )
    )
    if filename != '':
        fileout_path = manager.reports_directory + filename
        fig.write_html(fileout_path)
    return fig

def cof_builder(all_stars, manager, filename='', flag=False):

    fig = go.Figure()
    fig.update_layout(width=1200, height=800)
    burn_line = np.linspace(0, 100, len(manager.all_dates_array))
    for b in range(len(burn_line)):
        burn_line[b] = np.round(burn_line[b],2)
    fig.add_trace(go.Scatter(
        x=manager.all_dates_array,
        y=burn_line,
        mode='lines',
        line=dict(color='black', width=2, dash='dash'),
        name="Even Burn Rate",
        hovertemplate= 'Date: %{x}' + '<br>% Complete: %{y}'
    ))

    lines = []
    cume_observe = np.zeros(len(manager.all_dates_array))
    max_value = 0
    for i in range(len(all_stars)):
        fig.add_trace(go.Scatter(
            x=manager.all_dates_array,
            y=all_stars[i].cume_observe_pct,
            mode='lines',
            line=dict(color=all_stars[i].star_color_rgb, width=2),
            name=all_stars[i].starname,
            hovertemplate= 'Date: %{x}' + '<br>% Complete: %{y}' + '<br># Obs Requested: ' + \
                str(all_stars[i].total_observations_requested) + '<br>'
        ))
        lines.append(str(all_stars[i].starname) + "," + str(np.round(all_stars[i].cume_observe_pct[-1],2)))
        # Compute the COF data for all stars in the given program and plot total
        cume_observe += all_stars[i].cume_observe
        max_value += all_stars[i].total_observations_requested

    cume_observe_pct = (cume_observe / max_value) * 100
    fig.add_trace(go.Scatter(
        x=manager.all_dates_array,
        y=cume_observe_pct,
        mode='lines',
        line=dict(color=all_stars[0].program_color_rgb, width=2),
        name="Total",
        hovertemplate= 'Date: %{x}' + '<br>% Complete: %{y}' + '<br># Obs Requested: ' + \
            str(max_value) + '<br>'
    ))

    fig.add_vrect(
            x0=manager.current_day,
            x1=manager.current_day,
            annotation_text="Today",
            line_dash="dash",
            fillcolor=None,
            line_width=2,
            line_color='black',
            annotation_position="bottom left"
        )
    fig.update_layout(
        title="Cumulative Observation Function (COF)",
        xaxis_title="Calendar Date",
        yaxis_title="Request % Complete",
        showlegend=True,
        legend=dict(
            x=0.98,              # far right
            y=0.05,              # bottom
            xanchor='right',  # anchor legend box to the right side
            yanchor='bottom', # anchor legend box to the bottom
            bgcolor='rgba(255,255,255,0.7)',  # semi-transparent white background for readability
            bordercolor='black',
            borderwidth=1,
            font=dict(size=labelsize-18)
        ),
        xaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize-4)
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize-4)
        ),
    )
    if filename != '':
        fig.write_html(manager.reports_directory + filename)
        if flag:
            file = open(manager.reports_directory + "admin/" + manager.current_day + "/completion_Report.txt", "w")
            for l in lines:
                file.write(l + "\n")
            file.close()
    return fig


def write_star_objects(savepath, data):
    with open(savepath + "star_objects.pkl", "wb") as f:
        pickle.dump(data, f)
    print("Saved to " + savepath + "star_objects.pkl")

def read_star_objects(savepath):
    with open(savepath + "star_objects.pkl", "rb") as f:
        loaded_obj_list = pickle.load(f)
    return loaded_obj_list

def build_plot_file(manager):
    stars_in_program, programs_as_stars, nulltime = process_stars(manager)
    save_path = manager.reports_directory + "admin/" + manager.current_day + '/'
    os.makedirs(save_path, exist_ok = True)
    write_star_objects(save_path, [stars_in_program, programs_as_stars, nulltime])

def run_plot_suite(config_file):
    manager = mn.data_admin(config_file)
    manager.run_admin()
    build_plot_file(manager)

    print("Generating the full suite of plots. All under admin.")
    data = read_star_objects(manager.reports_directory + "admin/" + manager.current_day + '/')

    # build global plots
    all_stars_list = [obj for obj_list in data[0].values() for obj in obj_list]
    cof_builder(all_stars_list, manager, 'admin/' + manager.current_day + '/all_stars_COF.html')
    cof_builder(list(data[1].values()), manager, 'admin/' + manager.current_day + '/all_programs_COF.html', flag=True)
    generate_birds_eye(manager, data[2], all_stars_list, 'admin/' + manager.current_day + '/all_stars_birdseye.html')
    # generate_single_star_maps(manager, 'your star name here')
    # # build plots for each program
    # program_names = data[0].keys()
    # for p in program_names:
    #     cof_builder(data[0][p], manager, 'admin/'+ manager.current_day + "/" + p + "_COF.html")
    #     generate_birds_eye(manager, data[2], [data[1][p]], 'admin/'+ manager.current_day + "/" + p +"_birdseye.html")

def generate_single_star_maps(manager, starname):

    is_altaz, is_moon, is_night, is_inter, is_future, is_alloc = mp.mod_produce_ultimate_map(manager, starname)
    all_maps = [is_altaz, is_alloc, is_night, is_moon, is_inter, is_future]
    mapnames = ['is_altaz', 'is_alloc', 'is_night', 'is_moon', 'is_inter', 'is_future', ]

    forecast = pd.read_csv(manager.future_forecast, header=None).astype(str)
    forecastT = np.array(forecast)
    starmap = (forecast == starname).astype(int).to_numpy()
    starmap =  np.array(pad_rows_top(starmap, manager.semester_length)).T

    colors = sns.color_palette("deep", len(mapnames) + 1)
    rgb_strings = [f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})" for r, g, b in colors]
    fig = go.Figure()

    for i in range(len(all_maps)):
        # import pdb; pdb.set_trace()
        int_array = all_maps[i][0].astype(int).T
        int_array = 1 - int_array
        fig.add_trace(go.Heatmap(
            z=int_array,
            colorscale=[[0, 'rgba(0,0,0,0)'], [1, 'rgba(211, 211, 211, 1.0)']],
            zmin=0, zmax=1,
            opacity=0.5,
            showscale=False,
            name=mapnames[i],
            showlegend=True,
        ))

    fig.add_trace(go.Heatmap(
        z=starmap,
        colorscale=[[0, 'rgba(0,0,0,0)'], [1, 'rgba(255, 0, 0, 1.0)']],
        zmin=0, zmax=1,
        opacity=1.0,
        showscale=False,
        name=starname,
        showlegend=True,
    ))

    # Add vertical grid lines every slot (x)
    for x in np.arange(0.5, starmap.shape[1], 1):
        fig.add_shape(
            type="line",
            x0=x, x1=x,
            y0=0, y1=starmap.shape[0] - 1,
            line=dict(color="lightgray", width=1),
            layer="below"
        )

    # Add horizontal grid lines every slot (y)
    for y in np.arange(0.5, starmap.shape[0], 1):
        fig.add_shape(
            type="line",
            x0=0, x1=starmap.shape[1] - 1,
            y0=y, y1=y,
            line=dict(color="lightgray", width=1),
            layer="below"
        )

    fig.update_layout(
        title="Birds Eye View Semester Schedule",
        yaxis_title="Slot in Night",
        xaxis_title="Night in Semester",
        xaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=np.arange(0, manager.semester_length, 21),
            tickmode='array',
            showgrid=False,  # Turn off native grid
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=np.arange(0, manager.n_slots_in_night, 12),
            tickmode='array',
            showgrid=False,  # Turn off native grid
        ),
        template="plotly_white",
        showlegend=True,
        legend=dict(
            font=dict(size=labelsize-10)
        )
    )
    if filename != '':
        fileout_path = manager.reports_directory + filename
        fig.write_html(fileout_path)
    return fig
