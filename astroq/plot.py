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
import plotly.io as pio
from astropy.time import Time
from astropy.time import TimeDelta

from collections import defaultdict
import pickle

import astroq.management as mn
import astroq.history as hs
import astroq.access as ac
import astroq.io as io
import ttp.plotting as plotting

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

    def get_map(self, manager):
        """
        Build the starmap for this star using the new schedule format (future_forecast DataFrame with columns r, d, s).
        Only set starmap[d, s] = 1 if sched['r'] == self.starname.
        """
        forecast_df = pd.read_csv(manager.future_forecast)
        n_nights = manager.semester_length
        n_slots = int((24 * 60) / manager.slot_size)
        starmap = np.zeros((n_nights, n_slots), dtype=int)
        for _, sched in forecast_df.iterrows():
            if sched['r'] == self.starname:
                d = int(sched['d'])
                s = int(sched['s'])
                starmap[d, s] = 1
        self.starmap = starmap.T
        # assert np.shape(starmap) == (manager.semester_length, int((24 * 60) / manager.slot_size)), np.shape(starmap)
        # assert np.shape(starmap) == (10, 4), np.shape(starmap)

def process_stars(manager):

    # Create a starmap of the times when we cannot observe due to twilight and allocation constraints
    # Used in the birdseye view plot to blackout the unavailble squares

    access = ac.Access().produce_ultimate_map(manager, manager.requests_frame)
    nulltime = access['is_alloc'][0]
    nulltime = 1 - nulltime
    # nulltime = pad_rows_top(nulltime, manager.semester_length)
    nulltime = np.array(nulltime).T

    # Previously, there was a unique call to star names, every row of the request frame should
    # unique already
    starnames = manager.requests_frame['starname'].unique()
    programs = manager.requests_frame['program_code'].unique()

    # Make colors consistent for all stars in each program
    colors = sns.color_palette("deep", len(programs))
    rgb_strings = [f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})" for r, g, b in colors]
    program_colors_rgb_vals = dict(zip(programs, rgb_strings))

    all_stars = []
    i = 0 
    for i, row in manager.requests_frame.iterrows():
        # Create a StarPlotter object for each request, fill and compute relavant information
        newstar = StarPlotter(row['starname'])
        newstar.get_map(manager)
        newstar.get_stats(manager.requests_frame, manager.slot_size)
        print(list(manager.past_history.keys()))
        if newstar.starname in list(manager.past_history.keys()):
            newstar.observations_past = manager.past_history[newstar.starname].n_obs_on_nights
        else:
            newstar.observations_past = {}
        newstar.get_future(manager.future_forecast, manager.all_dates_array)

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

def pad_rows_top(arr, target_rows):
    current_rows, cols = arr.shape
    missing_rows = target_rows - current_rows
    if missing_rows < 0:
        raise ValueError("Array already has more rows than target.")
    padding = np.ones((missing_rows, cols), dtype=arr.dtype)
    padded_arr = np.vstack((padding, arr))
    return padded_arr

gray = 'rgb(210,210,210)'
clear = 'rgba(0,0,0,0)'
clear = 'rgba(255,255,255,1)'

def generate_birds_eye(manager, availablity, all_stars, filename=''):

    fig = go.Figure()
    # fig.update_layout(width=1200, height=800, plot_bgcolor=clear, paper_bgcolor=clear)
    fig.update_layout(width=800, height=600, plot_bgcolor=clear, paper_bgcolor=clear)

    if len(all_stars) > 1 or all_stars[0].allow_mapview == False:
        fig.add_trace(go.Heatmap(
            z=availablity,
            colorscale=[[0, 'rgba(0,0,0,0)'], [1, gray]],
            zmin=0, zmax=1,
            opacity=1.0,
            showscale=False,
            name="Not On Sky",
            showlegend=False,
        ))
    else:
        colors = sns.color_palette("deep", len(all_stars[0].maps_names) + 1)
        rgb_strings = [f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})" for r, g, b in colors]
        for m in range(len(all_stars[0].maps_names)):
            fig.add_trace(go.Heatmap(
                z=all_stars[0].maps[m][0].astype(int).T,
                colorscale=[[0, 'rgba(0,0,0,0)'], [1, rgb_strings[m]]],
                zmin=0, zmax=1,
                opacity=1.0,
                showscale=False,
                name=all_stars[0].maps_names[m],
                showlegend=True,
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
            showlegend=True,
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

    # X-axis: ticks every 23 days, plus the last day
    x_tick_step = 23
    x_tickvals = list(range(0, manager.semester_length, x_tick_step))
    if (manager.semester_length - 1) not in x_tickvals:
        x_tickvals.append(manager.semester_length - 1)
    x_ticktext = [str(val + 1) for val in x_tickvals]

    # Y-axis: ticks every 2 hours, using slot_size
    n_slots = int(24 * 60 // manager.slot_size)
    slots_per_2hr = int(2 * 60 // manager.slot_size)
    y_tickvals = list(range(0, n_slots, slots_per_2hr))
    y_ticktext = []
    for slot in y_tickvals:
        total_minutes = slot * manager.slot_size
        hours = total_minutes // 60
        minutes = total_minutes % 60
        y_ticktext.append(f"{hours:02.0f}:{minutes:02.0f}")

    fig.update_layout(
        yaxis_title="Slot in Night",
        xaxis_title="Night in Semester",
        xaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=x_tickvals,
            ticktext=x_ticktext,
            tickmode='array',
            showgrid=False,
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=y_tickvals,
            ticktext=y_ticktext,
            tickmode='array',
            showgrid=False,
        ),
        template="plotly_white",
        showlegend=True,
        legend=dict(
            font=dict(size=labelsize-10)
        )
    )
    # import pdb; pdb.set_trace()
    if filename != '':

        fileout_path = manager.reports_directory + 'birdseye_static.png'

        # If it's a file, get its directory
        dir_path = fileout_path if os.path.isdir(fileout_path) else os.path.dirname(fileout_path)
        os.makedirs(dir_path, exist_ok=True) # Make sure the directory exists

        fig.write_html(fileout_path)
        pio.write_image(
            fig,
            manager.reports_directory + "birdseye_static.png",
            format="png",
            width=1200,
            height=800,
            scale=1,
        )
    return fig

def cof_builder(all_stars, manager, filename='', flag=False):

    fig = go.Figure()
    fig.update_layout(
        autosize=True,
        margin=dict(l=40, r=40, t=40, b=40),
        plot_bgcolor=gray, paper_bgcolor=clear
    )
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
        # title="Cumulative Observation Function (COF)",
        xaxis_title="Calendar Date",
        yaxis_title="Request % Complete",
        showlegend=True,
        legend=dict(
            x=0.98,
            y=0.05,
            xanchor='right',
            yanchor='bottom',
            bgcolor='rgba(255,255,255,0.7)',
            bordercolor='black',
            borderwidth=1,
            font=dict(size=labelsize-18)
        ),
        xaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize-4),
            showgrid=False,
            zeroline=False
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize-4),
            showgrid=False,
            zeroline=False
        ),
    )
    # import pdb; pdb.set_trace()
    if filename != '':

        fileout_path = manager.reports_directory + 'cof_static.png'

        # If it's a file, get its directory
        dir_path = fileout_path if os.path.isdir(fileout_path) else os.path.dirname(fileout_path)
        os.makedirs(dir_path, exist_ok=True) # Make sure the directory exists

        fig.write_html(fileout_path)
        if flag:
            file = open(manager.reports_directory + 'completion_Report.txt', "w")
            for l in lines:
                file.write(l + "\n")
            print("JUST wrote figure to {}".format(manager.reports_directory, filename))
            file.close()

        pio.write_image(
            fig,
            manager.reports_directory + "cof_static.png",
            format="png",
            width=1200,
            height=800,
            scale=1,
        )
    return fig

def generate_single_star_maps(manager, starname):

    access = ac.mod_produce_ultimate_map(manager, starname)
    all_maps = [is_altaz, is_alloc, is_night, is_moon, is_inter, is_future]
    mapnames = ['is_altaz', 'is_alloc', 'is_custom', 'is_moon', 'is_inter', 'is_future', ]

    forecast = pd.read_csv(manager.future_forecast)  # columns: id, d, s
    n_nights = manager.semester_length
    n_slots = int((24 * 60) / manager.slot_size)  
    starmap = np.zeros((n_nights, n_slots), dtype=int)
    star_rows = forecast[forecast['r'] == starname]
    for _, sched in star_rows.iterrows():
        d = int(sched['d'])
        s = int(sched['s'])
        starmap[d, s] = 1

    colors = sns.color_palette("deep", len(mapnames) + 1)
    rgb_strings = [f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})" for r, g, b in colors]
    fig = go.Figure()

    for i in range(len(all_maps)):
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
            showgrid=False,
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=np.arange(0, manager.n_slots_in_night, 12),
            tickmode='array',
            showgrid=False,
        ),
        template="plotly_white",
        showlegend=True,
        legend=dict(
            font=dict(size=labelsize-10)
        )
    )
    # if filename != '':
    #     fileout_path = manager.reports_directory + filename
    #     fig.write_html(fileout_path)
    return fig

def write_star_objects(savepath, data):
    with open(savepath + "star_objects.pkl", "wb") as f:
        pickle.dump(data, f)
    print("Plotting objects saved to " + savepath + "star_objects.pkl")

def read_star_objects(savepath):
    with open(savepath, "rb") as f:
        loaded_obj_list = pickle.load(f)
    return loaded_obj_list

def build_plot_file(manager):
    stars_in_program, programs_as_stars, nulltime = process_stars(manager)
    save_path = manager.reports_directory
    os.makedirs(save_path, exist_ok = True)
    write_star_objects(save_path, [stars_in_program, programs_as_stars, nulltime])

def build_semester_webapp_pkl(config_file):
    manager = mn.data_admin(config_file)
    manager.run_admin()
    build_plot_file(manager)

def build_static_plots(config_file):
    manager = mn.data_admin(config_file)
    manager.run_admin()
    data = read_star_objects(manager.reports_directory + "star_objects.pkl")
    all_stars_list = [obj for obj_list in data[0].values() for obj in obj_list]
    cof_builder(all_stars_list, manager, 'all_stars_COF.html')
    cof_builder(list(data[1].values()), manager, 'all_programs_COF.html', flag=True)
    generate_birds_eye(manager, data[2], all_stars_list, 'all_stars_birdseye.html')
    # generate_single_star_maps(manager, manager.requests_frame['starname'][0])

def save_interactive_observing_plan(observing_plan):
    from pathlib import Path

    # Generate the core HTML table (just the <table>...</table>)
    observing_plan_html = observing_plan.to_html(
        index=False,
        border=0,
        classes='display',
        escape=False  # Set to True if you have HTML tags in data you want escaped
    )

    # Full HTML page with DataTables integration
    full_html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Observing Plan</title>

    <!-- DataTables CSS -->
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">

    <!-- jQuery and DataTables JS -->
    <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

    <style>
        body {{
            font-family: Arial, sans-serif;
            padding: 2em;
        }}
        table {{
            width: 100%;
        }}
    </style>
</head>
<body>

<h2>Observing Plan</h2>

{observing_plan_html}

<script>
    $(document).ready(function () {{
        $('table.display').DataTable({{
            paging: true,
            searching: true,
            responsive: true,
            pageLength: 100,
            order: [[11, 'asc']]
        }});
    }});
</script>

</body>
</html>
"""

    # Write to file
    # Path(output_path).write_text(full_html, encoding='utf-8')
    return full_html
