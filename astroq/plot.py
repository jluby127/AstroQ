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
import base64
from io import BytesIO

# Third-party imports
import numpy as np
import pandas as pd
import seaborn as sns
import astropy as apy
import astropy.units as u
from astropy.coordinates import SkyCoord
import astroplan as apl
import imageio.v3 as iio
import matplotlib
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from matplotlib.figure import Figure
from plotly.subplots import make_subplots
from astropy.time import Time, TimeDelta

# Local imports
import astroq.access as ac
import astroq.io as io_mine
import astroq.nplan
DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),'data')

# Configure matplotlib for headless rendering
matplotlib.use("Agg")

# just used for color reproducibility
np.random.seed(24)

# Global variables from dynamic.py
gray = 'rgb(210,210,210)'
clear = 'rgba(255,255,255,1)'
labelsize = 38
slew_overhead = 180.
readout_overhead = 45.
hours_per_night = 12.

class StarPlotter(object):
    """
        Define the StarPlotter class, which contains all information about a single star and its program.
        from which we can easily produce plots dynamically.
    """

    def __init__(self, unique_id):
        self.unique_id = unique_id

    def get_stats(self, requests_frame, slot_size):
        """
        Grab the observational stategy information for a given star

        Args:


        Returns:
            expected_nobs_per_night (int): how many exposures we expect to take
            total_observations_requested (int): sum of observational strategie values
            exposure_time (int): exposure time of single shot
            slots_per_night (int): number of slots required to complete all exposures in a night
            program (str): the program code
        """
        index = requests_frame.loc[requests_frame['unique_id'] == self.unique_id].index
        self.starname = requests_frame['starname'][index].values[0]
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
        self.expected_nobs_per_night = self.n_exp * self.n_intra_max
        self.total_observations_requested = self.expected_nobs_per_night * self.n_inter_max
        self.total_requested_seconds = self.n_exp*self.n_intra_max*self.n_inter_max*self.exptime + readout_overhead*(self.n_exp-1) + slew_overhead*self.n_intra_max*self.n_inter_max
        self.total_requested_hours = self.total_requested_seconds / 3600
        self.total_requested_nights = self.total_requested_hours / hours_per_night   

    def get_past(self, past):
        """
        Gather the information about a star's past observation history this semester.
        Args:
            past (DataFrame): DataFrame of past observations, must have 'target' and 'timestamp' columns.
        Returns:
            None. Sets self.observations_past as a dict: {date: n_obs}
        """
        # Filter to this star
        star_obs_past = past[past['target'] == str(self.unique_id)]
        # Parse date from timestamp and group by date
        star_obs_past = star_obs_past.copy()
        star_obs_past['date'] = star_obs_past['timestamp'].str[:10]
        observations_past = star_obs_past.groupby('date').size().to_dict()
        self.observations_past = observations_past

    def get_future(self, forecast_file, all_dates_array):
        """
        Process a DataFrame with columns ['r', 'd', 's'] to gather the future schedule for this star.

        Args:
            forecast_df (pd.DataFrame): DataFrame with columns ['r', 'd', 's']
            all_dates_array (list): List of all dates in the semester, indexed by 'd'
        """
        forecast_df = pd.read_csv(forecast_file)
        forecast_df['r'] = forecast_df['r'].astype(str)

        # Only keep rows for this star
        star_rows = forecast_df[forecast_df['r'] == str(self.unique_id)]
        # Count number of slots scheduled per night (d)
        observations_future = {}
        for d, group in star_rows.groupby('d'):
            # d may be int or str; ensure it's int for indexing
            date = all_dates_array[int(d)]
            n_slots = len(group)
            # Optionally, normalize by slots_per_visit or n_intra_max as before
            observations_future[date] = n_slots*self.n_exp
        self.observations_future = observations_future

    def get_map(self, semester_planner):
        """
        Build the starmap for this star using the new schedule format (future_forecast DataFrame with columns r, d, s).
        Only set starmap[d, s] = 1 if sched['r'] == self.unique_id.
        """
        forecast_df = pd.read_csv(semester_planner.output_directory + semester_planner.future_forecast)
        n_nights = semester_planner.semester_length
        n_slots = int((24 * 60) / semester_planner.slot_size)
        starmap = np.zeros((n_nights, n_slots), dtype=int)
        for _, sched in forecast_df.iterrows():
            if str(sched['r']) == str(self.unique_id):
                d = int(sched['d'])
                s = int(sched['s'])
                starmap[d, s] = 1
                reserve_slots = semester_planner.slots_needed_for_exposure_dict[str(self.unique_id)]
                for r in range(1, reserve_slots):
                    starmap[d, s+r] = 1
        self.starmap = starmap.T

def process_stars(semester_planner):

    # Create a starmap of the times when we cannot observe due to twilight and allocation constraints
    # Used in the birdseye view plot to blackout the unavailble squares

    # Use the stored access record from the semester planner instead of recomputing
    access = semester_planner.access_record
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
        newstar = StarPlotter(row['unique_id'])
        newstar.get_map(semester_planner)
        newstar.get_stats(semester_planner.requests_frame, semester_planner.slot_size)
        if newstar.unique_id in list(semester_planner.past_history.keys()):
            newstar.observations_past = semester_planner.past_history[newstar.unique_id].n_visits_on_nights
        else:
            newstar.observations_past = {}
        newstar.get_future(semester_planner.output_directory + semester_planner.future_forecast, semester_planner.all_dates_array)

        # Create COF arrays for each request
        combined_set = set(list(newstar.observations_past.keys()) + list(newstar.observations_future.keys()))
        newstar.dates_observe = [newstar.observations_past[date] if date in newstar.observations_past.keys() else (newstar.n_intra_max*newstar.n_exp if date in combined_set else 0) for date in semester_planner.all_dates_array]
        # don't assume that all future observations forecast for getting all desired n_intra_max
        # for b in range(len(newstar.dates_observe)):
        #     if semester_planner.all_dates_array[b] in list(newstar.observations_future.keys()):
        #         index = list(newstar.observations_future.keys()).index(semester_planner.all_dates_array[b])
        #         newstar.dates_observe[b] *= list(newstar.observations_future.values())[index]
        newstar.cume_observe = np.cumsum(newstar.dates_observe)
        newstar.cume_observe_pct = np.round((np.cumsum(newstar.dates_observe)/newstar.total_observations_requested)*100.,3)

        # Create consistent colors across programs, and random colors for each star within programs
        newstar.program_color_rgb = program_colors_rgb_vals[newstar.program]
        # Ensure rgb_strings has at least one element before random selection
        if len(rgb_strings) > 1:
            newstar.star_color_rgb = rgb_strings[np.random.randint(0, len(rgb_strings)-1)]
        else:
            newstar.star_color_rgb = rgb_strings[0]
        newstar.draw_lines = False
        newstar.maps_names = ['is_alloc', 'is_custom', 'is_altaz', 'is_moon', 'is_inter', 'is_future', 'is_clear', 'is_observable_now']
        # Find the target index for this star in the access record
        target_idx = np.where(semester_planner.requests_frame['unique_id'] == newstar.unique_id)[0][0]
        # Extract the 2D slice for this specific target from each 3D map
        newstar.maps = {name: access[name][target_idx] for name in newstar.maps_names}
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
        programmatic_star.starname = all_stars[prog_indices[0]].program
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

def get_cof(semester_planner, all_stars):
    '''
    Return the html string for a plotly figure showing the COF for a selection of stars

    all_stars must be an array, even if only 1 element long. It is an array of StarPlotter objects, as defined in astroq.plot
    '''

    fig = go.Figure()
    fig.update_layout(plot_bgcolor=gray, paper_bgcolor=clear) #autosize=True,margin=dict(l=40, r=40, t=40, b=40),
    burn_line = np.linspace(0, 100, len(semester_planner.all_dates_array))
    for b in range(len(burn_line)):
        burn_line[b] = np.round(burn_line[b],2)
    fig.add_trace(go.Scatter(
        x=semester_planner.all_dates_array,
        y=burn_line,
        mode='lines',
        line=dict(color='black', width=2, dash='dash'),
        name="Even Burn Rate",
        hovertemplate= 'Date: %{x}' + '<br>% Complete: %{y}'
    ))

    lines = []
    cume_observe = np.zeros(len(semester_planner.all_dates_array))
    max_value = 0
    
    # First, compute the total COF data for all stars
    for i in range(len(all_stars)):
        cume_observe += all_stars[i].cume_observe
        max_value += all_stars[i].total_observations_requested
    
    cume_observe_pct = (cume_observe / max_value) * 100
    
    # Add the Total trace first (so it appears below other traces)
    fig.add_trace(go.Scatter(
        x=semester_planner.all_dates_array,
        y=cume_observe_pct,
        mode='lines',
        line=dict(color=all_stars[0].program_color_rgb, width=2),
        name="Total",
        hovertemplate= 'Date: %{x}' + '<br>% Complete: %{y}' + '<br># Obs Requested: ' + \
            str(max_value) + '<br>'
    ))
    
    # Then add individual star traces (so they appear above the Total trace)
    for i in range(len(all_stars)):
        fig.add_trace(go.Scatter(
            x=semester_planner.all_dates_array,
            y=all_stars[i].cume_observe_pct,
            mode='lines',
            line=dict(color=all_stars[i].star_color_rgb, width=2),
            name=all_stars[i].starname,
            hovertemplate= 'Date: %{x}' + '<br>% Complete: %{y}' + '<br># Obs Requested: ' + \
                str(all_stars[i].total_observations_requested) + '<br>'
        ))
        lines.append(str(all_stars[i].starname) + "," + str(np.round(all_stars[i].cume_observe_pct[-1],2)))

    fig.add_vrect(
            x0=semester_planner.current_day,
            x1=semester_planner.current_day,
            annotation_text="Today",
            line_dash="dash",
            fillcolor=None,
            line_width=2,
            line_color='black',
            annotation_position="bottom left"
        )
    
    # Calculate legend height based on number of traces
    num_traces = len(all_stars) + 2  # +2 for "Even Burn Rate" and "Total"
    legend_height = min(300, max(150, num_traces * 25))  # Between 150-300px, 25px per trace
    
    fig.update_layout(
        width=1400,
        height=1000,
        xaxis_title="Calendar Date",
        yaxis_title="Request % Complete",
        showlegend=True,
        legend=dict(
            orientation="h",
            x=0.5,
            y=-0.15,  # Position below plot
            xanchor="center",
            yanchor="top",
            bgcolor='rgba(255,255,255,0.7)',
            bordercolor='black',
            borderwidth=1,
            font=dict(size=labelsize-18),
            # Standardize legend size
            itemsizing='constant',  # All legend items same size
            itemwidth=30,  # Fixed width for legend items
            # Make legend more compact
            groupclick="toggleitem",  # Click group to toggle all items
            # Standardize legend dimensions
            tracegroupgap=5,  # Gap between trace groups
            traceorder="normal"  # Keep order as traces were added
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
        margin=dict(b=200, t=50)  # Bottom margin for legend below, minimal top margin
    )
    return fig

def get_birdseye(semester_planner, availablity, all_stars):
    '''
    Return the html string for a plotly figure showing the COF for a selection of stars

    all_stars must be an array, even if only 1 element long. It is an array of StarPlotter objects, as defined in astroq.plot
    availability is a 2D array of N_slots by N_nights, binary 1/0, it is the intersection of is_alloc and is_night
    '''

    fig = go.Figure()
    # fig.update_layout(width=1200, height=800, plot_bgcolor=clear, paper_bgcolor=clear)
    fig.update_layout(plot_bgcolor=clear, paper_bgcolor=clear)

    # when multiple StarPlotter obects are submitted or a programmatic StarPlotter object,
    # show the grayed out slots from the intersection of is_alloc and is_night
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
    # when just one StarPlotter object is submitted, show the overlay of all maps
    else:
        colors = sns.color_palette("deep", len(all_stars[0].maps_names) + 1)
        rgb_strings = [f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})" for r, g, b in colors]
        for m in range(len(all_stars[0].maps_names)):
            # Skip the is_observable_now map
            if all_stars[0].maps_names[m] == 'is_observable_now':
                continue
            map_name = all_stars[0].maps_names[m]
            z_data = 1-all_stars[0].maps[map_name].astype(int).T  # Invert all other maps
            
            fig.add_trace(go.Heatmap(
                z=z_data,
                colorscale=[[0, 'rgba(0,0,0,0)'], [1, gray]],
                zmin=0, zmax=1,
                opacity=1.0,
                showscale=False,
                name=all_stars[0].maps_names[m],
                showlegend=True,
            ))

    for i in range(len(all_stars)):

        fig.add_trace(go.Heatmap(
            z=all_stars[i].starmap,
            colorscale=[[0, 'rgba(0,0,0,0)'], [1, all_stars[i].star_color_rgb]],
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
                line=dict(color=all_stars[i].star_color_rgb, width=2),
                marker=dict(size=6, color=all_stars[i].starcolor_rgb),
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
    
    # Add vertical dashed line denoting "today"
    fig.add_vrect(
        x0=semester_planner.today_starting_night-1, #The minus one is just for aesthetic purposes.
        x1=semester_planner.today_starting_night-1,
        annotation_text="Today",
        line_dash="dash",
        fillcolor=None,
        line_width=2,
        line_color='black',
        annotation_position="bottom left"
    )
     # X-axis: ticks every 23 days, plus the last day
    x_tick_step = 23
    x_tickvals = list(range(0, semester_planner.semester_length, x_tick_step))
    if (semester_planner.semester_length - 1) not in x_tickvals:
        x_tickvals.append(semester_planner.semester_length - 1)
    x_ticktext = [str(val + 1) for val in x_tickvals]

    # Y-axis: ticks every 2 hours, using slot_size
    n_slots = int(24 * 60 // semester_planner.slot_size)
    slots_per_2hr = int(2 * 60 // semester_planner.slot_size)
    y_tickvals = list(range(0, n_slots, slots_per_2hr))
    y_ticktext = []
    for slot in y_tickvals:
        total_minutes = slot * semester_planner.slot_size
        hours = total_minutes // 60
        minutes = total_minutes % 60
        y_ticktext.append(f"{hours:02.0f}:{minutes:02.0f}")

    # Calculate legend height based on number of traces
    num_traces = len(all_stars) + (1 if len(all_stars) > 1 or all_stars[0].allow_mapview == False else len([m for m in all_stars[0].maps_names if m != 'is_observable_now']))
    legend_height = min(300, max(150, num_traces * 25))  # Between 150-300px, 25px per trace
    
    fig.update_layout(
        width=1400,
        height=1000,
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
            orientation="h",
            x=0.5,
            y=-0.15,  # Position below plot
            xanchor="center",
            yanchor="top",
            font=dict(size=labelsize-18),
            bgcolor='rgba(255,255,255,0.7)',
            bordercolor='black',
            borderwidth=1,
            # Standardize legend size
            itemsizing='constant',  # All legend items same size
            itemwidth=30,  # Fixed width for legend items
            # Make legend more compact
            groupclick="toggleitem",  # Click group to toggle all items
            # Standardize legend dimensions
            tracegroupgap=5,  # Gap between trace groups
            traceorder="normal"  # Keep order as traces were added
        ),
        margin=dict(b=200, t=50)  # Bottom margin for legend below, minimal top margin
    )
    return fig

def get_tau_inter_line(semester_planner, all_stars, use_program_colors=False):
    """
    Create an interactive scatter plot grouped by star name, using provided colors for each point.
    Points from the same star will share a legend entry and color.

    Parameters:
        semester_planner: the semester planner object
        all_stars (list): array of StarPlotter objects
        use_program_colors (bool): If True, use program_color_rgb; if False, use star_color_rgb (default: False)

    Returns:
        plotly.graph_objects.Figure
    """

    request_tau_inter = []
    onsky_tau_inter = []
    starnames = []
    programs = []
    colors = []
    for starobj in all_stars:
        onsky_diffs = list(np.diff(np.where(np.diff(starobj.cume_observe) > 0)[0]))
        onsky_tau_inter.extend(onsky_diffs)
        request_tau_inter.extend([starobj.tau_inter] * len(onsky_diffs))
        starnames.extend([starobj.starname] * len(onsky_diffs))
        programs.extend([starobj.program] * len(onsky_diffs))
        # Choose color based on flag
        if use_program_colors:
            colors.extend([starobj.program_color_rgb] * len(onsky_diffs))
        else:
            colors.extend([starobj.star_color_rgb] * len(onsky_diffs))

    all_request_tau_inters = np.array(request_tau_inter)
    all_onsky_tau_inters = np.array(onsky_tau_inter)
    all_starnames = np.array(starnames)
    all_programs = np.array(programs)
    all_colors = np.array(colors)

    fig = go.Figure()

    # Build map from program to point indices
    program_to_indices = {}
    for i, prog in enumerate(all_programs):
        program_to_indices.setdefault(prog, []).append(i)

    #Create one trace per star (grouped by starname)
    maxyvals = []
    # Build map from starname to point indices
    starname_to_indices = {}
    for i, starname in enumerate(all_starnames):
        starname_to_indices.setdefault(starname, []).append(i)
    
    for starname, indices in starname_to_indices.items():
        x_vals = [all_request_tau_inters[i] for i in indices]
        y_vals = [all_onsky_tau_inters[i] for i in indices]
        text_vals = [f"{all_starnames[i]} in {all_programs[i]}" for i in indices]
        color_vals = [all_colors[i] for i in indices]
        maxyvals.append(max(y_vals))
        fig.add_trace(go.Scatter(
            x=x_vals,
            y=y_vals,
            mode='markers',
            name=starname,  # Use star name for legend
            marker=dict(size=10, color=color_vals),
            text=text_vals,
            hovertemplate="%{text}<br>X: %{x}<br>Y: %{y}<extra></extra>"
        ))

    # Add 1-to-1 line
    min_val = 0
    max_val = max(maxyvals)
    fig.add_trace(go.Scatter(
        x=[min_val, max_val],
        y=[min_val, max_val],
        mode='lines',
        line=dict(color='black', dash='dash'),
        name='1-to-1 line',
        showlegend=True
    ))

    fig.update_layout(
        width=1200,
        height=800,
        xaxis_title="Requested Minimum Inter-Night Cadence",
        yaxis_title="On Sky Inter-Night Cadence",
        template='plotly_white',
        xaxis=dict(
            type="log",
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize-4),
            showgrid=True,
            gridcolor='lightgray',
            gridwidth=0.5,
            tickmode='array',
            tickvals=[1, 10, 100],
            ticktext=['1', '10', '100'],
            range=[np.log10(0.5), np.log10(180)]  # Set range from 0.5 to 180 in log scale
        ),
        yaxis=dict(
            type="log",
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize-4),
            showgrid=True,
            gridcolor='lightgray',
            gridwidth=0.5,
            tickmode='array',
            tickvals=[1, 10, 100],
            ticktext=['1', '10', '100'],
            range=[np.log10(0.5), np.log10(180)]  # Set range from 0.5 to 180 in log scale
        )
    )
    return fig

def get_timepie(semester_planner, all_stars, use_program_colors=False):
    """
    Create an pie chart of the time used vs forecasted vs available

    Parameters:
        semester_planner: the semester planner object
        all_stars (list): array of StarPlotter objects
        use_program_colors (bool): If True, use program_color_rgb; if False, use star_color_rgb (default: False)

    Returns:
        plotly.graph_objects.Figure
    """
    programmatics = pd.read_csv(semester_planner.semester_directory + 'programmatics.csv')

    # Accumulate total times across all stars
    total_past = 0
    total_future = 0
    total_incomplete = 0
    total_requested_hours = 0
    
    programs_used = []
    for starobj in all_stars:
        total_past += sum(starobj.observations_past.values()) * starobj.exptime + len(starobj.observations_past.keys()) * slew_overhead
        total_future += sum(starobj.observations_future.values()) * starobj.exptime + len(starobj.observations_future.keys()) * slew_overhead
        total_requested_hours += starobj.total_requested_hours
        programs_used.append(starobj.program)
    
    # Convert to hours for better readability
    total_past_hours = total_past / 3600
    total_future_hours = total_future / 3600
    total_incomplete_hours = total_requested_hours - total_past_hours - total_future_hours

    if len(programs_used) > 1:
        total_allocated_hours = programmatics[programmatics['program'].isin(programs_used)]['hours'].sum()
        total_allocated_nights = programmatics[programmatics['program'].isin(programs_used)]['nights'].sum()
    else:
        total_allocated_hours = programmatics[programmatics['program'] == programs_used[0]]['hours'].sum()
        total_allocated_nights = programmatics[programmatics['program'] == programs_used[0]]['nights'].sum()

    # Create pie chart data
    labels = ['Completed', 'Scheduled', 'Incomplete', "Unused"]
    values = [total_past_hours, total_future_hours, total_incomplete_hours, total_allocated_hours-total_requested_hours]
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#000000']  # Blue, Purple, Orange, Black
    
    # Create the pie chart
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        marker=dict(colors=colors),
        textinfo='label+percent+value',
        texttemplate='%{label}<br>%{value:.1f} hrs<br>(%{percent})',
        hovertemplate='<b>%{label}</b><br>%{value:.2f} hours<br>%{percent}<extra></extra>'
    )])
    
    # Adjust margin if there's a warning to display
    top_margin = 150 if total_requested_hours > total_allocated_hours else 100
    
    fig.update_layout(
        title_text=f'<b>Total Requested:</b> {total_requested_hours:.1f} hours ≈ {total_requested_hours/hours_per_night:.1f} nights<br><b>Total Allocated:</b> {total_allocated_hours:.1f} hours = {total_allocated_nights:.1f} nights',
        template='plotly_white',
        showlegend=False,
        height=600,
        width=700,
        margin=dict(t=top_margin, b=50, l=50, r=50)
    )
    
    # Add warning annotation if requested time exceeds allocated time
    if total_requested_hours > total_allocated_hours:
        fig.add_annotation(
            text='<b>You have requested more time than you are allocated.</b>',
            xref='paper', yref='paper',
            x=0.5, y=1.05,
            showarrow=False,
            font=dict(size=18, color='red'),
            xanchor='center',
            yanchor='middle'
        )
    
    return fig


def compute_seasonality(semester_planner, starnames, ras, decs):
    """
    Compute seasonality using Access object's produce_ultimate_map function

    Args:
        semester_planner (SemesterPlanner): the semester planner object containing configuration
        starnames (list): list of star names
        ras (array): right ascension values in degrees
        decs (array): declination values in degrees
    Returns:
        available_nights_onsky (list): number of observable nights for each target

    """
    # Create a temporary requests frame from the input parameters
    temp_requests_frame = pd.DataFrame({
        'starname': starnames,
        'unique_id': starnames,
        'ra': ras,
        'dec': decs,
        'exptime': [300] * len(starnames),  # Default values
        'n_exp': [1] * len(starnames),
        'n_intra_max': [1] * len(starnames),
        'n_intra_min': [1] * len(starnames),
        'n_inter_max': [1] * len(starnames),
        'tau_inter': [1] * len(starnames),
        'tau_intra': [1] * len(starnames),
        'minimum_elevation': [30.] * len(starnames),
        'minimum_moon_separation': [30.] * len(starnames)
    })
    
    # Build or get the twilight allocation file
    twilight_allocation_file = ac.build_twilight_allocation_file(semester_planner)
    
    # Temporarily override the allocation file path in the access object
    original_allocation_file = semester_planner.access_obj.allocation_file
    semester_planner.access_obj.allocation_file = twilight_allocation_file
   
    # Create dummy allocation for if the try statement fails.
    is_alloc = np.ones((len(starnames), semester_planner.semester_length, semester_planner.n_slots_in_night), dtype=bool)
    try:
        # Use Access object to produce the ultimate map with our custom requests frame
        access_record = semester_planner.access_obj.produce_ultimate_map(temp_requests_frame, running_backup_stars=True)
        is_alloc = access_record.is_alloc
    finally:
        # Restore the original allocation file path
        semester_planner.access_obj.allocation_file = original_allocation_file
    
    # Extract is_altaz and is_moon arrays
    is_altaz = access_record.is_altaz
    is_moon = access_record.is_moon
    
    ntargets = len(starnames)
    nnights = semester_planner.n_nights_in_semester
    nslots = semester_planner.n_slots_in_night
    
    # Create the combined observability mask
    is_observable_now = np.logical_and.reduce([
        is_altaz,
        is_moon,
        is_alloc
    ])

    # specify indeces of 3D observability array
    itarget, inight, islot = np.mgrid[:ntargets,:nnights,:nslots]

    # define flat table to access maps
    df = pd.DataFrame(
        {'itarget':itarget.flatten(),
         'inight':inight.flatten(),
         'islot':islot.flatten()}
    )
    available_nights_onsky = []
    for itarget in range(ntargets):
        onskycount = 0
        for inight in range(nnights):
            temp = list(islot[itarget,inight,is_observable_now[itarget,inight,:]])
            if len(temp) > 0:
                onskycount += 1

        available_nights_onsky.append(onskycount)

    return available_nights_onsky

def get_football(semester_planner, all_stars, use_program_colors=False):
    """
    "Football plot"
    Interactive sky map with static heatmap background and interactive star points.

    Parameters:
        semester_planner: the semester planner object
        all_stars (list): array of StarPlotter objects
        use_program_colors (bool): If True, use program_color_rgb; if False, use star_color_rgb (default: False)

    Returns:
        plotly.graph_objects.Figure
    """

    starnames = [all_stars[r].starname for r in range(len(all_stars))]
    programs = [all_stars[r].program for r in range(len(all_stars))]
    ras = [all_stars[r].ra for r in range(len(all_stars))]
    decs = [all_stars[r].dec for r in range(len(all_stars))]
    # Choose color based on flag
    if use_program_colors:
        colors = [all_stars[r].program_color_rgb for r in range(len(all_stars))]
    else:
        colors = [all_stars[r].star_color_rgb for r in range(len(all_stars))]
    program_frame = pd.DataFrame({"starname":starnames, "program_code":programs, "color":colors, "ra":ras, "dec":decs})

    n_ra = 90
    ras = np.linspace(0,360,n_ra)
    n_dec = 90
    start_deg = -90
    stop_deg = 90
    # Split number of points proportionally between negative and positive ranges
    neg_points = int(n_dec * (0 - start_deg) / (stop_deg - start_deg))  # points from -30 to 0
    pos_points = n_dec - neg_points  # points from 0 to 90
    # Negative part: from -30 to 0 deg
    neg_deg = np.linspace(start_deg, 0, neg_points, endpoint=False)
    neg_cos = np.cos(np.radians(neg_deg))
    # Positive part: from 0 to 90 deg
    pos_deg = np.linspace(0, stop_deg, pos_points)
    pos_cos = np.cos(np.radians(pos_deg))
    cos_vals = np.concatenate([neg_cos, pos_cos])
    angles_rad = np.arccos(cos_vals)
    # Assign negative sign to angles that came from the negative part
    angles_rad[:neg_points] *= -1
    decs = np.degrees(angles_rad)
    RA_grid, DEC_grid = np.meshgrid(ras, decs)

    grid_stars = {
        'starname': ['noname_' + str(i) for i in range(n_dec*n_ra)],
        'program_code': ['noprog_' + str(i) for i in range(n_dec*n_ra)],
        'ra': RA_grid.flatten(),
        'dec': DEC_grid.flatten(),
        'exptime': [300]*(n_dec*n_ra),
        'n_exp': [1]*(n_dec*n_ra),
        'n_intra_max': [1]*(n_dec*n_ra),
        'n_intra_min': [1]*(n_dec*n_ra),
        'n_inter_max': [1]*(n_dec*n_ra),
        'tau_inter': [1]*(n_dec*n_ra),
        'tau_intra': [1]*(n_dec*n_ra)
    }
    grid_frame = pd.DataFrame(grid_stars)

    available_nights_onsky_requests = compute_seasonality(semester_planner, program_frame['starname'], program_frame['ra'], program_frame['dec'])
    
    # Check if cached sky availability data exists for the background grid. 
    semester = semester_planner.semester_start_date[:4] + semester_planner.semester_letter
    cache_file = f"{DATADIR}/{semester}_sky_availability.csv"
    
    if os.path.exists(cache_file):
        grid_frame = pd.read_csv(cache_file)
    else:
        grid_frame['nights_observable'] = compute_seasonality(semester_planner, grid_frame['starname'], grid_frame['ra'], grid_frame['dec'])
        cache_df = pd.DataFrame({
            'starname': grid_frame['starname'],
            'ra': grid_frame['ra'],
            'dec': grid_frame['dec'],
            'nights_observable': grid_frame['nights_observable']
        })
        cache_df.to_csv(cache_file, index=False)

    from scipy.interpolate import griddata
    NIGHTS_grid = griddata(
    points=(grid_frame.ra, grid_frame.dec),
    values=grid_frame.nights_observable,
    xi=(RA_grid, DEC_grid),
    method='linear'
    )

    # Step 1: Generate static heatmap image with matplotlib
    RA_shifted = np.radians(RA_grid - 180)
    DEC_rad = np.radians(DEC_grid)

    fig_mpl, ax = plt.subplots(subplot_kw={'projection': 'mollweide'}, figsize=(10, 5))
    im = ax.pcolormesh(RA_shifted, DEC_rad, NIGHTS_grid, cmap='gray', shading='nearest', vmin=70, vmax=184)
    ax.axis('off')

    # Save to buffer
    buf = BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', pad_inches=0, dpi=150)
    plt.close()
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode()

    # Step 2: Create Plotly figure with static background image
    fig = go.Figure()

    fig.add_layout_image(
        dict(
            source=f'data:image/png;base64,{img_base64}',
            xref="paper", yref="paper",
            x=0, y=1,
            sizex=1, sizey=1,
            xanchor="left", yanchor="top",
            sizing="stretch",
            layer="below",
            opacity=1
        )
    )

    # Step 3: Add dummy contour for colorbar
    fig.add_trace(go.Contour(
        z=NIGHTS_grid,
        x=RA_grid[0] - 180,
        y=DEC_grid[:, 0],
        showscale=True,
        colorscale='gray',
        contours=dict(start=70, end=184, size=10),
        opacity=0,  # Hide contour but keep colorbar
        colorbar=dict(
            title='Observable<br>Nights',
            titleside='top',
            x=-0.15,  # Place on left of plot
            len=0.75,
            thickness=15
        )
    ))

    # Step 4: Add interactive points grouped by program
    if not program_frame.empty:
        grouped = program_frame.groupby('program_code')
        for program, group in grouped:
            group.reset_index(inplace=True, drop=True)
            hover = [f"{name} in {program}" for name in group['starname']]
            color = group['color'].tolist()  # Use individual star colors

            if len(all_stars)==1:
                size=20
                marker='star'
            else:
                size=6
                marker='star'
            fig.add_trace(go.Scattergeo(
                lon=group['ra'] - 180,
                lat=group['dec'],
                mode='markers',
                name=program,
                marker=dict(symbol=marker, size=size, color=color, opacity=1),
                text=hover,
                hovertemplate="%{text}<br>RA: %{lon:.2f}°, Dec: %{lat:.2f}°<extra></extra>"
            ))

    fig.update_layout(shapes=[
        dict(
            type="circle",
            xref="paper", yref="paper",
            x0=0.0, y0=0.0, x1=1., y1=1.,
            line=dict(color="black", width=2)
        )
    ])

    # Step 5: Layout
    fig.update_layout(
        geo=dict(
            projection_type='mollweide',
            showland=False,
            showcoastlines=False,
            showframe=False,
            bgcolor='rgba(0,0,0,0)',
            lonaxis=dict(showgrid=False),
            lataxis=dict(showgrid=False),
            ),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        template='none',
        width=1000,
        height=600,
        xaxis=dict(showgrid=False, visible=True),
        yaxis=dict(showgrid=False, visible=True),
        annotations=[
            dict(
                text="RA (deg)",  # X-axis label
                x=0.5,
                y=-0.10,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=14)
            ),
            dict(
                text="Dec (deg)",  # Y-axis label
                x=-0.07,
                y=0.5,
                xref="paper",
                yref="paper",
                showarrow=False,
                textangle=-90,
                font=dict(size=14)
            )
        ]
    )
    return fig


def get_request_frame(semester_planner, all_stars):
    """
    Get a filtered request frame containing only the stars in all_stars.
    
    Args:
        semester_planner: the semester planner object
        all_stars (list): array of StarPlotter objects
        
    Returns:
        pd.DataFrame: filtered request frame with only the specified stars
    """
    # Extract starnames from the StarPlotter objects
    starids = [star.unique_id for star in all_stars]
    
    # Filter the request frame to only include the specified stars
    filtered_frame = semester_planner.requests_frame[
        semester_planner.requests_frame['unique_id'].isin(starids)
    ].copy()
    
    return filtered_frame

def get_ladder(data):
    """Create an interactive plot which illustrates the solution.

    Args:
        data: TTP data containing the schedule information
    """

    orderData = data[0].plotly

    # reverse the order so that the plot flows from top to bottom with time
    orderData = pd.DataFrame.from_dict(orderData)
    orderData = orderData.iloc[::-1]
    orderData.reset_index(inplace=True)

    # Each priority gets a different color. Make sure that each priority is actually included here or the plot will break. Recall bigger numbers are higher priorities.
    colordict = {'10':'red',
                 '9':'tomato',
                 '8':'darkorange',
                 '9':'sandybrown',
                 '7':'gold',
                 '6':'olive',
                 '5':'green',
                 '4':'cyan',
                 '3':'darkviolet',
                 '2':'magenta',
                '1':'blue'}

    # build the outline of the plot, add dummy points that are not displyed within the x/y limits so as to fill in the legend
    fig = px.scatter(orderData, x='Minutes the from Start of the Night', y="human_starname", hover_data=['First Available', 'Last Available', 'Exposure Time (min)', "N_shots", "Total Exp Time (min)"] ,title='Night Plan', width=800, height=1000) #color='Program'
    fig.add_shape(type="rect", x0=-100, x1=-80, y0=-0.5, y1=0.5, fillcolor='red', showlegend=True, name='Expose P1')
    fig.add_shape(type="rect", x0=-100, x1=-80, y0=-0.5, y1=0.5, fillcolor='blue', showlegend=True, name='Expose P3')
    fig.add_shape(type="rect", x0=-100, x1=-80, y0=-0.5, y1=0.5, fillcolor='lime', opacity=0.3, showlegend=True, name='Accessible')

    new_already_processed = []
    ifixer = 0 # for multi-visit targets, it throws off the one row per target plotting...this fixes it
    for i in range(len(orderData['Starname'])):
        if orderData['Starname'][i] not in new_already_processed:
            # find all the times in the night when the star is being visited
            indices = [k for k in range(len(orderData['Starname'])) if orderData['Starname'][k] == orderData['Starname'][i]]
            for j in range(len(indices)):
                fig.add_shape(type="rect", x0=orderData['Start Exposure'][indices[j]], x1=orderData['Start Exposure'][indices[j]] + orderData["Total Exp Time (min)"][indices[j]], y0=i+ifixer-0.5, y1=i+ifixer+0.5, fillcolor=colordict[str(orderData['Priority'][indices[j]])])
                if j == 0:
                    # only do this once, otherwise the green bar gets discolored compared to other rows
                    fig.add_shape(type="rect", x0=orderData['First Available'][indices[j]], x1=orderData['Last Available'][indices[j]], y0=i+ifixer-0.5, y1=i+ifixer+0.5, fillcolor='lime', opacity=0.3, showlegend=False)
            new_already_processed.append(orderData['Starname'][i])
        else:
            # if we already did this star, it is a multi-visit star and we need to adjust the row counter for plotting purposes
            ifixer -= 1

    fig.update_layout(xaxis_range=[0,orderData['Start Exposure'][0] + orderData["Total Exp Time (min)"][0]])
    return fig


def createTelSlewPath(stamps, changes, pointings, animationStep=120):
    '''
    Correctly assign each frame of the animation to the telescope pointing at that time

    stamps (list of zeros) - the list where each element represents a frame of the animation. We manipulate and return this at the end.
    changes (list) - the times at which the telescope pointing changes (in order of the slew path)
    poitings (list) - the astropy target objects of for the stars to be observed, in order of the slew path
    animationStep (int) - the time, in seconds, between frames

    return
        stamps - now a list where element holds the pointing of the telescope (aka the star object) at that frame

    '''
    # determine how many minutes each frame of the animation represents
    minPerStep = int(animationStep/60)
    mins = int(60/minPerStep)

    # offset the timestamps of the observations to a zero point
    changes = (changes - changes[0])*24*mins
    for c in range(len(changes)):
        changes[c] = int(changes[c])

    # determine telescope pointing at each frame
    for i in range(len(changes)-1):
        for j in range(len(stamps)):
            if j >= changes[i] and j < changes[i+1]:
                stamps[j] = pointings[i]

    # Add edge cases of the first and last telescope pointing
    if len(stamps) > 0:
        k = 0
        while k < len(stamps) and stamps[k] == 0:
            stamps[k] = pointings[0]
            k += 1
        l = len(stamps)-1
        while l >= 0 and stamps[l] == 0:
            stamps[l] = pointings[-1]
            l -= 1

    return stamps

def get_slew_animation(data, animationStep=120):
    """Create a list of matplotlib figures showing telescope slew path during observations.

    Args:
        data: TTP data containing the schedule information
        animationStep (int): the time, in seconds, between animation still frames. Default to 120s.
        
    Returns:
        list: List of matplotlib Figure objects representing each frame of the animation
    """

    model = data[0]

    # Set up animation times
    t = np.arange(model.nightstarts.jd, model.nightends.jd, TimeDelta(animationStep, format='sec').jd)
    t = Time(t, format='jd')

    # Get list of astropy target objects in scheduled order
    names = [s for s in model.schedule['Starname']]
    list_targets = []
    for n in range(len(model.schedule['Starname'])):
        for s in model.stars:
            if s.name == model.schedule['Starname'][n]:
                list_targets.append(s.target)

    # Compute alt/az of each target at each time
    AZ = model.observatory.observer.altaz(t, list_targets, grid_times_targets=True)
    alt = np.round(AZ.az.rad, 2).T.tolist()
    az = (90 - np.round(AZ.alt.deg, 2)).T.tolist()

    # Telescope slew path
    stamps = [0] * len(t)
    slewPath = createTelSlewPath(stamps, model.schedule['Time'], list_targets)
    AZ1 = model.observatory.observer.altaz(t, slewPath, grid_times_targets=False)
    tel_az = np.round(AZ1.az.rad, 2)
    tel_zen = 90 - np.round(AZ1.alt.deg, 2)

    # Set up plotting params
    plotlowlim = 80
    theta = np.arange(model.observatory.deckAzLim1 / 180, model.observatory.deckAzLim2 / 180, 1. / 180) * np.pi
    allaround = np.arange(0., 360 / 180, 1. / 180) * np.pi

    # Create list of matplotlib figures
    figures = []
    for i in range(len(t)):
        fig = Figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='polar')
        ax.set_ylim(0, plotlowlim)
        ax.set_title(t.isot[i])
        ax.set_yticklabels([])
        ax.fill_between(theta, 90 - model.observatory.deckAltLim, plotlowlim, color='red', alpha=.7)
        ax.fill_between(allaround, 90 - model.observatory.vigLim, plotlowlim, color='red', alpha=.7)
        ax.fill_between(allaround, 90 - model.observatory.zenLim, 0, color='red', alpha=.7)
        ax.set_theta_zero_location('N')
        ax.set_facecolor('black')

        # Determine which stars have been observed by this time
        wasObserved = np.array(model.schedule['Time']) <= float(t[i].jd)
        observed_list = np.array(model.schedule['Starname'])[wasObserved]

        # Plot stars
        for j in range(len(model.schedule['Starname'])):
            color = 'orange' if names[j] in observed_list else 'white'
            ax.scatter(alt[i][j], az[i][j], color=color, marker='*')

        # Draw telescope path
        ax.plot(tel_az[:i], tel_zen[:i], color='orange')
        
        figures.append(fig)

    return figures

def get_script_plan(night_planner):
    """Generate script plan DataFrame from semester planner and night planner objects.
    
    This function reads the request_selected.csv file from the semester planner's output directory,
    merges it with the night planner's solution data, and returns a properly formatted DataFrame
    with the same column structure as the original get_script_plan function.
    
    Args:
        night_planner: NightPlanner object containing solution attribute
        
    Returns:
        pd.DataFrame: Formatted observing plan DataFrame
    """
    
    # Read the request_selected.csv file from the semester planner's output directory
    request_selected_path = os.path.join(night_planner.output_directory, 'request_selected.csv')
    
    if not os.path.exists(request_selected_path):
        raise FileNotFoundError(f"request_selected.csv not found at {request_selected_path}")
    
    # Read the request_selected.csv file
    request_selected_df = pd.read_csv(request_selected_path)
    solution = night_planner.solution[0]  # First index as specified
    
    # Extract the schedule from the solution and convert to a DataFrame
    solution_schedule = solution.plotly
    solution_df = pd.DataFrame(solution_schedule)
    
    # Merge the solution DataFrame with the request_selected dataframe
    # Use starname as the key for merging
    merged_df = pd.merge(request_selected_df, solution_df, 
                         left_on='unique_id', right_on='Starname', 
                         how='inner')
                             
    # Select and reorder only the specific columns requested
    desired_columns = [
        'Start Exposure', 'unique_id', 'starname', 'program_code', 'ra', 'dec', 
        'exptime', 'n_exp', 'n_intra_max', 'tau_intra', 'weather_band_1', 'weather_band_2', 'weather_band_3', 'teff', 
        'jmag', 'gmag', 'epoch', 'gaia_id', 'First Available', 'Last Available'
    ]
    
    # Keep only the columns that exist in the merged dataframe
    available_columns = [col for col in desired_columns if col in merged_df.columns]
    
    # Reorder columns to match the desired structure
    final_df = merged_df[available_columns].copy()
    
    # Round numeric fields to appropriate decimal places
    if 'ra' in final_df.columns:
        # Ensure ra is numeric before rounding, handle 'None' strings
        final_df['ra'] = final_df['ra'].replace('None', pd.NA)
        final_df['ra'] = pd.to_numeric(final_df['ra'], errors='coerce').round(3)
    
    if 'dec' in final_df.columns:
        # Ensure dec is numeric before rounding, handle 'None' strings
        final_df['dec'] = final_df['dec'].replace('None', pd.NA)
        final_df['dec'] = pd.to_numeric(final_df['dec'], errors='coerce').round(3)
    
    if 'jmag' in final_df.columns:
        # Ensure jmag is numeric before rounding, handle 'None' strings
        final_df['jmag'] = final_df['jmag'].replace('None', pd.NA)
        final_df['jmag'] = pd.to_numeric(final_df['jmag'], errors='coerce').round(1)
    
    if 'gmag' in final_df.columns:
        # Ensure gmag is numeric before rounding, handle 'None' strings
        final_df['gmag'] = final_df['gmag'].replace('None', pd.NA)
        final_df['gmag'] = pd.to_numeric(final_df['gmag'], errors='coerce').round(1)
    
    # Convert time fields from "minutes from start of night" to HST timestamps
    try:
        # Get the night start time from the night planner
        from astroq.nplan import get_nightly_times_from_allocation
        from astropy.time import TimeDelta
        
        night_start_time, _ = get_nightly_times_from_allocation(
            night_planner.allocation_file, 
            night_planner.current_day
        )
        
        # Convert the time columns to HST timestamps
        if 'Start Exposure' in final_df.columns:
            final_df['Start Exposure'] = final_df['Start Exposure'].apply(
                lambda x: str(TimeDelta(x * 60, format='sec') + night_start_time)[11:16] if pd.notna(x) else ''
            )
        
        if 'First Available' in final_df.columns:
            final_df['First Available'] = final_df['First Available'].apply(
                lambda x: str(TimeDelta(x * 60, format='sec') + night_start_time)[11:16] if pd.notna(x) else ''
            )
        
        if 'Last Available' in final_df.columns:
            final_df['Last Available'] = final_df['Last Available'].apply(
                lambda x: str(TimeDelta(x * 60, format='sec') + night_start_time)[11:16] if pd.notna(x) else ''
            )
            
    except Exception as e:
        print(f"Warning: Could not convert time fields to HST timestamps: {e}")
        print("Time fields will remain as minutes from start of night")
    
    # Handle missing values and 'None' strings
    final_df = final_df.replace(['', 'NoGaiaName', 'None'], pd.NA)
    
    return final_df

def plot_path_2D_interactive(data):
    """Create an interactive Plotly plot showing telescope azimuth and altitude paths with UTC times and white background."""

    model = data[0]

    names = list(model.plotly['human_starname'])
    times = model.times
    az_path = model.az_path
    alt_path = model.alt_path
    wrap = model.observatory.wrapLimitAngle

    # Use times as floats (e.g., minutes from start of night or JD)
    # If times are astropy Time objects, convert to minutes from start of night
    # Here, we'll use JD as the linear x-axis, but you can adjust to minutes if you have a reference
    obs_time = np.array([t.jd for t in times])
    # For tick labels and hover, format as HH:MM
    time_labels = [Time(t, format='jd').isot[11:16] for t in obs_time]

    # Prepare custom hover text with HH:MM time
    hover_texts = [f"Time: {label}<br>Target: {name}" for label, name in zip(time_labels, names)]

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        subplot_titles=("Azimuth Path", "Elevation Path"),
        vertical_spacing=0.1
    )

    # Azimuth plot
    fig.add_trace(go.Scatter(
        x=obs_time, y=az_path,
        mode='lines+markers',
        marker=dict(color='indigo'),
        name='Azimuth',
        text=hover_texts,
        hovertemplate='%{text}<br>Az: %{y}'
    ), row=1, col=1)

    # Elevation plot
    fig.add_trace(go.Scatter(
        x=obs_time, y=alt_path,
        mode='lines+markers',
        marker=dict(color='seagreen'),
        name='Elevation',
        text=hover_texts,
        hovertemplate='%{text}<br>Alt: %{y}'
    ), row=2, col=1)

    # Optional: wrap line on azimuth
    if wrap is not None:
        fig.add_shape(
            type="line",
            x0=obs_time[0], x1=obs_time[-1],
            y0=wrap, y1=wrap,
            line=dict(color="red", dash="dash"),
            row=1, col=1
        )
        fig.add_annotation(
            x=obs_time[-1], y=wrap,
            text=f"Wrap = {wrap}",
            showarrow=False,
            font=dict(color="red"),
            row=1, col=1
        )

    # Highlight observed intervals (use obs_time for x0/x1)
    for i in range(0, len(obs_time)-1, 2):
        fig.add_vrect(
            x0=obs_time[i], x1=obs_time[i+1],
            fillcolor="orange", opacity=0.2,
            layer="below", line_width=0,
            row=1, col=1
        )
        fig.add_vrect(
            x0=obs_time[i], x1=obs_time[i+1],
            fillcolor="orange", opacity=0.2,
            layer="below", line_width=0,
            row=2, col=1
        )

    # Set x-axis tick labels as HH:MM with evenly spaced grid
    # Create evenly spaced time ticks (e.g., every hour or every 30 minutes)
    time_span = obs_time[-1] - obs_time[0]
    # Determine appropriate interval based on time span
    if time_span < 0.1:  # Less than ~2.4 hours
        interval_hours = 0.5  # 30 minutes
    elif time_span < 0.3:  # Less than ~7 hours  
        interval_hours = 1.0  # 1 hour
    else:
        interval_hours = 2.0  # 2 hours
    
    # Convert interval to JD units (1 hour = 1/24 JD)
    interval_jd = interval_hours / 24
    
    # Create evenly spaced tick positions
    start_time = obs_time[0]
    end_time = obs_time[-1]
    num_ticks = int((end_time - start_time) / interval_jd) + 2
    tick_positions = np.linspace(start_time, end_time, num_ticks)
    
    # Convert tick positions to time labels
    tick_labels = [Time(t, format='jd').isot[11:16] for t in tick_positions]
    
    fig.update_xaxes(
        tickmode='array',
        tickvals=tick_positions,
        ticktext=tick_labels,
        title_text='Time (UTC)',
        row=2, col=1
    )

    fig.update_layout(
        height=600,
        width=1000,
        yaxis_title="Azimuth (deg)",
        yaxis2_title="Altitude (deg)",
        template="plotly_white"
    )
    return fig

def dataframe_to_html(dataframe, sort_column=2, page_size=10, table_id='request-table'):
    """
    Convert a pandas dataframe into an HTML string for rendering
    on the webapp pages.
    
    Args:
        dataframe (pd.DataFrame): The dataframe to convert
        sort_column (int): Column index to sort by (default: 2 for starname)
        page_size (int): Default number of rows per page (default: 25)
        table_id (str): Unique ID for the table (default: 'request-table')
        
    Returns:
        str: HTML string with table and DataTables initialization
    """
    # Convert DataFrame to HTML table with unique ID
    table_html = dataframe.to_html(
        classes='table table-striped table-hover', 
        index=False, 
        escape=False, 
        table_id=table_id
    )
    
    # Custom CSS for beautiful table styling
    custom_css = f"""
    <style>
    #{table_id} {{
        border-collapse: separate !important;
        border-spacing: 0 !important;
        border-radius: 8px !important;
        overflow: hidden !important;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1) !important;
        margin: 20px 0 !important;
    }}
    
    #{table_id} thead th {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        font-weight: 600 !important;
        padding: 15px 12px !important;
        border: none !important;
        text-align: center !important;
        font-size: 14px !important;
        text-transform: uppercase !important;
        letter-spacing: 0.5px !important;
    }}
    
    /* Additional DataTables header styling to ensure purple background */
    #{table_id} .dataTables_scrollHead thead th {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        font-weight: 600 !important;
        padding: 15px 12px !important;
        border: none !important;
        text-align: left !important;
        font-size: 14px !important;
        text-transform: uppercase !important;
        letter-spacing: 0.5px !important;
    }}
    
    #{table_id} .dataTables_wrapper .dataTables_scrollHead thead th {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        font-weight: 600 !important;
        padding: 15px 12px !important;
        border: none !important;
        text-align: left !important;
        font-size: 14px !important;
        text-transform: uppercase !important;
        letter-spacing: 0.5px !important;
    }}
    
    /* Target any header cells that might be created by DataTables */
    #{table_id} th, #{table_id} .dataTables_scrollHead th {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        font-weight: 600 !important;
        padding: 15px 12px !important;
        border: none !important;
        text-align: left !important;
        font-size: 14px !important;
        text-transform: uppercase !important;
        letter-spacing: 0.5px !important;
    }}
    
    /* More aggressive DataTables header targeting */
    #{table_id} .dataTables_wrapper thead th {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        font-weight: 600 !important;
        padding: 15px 12px !important;
        border: none !important;
        text-align: left !important;
        font-size: 14px !important;
        text-transform: uppercase !important;
        letter-spacing: 0.5px !important;
    }}
    
    /* Force all header cells to have purple background */
    #{table_id} thead th, #{table_id} th, #{table_id} .dataTables_scrollHead th, #{table_id} .dataTables_wrapper th {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        font-weight: 600 !important;
        padding: 15px 12px !important;
        border: none !important;
        text-align: left !important;
        font-size: 14px !important;
        text-transform: uppercase !important;
        letter-spacing: 0.5px !important;
    }}
    
    /* Override any DataTables default styling */
    #{table_id} .dataTables_wrapper .dataTables_scrollHead thead th,
    #{table_id} .dataTables_wrapper .dataTables_scrollHead thead td {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        font-weight: 600 !important;
        padding: 15px 12px !important;
        border: none !important;
        text-align: left !important;
        font-size: 14px !important;
        text-transform: uppercase !important;
        letter-spacing: 0.5px !important;
    }}
    
    /* Nuclear option - target everything with maximum specificity */
    #{table_id} .dataTables_wrapper .dataTables_scrollHead thead th,
    #{table_id} .dataTables_wrapper .dataTables_scrollHead thead td,
    #{table_id} .dataTables_wrapper thead th,
    #{table_id} .dataTables_wrapper thead td,
    #{table_id} thead th,
    #{table_id} thead td {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        font-weight: 600 !important;
        padding: 15px 12px !important;
        border: none !important;
        text-align: left !important;
        font-size: 14px !important;
        text-transform: uppercase !important;
        letter-spacing: 0.5px !important;
    }}
    
    /* Force override any inline styles that DataTables might add */
    #{table_id} thead th[style*="background"],
    #{table_id} .dataTables_wrapper thead th[style*="background"] {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        font-weight: 600 !important;
        padding: 15px 12px !important;
        border: none !important;
        text-align: left !important;
        font-size: 14px !important;
        text-transform: uppercase !important;
        letter-spacing: 0.5px !important;
    }}
    
    #{table_id} tbody tr {{
        transition: all 0.2s ease !important;
    }}
    
    #{table_id} tbody tr:nth-child(even) {{
        background-color: #f8f9fa !important;
    }}
    
    #{table_id} tbody tr:nth-child(odd) {{
        background-color: white !important;
    }}
    
    #{table_id} tbody tr:hover {{
        background-color: #e3f2fd !important;
        transform: translateY(-1px) !important;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1) !important;
    }}
    
    #{table_id} tbody td {{
        padding: 12px !important;
        border: none !important;
        border-bottom: 1px solid #e9ecef !important;
        font-size: 14px !important;
        color: #495057 !important;
        text-align: center !important;
    }}
    
    #{table_id} tbody tr:last-child td {{
        border-bottom: none !important;
    }}
    
    /* Column widths are now controlled by DataTables columnDefs for proper alignment */
    
    /* DataTables controls styling */
    .dataTables_length select {{
        border: 1px solid #ddd !important;
        border-radius: 4px !important;
        padding: 4px 8px !important;
        background: white !important;
    }}
    
    .dataTables_filter input {{
        border: 1px solid #ddd !important;
        border-radius: 4px !important;
        padding: 6px 12px !important;
        background: white !important;
    }}
    
    .dt-buttons .dt-button {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        border: none !important;
        border-radius: 4px !important;
        padding: 8px 16px !important;
        margin: 2px !important;
        font-size: 12px !important;
        transition: all 0.2s ease !important;
    }}
    
    .dt-buttons .dt-button:hover {{
        background: linear-gradient(135deg, #5a6fd8 0%, #6a4190 100%) !important;
        transform: translateY(-1px) !important;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2) !important;
    }}
    
    .dataTables_info {{
        color: #6c757d !important;
        font-size: 14px !important;
        margin-top: 10px !important;
    }}
    
    .dataTables_paginate .paginate_button {{
        border: 1px solid #ddd !important;
        border-radius: 4px !important;
        padding: 6px 12px !important;
        margin: 2px !important;
        background: white !important;
        color: #495057 !important;
        transition: all 0.2s ease !important;
    }}
    
    .dataTables_paginate .paginate_button:hover {{
        background: #e9ecef !important;
        border-color: #adb5bd !important;
    }}
    
    .dataTables_paginate .paginate_button.current {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        border-color: #667eea !important;
    }}
    </style>
    """
    
    # DataTables initialization script - destroy existing instance first
    init_script = f"""
    <script>
    $(document).ready(function() {{
        // Destroy existing DataTable if it exists
        if ($.fn.DataTable.isDataTable('#{table_id}')) {{
            $('#{table_id}').DataTable().destroy();
        }}
        
        // Initialize new DataTable
        var table = $('#{table_id}').DataTable({{
            pageLength: {page_size},
            order: [[{sort_column}, 'asc']],
            dom: 'lBfrtip',
            buttons: ['copy', 'csv', 'excel', 'print'],
            scrollX: false,  // Disable horizontal scrolling to prevent header misalignment
            responsive: false,  // Disable responsive features that can cause header issues
            lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
            tableLayout: 'fixed',  // Force fixed table layout for consistent column widths
            columnDefs: [
                {{ targets: 0, width: '60px' }},   // Start Exposure
                {{ targets: 1, width: '260px' }},  // unique_id
                {{ targets: 2, width: '260px' }},  // starname
                {{ targets: 3, width: '120px' }},  // program_code
                {{ targets: 4, width: '100px' }},  // ra
                {{ targets: 5, width: '100px' }},  // dec
                {{ targets: 6, width: '50px' }},   // exptime
                {{ targets: 7, width: '50px' }},   // n_exp
                {{ targets: 8, width: '50px' }},   // n_intra_max
                {{ targets: 9, width: '50px' }},   // tau_intra
                {{ targets: 10, width: '50px' }},  // weather_band_1
                {{ targets: 11, width: '50px' }},  // weather_band_2
                {{ targets: 12, width: '50px' }},  // weather_band_3
                {{ targets: 13, width: '50px' }},  // teff
                {{ targets: 14, width: '50px' }},  // jmag
                {{ targets: 15, width: '50px' }},  // gmag
                {{ targets: 16, width: '50px' }},  // epoch
                {{ targets: 17, width: '260px' }}, // gaia_id
                {{ targets: 18, width: '60px' }},  // First Available
                {{ targets: 19, width: '60px' }}   // Last Available
            ],
            initComplete: function() {{
                // Simple styling after DataTables is initialized
                $('#{table_id} thead th').css({{
                    'background': 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                    'color': 'white',
                    'font-weight': '600',
                    'padding': '15px 12px',
                    'border': 'none',
                    'text-align': 'center',
                    'font-size': '14px',
                    'text-transform': 'uppercase',
                    'letter-spacing': '0.5px'
                }});
            }}
        }});
    }});
    </script>
    """
    
    return custom_css + table_html + init_script

