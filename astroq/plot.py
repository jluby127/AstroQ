"""
Module for constructing the standard AstroQ plots. All plots are returned as html strings. 
From there, they can be used as is or saved as png files.
"""

# Standard library imports
from collections import defaultdict
from datetime import datetime, timedelta
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
        Define the StarPlotter class, which contains all information about a single request 
        which is used for standardizing the plot inputs.
    """

    def __init__(self, unique_id):
        """
        Initialize the StarPlotter object.

        Args:
            unique_id (str): the unique id of the request, from the request.csv file.

        Returns:
            None
        """
        self.unique_id = unique_id

    def get_stats(self, row, slot_size):
        """
        Grab the observational stategy information for a given star from the requests.csv file.

        Args:
            row (pd.Series): A row from the requests.csv file as a DataFrame 
            slot_size (int): The slot size in minutes

        Returns:
            expected_nobs_per_night (int): how many exposures we expect to take
            total_observations_requested (int): sum of observational strategie values
            exposure_time (int): exposure time of single shot
            slots_per_night (int): number of slots required to complete all exposures in a night
            program (str): the program code
        """
        # Access row data directly instead of filtering the entire DataFrame (PERFORMANCE OPTIMIZATION)
        self.starname = row['starname']
        self.inactive = row['inactive']
        self.ra = float(row['ra'])
        self.dec = float(row['dec'])
        self.program = str(row['program_code'])
        self.exptime = int(row['exptime'])
        self.n_exp = int(row['n_exp'])
        self.n_intra_max = int(row['n_intra_max'])
        self.n_intra_min = int(row['n_intra_min'])
        self.tau_intra = int(row['tau_intra'])
        self.n_inter_max = int(row['n_inter_max'])
        self.tau_inter = int(row['tau_inter'])
        self.total_observations_requested = self.n_exp * self.n_intra_max * self.n_inter_max
        self.total_requested_seconds = self.total_observations_requested*self.exptime + readout_overhead*(self.n_exp-1)* self.n_inter_max + slew_overhead*self.n_intra_max*self.n_inter_max
        self.total_requested_hours = self.total_requested_seconds / 3600
        self.total_requested_nights = self.total_requested_hours / hours_per_night   


    def get_past(self, past):
        """
        Gather the information about a star's past observation history this semester in standard format.

        Args:
            past (DataFrame): A DataFrame version of past.csv.
        Returns:
            None. 
        """
        # Filter to this star
        star_obs_past = past[past['target'] == str(self.unique_id)]
        # Parse date from timestamp and group by date
        star_obs_past = star_obs_past.copy()
        star_obs_past['date'] = star_obs_past['timestamp'].str[:10]
        observations_past = star_obs_past.groupby('date').size().to_dict()
        self.observations_past = observations_past

    def get_future(self, forecast_df, all_dates_array):
        """
        Gather the star's future schedule out of the semester_planner solution from semester_plan.csv.

        Args:
            forecast_df (pd.DataFrame): Pre-loaded forecast DataFrame with minimum columns ['r', 'd', 's']
            all_dates_array (list): List of all dates in the semester, indexed by 'd'
        
        Returns:
            None
        """
        # Only keep rows for this star
        star_rows = forecast_df[forecast_df['r'] == str(self.unique_id)]
        # Count number of slots scheduled per night (d)
        observations_future = {}
        for d, group in star_rows.groupby('d'):
            # d may be int or str; ensure it's int for indexing
            date = all_dates_array[int(d)]
            n_slots = len(group) # this is the number of starting slots in given to this target in this night
            observations_future[date] = n_slots # no need to multiply by nexp here because we do it later in timebar; so that COF has right values.
        self.observations_future = observations_future

    def get_map(self, semester_planner, forecast_df):
        """
        Build the 2D d/s matrix starmap for teh given star using semester_plan.csv.
        Only set starmap[d, s] = 1 if sched['r'] == self.unique_id.
        
        Args:
            semester_planner: The semester planner object
            forecast_df (pd.DataFrame): Pre-loaded forecast DataFrame with columns ['r', 'd', 's']
        """
        n_nights = semester_planner.semester_length
        n_slots = int((24 * 60) / semester_planner.slot_size)
        starmap = np.zeros((n_nights, n_slots), dtype=int)
        
        # Filter to only this star's rows
        star_forecast = forecast_df[forecast_df['r'] == str(self.unique_id)]
        
        if len(star_forecast) > 0:
            # Vectorized approach: extract d,s values as numpy arrays and set all at once (PERFORMANCE OPTIMIZATION)
            d_values = star_forecast['d'].values.astype(int)
            s_values = star_forecast['s'].values.astype(int)
            
            # Set the primary slots
            starmap[d_values, s_values] = 1
            
            # Set the reserve slots
            reserve_slots = semester_planner.slots_needed_for_exposure_dict[str(self.unique_id)]
            for r in range(1, reserve_slots):
                starmap[d_values, s_values + r] = 1
        
        self.starmap = starmap.T

def process_stars(semester_planner):
    """
    Construct the StarPlotter objects for all the stars in the semester planner.

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py

    Returns:
        program_dict (dict): a dictionary of program names and their corresponding StarPlotter objects
        programs_as_stars (dict): a dictionary of program names and their corresponding StarPlotter objects
        nulltime (array): a 2D array of N_slots by N_nights, binary 1/0, it is the intersection of is_alloc and is_night
    """

    # Create a starmap of the times when we cannot observe due to twilight and allocation constraints
    # Used in the birdseye view plot to blackout the unavailble squares

    # Use the stored access record from the semester planner instead of recomputing
    access = semester_planner.access_record
    nulltime = access['is_alloc'][0]
    nulltime = 1 - nulltime
    nulltime = np.array(nulltime).T

    # Read forecast CSV once instead of once per star (PERFORMANCE OPTIMIZATION)
    forecast_df = semester_planner.serialized_schedule # pd.read_csv(semester_planner.output_directory + semester_planner.future_forecast)
    forecast_df['r'] = forecast_df['r'].astype(str)  # Convert to string once

    # Previously, there was a unique call to star names, every row of the request frame will be unique already when we switch to "id"
    starnames = semester_planner.requests_frame_all['starname'].unique()
    programs = semester_planner.requests_frame_all['program_code'].unique()

    # Make colors consistent for all stars in each program
    colors = sns.color_palette("deep", len(programs))
    rgb_strings = [f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})" for r, g, b in colors]
    program_colors_rgb_vals = dict(zip(programs, rgb_strings))

    all_stars = []
    i = 0 
    for i, row in semester_planner.requests_frame_all.iterrows():
        # Create a StarPlotter object for each request, fill and compute relavant information
        newstar = StarPlotter(row['unique_id'])
        newstar.get_map(semester_planner, forecast_df)
        newstar.get_stats(row, semester_planner.slot_size)
        if newstar.unique_id in list(semester_planner.past_history.keys()):
            newstar.observations_past = semester_planner.past_history[newstar.unique_id].n_visits_on_nights
            newstar.observations_past_exposures = semester_planner.past_history[newstar.unique_id].n_obs_on_nights
        else:
            newstar.observations_past = {}
            newstar.observations_past_exposures = {}
        newstar.get_future(forecast_df, semester_planner.all_dates_array)

        # Create COF arrays for each request
        combined_set = set(list(newstar.observations_past.keys()) + list(newstar.observations_future.keys()))
        # For inactive stars, only include past observations; for active stars, include both past and future
        if newstar.inactive == False:
            newstar.dates_observe = [newstar.observations_past[date] if date in newstar.observations_past.keys() else (newstar.observations_future[date]*newstar.n_exp if date in combined_set else 0) for date in semester_planner.all_dates_array]
            newstar.dates_observe_time = [(newstar.observations_past_exposures[date]*newstar.exptime + readout_overhead*(newstar.observations_past[date]-1) + slew_overhead*(newstar.observations_past[date]-1)) / 3600 if date in newstar.observations_past_exposures.keys() else ((newstar.observations_future[date]*newstar.n_exp*newstar.exptime + readout_overhead*(newstar.n_exp-1)*newstar.observations_future[date] + slew_overhead*newstar.observations_future[date]) / 3600 if date in combined_set else 0) for date in semester_planner.all_dates_array]
        else:
            # For inactive stars, only show past observations
            newstar.dates_observe = [newstar.observations_past[date] if date in newstar.observations_past.keys() else 0 for date in semester_planner.all_dates_array]
            newstar.dates_observe_time = [(newstar.observations_past_exposures[date]*newstar.exptime + readout_overhead*(newstar.observations_past[date]-1) + slew_overhead*(newstar.observations_past[date]-1)) / 3600 if date in newstar.observations_past_exposures.keys() else 0 for date in semester_planner.all_dates_array]

        newstar.cume_observe = np.cumsum(newstar.dates_observe)
        newstar.cume_observe_time = np.cumsum(newstar.dates_observe_time)  # in hours

        if newstar.inactive:
            newstar.total_observations_requested = np.max(newstar.cume_observe)
            newstar.total_requested_seconds =newstar.total_observations_requested*newstar.exptime + slew_overhead*newstar.total_observations_requested
            newstar.total_requested_hours = newstar.total_requested_seconds / 3600
            newstar.total_requested_nights = newstar.total_requested_hours / hours_per_night   

        # Handle division by zero for inactive stars (total_observations_requested = 0)
        if newstar.total_observations_requested > 0:
            newstar.cume_observe_pct = np.round((np.cumsum(newstar.dates_observe)/newstar.total_observations_requested)*100.,3)
        else:
            # For inactive stars, show percentage based on total past observations if any exist
            total_past_obs = sum(newstar.observations_past.values()) if newstar.observations_past else 0
            if total_past_obs > 0:
                newstar.cume_observe_pct = np.round((np.cumsum(newstar.dates_observe)/total_past_obs)*100.,3)
            else:
                newstar.cume_observe_pct = np.zeros(len(semester_planner.all_dates_array))

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
        # For inactive targets, they won't be in requests_frame, so create zero maps
        try:
            target_idx = np.where(semester_planner.requests_frame['unique_id'] == newstar.unique_id)[0][0]
            # Extract the 2D slice for this specific target from each 3D map
            newstar.maps = {name: access[name][target_idx] for name in newstar.maps_names}
            newstar.allow_mapview = True
        except (IndexError, KeyError):
            # Target is inactive (not in access record) - create zero maps with appropriate shape
            n_nights = semester_planner.semester_length
            n_slots = int((24 * 60) / semester_planner.slot_size)
            newstar.maps = {name: np.zeros((n_nights, n_slots), dtype=bool) for name in newstar.maps_names}
            newstar.allow_mapview = False

        all_stars.append(newstar)
        i += 1

    # Now create StarPlotter objects for each program, as it were one star.
    # These will not have all the attributes, but we only need these for the admin COF plot
    # These StarPlotter objects cannot be used to create a birdseye plot, they don't have all attributes
    programmatics = pd.read_csv(os.path.join(semester_planner.semester_directory, 'programs.csv'))

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
        programmatic_star.cume_observe_pct = np.round(summed_cumulative / max_value * 100, 2)

        # Compute the cumulative observe time for all stars in the given program
        cume_observe_time = [all_stars[k].cume_observe_time for k in prog_indices]
        stars_stacked_time = np.vstack(cume_observe_time)
        summed_cumulative_time = np.sum(stars_stacked_time, axis=0)
        total_requested_prog = np.sum([all_stars[k].total_requested_hours for k in prog_indices])
        allocated = programmatics[programmatics['program'] == unique_programs[i]]['hours'].sum()
        # Use requested as divisor when requested < allocated, else allocated
        max_value_time = min(total_requested_prog, allocated)
        # summed_cumulative_time and max_value_time are both in hours
        if max_value_time > 0:
            programmatic_star.cume_observe_time_pct = np.round(summed_cumulative_time / max_value_time * 100, 2)
        else:
            programmatic_star.cume_observe_time_pct = np.zeros(len(semester_planner.all_dates_array))
        programmatic_star.cume_observe_time = summed_cumulative_time  # in hours

        # Handle division by zero for programs with only inactive stars
        if max_value > 0:
            programmatic_star.cume_observe_pct = np.round(summed_cumulative / max_value * 100, 2)
        else:
            # For inactive-only programs, use total past observations as denominator
            total_past_obs = sum(sum(all_stars[k].observations_past.values()) if all_stars[k].observations_past else 0 for k in prog_indices)
            if total_past_obs > 0:
                programmatic_star.cume_observe_pct = summed_cumulative / total_past_obs * 100
            else:
                programmatic_star.cume_observe_pct = np.zeros(len(semester_planner.all_dates_array))

        # Compute sum of starmaps
        super_map = np.zeros(np.shape(all_stars[prog_indices[0]].starmap))
        for m in range(len(prog_indices)):
            super_map += all_stars[prog_indices[m]].starmap
        programmatic_star.starmap = super_map
    
        # Aggregate observations_past for the program
        combined_past = {}
        for k in prog_indices:
            for date, count in all_stars[k].observations_past.items():
                combined_past[date] = combined_past.get(date, 0) + count
        programmatic_star.observations_past = combined_past

        programmatic_star.total_observations_requested = np.sum([all_stars[k].total_observations_requested for k in prog_indices])
        programmatic_star.total_requested_hours = np.sum([all_stars[k].total_requested_hours for k in prog_indices])
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

def get_cof(semester_planner, all_stars, use_time=False):
    '''
    Produce a plotly figure showing the Cumulative Observability Function (COF) for a selection of stars

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        all_stars (array): a array of StarPlotter objects
        use_time (bool): if True, use the cumulative observe time percentage instead of the cumulative observe percentage

    Returns:
        fig (plotly figure): a plotly figure showing the COF for a selection of stars
    '''

    fig = go.Figure()
    fig.update_layout(plot_bgcolor=gray, paper_bgcolor=clear) #autosize=True,margin=dict(l=40, r=40, t=40, b=40),
    
    # Convert calendar dates to night indices (0, 1, 2, ...)
    night_indices = np.arange(len(semester_planner.all_dates_array))
    
    burn_line = np.linspace(0, 100, len(semester_planner.all_dates_array))
    burn_line = np.round(burn_line, 2)

    # Add "Even Burn Rate" line as a shape so it's always visible and can't be toggled
    # Use add_shape to create a line that spans the entire plot
    fig.add_shape(
        type="line",
        x0=night_indices[0],
        y0=burn_line[0],
        x1=night_indices[-1],
        y1=burn_line[-1],
        line=dict(color='black', width=2, dash='dash'),
        layer='below',  # Draw below traces so it doesn't obscure data
    )
    
    # Add an invisible trace just for the legend entry (so users know what the line represents)
    # This trace will be visible in legend but clicking it won't hide the actual line
    fig.add_trace(go.Scatter(
        x=[None],  # No actual data points
        y=[None],
        mode='lines',
        line=dict(color='black', width=2, dash='dash'),
        name="Even Burn Rate",
        showlegend=True,
        hoverinfo='skip',  # Don't show hover for this dummy trace
    ))
    lines = []
    if use_time is False:
        cume_observe = np.zeros(len(semester_planner.all_dates_array))
        max_value = 0
        cume_observe = np.sum([star.cume_observe for star in all_stars], axis=0)
        max_value = sum(star.total_observations_requested for star in all_stars)
        # Handle division by zero: if all stars are inactive, use total past observations as denominator
        if max_value > 0:
            cume_observe_pct = np.round((cume_observe / max_value) * 100, 2)
        else:
            # For inactive-only programs, calculate total past observations
            total_past_obs = sum(sum(star.observations_past.values()) if star.observations_past else 0 for star in all_stars)
            if total_past_obs > 0:
                cume_observe_pct = (cume_observe / total_past_obs) * 100
            else:
                cume_observe_pct = np.zeros(len(semester_planner.all_dates_array))

        # Add the Total trace first (so it appears below other traces)
        fig.add_trace(go.Scatter(
            x=night_indices,
            y=cume_observe_pct,
            mode='lines',
            line=dict(color=all_stars[0].program_color_rgb, width=2),
            name="Total",
            hovertemplate= 'Night: %{x}' + '<br>Date: ' + '%{customdata}' + '<br>% Complete: %{y}' + '<br># Obs Requested: ' + \
                str(max_value) + '<br>',
            customdata=semester_planner.all_dates_array
        ))
    else:
        # use_time=True: normalize by program hours from programs.csv
        programmatics_cof = pd.read_csv(os.path.join(semester_planner.semester_directory, 'programs.csv'))
        programs_in_stars = set(getattr(s, 'program', getattr(s, 'starname', None)) for s in all_stars)
        programs_in_stars = {p for p in programs_in_stars if p is not None}
        summed_cume_time = np.sum([getattr(s, 'cume_observe_time', np.zeros(len(semester_planner.all_dates_array))) for s in all_stars], axis=0)
        total_program_hours = programmatics_cof[programmatics_cof['program'].isin(programs_in_stars)]['hours'].sum()
        # summed_cume_time and total_program_hours are both in hours
        if total_program_hours > 0:
            cume_time_pct = np.round(summed_cume_time / total_program_hours * 100, 2)
        else:
            cume_time_pct = np.zeros(len(semester_planner.all_dates_array))

        # Add the Total trace (time-based)
        # Build program label for hover: when multiple programs, show "All programs"; when one, show its name
        if len(programs_in_stars) == 1:
            total_trace_label = '<b>' + list(programs_in_stars)[0] + '</b> (Total)<br>'
        else:
            total_trace_label = '<b>All programs (Total)</b><br>'
        fig.add_trace(go.Scatter(
            x=night_indices,
            y=cume_time_pct,
            mode='lines',
            line=dict(color=all_stars[0].program_color_rgb, width=2),
            name="Total",
            hovertemplate= total_trace_label + 'Night: %{x}' + '<br>Date: ' + '%{customdata}' + '<br>Time % Complete: %{y}' + '<br>Total program time: ' + f'{total_program_hours:.1f} hours<br>' + '<extra></extra>',
            customdata=semester_planner.all_dates_array
        ))

    # Then add individual star traces (so they appear above the Total trace)
    for i in range(len(all_stars)):
        if use_time:
            y_vals = getattr(all_stars[i], 'cume_observe_time_pct', None)
            prog_for_star = getattr(all_stars[i], 'program', all_stars[i].starname)
            total_prog_hours = programmatics_cof.loc[programmatics_cof['program'] == prog_for_star, 'hours'].iloc[0] if prog_for_star in programmatics_cof['program'].values else 0.0
            if y_vals is None:
                # Individual stars: compute from cume_observe_time (hours) / program hours
                y_vals = np.round(all_stars[i].cume_observe_time / total_prog_hours * 100, 2) if total_prog_hours > 0 else np.zeros(len(semester_planner.all_dates_array))
            hovertemplate = '<b>' + str(prog_for_star) + '</b><br>Night: %{x}' + '<br>Date: ' + '%{customdata}' + '<br>Time % Complete: %{y}<br>Total program time: ' + f'{total_prog_hours:.1f} hours<br>' + '<extra></extra>'
        else:
            y_vals = all_stars[i].cume_observe_pct
            hovertemplate = 'Night: %{x}' + '<br>Date: ' + '%{customdata}' + '<br>% Complete: %{y}' + '<br># Obs Requested: ' + str(all_stars[i].total_observations_requested) + '<br>'

        fig.add_trace(go.Scatter(
            x=night_indices,
            y=y_vals,
            mode='lines',
            line=dict(color=all_stars[i].star_color_rgb, width=2),
            name=all_stars[i].starname,
            hovertemplate=hovertemplate,
            customdata=semester_planner.all_dates_array
        ))
        last_pct = float(np.round(y_vals[-1], 2)) if len(y_vals) else 0
        lines.append(str(all_stars[i].starname) + "," + str(last_pct))

    # Find the night index for "today" (current_day)
    try:
        today_night_index = semester_planner.all_dates_array.index(semester_planner.current_day)
    except (ValueError, AttributeError):
        # Fallback to today_starting_night if available, otherwise use 0
        today_night_index = getattr(semester_planner, 'today_starting_night', 0) - 1

    fig.add_vrect(
            x0=today_night_index,
            x1=today_night_index,
            annotation_text="Today",
            line_dash="dash",
            fillcolor=None,
            line_width=2,
            line_color='black',
            annotation_position="bottom left"
        )
    
    # X-axis: ticks every 23 days, plus the last day (matching birdseye)
    x_tick_step = 23
    x_tickvals = list(range(0, semester_planner.semester_length, x_tick_step))
    if (semester_planner.semester_length - 1) not in x_tickvals:
        x_tickvals.append(semester_planner.semester_length - 1)
    x_ticktext = [str(val + 1) for val in x_tickvals]  # Night indices (1-indexed for display, matching birdseye)
    
    # Create calendar date labels for secondary x-axis (top axis)
    # Format dates as "Feb<br>01" (month and day on separate lines)
    from datetime import datetime
    x_ticktext_dates = []
    for day_idx in x_tickvals:
        if day_idx < len(semester_planner.all_dates_array):
            date_str = semester_planner.all_dates_array[day_idx]
            # Parse date and format as "Feb<br>01" using HTML break tag
            date_obj = datetime.strptime(date_str, '%Y-%m-%d')
            month = date_obj.strftime('%b')
            day = date_obj.strftime('%d')
            x_ticktext_dates.append(f'{month}<br>{day}')
        else:
            x_ticktext_dates.append('')
    
    # Calculate legend height based on number of traces
    num_traces = len(all_stars) + 2  # +2 for "Even Burn Rate" and "Total"
    legend_height = min(300, max(150, num_traces * 25))  # Between 150-300px, 25px per trace
    
    yaxis_title = "Time % Complete (vs program hours)" if use_time else "Request % Complete"
    fig.update_layout(
        width=1400,
        height=1000,
        xaxis_title="Night in Semester",
        yaxis_title=yaxis_title,
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
            tickvals=x_tickvals,
            ticktext=x_ticktext,
            tickmode='array',
            showgrid=False,
            zeroline=False,
            anchor='y',
            side='bottom',
            range=[0, semester_planner.semester_length - 1],  # Explicitly set range
        ),
        xaxis2=dict(
            title='',
            tickvals=x_tickvals,
            ticktext=x_ticktext_dates,
            tickmode='array',
            showgrid=False,
            side='top',
            overlaying='x',
            tickfont=dict(size=labelsize - 6),
            showticklabels=True,
            range=[0, semester_planner.semester_length - 1],  # Match primary x-axis range
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize-4),
            showgrid=False,
            zeroline=False
        ),
        margin=dict(b=200, t=100)  # Bottom margin for legend below, top margin for date labels
    )
    
    # Add an invisible trace AFTER layout to force the secondary x-axis to appear
    # This trace must be associated with xaxis='x2' to make the secondary axis visible
    fig.add_trace(go.Scatter(
        x=[0, len(semester_planner.all_dates_array) - 1],
        y=[100, 100],  # Position at top of y-axis range
        mode='markers',
        marker=dict(size=0.01, opacity=0),
        showlegend=False,
        hoverinfo='skip',
        xaxis='x2',
        name='',  # Empty name to prevent legend entry
    ))
    
    # Explicitly hide any trace with xaxis='x2' or empty name from the legend
    for trace in fig.data:
        if hasattr(trace, 'xaxis') and str(trace.xaxis) == 'x2':
            trace.update(showlegend=False)
        if hasattr(trace, 'name') and (trace.name == '' or trace.name is None):
            trace.update(showlegend=False)
    
    return fig

def get_birdseye(semester_planner, availablity, all_stars):
    '''
    Produce the plotly figure showing the day/slot matrix intersection for a selection of stars

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        availability (array): a 2D array of N_slots by N_nights, binary 1/0, it is the intersection of is_alloc and is_night
        all_stars (array): a array of StarPlotter objects

    Returns:
        fig (plotly figure): a plotly figure showing the day/slot matrix intersection for a selection of stars
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

    add_grid_lines = False # this takes a long time to plot. Might not be necessary/worth it. 
    if add_grid_lines:
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
    
    # Create calendar date labels for secondary x-axis (top axis)
    # Format dates as "Jan<br>15" or "Aug<br>12" (month and day on separate lines)
    from datetime import datetime
    x_ticktext_dates = []
    for day_idx in x_tickvals:
        if day_idx < len(semester_planner.all_dates_array):
            date_str = semester_planner.all_dates_array[day_idx]
            # Parse date and format as "Jan<br>15" or "Aug<br>12" using HTML break tag
            date_obj = datetime.strptime(date_str, '%Y-%m-%d')
            month = date_obj.strftime('%b')
            day = date_obj.strftime('%d')
            x_ticktext_dates.append(f'{month}<br>{day}')
        else:
            x_ticktext_dates.append('')

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
    
    # Add an invisible trace to force the secondary x-axis to appear
    # This trace must be associated with xaxis='x2' to make the secondary axis visible
    n_slots = int(24 * 60 // semester_planner.slot_size)
    fig.add_trace(go.Scatter(
        x=[0, len(semester_planner.all_dates_array) - 1],
        y=[n_slots + 1, n_slots + 1],  # Position just above visible area
        mode='markers',
        marker=dict(size=0.01, opacity=0),
        showlegend=False,
        legendgroup=None,
        hoverinfo='skip',
        xaxis='x2',
        name='',  # Empty name to prevent legend entry
    ))
    
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
            anchor='y',
            side='bottom',
            range=[0, semester_planner.semester_length - 1],  # Explicitly set range
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
        xaxis2=dict(
            title='',
            tickvals=x_tickvals,
            ticktext=x_ticktext_dates,
            tickmode='array',
            showgrid=False,
            side='top',
            overlaying='x',
            tickfont=dict(size=labelsize - 6),
            showticklabels=True,
            range=[0, semester_planner.semester_length - 1],  # Match primary x-axis range
        ),
        margin=dict(b=200, t=100)  # Bottom margin for legend below, top margin for date labels
    )
    return fig

def get_tau_inter_line(semester_planner, all_stars, use_program_colors=False):
    """
    Produce a plotly figure showing requested vs on sky inter-night cadences, grouped by star name.

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        all_stars (array): a array of StarPlotter objects
        use_program_colors (bool): If True, use program_color_rgb; if False, use star_color_rgb (default: False)

    Returns:
        fig (plotly figure): a plotly figure showing requested vs on sky inter-night cadences, grouped by star name.
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
        idx_array = np.array(indices)
        x_vals = all_request_tau_inters[idx_array]
        y_vals = all_onsky_tau_inters[idx_array]
        text_vals = [f"{all_starnames[i]} in {all_programs[i]}" for i in indices]
        color_vals = all_colors[idx_array].tolist()  # Convert to list for Plotly
        maxyvals.append(np.max(y_vals))
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
    if maxyvals == []:
        max_val = 0
    else:
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
        width=1400,
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

def get_rawobs(semester_planner, all_stars, use_program_colors=False):
    '''
    Produce a plotly figure showing a scatter plot of observation counts for each star.
    X-axis: total requested observations
    Y-axis: sum of past and scheduled observations
    Each point represents one StarPlotter object.
    
    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        all_stars (array): an array of StarPlotter objects
        use_program_colors (bool): If True, use program_color_rgb; if False, use star_color_rgb (default: False)
    
    Returns:
        fig (plotly figure): a plotly figure showing observation counts as a scatter plot
    '''
    
    fig = go.Figure()
    fig.update_layout(plot_bgcolor=clear, paper_bgcolor=clear)
    
    # Prepare data for each star
    starnames = []
    total_requested = []
    past_obs = []
    future_obs = []
    total_completed = []  # past + scheduled
    pct_complete = []
    star_colors = []
    
    for star in all_stars:
        starnames.append(star.starname)
        total = star.total_observations_requested
        
        # Sum past observations
        past_total = sum(star.observations_past.values()) if star.observations_past else 0
        
        # Sum future observations
        future_total = sum(star.observations_future.values()) if star.observations_future else 0
        
        total_completed_val = past_total + future_total
        
        total_requested.append(total)
        past_obs.append(past_total)
        future_obs.append(future_total)
        total_completed.append(total_completed_val)
        
        # Choose color based on flag
        if use_program_colors:
            star_colors.append(star.program_color_rgb)
        else:
            star_colors.append(star.star_color_rgb)
        
        # Calculate percentage complete
        if total > 0:
            pct_complete.append((total_completed_val / total) * 100)
        else:
            pct_complete.append(0)
    
    # Create one trace per star so they can be toggled on/off in legend
    for i, star in enumerate(all_stars):
        fig.add_trace(go.Scatter(
            x=[total_requested[i]],
            y=[total_completed[i]],
            mode='markers',
            marker=dict(
                size=10,
                color=star_colors[i],  # Use each star's individual color
                opacity=0.7,
            ),
            name=starnames[i],  # Star name for legend (allows toggling)
            text=[starnames[i]],  # Star name for hover
            hovertemplate='<b>%{text}</b><br>' +
                          'Total Requested: %{x}<br>' +
                          'Past: %{customdata[0]}<br>' +
                          'Scheduled: %{customdata[1]}<br>' +
                          'Total (Past + Scheduled): %{y}<br>' +
                          '% Complete: %{customdata[2]:.1f}%<extra></extra>',
            customdata=[[past_obs[i], future_obs[i], pct_complete[i]]],
        ))
    
    # Add diagonal lines for reference (y = x for 100% complete, y = 0.5x for 50% complete)
    # For log scale, we need to use log values
    min_val = min(min(total_requested) if total_requested else 1, min(total_completed) if total_completed else 1)
    max_val = max(max(total_requested) if total_requested else 1, max(total_completed) if total_completed else 1)
    # Ensure min_val is at least 1 for log scale
    if min_val < 1:
        min_val = 1
    
    # Add 100% complete reference line (y = x) - solid black line
    fig.add_trace(go.Scatter(
        x=[min_val, max_val],
        y=[min_val, max_val],
        mode='lines',
        line=dict(color='black', width=1, dash='solid'),
        name='100% Complete',
        showlegend=False,  # Hide reference line from legend
        hovertemplate='100% Complete Reference Line<extra></extra>',
    ))
    
    # Add 50% complete reference line (y = 0.5x)
    fig.add_trace(go.Scatter(
        x=[min_val, max_val],
        y=[min_val * 0.5, max_val * 0.5],
        mode='lines',
        line=dict(color='gray', width=1, dash='dash'),
        name='50% Complete',
        showlegend=False,  # Hide reference line from legend
        hovertemplate='50% Complete Reference Line<extra></extra>',
    ))
    
    # Add annotation at the top explaining the reference lines
    fig.add_annotation(
        x=0.5,  # Center horizontally
        y=1.02,  # Just above the plot
        xref='paper',
        yref='paper',
        text="solid = 1:1<br>dashed = 1:2",
        showarrow=False,
        font=dict(size=labelsize-8, color='black'),
        align='center',
    )
    
    fig.update_layout(
        width=1400,
        height=800,
        xaxis_title="Total Requested Observations",
        yaxis_title="Total Observations (Past + Scheduled)",
        template='plotly_white',
        showlegend=True,  # Show legend so stars can be toggled on/off
        xaxis=dict(
            type="log",  # Log scale for x-axis
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize-4),
            showgrid=True,
            gridcolor='lightgray',
            minor=dict(
                showgrid=False,  # Hide minor grid lines
                ticks="",  # Hide minor tick marks
            ),
            dtick=1,  # Major ticks at powers of 10
        ),
        yaxis=dict(
            type="log",  # Log scale for y-axis
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize-4),
            showgrid=True,
            gridcolor='lightgray',
            minor=dict(
                showgrid=False,  # Hide minor grid lines
                ticks="",  # Hide minor tick marks
            ),
            dtick=1,  # Major ticks at powers of 10
        ),
        margin=dict(b=100, t=50),
    )
    
    return fig

def get_timebar(semester_planner, all_stars, use_program_colors=False, prevent_negative=False):
    """
    Create a horizontal bar chart of the time used vs forecasted vs available

    Parameters:
        semester_planner: the semester planner object
        all_stars (list): array of StarPlotter objects
        use_program_colors (bool): If True, use program_color_rgb; if False, use star_color_rgb (default: False)
        prevent_negative (bool): If True, set Incomplete and Not used categories to zero if they are negative (default: True)

    Returns:
        fig (plotly figure): a plotly figure showing the time used vs forecasted vs available as a horizontal bar chart
    """
    programmatics = pd.read_csv(os.path.join(semester_planner.semester_directory, 'programs.csv'))

    # Accumulate total times across all stars
    total_past = 0
    total_future = 0
    total_incomplete = 0
    total_requested_hours = 0
    
    programs_used = []
    for starobj in all_stars:
        # Past: day-by-day sum of (exposure time) + (readout) + (slew) per visit
        # Per date, visits = observations_past[date]: exposure = exptime * n_exp * visits; readout = readout_overhead * (n_exp - 1) * visits; slew = slew_overhead * visits
        for visits in starobj.observations_past.values():
            total_past += visits * (starobj.exptime * starobj.n_exp + readout_overhead * (starobj.n_exp - 1) + slew_overhead)
        # Future: same day-by-day formula (a) exposures*visits, (b) readout*(n_exp-1)*visits, (c) slew*visits
        for visits in starobj.observations_future.values():
            total_future += visits * (starobj.exptime * starobj.n_exp + readout_overhead * (starobj.n_exp - 1) + slew_overhead)
        total_requested_hours += starobj.total_requested_hours
        programs_used.append(starobj.program)
    
    # Convert to hours for better readability
    total_past_hours = total_past / 3600
    total_future_hours = total_future / 3600
    total_incomplete_hours = total_requested_hours - total_past_hours - total_future_hours

    if len(programs_used) > 1:
        program_rows = programmatics[programmatics['program'].isin(programs_used)]
        total_allocated_hours = program_rows['hours'].sum()
        total_allocated_nights = program_rows['nights'].sum()
    else:
        program_rows = programmatics[programmatics['program'] == programs_used[0]]
        total_allocated_hours = program_rows['hours'].sum()
        total_allocated_nights = program_rows['nights'].sum()

    # Calculate unused hours
    unused_hours = total_allocated_hours - total_future_hours - total_past_hours
    
    # Apply negative value prevention if enabled
    if prevent_negative:
        total_incomplete_hours = max(0, total_incomplete_hours)
        unused_hours = max(0, unused_hours)

    # Create bar chart data
    # Reverse order so bars appear top to bottom: Requested, Completed, Scheduled, Incomplete, Not used, Sum
    # Labels include descriptions for clarity
    labels = [
        "<b>Unused Time</b><br>(allocation - past - future)<br>If you have positive unused time, <br>consider adding or changing requests", 
        '<b>Incomplete Time</b><br>(requested - past - future)<br>If you have incomplete time, <br>some of your requests are infeasible <br> consider changing them, <br> i.e. cadence or redistributing', 
        '<b>Future Scheduled Time</b>', 
        '<b>Past Completed Time</b>', 
        '<b>Requested Time</b>'
    ]
    sum_hours = total_past_hours + total_future_hours + total_incomplete_hours + unused_hours
    values = [unused_hours, total_incomplete_hours, total_future_hours, total_past_hours, total_requested_hours]
    colors = ['#FF0000', '#F18F01', '#A23B72', '#2E86AB', '#00FF00']  # Red, Orange, Purple, Blue, Green
    
    # Create the horizontal bar chart
    # Calculate percentages based on total allocated hours for all bars
    text_labels = []
    for i, (label, val) in enumerate(zip(labels, values)):
        # Calculate percentage relative to total allocated hours
        pct = (val / total_allocated_hours * 100) if total_allocated_hours > 0 else 0
        text_labels.append(f'{val:.1f} hrs ({pct:.1f}%)')
    
    fig = go.Figure(data=[go.Bar(
        x=values,
        y=labels,
        orientation='h',
        marker=dict(color=colors),
        text=text_labels,
        textposition='auto',
        hovertemplate='<b>%{y}</b><br>%{x:.2f} hours<br><extra></extra>',
    )])
    
    # Adjust margin if there's a warning to display
    top_margin = 180 if total_requested_hours > total_allocated_hours else 130
    
    fig.update_layout(
        title_text=f'<b>Total Requested:</b> {total_requested_hours:.1f} hours ≈ {total_requested_hours/hours_per_night:.1f} nights<br><b>Total Allocated:</b> {total_allocated_hours:.1f} hours = {total_allocated_nights:.1f} nights ----> w/ losses = {total_allocated_nights*0.75:.1f} nights <br>Requested time is measured in hours. Allocated time is measured in nights. Conversion is 12 hours per night.<br>All bars include exposure times and standard overheads.',
        template='plotly_white',
        showlegend=False,
        height=710,  # Increased height for more vertical spacing between labels
        width=1400,
        margin=dict(t=top_margin, b=50, l=200, r=50),
        bargap=0.2,
        xaxis=dict(
            title='Hours',
            titlefont=dict(size=14),
            tickfont=dict(size=12)
        ),
        yaxis=dict(
            title='',
            titlefont=dict(size=14),
            tickfont=dict(size=11)  
        )
    )
    
    # Add black vertical dashed line at total_allocated_hours
    fig.add_shape(
        type="line",
        x0=total_allocated_hours,
        x1=total_allocated_hours,
        y0=-0.5,
        y1=len(labels) - 0.5,
        line=dict(color="black", width=2, dash="dash"),
        xref="x",
        yref="y"
    )
    
    # Add gray vertical dashed line for weather loss factor
    weather_loss_factor = 0.2
    fig.add_shape(
        type="line",
        x0=total_allocated_hours - total_allocated_hours * weather_loss_factor,
        x1=total_allocated_hours - total_allocated_hours * weather_loss_factor,
        y0=-0.5,
        y1=len(labels) - 0.5,
        line=dict(color="gray", width=2, dash="dash"),
        xref="x",
        yref="y"
    )
    
    # Add gray vertical dashed line at total_allocated_hours * throttle_grace
    grace_factor = semester_planner.throttle_grace
    fig.add_shape(
        type="line",
        x0=total_allocated_hours * grace_factor,
        x1=total_allocated_hours * grace_factor,
        y0=-0.5,
        y1=len(labels) - 0.5,
        line=dict(color="gray", width=2, dash="dash"),
        xref="x",
        yref="y"
        )

    # Add invisible scatter trace for hover text on the allocated time line
    # Use the same categorical labels as the bar chart to avoid numeric y-axis ticks
    fig.add_trace(go.Scatter(
        x=[total_allocated_hours] * len(labels),
        y=labels,  # Use categorical labels instead of numeric positions
        mode='markers',
        marker=dict(size=20, opacity=0),  # Invisible but hoverable markers
        hovertemplate=f'<b>Allocated Time</b><br>{total_allocated_hours:.2f} hours<br>This line represents the total allocated time for your program<extra></extra>',
        hoverlabel=dict(bgcolor='black', font_color='white'),
        showlegend=False
    ))
    
    # Add invisible scatter trace for hover text on the weather loss factor line
    weather_loss_value = total_allocated_hours - total_allocated_hours * weather_loss_factor
    fig.add_trace(go.Scatter(
        x=[weather_loss_value] * len(labels),
        y=labels,  # Use categorical labels instead of numeric positions
        mode='markers',
        marker=dict(size=20, opacity=0),  # Invisible but hoverable markers
        hovertemplate=f'<b>Weather Loss Factor</b><br>{weather_loss_value:.2f} hours<br>Allocated time minus {weather_loss_factor*100:.0f}% weather loss<br>This is only a first order estimate based on historical losses.<extra></extra>',
        hoverlabel=dict(bgcolor='gray', font_color='white'),
        showlegend=False
    ))
    
    # Add invisible scatter trace for hover text on the throttle grace line
    grace_value = total_allocated_hours * grace_factor
    fig.add_trace(go.Scatter(
        x=[grace_value] * len(labels),
        y=labels,  # Use categorical labels instead of numeric positions
        mode='markers',
        marker=dict(size=20, opacity=0),  # Invisible but hoverable markers
        hovertemplate=f'<b>Maximum Schedulable Time</b><br>{grace_value:.2f} hours<br>We allow for over-filled requests by a factor of up to {grace_factor:.2f} your allocation<br>Algorithmically, you are forbidden from getting more time than this.<extra></extra>',
        hoverlabel=dict(bgcolor='gray', font_color='white'),
        showlegend=False
    ))
    
    # Add warning annotation if requested time exceeds allocated time
    if total_requested_hours > total_allocated_hours*1.1:
        fig.add_annotation(
            text='<b>You have requested more time than you are allocated.</b>',
            xref='paper', yref='paper',
            x=0.5, y=1.35,
            showarrow=False,
            font=dict(size=18, color='red'),
            xanchor='center',
            yanchor='middle'
        )
    
    return fig


def get_timebar_by_program(semester_planner, programs_dict, prevent_negative=False):
    """
    Create a grid of horizontal bar charts showing time breakdown for each program individually
    
    Each program displays 5 bars: Unused, Incomplete, Future Scheduled, Past Completed, and Requested.
    A dashed vertical line represents their total allocated time.
    Programs are arranged in a grid with 3 columns.
    All bars use the same scale for easy comparison across programs.

    Parameters:
        semester_planner: the semester planner object
        programs_dict (dict): dictionary mapping program codes to lists of StarPlotter objects (e.g., data_astroq[0])
        prevent_negative (bool): If True, set Incomplete and Not used categories to zero if they are negative (default: False)

    Returns:
        fig (plotly figure): a plotly figure showing time breakdown per program as a grid of horizontal bar charts
    """
    programmatics = pd.read_csv(os.path.join(semester_planner.semester_directory, 'programs.csv'))
    
    # Get all programs from programs.csv
    all_programs_in_csv = set(programmatics['program'].unique())
    programs_with_requests = set(programs_dict.keys())
    
    # Find programs in CSV that don't have any requests
    programs_without_requests = all_programs_in_csv - programs_with_requests
    
    # Combine all programs: those with requests and those without
    all_program_codes = sorted(list(programs_with_requests) + list(programs_without_requests))
    
    # Store data for each program
    program_data = {}
    max_x_value = 0  # Track maximum x value for consistent scaling
    
    # Process programs with requests
    for program_code in sorted(programs_with_requests):
        program_stars = programs_dict[program_code]
        
        # Calculate times for this program (same logic as get_timebar)
        total_past = 0
        total_future = 0
        total_requested_hours = 0
        
        for starobj in program_stars:
            # Past: day-by-day sum of (exposure) + (readout) + (slew) per visit; per date: visits * (exptime*n_exp + readout*(n_exp-1) + slew)
            for visits in starobj.observations_past.values():
                total_past += visits * (starobj.exptime * starobj.n_exp + readout_overhead * (starobj.n_exp - 1) + slew_overhead)
            # Future: same day-by-day formula
            for visits in starobj.observations_future.values():
                total_future += visits * (starobj.exptime * starobj.n_exp + readout_overhead * (starobj.n_exp - 1) + slew_overhead)
            total_requested_hours += starobj.total_requested_hours
        
        # Convert to hours
        total_past_hours = total_past / 3600
        total_future_hours = total_future / 3600
        total_incomplete_hours = total_requested_hours - total_past_hours - total_future_hours
        
        # Get allocated hours for this program
        program_row = programmatics[programmatics['program'] == program_code]
        if len(program_row) > 0:
            total_allocated_hours = program_row['hours'].sum()
        else:
            total_allocated_hours = 0
        
        # Calculate unused hours
        unused_hours = total_allocated_hours - total_future_hours - total_past_hours
        
        # Apply negative value prevention if enabled
        if prevent_negative:
            total_incomplete_hours = max(0, total_incomplete_hours)
            unused_hours = max(0, unused_hours)
        
        program_data[program_code] = {
            'unused': unused_hours,
            'incomplete': total_incomplete_hours,
            'future': total_future_hours,
            'past': total_past_hours,
            'requested': total_requested_hours,
            'allocated': total_allocated_hours
        }
        
        # Update max value for scaling
        max_x_value = max(max_x_value, total_requested_hours, total_allocated_hours, 
                         unused_hours, total_incomplete_hours, total_future_hours, total_past_hours)
    
    # Process programs without requests (all bars = 0, but show allocated time)
    for program_code in sorted(programs_without_requests):
        # Get allocated hours for this program from programs.csv
        program_row = programmatics[programmatics['program'] == program_code]
        if len(program_row) > 0:
            total_allocated_hours = program_row['hours'].sum()
        else:
            total_allocated_hours = 0
        
        # All values are zero for programs with no requests
        program_data[program_code] = {
            'unused': total_allocated_hours,  # All allocated time is unused
            'incomplete': 0,
            'future': 0,
            'past': 0,
            'requested': 0,
            'allocated': total_allocated_hours
        }
        
        # Update max value for scaling
        max_x_value = max(max_x_value, total_allocated_hours)
    
    # Calculate grid dimensions: 3 columns, as many rows as needed
    num_programs = len(all_program_codes)
    num_cols = 3
    num_rows = (num_programs + num_cols - 1) // num_cols  # Ceiling division
    
    # Create subplots grid
    fig = make_subplots(
        rows=num_rows,
        cols=num_cols,
        subplot_titles=[f"<b>{prog}</b>" for prog in all_program_codes],
        horizontal_spacing=0.15,
        vertical_spacing=0.12
    )
    
    # Colors in display order: Red, Orange, Purple, Blue, Green
    display_colors = ['#FF0000', '#F18F01', '#A23B72', '#2E86AB', '#00FF00']
    category_names = ['Unused', 'Incomplete', 'Future Scheduled', 'Past Completed', 'Requested']
    
    # Add bars for each program in its own subplot
    for idx, program_code in enumerate(all_program_codes):
        data = program_data[program_code]
        
        # Calculate row and column position (1-indexed)
        row = (idx // num_cols) + 1
        col = (idx % num_cols) + 1
        
        # Prepare bar data for this program
        program_values = [data['unused'], data['incomplete'], data['future'], data['past'], data['requested']]
        
        # Add bars to this subplot
        fig.add_trace(
            go.Bar(
                x=program_values,
                y=category_names,
                orientation='h',
                marker=dict(color=display_colors),
                text=[f'{v:.1f}' if v > 0 else '' for v in program_values],
                textposition='auto',
                hovertemplate=f'<b>{program_code}</b><br>%{{y}}<br>%{{x:.2f}} hours<extra></extra>',
                showlegend=False
            ),
            row=row,
            col=col
        )
        
        # Add vertical dashed line for allocated time
        allocated = data['allocated']
        # For subplots, determine the correct axis reference
        # In make_subplots, axes are numbered: x, x2, x3, ... and y, y2, y3, ...
        if idx == 0:
            xref, yref = "x", "y"
        else:
            xref, yref = f"x{idx+1}", f"y{idx+1}"
        
        fig.add_shape(
            type="line",
            x0=allocated,
            x1=allocated,
            y0=-0.5,
            y1=4.5,
            line=dict(color="black", width=2, dash="dash"),
            xref=xref,
            yref=yref
        )

        # Add gray vertical dashed line at allocated * throttle_grace
        weather_loss_factor = 0.2
        fig.add_shape(
            type="line",
            x0=allocated - allocated * weather_loss_factor,
            x1=allocated - allocated * weather_loss_factor,
            y0=-0.5,
            y1=4.5,
            line=dict(color="gray", width=2, dash="dash"),
            xref=xref,
            yref=yref
        )
        
        # Add gray vertical dashed line at allocated * throttle_grace
        grace_factor = semester_planner.throttle_grace
        fig.add_shape(
            type="line",
            x0=allocated * grace_factor,
            x1=allocated * grace_factor,
            y0=-0.5,
            y1=4.5,
            line=dict(color="gray", width=2, dash="dash"),
            xref=xref,
            yref=yref
        )
        
        # Add invisible scatter for hover on allocated line
        fig.add_trace(
            go.Scatter(
                x=[allocated],
                y=[category_names[2]],  # Middle bar (Future Scheduled)
                mode='markers',
                marker=dict(size=15, opacity=0),
                hovertemplate=f'<b>{program_code} Allocated Time</b><br>{allocated:.2f} hours<br>Total allocated time for this program<extra></extra>',
                hoverlabel=dict(bgcolor='black', font_color='white'),
                showlegend=False
            ),
            row=row,
            col=col
        )
        
        # Add invisible scatter for hover on weather loss line
        weather_loss_value = allocated - allocated * weather_loss_factor
        fig.add_trace(
            go.Scatter(
                x=[weather_loss_value],
                y=[category_names[2]],  # Middle bar (Future Scheduled)
                mode='markers',
                marker=dict(size=15, opacity=0),
                hovertemplate=f'<b>{program_code} Weather Loss Factor</b><br>{weather_loss_value:.2f} hours<br>Allocated time minus {weather_loss_factor*100:.0f}% weather loss<extra></extra>',
                hoverlabel=dict(bgcolor='gray', font_color='white'),
                showlegend=False
            ),
            row=row,
            col=col
        )
        
        # Add invisible scatter for hover on throttle grace line
        grace_value = allocated * grace_factor
        fig.add_trace(
            go.Scatter(
                x=[grace_value],
                y=[category_names[2]],  # Middle bar (Future Scheduled)
                mode='markers',
                marker=dict(size=15, opacity=0),
                hovertemplate=f'<b>{program_code} Throttle Grace</b><br>{grace_value:.2f} hours<br>Allocated time times throttle grace factor ({grace_factor:.2f})<extra></extra>',
                hoverlabel=dict(bgcolor='gray', font_color='white'),
                showlegend=False
            ),
            row=row,
            col=col
        )
        
        # Update x-axis for this subplot (scaled to this program's data)
        # Include allocated*grace and weather loss so the gray lines are visible when they exceed the bars
        weather_loss_value = allocated - allocated * weather_loss_factor
        program_max = max(data['unused'], data['incomplete'], data['future'], 
                         data['past'], data['requested'], data['allocated'],
                         allocated * grace_factor, weather_loss_value)
        program_max = max(program_max, 1.0)  # Ensure at least 1.0 to avoid empty scale
        
        fig.update_xaxes(
            title='Hours',
            range=[0, program_max * 1.1],
            row=row,
            col=col
        )
        
        # Update y-axis for this subplot (no labels)
        fig.update_yaxes(
            title='',
            showticklabels=False,
            row=row,
            col=col
        )
    
    # Update overall layout
    fig.update_layout(
        title_text="<b>Time Breakdown by Program</b><br>Each program shows 5 bars (top to bottom): Requested (green), Past Completed (blue), Future Scheduled (purple), Incomplete (orange), Unused (red)<br>Dashed vertical line represents total allocated time. Note each grid is on its own scaling.",
        template='plotly_white',
        showlegend=False,
        height=max(600, num_rows * 250),
        width=1400,
        margin=dict(t=150, b=50, l=50, r=50)
        )
    
    return fig


def compute_seasonality(semester_planner, starnames, ras, decs):
    """
    Compute the number of days a RA/Dec point is observable in the semester using Access object

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
    
    # Temporarily override the allocation file path and request frame in the access object
    original_allocation_file = semester_planner.access_obj.allocation_file
    original_request_frame = semester_planner.access_obj.request_frame
    original_targets = semester_planner.access_obj.targets
    original_ntargets = semester_planner.access_obj.ntargets
    
    semester_planner.access_obj.allocation_file = twilight_allocation_file
    semester_planner.access_obj.request_frame = temp_requests_frame
    # Recompute targets and ntargets for the new request frame
    coords = SkyCoord(temp_requests_frame.ra * u.deg, temp_requests_frame.dec * u.deg, frame='icrs')
    semester_planner.access_obj.targets = apl.FixedTarget(name=temp_requests_frame.unique_id, coord=coords)
    semester_planner.access_obj.ntargets = len(temp_requests_frame)
   
    # Create dummy allocation for if the try statement fails.
    is_alloc = np.ones((len(starnames), semester_planner.semester_length, semester_planner.n_slots_in_night), dtype=bool)
    try:
        # Use Access object to produce the ultimate map with our custom requests frame
        access_record = semester_planner.access_obj.produce_ultimate_map(running_backup_stars=True)
        is_alloc = access_record.is_alloc
    finally:
        # Restore the original allocation file path and request frame
        semester_planner.access_obj.allocation_file = original_allocation_file
        semester_planner.access_obj.request_frame = original_request_frame
        semester_planner.access_obj.targets = original_targets
        semester_planner.access_obj.ntargets = original_ntargets
    
    # Extract is_altaz and is_moon arrays
    is_altaz = access_record.is_altaz
    is_moon = access_record.is_moon
    
    ntargets = len(starnames)
    nnights = semester_planner.semester_length
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
    Produce a plotly figure showing the sky map of the locations of requests with a static heatmap background and interactive star points.

    Parameters:
        semester_planner: the semester planner object
        all_stars (list): array of StarPlotter objects
        use_program_colors (bool): If True, use program_color_rgb; if False, use star_color_rgb (default: False)

    Returns:
        fig (plotly figure): a plotly figure showing the sky map of the locations of requests with a static heatmap background and interactive star points.
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

    n_points = n_dec * n_ra
    grid_stars = {
        'starname': [f'noname_{i}' for i in range(n_points)],
        'program_code': [f'noprog_{i}' for i in range(n_points)],
        'ra': RA_grid.flatten(),
        'dec': DEC_grid.flatten(),
        'exptime': np.full(n_points, 300),
        'n_exp': np.full(n_points, 1),
        'n_intra_max': np.full(n_points, 1),
        'n_intra_min': np.full(n_points, 1),
        'n_inter_max': np.full(n_points, 1),
        'tau_inter': np.full(n_points, 1),
        'tau_intra': np.full(n_points, 1)
    }
    grid_frame = pd.DataFrame(grid_stars)
    
    # Check if cached grid data AND background image exist
    semester = semester_planner.semester_start_date[:4] + semester_planner.semester_letter
    cache_grids_file = f"{DATADIR}/{semester}_sky_grids.npz"
    cache_image_file = f"{DATADIR}/{semester}_sky_availability_image.txt"
    
    # Try to load cached grid arrays (RA_grid, DEC_grid, NIGHTS_grid)
    if os.path.exists(cache_grids_file):
        # Load pre-computed grids from cache (FAST - skips griddata interpolation)
        cached_data = np.load(cache_grids_file)
        RA_grid = cached_data['RA_grid']
        DEC_grid = cached_data['DEC_grid']
        NIGHTS_grid = cached_data['NIGHTS_grid']
    else:
        # Need to compute the grids from scratch
        # First compute seasonality for the grid points
        grid_frame['nights_observable'] = compute_seasonality(semester_planner, grid_frame['starname'], grid_frame['ra'], grid_frame['dec'])
        
        # Perform griddata interpolation
        from scipy.interpolate import griddata
        NIGHTS_grid = griddata(
            points=(grid_frame.ra, grid_frame.dec),
            values=grid_frame.nights_observable,
            xi=(RA_grid, DEC_grid),
            method='linear'
        )
        
        # Cache the grid arrays for next time
        np.savez(cache_grids_file, RA_grid=RA_grid, DEC_grid=DEC_grid, NIGHTS_grid=NIGHTS_grid)
    
    # Try to load cached image (fastest path - skips matplotlib rendering)
    if os.path.exists(cache_image_file):
        with open(cache_image_file, 'r') as f:
            img_base64 = f.read()
    else:
        # Need to generate and cache the matplotlib image
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
        
        # Cache the base64 image for next time
        with open(cache_image_file, 'w') as f:
            f.write(img_base64)

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
                size=10
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
        width=1400,
        height=800,
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
        filtered_frame (pd.DataFrame): filtered request frame with only the specified stars
    """
    # Extract starnames from the StarPlotter objects
    starids = [star.unique_id for star in all_stars]

    # Filter the request frame to only include the specified stars
    filtered_frame = semester_planner.requests_frame_all[
        semester_planner.requests_frame_all['unique_id'].isin(starids)
    ].copy()

    return filtered_frame


def add_star_links(request_df, semester_code, date, band):
    """
    Convert starname column to links: /semester/date/band/program_code/starname

    Args:
        request_df (pd.DataFrame): request frame with starname and program_code columns
        semester_code (str): e.g. 2025B
        date (str): e.g. 2025-01-15
        band (str): e.g. band1

    Returns:
        request_df (pd.DataFrame): df with starname as HTML links
    """
    from urllib.parse import quote
    if 'program_code' not in request_df.columns or 'starname' not in request_df.columns:
        return request_df.copy()
    df = request_df.copy()
    df['starname'] = df.apply(
        lambda row: f'<a href="/{semester_code}/{date}/{band}/{quote(str(row["program_code"]))}/{quote(str(row["starname"]))}">{row["starname"]}</a>',
        axis=1
    )
    return df

def get_ladder(data, tonight_start_time):
    """Produce a plotly figure which illustrates the night plan solution.

    Args:
        data (obj): a TTP data object containing the schedule information

    Returns:
        fig (plotly figure): a plotly figure illustrating the night plan solution.
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
    fig = px.scatter(orderData, x='Minutes the from Start of the Night', y='human_starname', hover_data=['First Available', 'Last Available', 'Exposure Time (min)', "N_shots", "Total Exp Time (min)", 'UTC Start Time'] ,title='Night Plan', width=800, height=1000) #color='Program'
    # Hide the y-axis label
    fig.update_layout(yaxis_title='')
    fig.add_shape(type="rect", x0=-100, x1=-80, y0=-0.5, y1=0.5, fillcolor='red', showlegend=True, name='Exposure')
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

    # Get the x-axis range
    x_min = 0
    if len(orderData) > 0:
        # Calculate the maximum end time (start + duration) across all observations
        end_times = orderData['Start Exposure'] + orderData["Total Exp Time (min)"]
        x_max = end_times.max()
    else:
        x_max = 600
    fig.update_layout(xaxis_range=[x_min, x_max])
    # Add secondary x-axis with UTC time
    start_time = tonight_start_time.to_datetime()
    # Create tick positions (every 60 minutes or so, adjust as needed)
    tick_interval = 60  # minutes
    tick_positions = list(range(0, int(x_max) + tick_interval, tick_interval))
    tick_labels = [(start_time + timedelta(minutes=pos)).strftime('%H:%M') for pos in tick_positions]
    # Add secondary x-axis
    # Add an invisible trace to force the secondary axis to appear
    fig.add_trace(go.Scatter(
        x=[x_min, x_max],
        y=["Starname","Starname"],  # Place just below the visible range
        mode='markers',
        marker=dict(size=0.1, opacity=0),
        showlegend=False,
        hoverinfo='skip',
        xaxis='x2'
    ))
    # Create the secondary x-axis configuration
    fig.update_layout(
        xaxis2=dict(
            title=dict(text='UTC Time', standoff=0),
            overlaying='x',
            side='top',
            range=[x_min, x_max],
            tickmode='array',
            tickvals=tick_positions,
            ticktext=tick_labels,
            showgrid=False,
            showline=True,
            mirror=True
        )
    )
    
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

def get_slew_animation_plotly(data, request_selected_path, animationStep=120):
    """Create a Plotly animated polar plot showing telescope slew path during observations.

    Args:
        data: TTP data containing the schedule information
        request_selected_path: Path to request_selected.csv file
        animationStep (int): the time, in seconds, between animation frames. Default to 120s.
        
    Returns:
        fig (plotly figure): an interactive animated figure with play/pause controls
    """

    model = data[0]

    # Read the request_selected.csv file
    request_selected_df = pd.read_csv(request_selected_path)

    # Set up animation times
    t = np.arange(model.nightstarts.jd, model.nightends.jd, TimeDelta(animationStep, format='sec').jd)
    t = Time(t, format='jd')

    # Get list of astropy target objects in scheduled order (OPTIMIZED - use dict lookup)
    star_dict = {s.name: s.target for s in model.stars}
    names = list(model.schedule['Starname'])
    list_targets = [star_dict[name] for name in names]

    # Compute alt/az of each target at each time
    AZ = model.observatory.observer.altaz(t, list_targets, grid_times_targets=True)
    alt = np.round(AZ.az.rad, 2)
    az = 90 - np.round(AZ.alt.deg, 2)

    # Telescope slew path
    stamps = [0] * len(t)
    slewPath = createTelSlewPath(stamps, model.schedule['Time'], list_targets)
    AZ1 = model.observatory.observer.altaz(t, slewPath, grid_times_targets=False)
    tel_az = np.round(AZ1.az.rad, 2)
    tel_zen = 90 - np.round(AZ1.alt.deg, 2)

    # Pre-compute arrays
    schedule_times = np.array(model.schedule['Time'])
    schedule_names = np.array(model.schedule['Starname'])
    names_array = np.array(names)
    
    # Cross-match unique_id (names_array) with starname from request_selected_df
    # Create a dictionary for fast lookup
    unique_id_to_starname = dict(zip(request_selected_df['unique_id'].astype(str), 
                                      request_selected_df['starname']))
    # Map each unique_id in names_array to its human-readable starname
    human_starname_array = np.array([unique_id_to_starname.get(str(uid), str(uid)) 
                                      for uid in names_array])

    # Create telescope limit zones (red areas) - matching matplotlib version exactly
    # In polar plot: r represents zenith distance (90 - altitude), where 0 is zenith and 90 is horizon
    plotlowlim = 90  # Plotting limit to match matplotlib version
    
    # Deck limit zone (specific azimuth range)
    theta_deck = np.linspace(model.observatory.deckAzLim1, model.observatory.deckAzLim2, 100)
    r_deck_lower = np.full(100, 90 - model.observatory.deckAltLim)
    r_deck_upper = np.full(100, plotlowlim - model.observatory.vigLim)
    
    # Vignetting limit (all around)
    theta_all = np.linspace(0, 360, 100)
    r_vig_lower = np.full(100, 90 - model.observatory.vigLim)
    r_vig_upper = np.full(100, plotlowlim)
    
    # Zenith limit (all around, near center)
    r_zen_lower = np.full(100, 90 - model.observatory.zenLim)
    r_zen_upper = np.full(100, 0)

    # Create frames for animation
    frames = []
    for i in range(len(t)):
        # Determine which stars have been observed
        wasObserved = schedule_times <= float(t[i].jd)
        observed_list = schedule_names[wasObserved]
        is_observed = np.isin(names_array, observed_list)
        
        # Create frame data
        frame_data = [
            # Telescope limit zones (red areas) - using fill='toself' with closed paths
            # Deck limit
            go.Scatterpolar(
                r=np.concatenate([r_deck_lower, r_deck_upper[::-1], [r_deck_lower[0]]]),
                theta=np.concatenate([theta_deck, theta_deck[::-1], [theta_deck[0]]]),
                fill='toself',
                fillcolor='rgba(255, 0, 0, 0.7)',
                line=dict(color='rgba(255, 0, 0, 0)'),
                showlegend=(i==0),
                name='Deck Limit',
                hoverinfo='skip'
            ),
            # Vignetting limit
            go.Scatterpolar(
                r=np.concatenate([r_vig_lower, r_vig_upper[::-1], [r_vig_lower[0]]]),
                theta=np.concatenate([theta_all, theta_all[::-1], [theta_all[0]]]),
                fill='toself',
                fillcolor='rgba(255, 0, 0, 0.7)',
                line=dict(color='rgba(255, 0, 0, 0)'),
                showlegend=(i==0),
                name='Vignetting Limit',
                hoverinfo='skip'
            ),
            # Zenith limit
            go.Scatterpolar(
                r=np.concatenate([r_zen_lower, r_zen_upper[::-1], [r_zen_lower[0]]]),
                theta=np.concatenate([theta_all, theta_all[::-1], [theta_all[0]]]),
                fill='toself',
                fillcolor='rgba(255, 0, 0, 0.7)',
                line=dict(color='rgba(255, 0, 0, 0)'),
                showlegend=(i==0),
                name='Zenith Limit',
                hoverinfo='skip'
            ),
            # Stars - observed
            go.Scatterpolar(
                r=az[:, i][is_observed],
                theta=np.degrees(alt[:, i][is_observed]),
                mode='markers',
                marker=dict(size=10, color='orange', symbol='star'),
                name='Observed',
                showlegend=(i==0),
                text=human_starname_array[is_observed],
                hovertemplate='<b>%{text}</b><br>Az: %{theta:.1f}°<br>ZD: %{r:.1f}°<extra></extra>'
            ),
            # Stars - not observed
            go.Scatterpolar(
                r=az[:, i][~is_observed],
                theta=np.degrees(alt[:, i][~is_observed]),
                mode='markers',
                marker=dict(size=10, color='white', symbol='star'),
                name='Scheduled',
                showlegend=(i==0),
                text=human_starname_array[~is_observed],
                hovertemplate='<b>%{text}</b><br>Az: %{theta:.1f}°<br>ZD: %{r:.1f}°<extra></extra>'
            ),
            # Telescope path
            go.Scatterpolar(
                r=tel_zen[:i+1] if i > 0 else tel_zen[:1],
                theta=np.degrees(tel_az[:i+1] if i > 0 else tel_az[:1]),
                mode='lines',
                line=dict(color='orange', width=2),
                name='Telescope Path',
                showlegend=(i==0)
            )
        ]
        
        frames.append(go.Frame(data=frame_data, name=str(i)))

    # Create initial figure with first frame
    fig = go.Figure(
        data=frames[0].data if frames else [],
        frames=frames
    )

    # Update layout for polar plot
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                range=[0, 90],  # 0=zenith (90° altitude), 90=horizon (0° altitude)
                showticklabels=False,  # Hide altitude degree labels
                ticks='',
                showline=False,  # Hide the radial axis line
                gridcolor='rgba(255, 255, 255, 0.2)',  # Soft white grid lines
                gridwidth=1
            ),
            angularaxis=dict(
                direction='counterclockwise',
                rotation=90,
                gridcolor='rgba(255, 255, 255, 0.2)',  # Soft white grid lines
                gridwidth=1,
                tickfont=dict(size=18, color='black'),  # Bigger azimuthal labels in black
                showticklabels=True
            ),
            bgcolor='black'
        ),
        # Add cardinal direction annotations
        annotations=[
            dict(text='<b>N</b>', x=0.495, y=1.1, xref='paper', yref='paper', 
                 showarrow=False, font=dict(size=22, color='black')),
            dict(text='<b>W</b>', x=1.0, y=0.5, xref='paper', yref='paper', 
                 showarrow=False, font=dict(size=22, color='black')),
            dict(text='<b>S</b>', x=0.495, y=-0.1, xref='paper', yref='paper', 
                 showarrow=False, font=dict(size=22, color='black')),
            dict(text='<b>E</b>', x=-0.0, y=0.5, xref='paper', yref='paper', 
                 showarrow=False, font=dict(size=22, color='black'))
        ],
        # Set default animation settings to match Play button
        transition={'duration': 0},
        updatemenus=[{
            'type': 'buttons',
            'showactive': False,
            'direction': 'left',
            'x': 0.35,
            'y': -0.2,
            'xanchor': 'left',
            'yanchor': 'bottom',
            'buttons': [
                {
                    'label': '  ▶ Play  ',
                    'method': 'animate',
                    'args': [None, {
                        'frame': {'duration': 100, 'redraw': True},
                        'fromcurrent': True,
                        'mode': 'immediate',
                        'transition': {'duration': 0}
                    }]
                },
                {
                    'label': '  ⏸ Pause  ',
                    'method': 'animate',
                    'args': [[None], {
                        'frame': {'duration': 0, 'redraw': False},
                        'mode': 'immediate',
                        'transition': {'duration': 0}
                    }]
                }
            ],
            'bgcolor': 'white',
            'bordercolor': 'black',
            'borderwidth': 2,
            'font': {'size': 16, 'color': 'black', 'family': 'Arial'},
        }],
        sliders=[{
            'active': 0,
            'yanchor': 'top',
            'y': -0.15,
            'xanchor': 'left',
            'currentvalue': {
                'prefix': 'Time: ',
                'visible': True,
                'xanchor': 'right',
                'font': {'size': 14, 'color': 'black'}
            },
            'pad': {'b': 10, 't': 50},
            'len': 0.9,
            'x': 0.1,
            'font': {'size': 12, 'color': 'black'},
            'steps': [
                {
                    'args': [[f.name], {
                        'frame': {'duration': 100, 'redraw': True},  # Match button duration
                        'mode': 'immediate',
                        'transition': {'duration': 0}
                    }],
                    'label': t[k].datetime.strftime('%H:%M'),  # HH:MM format only
                    'method': 'animate'
                }
                for k, f in enumerate(frames)
            ],
            'transition': {'duration': 100}  # Match frame duration for consistent speed
        }],
        width=800,
        height=800,
        title=dict(
            text='Telescope Slew Animation',
            font=dict(color='black', size=20)
        ),
        template='plotly_white',
        paper_bgcolor='white',
        plot_bgcolor='white',
        font=dict(color='black'),
        # Configure animation behavior
        hovermode='closest'
    )

    return fig

def get_script_plan(night_planner):
    """Generate script plan DataFrame from semester planner and night planner objects.
    
    This function reads the request_selected.csv file from the semester planner's output directory,
    merges it with the night planner's solution data, and returns a properly formatted DataFrame
    with the same column structure as the original get_script_plan function.
    
    Args:
        night_planner: NightPlanner object containing solution attribute
        
    Returns:
        final_df (pd.DataFrame): a formatted observing plan DataFrame
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
    # desired_columns = [
    #     'Start Exposure', 'unique_id', 'starname', 'program_code', 'ra', 'dec', 
    #     'exptime', 'n_exp', 'n_intra_max', 'tau_intra', 'weather_band_1', 'weather_band_2', 'weather_band_3', 'teff', 
    #     'jmag', 'gmag', 'epoch', 'gaia_id', 'First Available', 'Last Available'
    # ]
    desired_columns = [
         'First Available', 'Start Exposure', 'Last Available', 'unique_id', 'starname', 'program_code', 'ra', 'dec', 
        'exptime', 'n_exp', 'n_intra_max', 'tau_intra', 'jmag', 'gmag',]
    
    # Keep only the columns that exist in the merged dataframe
    available_columns = [col for col in desired_columns if col in merged_df.columns]
    
    # Reorder columns to match the desired structure
    final_df = merged_df[available_columns].copy()
    
    # Round numeric fields to appropriate decimal places
    if 'ra' in final_df.columns:
        # Ensure ra is numeric before rounding, handle 'None' strings
        final_df['ra'] = final_df['ra'].replace('None', pd.NA)
        final_df['ra'] = pd.to_numeric(final_df['ra'], errors='coerce').round(1)
    
    if 'dec' in final_df.columns:
        # Ensure dec is numeric before rounding, handle 'None' strings
        final_df['dec'] = final_df['dec'].replace('None', pd.NA)
        final_df['dec'] = pd.to_numeric(final_df['dec'], errors='coerce').round(1)
    
    if 'jmag' in final_df.columns:
        # Ensure jmag is numeric before rounding, handle 'None' strings
        final_df['jmag'] = final_df['jmag'].replace('None', pd.NA)
        final_df['jmag'] = pd.to_numeric(final_df['jmag'], errors='coerce').round(1)
    
    if 'gmag' in final_df.columns:
        # Ensure gmag is numeric before rounding, handle 'None' strings
        final_df['gmag'] = final_df['gmag'].replace('None', pd.NA)
        final_df['gmag'] = pd.to_numeric(final_df['gmag'], errors='coerce').round(1)

    # if 'teff' in final_df.columns:
    #     # Ensure teff is numeric before rounding, handle 'None' strings
    #     final_df['teff'] = final_df['teff'].replace('None', pd.NA)
    #     final_df['teff'] = pd.to_numeric(final_df['teff'], errors='coerce').round(0)
    
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
    
    # Ensure DataFrame is clean and properly structured for DataTables
    final_df = final_df.reset_index(drop=True)
    # Remove duplicate column names if any exist
    final_df = final_df.loc[:, ~final_df.columns.duplicated(keep='first')]
    # Fill NaN values with empty strings to ensure consistent structure
    final_df = final_df.fillna('')
    # Ensure all columns have consistent data types (convert objects to strings)
    for col in final_df.columns:
        if final_df[col].dtype == 'object':
            final_df[col] = final_df[col].astype(str).replace('nan', '').replace('None', '').replace('', '')
    
    return final_df

def plot_path_2D_interactive(data, night_start_time=None):
    """Create an interactive Plotly plot showing telescope azimuth and altitude paths with UTC times and white background.
    
    Args:
        data (list): a list containing the TTP model solution
        night_start_time: Astropy Time object representing the start of night (Minute 0) from allocation file

    Returns:
        fig (plotly figure): an interactive plot showing telescope azimuth and altitude paths with UTC times and white background.
    """

    model = data[0]

    names = list(model.plotly['human_starname'])
    times = model.times
    az_path = model.az_path
    alt_path = model.alt_path
    wrap = model.observatory.wrapLimitAngle
    
    # Use night_start_time as "Minute 0" reference
    if night_start_time is None:
        # Fallback: use first time from model
        night_start_jd = times[0].jd
    else:
        night_start_jd = night_start_time.jd

    # Use times from model for telescope path (these are waypoints from TTP solver)
    # These represent END of exposure times, so we need to subtract exposure duration
    obs_time = np.array([t.jd for t in times])
    
    # Adjust times to represent START of exposure instead of END
    # The times array has 2 points per observation (both at end of exposure)
    if hasattr(model, 'plotly') and 'Total Exp Time (min)' in model.plotly:
        total_exp_times = model.plotly['Total Exp Time (min)']
        
        # Subtract exposure time from each pair of points
        for i in range(len(total_exp_times)):
            idx1 = i * 2
            idx2 = i * 2 + 1
            if idx2 < len(obs_time):
                duration_days = total_exp_times[i] / 1440.0
                obs_time[idx1] -= duration_days
                obs_time[idx2] -= duration_days
    
    # The times/az/alt arrays have 2 points per observation (start and end)
    # but names only has one entry per observation. Expand names to match.
    if len(obs_time) == 2 * len(names):
        # Each name should appear twice (for start and end of observation)
        expanded_names = []
        for name in names:
            expanded_names.append(name)  # Start point
            expanded_names.append(name)  # End point
        names = expanded_names
    elif len(obs_time) != len(names):
        # Repeat names to match the length
        names = names * (len(obs_time) // len(names) + 1)
        names = names[:len(obs_time)]
    
    # Ensure all arrays are the same length
    min_len = min(len(obs_time), len(az_path), len(alt_path), len(names))
    obs_time = obs_time[:min_len]
    az_path = np.array(az_path[:min_len])
    alt_path = np.array(alt_path[:min_len])
    names = names[:min_len]
    
    # First, ensure all azimuth values are within 0-360 degrees using mod 360
    # This prevents values like -10° or 370° from appearing
    az_path = np.mod(az_path, 360)
    
    # Store original azimuth for hover text (0-360°)
    az_path_original = az_path.copy()
    
    # For display: values above 270° should be shown as negative (subtract 360)
    # e.g., 290° becomes -70°, 350° becomes -10°
    az_path_display = az_path.copy()
    az_path_display[az_path_display > 270] -= 360
    
    # For tick labels and hover, format as HH:MM
    time_labels = [Time(t, format='jd').isot[11:16] for t in obs_time]
    
    # Create hover text arrays using ORIGINAL azimuth (0-360°)
    hover_text_az = [f"Time: {time_labels[i]}<br>Target: {names[i]}<br>Az: {az_path_original[i]:.1f}°" 
                     for i in range(len(obs_time))]
    hover_text_alt = [f"Time: {time_labels[i]}<br>Target: {names[i]}<br>Alt: {alt_path[i]:.1f}°" 
                      for i in range(len(obs_time))]

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        subplot_titles=("Azimuth Path", "Elevation Path"),
        vertical_spacing=0.1
    )

    # Azimuth plot (use display values: >270° shown as negative)
    fig.add_trace(go.Scatter(
        x=obs_time, y=az_path_display,
        mode='lines+markers',
        marker=dict(color='indigo'),
        name='Azimuth',
        text=hover_text_az,
        hovertemplate='%{text}<extra></extra>'
    ), row=1, col=1)

    # Elevation plot
    fig.add_trace(go.Scatter(
        x=obs_time, y=alt_path,
        mode='lines+markers',
        marker=dict(color='seagreen'),
        name='Elevation',
        text=hover_text_alt,
        hovertemplate='%{text}<extra></extra>'
    ), row=2, col=1)

    # Add wrap limit line at the wrap position
    # Apply same conversion as display values: if wrap > 270°, subtract 360
    if wrap is not None:
        # Normalize wrap to 0-360 range
        wrap_normalized = wrap % 360
        
        # Apply same display conversion: if > 270°, show as negative
        wrap_display = wrap_normalized
        if wrap_display > 270:
            wrap_display -= 360
        
        # Draw wrap limit line at the display position
        fig.add_shape(
            type="line",
            x0=obs_time[0], x1=obs_time[-1],
            y0=wrap_display, y1=wrap_display,
            line=dict(color="red", dash="dash", width=2),
            row=1, col=1
        )
        fig.add_annotation(
            x=obs_time[-1], y=wrap_display,
            text=f"Wrap = {wrap_normalized}°",
            showarrow=False,
            font=dict(color="red", size=10),
            row=1, col=1
        )

    # Highlight observed intervals using Start/Stop Exposure times
    # Shade from Start Exposure to Stop Exposure (both in minutes from night start)
    if hasattr(model, 'plotly') and 'Start Exposure' in model.plotly and 'Stop Exposure' in model.plotly:
        start_exposures = model.plotly['Start Exposure']  # Minutes from start of night
        stop_exposures = model.plotly['Stop Exposure']    # Minutes from start of night
        
        for i, (start_min, stop_min) in enumerate(zip(start_exposures, stop_exposures)):
            # Convert minutes from night start to JD
            # 1 day = 1440 minutes, so minutes / 1440 = fraction of a day
            start_jd = night_start_jd + (start_min / 1440.0)
            stop_jd = night_start_jd + (stop_min / 1440.0)
            
            # Add yellow shaded region for this exposure
            fig.add_vrect(
                x0=start_jd, x1=stop_jd,
                fillcolor="yellow", opacity=0.3,
                layer="below", line_width=0,
                row=1, col=1
            )
            fig.add_vrect(
                x0=start_jd, x1=stop_jd,
                fillcolor="yellow", opacity=0.3,
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

    # Update y-axis for azimuth to always show consistent range with 270° at the top
    # Values >270° are displayed as negative (by subtracting 360)
    # Range goes from -95° to 275° with 5° buffer on both ends
    az_y_min = -95
    az_y_max = 275
    
    # Generate tick positions every 45 degrees from -90 to 270
    tick_interval = 45
    az_tick_positions = np.arange(-90, 271, tick_interval)  # -90, -45, 0, 45, 90, 135, 180, 225, 270
    
    # Create labels - convert negative angles to their 360° equivalents
    # -90° → 270°, -45° → 315°, etc.
    az_tick_labels = []
    for pos in az_tick_positions:
        if pos < 0:
            # Convert negative to 360° equivalent
            label = int(pos + 360)
        else:
            label = int(pos)
        az_tick_labels.append(f"{label}°")
    
    fig.update_yaxes(
        tickmode='array',
        tickvals=az_tick_positions,
        ticktext=az_tick_labels,
        range=[az_y_min, az_y_max],
        title_text="Azimuth (deg)",
        row=1, col=1
    )
    
    # Update y-axis for altitude to always show 0° to 90°
    fig.update_yaxes(
        range=[0, 90],
        title_text="Altitude (deg)",
        row=2, col=1
    )
    
    fig.update_layout(
        height=600,
        width=1000,
        template="plotly_white"
    )
    return fig

REQUEST_FRAME_COLUMNS = [
    'starname', 'unique_id', 'program_code', 'ra', 'dec', 'exptime', 'n_exp',
    'n_inter_max', 'tau_inter', 'n_intra_max', 'n_intra_min', 'tau_intra',
    'weather_band_1', 'weather_band_2', 'weather_band_3', 'inactive'
]
BOOLEAN_COLUMNS = {'weather_band_1': 'Band1', 'weather_band_2': 'Band2', 'weather_band_3': 'Band3', 'inactive': 'Inactive'}
REQUEST_FRAME_DISPLAY_NAMES = {
    'starname': 'Star', 'unique_id': 'ID', 'program_code': 'Program',
    'ra': 'RA', 'dec': 'Dec', 'exptime': 'ExpTime'
}
# Tooltips shown when hovering over column headers (add your custom text here)
REQUEST_FRAME_COLUMN_TOOLTIPS = {
    'Star': 'Name of the star',
    'ID': 'Keck OB database unique ID',
    'Program': 'Program Code',
    'RA': 'RA in decimal degrees',
    'Dec': 'Declination in decimal degrees',
    'ExpTime': 'Exposure time in seconds',
    'n_exp': 'Number of Exposures per Visit',
    'n_inter_max': 'Maximum number of unique nights to observe the star',
    'tau_inter': 'The minimum inter-night cadence between unique night observations',
    'n_intra_max': 'The desired number of visits to the star in each night it is observed',
    'n_intra_min': 'The accepted minimum number of visits to the star in each night it is observed',
    'tau_intra': 'The minimum intra-night cadence between visits within a night in hours',
    'Band1': 'Allowed to observe in Band1?',
    'Band2': 'Allowed to observe in Band2?',
    'Band3': 'Allowed to observe in Band3?',
    'Inactive': 'Is the star set to inactive?',
}


def request_frame_to_html(request_df, semester_code=None, date=None, band=None, table_id='request-table', page_size=25):
    """
    Convert a request frame (from request.csv) to HTML for admin/program/star pages.

    Displays only: starname, unique_id, program_code, ra, dec, exptime, n_exp,
    n_inter_max, tau_inter, n_intra_max, n_intra_min, tau_intra, Band1, Band2, Band3, Inactive.
    Boolean columns (weather bands, inactive) are shown as Y/N with transparent green/red.

    Args:
        request_df (pd.DataFrame): request frame, e.g. from get_request_frame
        semester_code (str, optional): for star links
        date (str, optional): for star links
        band (str, optional): for star links
        table_id (str): HTML table id
        page_size (int): rows per page

    Returns:
        str: HTML string with table and DataTables
    """
    import re
    from urllib.parse import quote
    df = request_df.copy()
    df = df.reset_index(drop=True)
    # Select only columns we want, in order; ignore any extra columns
    cols = [c for c in REQUEST_FRAME_COLUMNS if c in df.columns]
    df = df[cols].copy()
    df = df.fillna('')
    # Round RA and Dec to 2 decimals
    for coord in ('ra', 'dec'):
        if coord in df.columns:
            df[coord] = pd.to_numeric(df[coord], errors='coerce')
            df[coord] = df[coord].apply(lambda x: f'{x:.2f}' if pd.notna(x) else '')
    # Add star links if URL context provided
    if semester_code and date and band and 'program_code' in df.columns and 'starname' in df.columns:
        df['starname'] = df.apply(
            lambda row: f'<a href="/{semester_code}/{date}/{band}/{quote(str(row["program_code"]))}/{quote(str(row["starname"]))}">{row["starname"]}</a>',
            axis=1
        )
    # Convert boolean columns to Y/N with color
    green_bg = 'rgba(34, 139, 34, 0.25)'
    red_bg = 'rgba(220, 53, 69, 0.25)'

    def _is_true(val):
        if pd.isna(val) or val == '':
            return False
        s = str(val).lower()
        if s in ('true', '1', 'yes'):
            return True
        if s in ('false', '0', 'no'):
            return False
        try:
            return bool(float(val))
        except (ValueError, TypeError):
            return False

    for orig in BOOLEAN_COLUMNS:
        if orig not in df.columns:
            continue
        is_inactive = orig == 'inactive'
        def _cell(val, _inactive=is_inactive):
            truth = _is_true(val)
            if _inactive:
                y_n, bg = ('Y', red_bg) if truth else ('N', green_bg)
            else:
                y_n, bg = ('Y', green_bg) if truth else ('N', red_bg)
            return f'<span style="background:{bg};padding:2px 6px;border-radius:4px;">{y_n}</span>'
        df[orig] = df[orig].apply(lambda v: _cell(v))
    # Rename columns for display
    df = df.rename(columns={**BOOLEAN_COLUMNS, **REQUEST_FRAME_DISPLAY_NAMES})
    # Ensure object columns are strings
    for col in df.columns:
        if df[col].dtype == 'object':
            s = df[col].astype(str).replace('nan', '').replace('None', '')
            if col not in BOOLEAN_COLUMNS.values():  # Don't overwrite our HTML
                df[col] = s
    table_html = df.to_html(classes='table table-striped table-hover', index=False, escape=False, table_id=table_id)
    # Add data-tooltip to column headers (custom CSS tooltip, shows immediately)
    from html import escape
    tooltips = [REQUEST_FRAME_COLUMN_TOOLTIPS.get(col, '') for col in df.columns]
    def _add_th_tooltip(m):
        idx = _add_th_tooltip.idx
        _add_th_tooltip.idx += 1
        t = tooltips[idx] if idx < len(tooltips) else ''
        return f'<th data-tooltip="{escape(t)}">{m.group(1)}</th>' if t else m.group(0)
    _add_th_tooltip.idx = 0
    table_html = re.sub(r'<th>([^<]*)</th>', _add_th_tooltip, table_html, count=len(df.columns))
    # Add tfoot for column filter dropdowns
    tfoot_cells = ''.join(['<th></th>' for _ in df.columns])
    table_html = table_html.replace('</tbody>', '</tbody><tfoot><tr>' + tfoot_cells + '</tr></tfoot>')
    # Column widths = longest value in column (content-based)
    def visible_len(s):
        return len(re.sub(r'<[^>]+>', '', str(s)).strip())
    widths = {}
    band_cols = ('Band1', 'Band2', 'Band3', 'Inactive')
    no_padding_cols = ('n_inter_max', 'tau_inter', 'n_intra_max', 'n_intra_min', 'tau_intra')
    for col in df.columns:
        content_max = max((visible_len(c) for c in df[col]), default=0)
        header_len = len(str(col))
        pad = 3 if col == 'n_exp' else (0 if col in no_padding_cols else 2)
        ch_width = max(content_max, header_len, 1) + pad
        if col in band_cols:
            ch_width = max(ch_width, 3)  # Y/N box needs ~3ch with padding
        widths[col] = f'{ch_width}ch'
    column_defs = [f"{{ targets: {i}, width: '{widths[col]}' }}" for i, col in enumerate(df.columns)]
    column_defs_str = ',\n                '.join(column_defs)
    # CSS column widths (Star=1st, ID=2nd, etc.) - force narrow to override DataTables auto-sizing
    col_widths = [widths[col] for col in df.columns]
    col_css = ' '.join([f"#{table_id} th:nth-child({i+1}), #{table_id} td:nth-child({i+1}) {{ width: {w} !important; max-width: {w} !important; }}" for i, w in enumerate(col_widths)])
    custom_css = f"""
    <style>
    /* Override DataTables width:100% - table should shrink to fit column widths, not stretch to page */
    #{table_id} {{ width: auto !important; max-width: 100%; border-collapse: collapse; font-size: 21px; margin: 10px 0; table-layout: fixed !important; }}
    {col_css}
    #{table_id} thead th {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; font-weight: 600;
        padding: 3pt; text-align: center; font-size: 20px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; cursor: help; }}
    #{table_id} tbody td {{ padding: 3pt; text-align: center; font-size: 20px; border-bottom: 1px solid #e9ecef;
        overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
    #{table_id} tbody td:nth-child(7), #{table_id} thead th:nth-child(7) {{ padding: 3pt !important; }}
    #{table_id} tbody td:nth-child(8), #{table_id} tbody td:nth-child(9), #{table_id} tbody td:nth-child(10), #{table_id} tbody td:nth-child(11), #{table_id} tbody td:nth-child(12),
    #{table_id} thead th:nth-child(8), #{table_id} thead th:nth-child(9), #{table_id} thead th:nth-child(10), #{table_id} thead th:nth-child(11), #{table_id} thead th:nth-child(12) {{ padding: 0 !important; }}
    #{table_id} tbody td:nth-child(13), #{table_id} tbody td:nth-child(14), #{table_id} tbody td:nth-child(15), #{table_id} tbody td:nth-child(16),
    #{table_id} thead th:nth-child(13), #{table_id} thead th:nth-child(14), #{table_id} thead th:nth-child(15), #{table_id} thead th:nth-child(16) {{ padding: 3pt !important; }}
    #{table_id} tbody tr:nth-child(even) {{ background-color: #dee2e6 !important; }}
    #{table_id} tbody tr:nth-child(odd) {{ background-color: white !important; }}
    #{table_id} tbody tr:hover {{ background-color: #e3f2fd !important; }}
    #{table_id} tfoot th {{ padding: 8px 4px; background: #f1f3f5; border-top: 2px solid #dee2e6; }}
    #{table_id} tfoot .column-filter {{ width: 100%; padding: 4px 8px; font-size: 14px; border: 1px solid #ced4da; border-radius: 4px; box-sizing: border-box; }}
    /* Custom tooltip - bold, opaque, high visibility */
    #{table_id} thead th[data-tooltip] {{ position: relative; }}
    #{table_id} thead th[data-tooltip]:hover::after {{
        content: attr(data-tooltip);
        position: fixed;
        top: 140px;
        left: 50%;
        transform: translateX(-50%);
        padding: 16px 24px;
        background: #5a4d9e !important;
        color: white !important;
        font-size: 18px !important;
        font-weight: 700 !important;
        white-space: normal;
        max-width: 675px;
        width: 675px;
        text-align: center;
        border-radius: 8px;
        border: 3px solid #4a3d8e;
        box-shadow: 0 6px 24px rgba(0,0,0,0.5);
        z-index: 10000;
        pointer-events: none;
        opacity: 1 !important;
    }}
    </style>
    """
    # Column indices for numeric comparison: Star(0), ID(1), Program(2), RA(3), Dec(4), ExpTime(5), n_exp(6), n_inter_max(7), tau_inter(8), n_intra_max(9), n_intra_min(10), tau_intra(11), Band1-4(12-15)
    numeric_col_indices = [3, 4, 5, 6, 7, 8, 9, 10, 11]
    init_script = f"""
    <script>
    $(document).ready(function() {{
        if ($.fn.DataTable.isDataTable('#{table_id}')) {{ $('#{table_id}').DataTable().destroy(); }}
        var numericCols = {numeric_col_indices};
        var table = $('#{table_id}').DataTable({{
            autoWidth: false,
            pageLength: {page_size},
            order: [[0, 'asc']],
            dom: 'lBfrtip',
            buttons: ['copy', 'csv', 'excel', 'print'],
            lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
            columnDefs: [ {column_defs_str} ],
            initComplete: function() {{
                var api = this.api();
                var dt = this;
                $.fn.dataTable.ext.search.push(function(settings, data, dataIndex) {{
                    if (settings.nTable.id !== '{table_id}') return true;
                    var tbl = $('#' + '{table_id}');
                    for (var i = 0; i < data.length; i++) {{
                        var input = tbl.find('tfoot th').eq(i).find('input.column-filter');
                        var val = (input.val() || '').trim();
                        if (!val) continue;
                        var cellVal = data[i];
                        if (typeof cellVal === 'string' && cellVal.indexOf('<') >= 0) {{
                            cellVal = $('<div>').html(cellVal).text().trim();
                        }} else {{
                            cellVal = (cellVal || '').toString().trim();
                        }}
                        if (numericCols.indexOf(i) >= 0) {{
                            var match = val.match(/^(>=|<=|>|<|=|==)\s*(-?[\\d.]+)$/);
                            if (match) {{
                                var op = match[1];
                                var numVal = parseFloat(match[2]);
                                var cellNum = parseFloat(cellVal);
                                if (isNaN(cellNum)) return false;
                                switch(op) {{
                                    case '>': if (!(cellNum > numVal)) return false; break;
                                    case '<': if (!(cellNum < numVal)) return false; break;
                                    case '>=': if (!(cellNum >= numVal)) return false; break;
                                    case '<=': if (!(cellNum <= numVal)) return false; break;
                                    case '=':
                                    case '==': if (cellNum != numVal) return false; break;
                                }}
                            }} else {{
                                if (cellVal.toLowerCase().indexOf(val.toLowerCase()) < 0) return false;
                            }}
                        }} else {{
                            if (cellVal.toLowerCase().indexOf(val.toLowerCase()) < 0) return false;
                        }}
                    }}
                    return true;
                }});
                api.columns().every(function() {{
                    var column = this;
                    var input = $('<input type="text" class="column-filter" placeholder="Filter... (use > < >= <= for numbers)">')
                        .appendTo($(column.footer()).empty())
                        .on('keyup change', function() {{
                            api.draw();
                        }});
                }});
            }}
        }});
    }});
    </script>
    """
    return custom_css + table_html + init_script


NIGHTPLAN_COLUMNS = [
    'First Available', 'Start Exposure', 'Last Available', 'unique_id', 'starname',
    'program_code', 'ra', 'dec', 'exptime', 'n_exp', 'n_intra_max', 'tau_intra', 'jmag', 'gmag'
]
NIGHTPLAN_COLUMN_TOOLTIPS = {
    'First Available': 'First available time to observe (HH:MM). Use > < >= <= with HH:MM to filter.',
    'Start Exposure': 'Scheduled start time (HH:MM). Use > < >= <= with HH:MM to filter.',
    'Last Available': 'Last available time to observe (HH:MM). Use > < >= <= with HH:MM to filter.',
    'unique_id': 'Keck OB database unique ID',
    'starname': 'Name of the star',
    'program_code': 'Program Code',
    'ra': 'Right ascension in decimal degrees',
    'dec': 'Declination in decimal degrees',
    'exptime': 'Exposure time in seconds',
    'n_exp': 'Number of exposures per visit',
    'n_intra_max': 'Maximum intra-night visits',
    'tau_intra': 'Minimum intra-night cadence in hours',
    'jmag': 'J-band magnitude',
    'gmag': 'G-band magnitude',
}


def nightplan_table_to_html(script_df, table_id='script-table', page_size=100):
    """
    Convert nightplan script DataFrame to HTML with same styling as request_frame_to_html.

    Same colors, fonts, fontsize, filtering (partial match, numeric > < >= <=), hover tooltips.
    Displays: First Available, Start Exposure, Last Available, unique_id, starname, program_code,
    ra, dec, exptime, n_exp, n_intra_max, tau_intra, jmag, gmag.
    """
    import re
    from html import escape
    df = script_df.copy()
    df = df.reset_index(drop=True)
    cols = [c for c in NIGHTPLAN_COLUMNS if c in df.columns]
    df = df[cols].copy()
    df = df.fillna('')
    for col in df.columns:
        if df[col].dtype == 'object':
            df[col] = df[col].astype(str).replace('nan', '').replace('None', '')
    table_html = df.to_html(classes='table table-striped table-hover', index=False, escape=False, table_id=table_id)
    # Add tooltips to headers
    tooltips = [NIGHTPLAN_COLUMN_TOOLTIPS.get(col, '') for col in df.columns]
    def _add_th_tooltip(m):
        idx = _add_th_tooltip.idx
        _add_th_tooltip.idx += 1
        t = tooltips[idx] if idx < len(tooltips) else ''
        return f'<th data-tooltip="{escape(t)}">{m.group(1)}</th>' if t else m.group(0)
    _add_th_tooltip.idx = 0
    table_html = re.sub(r'<th>([^<]*)</th>', _add_th_tooltip, table_html, count=len(df.columns))
    tfoot_cells = ''.join(['<th></th>' for _ in df.columns])
    table_html = table_html.replace('</tbody>', '</tbody><tfoot><tr>' + tfoot_cells + '</tr></tfoot>')
    # Column widths
    def visible_len(s):
        return len(re.sub(r'<[^>]+>', '', str(s)).strip())
    widths = {}
    for col in df.columns:
        content_max = max((visible_len(c) for c in df[col]), default=0)
        header_len = len(str(col))
        ch_width = max(content_max, header_len, 1) + 2
        widths[col] = f'{ch_width}ch'
    column_defs = [f"{{ targets: {i}, width: '{widths[col]}' }}" for i, col in enumerate(df.columns)]
    column_defs_str = ',\n                '.join(column_defs)
    col_widths = [widths[col] for col in df.columns]
    col_css = ' '.join([f"#{table_id} th:nth-child({i+1}), #{table_id} td:nth-child({i+1}) {{ width: {w} !important; max-width: {w} !important; }}" for i, w in enumerate(col_widths)])
    custom_css = f"""
    <style>
    #{table_id} {{ width: auto !important; max-width: 100%; border-collapse: collapse; font-size: 21px; margin: 10px 0; table-layout: fixed !important; }}
    {col_css}
    #{table_id} thead th {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; font-weight: 600;
        padding: 3pt; text-align: center; font-size: 20px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; cursor: help; }}
    #{table_id} tbody td {{ padding: 3pt; text-align: center; font-size: 20px; border-bottom: 1px solid #e9ecef;
        overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
    #{table_id} tbody tr:nth-child(even) {{ background-color: #dee2e6 !important; }}
    #{table_id} tbody tr:nth-child(odd) {{ background-color: white !important; }}
    #{table_id} tbody tr:hover {{ background-color: #e3f2fd !important; }}
    #{table_id} tfoot th {{ padding: 8px 4px; background: #f1f3f5; border-top: 2px solid #dee2e6; }}
    #{table_id} tfoot .column-filter {{ width: 100%; padding: 4px 8px; font-size: 14px; border: 1px solid #ced4da; border-radius: 4px; box-sizing: border-box; }}
    #{table_id} thead th[data-tooltip] {{ position: relative; }}
    #{table_id} thead th[data-tooltip]:hover::after {{
        content: attr(data-tooltip);
        position: fixed;
        top: 140px;
        left: 50%;
        transform: translateX(-50%);
        padding: 16px 24px;
        background: #5a4d9e !important;
        color: white !important;
        font-size: 18px !important;
        font-weight: 700 !important;
        white-space: normal;
        max-width: 675px;
        width: 675px;
        text-align: center;
        border-radius: 8px;
        border: 3px solid #4a3d8e;
        box-shadow: 0 6px 24px rgba(0,0,0,0.5);
        z-index: 10000;
        pointer-events: none;
        opacity: 1 !important;
    }}
    </style>
    """
    # Numeric columns for > < >= <= : ra(6), dec(7), exptime(8), n_exp(9), n_intra_max(10), tau_intra(11), jmag(12), gmag(13)
    # Time columns (HH:MM): First Available(0), Start Exposure(1), Last Available(2)
    numeric_col_indices = [6, 7, 8, 9, 10, 11, 12, 13]
    time_col_indices = [0, 1, 2]
    init_script = f"""
    <script>
    function parseHHMM(s) {{
        var m = (s || '').match(/(\\d{{1,2}}):(\\d{{2}})/);
        return m ? parseInt(m[1], 10) * 60 + parseInt(m[2], 10) : NaN;
    }}
    $(document).ready(function() {{
        if ($.fn.DataTable.isDataTable('#{table_id}')) {{ $('#{table_id}').DataTable().destroy(); }}
        var numericCols = {numeric_col_indices};
        var timeCols = {time_col_indices};
        var table = $('#{table_id}').DataTable({{
            autoWidth: false,
            pageLength: {page_size},
            order: [[1, 'asc']],
            dom: 'lBfrtip',
            buttons: ['copy', 'csv', 'excel', 'print'],
            lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
            columnDefs: [ {column_defs_str} ],
            initComplete: function() {{
                var api = this.api();
                $.fn.dataTable.ext.search.push(function(settings, data, dataIndex) {{
                    if (settings.nTable.id !== '{table_id}') return true;
                    var tbl = $('#' + '{table_id}');
                    for (var i = 0; i < data.length; i++) {{
                        var input = tbl.find('tfoot th').eq(i).find('input.column-filter');
                        var val = (input.val() || '').trim();
                        if (!val) continue;
                        var cellVal = data[i];
                        if (typeof cellVal === 'string' && cellVal.indexOf('<') >= 0) {{
                            cellVal = $('<div>').html(cellVal).text().trim();
                        }} else {{
                            cellVal = (cellVal || '').toString().trim();
                        }}
                        if (timeCols.indexOf(i) >= 0) {{
                            var tMatch = val.match(/^(>=|<=|>|<|=|==)\s*(\\d{{1,2}}):(\\d{{2}})$/);
                            if (tMatch) {{
                                var op = tMatch[1];
                                var filterMins = parseInt(tMatch[2], 10) * 60 + parseInt(tMatch[3], 10);
                                var cellMins = parseHHMM(cellVal);
                                if (isNaN(cellMins)) return false;
                                switch(op) {{
                                    case '>': if (!(cellMins > filterMins)) return false; break;
                                    case '<': if (!(cellMins < filterMins)) return false; break;
                                    case '>=': if (!(cellMins >= filterMins)) return false; break;
                                    case '<=': if (!(cellMins <= filterMins)) return false; break;
                                    case '=':
                                    case '==': if (cellMins != filterMins) return false; break;
                                }}
                            }} else {{
                                if (cellVal.toLowerCase().indexOf(val.toLowerCase()) < 0) return false;
                            }}
                        }} else if (numericCols.indexOf(i) >= 0) {{
                            var match = val.match(/^(>=|<=|>|<|=|==)\s*(-?[\\d.]+)$/);
                            if (match) {{
                                var op = match[1];
                                var numVal = parseFloat(match[2]);
                                var cellNum = parseFloat(cellVal);
                                if (isNaN(cellNum)) return false;
                                switch(op) {{
                                    case '>': if (!(cellNum > numVal)) return false; break;
                                    case '<': if (!(cellNum < numVal)) return false; break;
                                    case '>=': if (!(cellNum >= numVal)) return false; break;
                                    case '<=': if (!(cellNum <= numVal)) return false; break;
                                    case '=':
                                    case '==': if (cellNum != numVal) return false; break;
                                }}
                            }} else {{
                                if (cellVal.toLowerCase().indexOf(val.toLowerCase()) < 0) return false;
                            }}
                        }} else {{
                            if (cellVal.toLowerCase().indexOf(val.toLowerCase()) < 0) return false;
                        }}
                    }}
                    return true;
                }});
                api.columns().every(function() {{
                    var column = this;
                    var input = $('<input type="text" class="column-filter" placeholder="Filter... (> < for HH:MM or numbers)">')
                        .appendTo($(column.footer()).empty())
                        .on('keyup change', function() {{ api.draw(); }});
                }});
            }}
        }});
    }});
    </script>
    """
    return custom_css + table_html + init_script


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
        table_html (str): HTML string with table and DataTables initialization
    """
    # Ensure DataFrame is clean and properly structured
    dataframe = dataframe.reset_index(drop=True)
    # Remove duplicate column names if any exist
    dataframe = dataframe.loc[:, ~dataframe.columns.duplicated(keep='first')]
    # Fill NaN values with empty strings
    dataframe = dataframe.fillna('')
    # Ensure all object columns are strings
    for col in dataframe.columns:
        if dataframe[col].dtype == 'object':
            dataframe[col] = dataframe[col].astype(str).replace('nan', '').replace('None', '')
    
    # Validate sort_column is within bounds
    num_columns = len(dataframe.columns)
    if sort_column >= num_columns:
        sort_column = 0  # Default to first column if out of bounds

    # if 'exptime' in dataframe.columns:
    #     # Ensure exptime is an integer, handle 'None' strings
    #     dataframe['exptime'] = dataframe['exptime'].replace('None', pd.NA)
    #     dataframe['exptime'] = pd.to_numeric(dataframe['exptime'], errors='coerce').fillna(0).astype(int)
    
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
    
    # Generate columnDefs dynamically based on actual number of columns
    num_columns = len(dataframe.columns)
    column_defs = []
    # Default width mapping for common column names
    width_map = {
        'First Available': '80px',
        'Start Exposure': '80px',
        'Last Available': '80px',
        'unique_id': '200px',
        'starname': '200px',
        'program_code': '120px',
        'ra': '100px',
        'dec': '100px',
        'exptime': '80px',
        'n_exp': '60px',
        'n_intra_max': '80px',
        'tau_intra': '80px',
        'jmag': '60px',
        'gmag': '60px'
    }
    
    for i, col in enumerate(dataframe.columns):
        width = width_map.get(col, '100px')  # Default width if not in map
        column_defs.append(f"{{ targets: {i}, width: '{width}' }}")
    
    column_defs_str = ',\n                '.join(column_defs)
    
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
            tableLayout: 'auto',  // Use auto layout for better column width handling
            columnDefs: [
                {column_defs_str}
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

