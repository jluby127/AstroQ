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
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")  # Ensures safe, headless rendering
from matplotlib.figure import Figure

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.io as pio
from astropy.time import Time
from astropy.time import TimeDelta
# import io
from io import BytesIO
import imageio.v3 as iio
import base64
from collections import defaultdict
import pickle

from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astropy.units as u
import astroplan as apl

import kpfcc.management as mn
import kpfcc.io as io_mine
import sys
sys.path.append("/Users/jack/Documents/github/ttp/ttp/")
import plotting

gray = 'rgb(210,210,210)'
clear = 'rgba(255,255,255,1)'
labelsize = 38

def get_cof(manager, all_stars):
    '''
    Return the html string for a plotly figure showing the COF for a selection of stars

    all_stars must be an array, even if only 1 element long. It is an array of StarPlotter objects, as defined in kpfcc.plot
    '''

    fig = go.Figure()
    fig.update_layout(plot_bgcolor=gray, paper_bgcolor=clear) #autosize=True,margin=dict(l=40, r=40, t=40, b=40),
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
    cof_html = pio.to_html(fig, full_html=True, include_plotlyjs='cdn')
    return cof_html

def get_birdseye(manager, availablity, all_stars):
    '''
    Return the html string for a plotly figure showing the COF for a selection of stars

    all_stars must be an array, even if only 1 element long. It is an array of StarPlotter objects, as defined in kpfcc.plot
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
    fig.update_layout(
        yaxis_title="Slot in Night",
        xaxis_title="Night in Semester",
        xaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=np.append(np.arange(0, manager.semester_length, 23), 183),
            ticktext=np.append(np.arange(0, manager.semester_length, 23), 184),
            tickmode='array',
            showgrid=False,
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=np.append(np.arange(0, manager.n_slots_in_night, 28), manager.n_slots_in_night-1),
            ticktext=[manager.n_slots_in_night - x for x in np.arange(0, manager.n_slots_in_night+1, 28)],
            tickmode='array',
            showgrid=False,
        ),
        template="plotly_white",
        showlegend=True,
        legend=dict(
            font=dict(size=labelsize-10)
        )
    )
    birdseye_html = pio.to_html(fig, full_html=True, include_plotlyjs='cdn')
    return birdseye_html

def get_tau_inter_line(manager, all_stars):
    """
    Create an interactive scatter plot grouped by program, using provided colors for each point.
    Points from the same program will share a legend entry and color.

    Parameters:
        x (list or array): X coordinates
        y (list or array): Y coordinates
        names (list of str): Labels for each point
        programs (list of str): Program name for each point
        colors (list of str): Color value for each point (e.g. hex code or named color)

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
        colors.extend([starobj.program_color_rgb] * len(onsky_diffs))

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

    # Create one trace per program
    for program, indices in program_to_indices.items():
        x_vals = [all_request_tau_inters[i] for i in indices]
        y_vals = [all_onsky_tau_inters[i] for i in indices]
        text_vals = [f"{all_starnames[i]} in {program}" for i in indices]
        color_vals = all_colors[indices[0]]  # Use the first point's color for the group

        fig.add_trace(go.Scatter(
            x=x_vals,
            y=y_vals,
            mode='markers',
            name=program,
            marker=dict(size=10, color=color_vals),
            text=text_vals,
            hovertemplate="%{text}<br>X: %{x}<br>Y: %{y}<extra></extra>"
        ))

    # Add 1-to-1 line
    min_val = min(min(x), min(y))
    max_val = max(max(x), max(y))
    fig.add_trace(go.Scatter(
        x=[min_val, max_val],
        y=[min_val, max_val],
        mode='lines',
        line=dict(color='black', dash='dash'),
        name='1-to-1 line',
        showlegend=True
    ))

    fig.update_layout(
        xaxis_title="Requested Minimum Inter-Night Cadence",
        yaxis_title="On Sky Inter-Night Cadence",
        template='plotly_white',
        height=600,
        width=800
    )
    tau_inter_line_html = pio.to_html(fig, full_html=True, include_plotlyjs='cdn')
    return tau_inter_line_html


def compute_seasonality(manager, starnames, ras, decs, semester='B', nslots=168, nnights=184, observatory='Keck Observatory', start_time='03:30', slot_size=5):
    """
    Combine all maps for a target to produce the final map

    Args:
        requests_frame (dataframe): the pandas dataframe containing request information
        running_backup_stars (bool): if true, then do not run the extra map of stepping back in time to account for the starting slot fitting into the night
    Returns:
        available_indices_for_request (dictionary): keys are the starnames and values are a 1D array
                                                  the indices where available_slots_for_request is 1.

    """
    # Prepatory work
    if semester=='A':
        start_date = '2025-02-01'
    elif semester=='B':
        start_date = '2025-08-01'
    else:
        print("invalid semester. Choose A or B.")
    date_formal = Time(start_date,format='iso',scale='utc')
    date = str(date_formal)[:10]
    ntargets = len(starnames)

    # Determine observability
    coords = apy.coordinates.SkyCoord(ras * u.deg, decs * u.deg, frame='icrs')
    targets = apl.FixedTarget(name=starnames, coord=coords)
    keck = apl.Observer.at_site(observatory)

    # Set up time grid for one night, first night of the semester
    start = date + "T" + start_time
    daily_start = Time(start, location=keck.location)
    daily_end = daily_start + TimeDelta(1.0, format='jd') # full day from start of first night
    tmp_slot_size = TimeDelta(5.0*u.min)
    t = Time(np.arange(daily_start.jd, daily_end.jd, tmp_slot_size.jd), format='jd',location=keck.location)
    t = t[np.argsort(t.sidereal_time('mean'))] # sort by lst

    # Compute base alt/az pattern, shape = (ntargets, nslots)
    coord0 = keck.altaz(t, targets, grid_times_targets=True)
    alt0 = coord0.alt.deg
    az0 = coord0.az.deg

    # 2D mask (n targets, n slots))
    is_altaz0 = np.ones_like(alt0, dtype=bool)
    is_altaz0 &= ~((5.3 < az0 ) & (az0 < 146.2) & (alt0 < 33.3)) # remove nasymth deck
    # remove min elevation for mid declination stars
    ismiddec = ((-30 < targets.dec.deg) & (targets.dec.deg < 75))
    fail = ismiddec[:,np.newaxis] & (alt0 < 30) # broadcast declination array
    is_altaz0 &= ~fail
    # all stars must be between 18 and 85 deg
    fail = (alt0 < 18) | (alt0 > 85)
    is_altaz0 &= ~fail
    # computing slot midpoint for all nights in semester 2D array (slots, nights)
    slotmidpoint0 = daily_start + (np.arange(nslots) + 0.5) *  slot_size * u.min
    # days = np.arange(manager.n_nights_in_semester) * u.day
    days = np.arange(nnights) * u.day
    slotmidpoint = (slotmidpoint0[np.newaxis,:] + days[:,np.newaxis])
    # 3D mask
    is_altaz = np.empty((ntargets, nnights, nslots),dtype=bool)

    # Pre-compute the sidereal times for interpolation
    x = t.sidereal_time('mean').value
    x_new = slotmidpoint.sidereal_time('mean').value
    idx = np.searchsorted(x, x_new, side='left')
    idx = np.clip(idx, 0, len(x)-1) # Handle edge cases
    is_altaz = is_altaz0[:,idx]

    # Compute moon accessibility
    is_moon = np.ones_like(is_altaz, dtype=bool)
    moon = apy.coordinates.get_moon(slotmidpoint[:,0] , keck.location)
    # Reshaping uses broadcasting to achieve a (ntarget, night) array
    ang_dist = apy.coordinates.angular_separation(
        targets.ra.reshape(-1,1), targets.dec.reshape(-1,1),
        moon.ra.reshape(1,-1), moon.dec.reshape(1,-1),
    ) # (ntargets)
    is_moon = is_moon & (ang_dist.to(u.deg) > 30*u.deg)[:, :, np.newaxis]

    # True if obseravtion occurs at night
    is_night = manager.twilight_map_remaining_2D.astype(bool) # shape = (nnights, nslots)
    is_night = np.ones_like(is_altaz, dtype=bool) & is_night[np.newaxis,:,:]

#     is_alloc = manager.allocation_map_2D.astype(bool) # shape = (nnights, nslots)
#     is_alloc = np.ones_like(is_altaz, dtype=bool) & is_alloc[np.newaxis,:,:] # shape = (ntargets, nnights, nslots)

    is_observable_now = np.logical_and.reduce([
        is_altaz,
        is_moon,
        is_night,
    ])

    # specify indeces of 3D observability array
    itarget, inight, islot = np.mgrid[:ntargets,:nnights,:nslots]

    # define flat table to access maps
    df = pd.DataFrame(
        {'itarget':itarget.flatten(),
         'inight':inight.flatten(),
         'islot':islot.flatten()}
    )
    # df['is_observable'] = is_observable_now.flatten()
    available_nights_onsky = []
    for itarget in range(ntargets):
        onskycount = 0
        for inight in range(nnights):
            temp = list(islot[itarget,inight,is_observable_now[itarget,inight,:]])
            if len(temp) > 0:
                onskycount += 1

        available_nights_onsky.append(onskycount)

    return available_nights_onsky


def interactive_sky_with_static_heatmap(manager, progcode):#RA_grid, DEC_grid, Z_grid, request_df):
    """
    Interactive sky map with static heatmap background and interactive star points.

    Parameters:
        RA_grid, DEC_grid, Z_grid: 2D arrays defining the heatmap in degrees.
        request_df (pd.DataFrame): Contains 'ra', 'dec', 'starname', 'program_code'.
        colors (dict): Optional color mapping by program_code.

    Returns:
        plotly.graph_objects.Figure
    """

    with open(manager.reports_directory + 'admin/' + manager.current_day + '/star_objects.pkl', 'rb') as f:
        data = pickle.load(f)

    all_stars = data[0][progcode]

    starnames = [all_stars[r].starname for r in range(len(all_stars))]
    programs = [all_stars[r].program for r in range(len(all_stars))]
    ras = [all_stars[r].ra for r in range(len(all_stars))]
    decs = [all_stars[r].dec for r in range(len(all_stars))]
    colors = [all_stars[r].program_color_rgb for r in range(len(all_stars))]
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

    available_nights_onsky_requests = compute_seasonality(manager, starnames, ras, decs, semester='B', nslots=168, nnights=184, observatory='Keck Observatory', start_time='03:30')
    grid_frame['nights_observable'] = compute_seasonality(manager, grid_frame['starname'], grid_frame['ra'], grid_frame['dec'], semester='B', nslots=168, nnights=184, observatory='Keck Observatory', start_time='03:30')

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
            title='Observable Nights',
            titleside='right',
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
            color = [group['color'][0]]*len(group)

            fig.add_trace(go.Scattergeo(
                lon=group['ra'] - 180,
                lat=group['dec'],
                mode='markers',
                name=program,
                marker=dict(size=6, color=color, opacity=1),
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
        height=600,
        width=1000
    )
    skymap_html = pio.to_html(fig, full_html=True, include_plotlyjs='cdn')
    return skymap_html





def get_ladder(manager, data):
    """Create an interactive plot which illustrates the solution.

    Args:
        orderData (dict): Formatted and ordered schedule from the
            model object
        current_day (string): Date of observations, to be included in the
            filename
        outputdir (str): Folder in which to save the html link to the plot
    """
    # with open(manager.reports_directory + 'observer/' + manager.current_day + '/ttp_data.pkl', 'rb') as f:
    #     data = pickle.load(f)
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
    fig = px.scatter(orderData, x='Minutes the from Start of the Night', y="Starname", hover_data=['First Available', 'Last Available', 'Exposure Time (min)', "N_shots", "Total Exp Time (min)"] ,title='Night Plan', width=800, height=1000) #color='Program'
    fig.add_shape(type="rect", x0=-100, x1=-80, y0=-0.5, y1=0.5, fillcolor='red', showlegend=True, name='Expose P1')
    fig.add_shape(type="rect", x0=-100, x1=-80, y0=-0.5, y1=0.5, fillcolor='blue', showlegend=True, name='Expose P3')
    fig.add_shape(type="rect", x0=-100, x1=-80, y0=-0.5, y1=0.5, fillcolor='lime', opacity=0.3, showlegend=True, name='Accessible')

    new_already_processed = []
    ifixer = 0 # for multi-visit targets, it throws off the one row per target plotting...this fixes it
    for i in range(len(orderData['Starname'])):
        if orderData['Starname'][i] not in new_already_processed:
            counter1 = list(orderData['Starname']).count(orderData['Starname'][i])
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
    ladder_html = pio.to_html(fig, full_html=True, include_plotlyjs='cdn')
    return ladder_html


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
    k = 0
    while stamps[k] == 0:
        stamps[k] = pointings[0]
        k += 1
    l = len(stamps)-1
    while stamps[l] == 0:
        stamps[l] = pointings[-1]
        l -= 1

    return stamps

def get_slew_animation(manager, data, animationStep=120):
    '''
    Produce the animation slew path GIF (in memory, no intermediate files).

    model (object) - the model object returned out from the TTP's model.py
    animationStep (int) - the time, in seconds, between animation still frames. Default to 120s.
    '''
    # with open(manager.reports_directory + 'observer/' + manager.current_day + '/ttp_data.pkl', 'rb') as f:
    #     data = pickle.load(f)
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

    # Create GIF frames in memory
    gif_frames = []
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
            # ax.scatter(alt[j][i], az[j][i], color=color, marker='*')
            ax.scatter(alt[i][j], az[i][j], color=color, marker='*')

        # Draw telescope path
        ax.plot(tel_az[:i], tel_zen[:i], color='orange')

        # Save frame to memory
        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=100)
        buf.seek(0)
        gif_frames.append(iio.imread(buf))
        buf.close()

    # Write GIF to memory
    gif_buf = BytesIO()
    iio.imwrite(gif_buf, gif_frames, format='gif', loop=0, duration=0.3)  # 0.3s per frame
    gif_buf.seek(0)

    # Encode in base64 for HTML embedding
    gif_base64 = base64.b64encode(gif_buf.read()).decode('utf-8')
    gif_buf.close()

    slew_animation_html = f'<img src="data:image/gif;base64,{gif_base64}" alt="Observing Animation"/>'
    return slew_animation_html

def get_script_plan(manager, data):

    # with open(manager.reports_directory + 'observer/' + manager.current_day + '/ttp_data.pkl', 'rb') as f:
    #     data = pickle.load(f)

    round_two_requests = [] # just for now
    nightly_start_stop_times = pd.read_csv(manager.nightly_start_stop_times_file)
    idx = nightly_start_stop_times[nightly_start_stop_times['Date'] == str(manager.current_day)].index[0]
    obs_start_time = Time(str(manager.current_day) + "T" + str(nightly_start_stop_times['Start'][idx]))
    lines = io_mine.write_starlist(manager.requests_frame, data[0].plotly, obs_start_time, data[0].extras,round_two_requests, str(manager.current_day), manager.reports_directory + 'observer/' + manager.current_day)
    observing_plan = pd.DataFrame([io_mine.parse_star_line(line) for line in lines])
    observing_plan = observing_plan[observing_plan.iloc[:, 0].str.strip() != '']
    observing_plan_html = save_interactive_observing_plan(observing_plan)
    return observing_plan_html

def plot_path_2D_interactive(manager, data):
    """Create an interactive Plotly plot showing telescope azimuth and altitude paths with UTC times and white background."""

    # with open(manager.reports_directory + 'observer/' + manager.current_day + '/ttp_data.pkl', 'rb') as f:
    #     data = pickle.load(f)
    model = data[0]

    names = list(model.schedule['Starname'])
    times = model.times
    az_path = model.az_path
    alt_path = model.alt_path
    wrap = model.observatory.wrapLimitAngle

    # Convert times to UTC string labels (HH:MM)
    obs_time = np.array([t.jd for t in times])
    time_labels = [Time(t, format='jd').isot[11:16] for t in obs_time]

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        subplot_titles=("Azimuth Path", "Elevation Path"),
        vertical_spacing=0.1
    )

    # Azimuth plot
    fig.add_trace(go.Scatter(
        x=time_labels, y=az_path,
        mode='lines+markers',
        marker=dict(color='indigo'),
        name='Azimuth',
        text=names,
        hovertemplate='Time: %{x}<br>Az: %{y}<br>Target: %{text}'
    ), row=1, col=1)

    # Elevation plot
    fig.add_trace(go.Scatter(
        x=time_labels, y=alt_path,
        mode='lines+markers',
        marker=dict(color='seagreen'),
        name='Elevation',
        text=names,
        hovertemplate='Time: %{x}<br>Alt: %{y}<br>Target: %{text}'
    ), row=2, col=1)

    # Optional: wrap line on azimuth
    if wrap is not None:
        fig.add_shape(
            type="line",
            x0=time_labels[0], x1=time_labels[-1],
            y0=wrap, y1=wrap,
            line=dict(color="red", dash="dash"),
            row=1, col=1
        )
        fig.add_annotation(
            x=time_labels[-1], y=wrap,
            text=f"Wrap = {wrap}",
            showarrow=False,
            font=dict(color="red"),
            row=1, col=1
        )

    # Highlight observed intervals
    for i in range(0, len(obs_time)-1, 2):
        fig.add_vrect(
            x0=time_labels[i], x1=time_labels[i+1],
            fillcolor="orange", opacity=0.2,
            layer="below", line_width=0,
            row=1, col=1
        )
        fig.add_vrect(
            x0=time_labels[i], x1=time_labels[i+1],
            fillcolor="orange", opacity=0.2,
            layer="below", line_width=0,
            row=2, col=1
        )

    fig.update_layout(
        height=600,
        width=1000,
        xaxis2_title="Time (UTC)",
        yaxis_title="Azimuth (deg)",
        yaxis2_title="Altitude (deg)",
        template="plotly_white"
    )
    slew_path_html = pio.to_html(fig, full_html=True, include_plotlyjs='cdn')
    return slew_path_html

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
