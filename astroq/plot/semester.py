"""Semester-scope Plotly figures."""

from astroq.plot._common import (  # noqa: F401
    cumulative_by_night,
    daily_visits_by_night,
    gray,
    clear,
    labelsize,
    hours_per_night,
    go,
    px,
    np,
    pd,
    os,
    sns,
    datetime,
    timedelta,
    html_escape,
    quote,
    base64,
    BytesIO,
    re,
    plt,
    make_subplots,
    u,
    SkyCoord,
    apl,
    ac,
    griddata,
    TimeDelta,
    _cof_pct_curve,
    _cof_group_for_star,
    _visit_denominator,
    _charged_hours_from_ps,
    _football_cache_dir,
    _visit_counts_by_date,
    programs_ledger_for_plot,
    _render_datatable,
    _TEMPLATE_ENV,
    _Path,
)
from astroq.plot.context import PlotData, PlotSelection, RequestView

def get_cof(plot_data, selection, use_time=False):
    """
    Produce a plotly figure showing the Cumulative Observability Function (COF) for a selection of stars

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        all_stars (array): a array of StarPlotter objects
        use_time (bool): if True, use the cumulative observe time percentage instead of the cumulative observe percentage

    Returns:
        fig (plotly figure): a plotly figure showing the COF for a selection of stars
    """

    all_stars = plot_data.views(selection)
    semester_planner = plot_data.semester_planner
    fig = go.Figure()
    fig.update_layout(
        plot_bgcolor=gray, paper_bgcolor=clear
    )  # autosize=True,margin=dict(l=40, r=40, t=40, b=40),

    # Convert calendar dates to night indices (0, 1, 2, ...)
    night_indices = np.arange(len(plot_data.semester_planner.access_obj.all_dates_array))

    burn_line = np.linspace(0, 100, len(plot_data.semester_planner.access_obj.all_dates_array))
    burn_line = np.round(burn_line, 2)

    # Add "Even Burn Rate" line as a shape so it's always visible and can't be toggled
    # Use add_shape to create a line that spans the entire plot
    fig.add_shape(
        type="line",
        x0=night_indices[0],
        y0=burn_line[0],
        x1=night_indices[-1],
        y1=burn_line[-1],
        line=dict(color="black", width=2, dash="dash"),
        layer="below",  # Draw below traces so it doesn't obscure data
    )

    # Add an invisible trace just for the legend entry (so users know what the line represents)
    # This trace will be visible in legend but clicking it won't hide the actual line
    fig.add_trace(
        go.Scatter(
            x=[None],  # No actual data points
            y=[None],
            mode="lines",
            line=dict(color="black", width=2, dash="dash"),
            name="Even Burn Rate",
            showlegend=True,
            hoverinfo="skip",  # Don't show hover for this dummy trace
        )
    )
    lines = []
    ps = plot_data.semester_planner.timeline
    n_nights = len(plot_data.semester_planner.access_obj.all_dates_array)
    is_programmatic = not getattr(all_stars[0], "allow_mapview", True)

    if is_programmatic:
        program_codes = {s.program for s in all_stars}
        total_sub = ps[ps["program_code"].isin(program_codes)]
    else:
        uids = [str(s.unique_id) for s in all_stars]
        total_sub = ps[ps["unique_id"].isin(uids)]

    if use_time is False:
        total_denom = sum(_visit_denominator(s) for s in all_stars)
        daily = total_sub.groupby(total_sub.index).size().reindex(
            range(n_nights), fill_value=0
        )
        cume_observe = daily.cumsum().to_numpy(dtype=float)
        if total_denom > 0:
            cume_observe_pct = np.round(cume_observe / total_denom * 100, 2)
        else:
            total_past = int(cume_observe[-1]) if len(cume_observe) else 0
            cume_observe_pct = (
                np.round(cume_observe / total_past * 100, 2)
                if total_past > 0
                else np.zeros(n_nights)
            )

        # Add the Total trace first (so it appears below other traces)
        fig.add_trace(
            go.Scatter(
                x=night_indices,
                y=cume_observe_pct,
                mode="lines",
                line=dict(color=all_stars[0].program_color_rgb, width=2),
                name="Total",
                hovertemplate="Night: %{x}"
                + "<br>Date: "
                + "%{customdata}"
                + "<br>% Complete: %{y}"
                + "<br># Visits Requested: "
                + str(int(total_denom))
                + "<br>",
                customdata=plot_data.semester_planner.access_obj.all_dates_array,
            )
        )
    else:
        # use_time=True: normalize by awarded program hours from programs.csv
        programmatics_cof = pd.read_csv(
            os.path.join(plot_data.semester_planner.config.get("global", "workdir"), "programs.csv")
        )
        programs_in_stars = {s.program for s in all_stars}
        total_program_hours = programmatics_cof.loc[
            programmatics_cof["program"].isin(programs_in_stars), "hours"
        ].sum()
        slot_size = plot_data.semester_planner.config.getfloat("semester", "slot_size")
        slots_per_hour = 60 / slot_size
        daily_slots = total_sub.groupby(total_sub.index)["t_visit_slots"].sum().reindex(
            range(n_nights), fill_value=0
        )
        cume_hours = daily_slots.cumsum().to_numpy(dtype=float) / slots_per_hour
        if total_program_hours > 0:
            cume_time_pct = np.round(cume_hours / total_program_hours * 100, 2)
        else:
            cume_time_pct = np.zeros(n_nights)

        # Build program label for hover: when multiple programs, show "All programs"; when one, show its name
        if len(programs_in_stars) == 1:
            total_trace_label = "<b>" + list(programs_in_stars)[0] + "</b> (Total)<br>"
        else:
            total_trace_label = "<b>All programs (Total)</b><br>"
        fig.add_trace(
            go.Scatter(
                x=night_indices,
                y=cume_time_pct,
                mode="lines",
                line=dict(color=all_stars[0].program_color_rgb, width=2),
                name="Total",
                hovertemplate=total_trace_label
                + "Night: %{x}"
                + "<br>Date: "
                + "%{customdata}"
                + "<br>Time charged (% of awarded hours): %{y}"
                + "<br>Total program time: "
                + f"{total_program_hours:.1f} hours<br>"
                + "<extra></extra>",
                customdata=plot_data.semester_planner.access_obj.all_dates_array,
            )
        )

    programmatics_cof = None
    if use_time:
        programmatics_cof = pd.read_csv(
            os.path.join(plot_data.semester_planner.config.get("global", "workdir"), "programs.csv")
        )

    # Then add individual star traces (so they appear above the Total trace)
    for i in range(len(all_stars)):
        group_col, group_val = _cof_group_for_star(all_stars[i])
        if use_time:
            prog_for_star = all_stars[i].program
            total_prog_hours = programmatics_cof.loc[
                programmatics_cof["program"] == prog_for_star, "hours"
            ].sum()
            y_vals = _cof_pct_curve(plot_data.semester_planner,
                ps,
                n_nights,
                group_col=group_col,
                group_val=group_val,
                use_time=True,
                denominator=total_prog_hours,
            )
            hovertemplate = (
                "<b>"
                + str(prog_for_star)
                + "</b><br>Night: %{x}"
                + "<br>Date: "
                + "%{customdata}"
                + "<br>Time charged (% of awarded hours): %{y}<br>Total program time: "
                + f"{total_prog_hours:.1f} hours<br>"
                + "<extra></extra>"
            )
        else:
            y_vals = _cof_pct_curve(plot_data.semester_planner,
                ps,
                n_nights,
                group_col=group_col,
                group_val=group_val,
                use_time=False,
                denominator=_visit_denominator(all_stars[i]),
            )
            hovertemplate = (
                "Night: %{x}"
                + "<br>Date: "
                + "%{customdata}"
                + "<br>% Complete: %{y}"
                + "<br># Visits Requested: "
                + str(_visit_denominator(all_stars[i]))
                + "<br>"
            )

        fig.add_trace(
            go.Scatter(
                x=night_indices,
                y=y_vals,
                mode="lines",
                line=dict(color=all_stars[i].star_color_rgb, width=2),
                name=all_stars[i].target,
                hovertemplate=hovertemplate,
                customdata=plot_data.semester_planner.access_obj.all_dates_array,
            )
        )
        last_pct = float(np.round(y_vals[-1], 2)) if len(y_vals) else 0
        lines.append(str(all_stars[i].target) + "," + str(last_pct))

    today_night_index = plot_data.semester_planner.access_obj.current_night_index

    fig.add_vrect(
        x0=today_night_index,
        x1=today_night_index,
        annotation_text="Today",
        line_dash="dash",
        fillcolor=None,
        line_width=2,
        line_color="black",
        annotation_position="bottom left",
    )

    # X-axis: ticks every 23 days, plus the last day (matching birdseye)
    x_tick_step = 23
    x_tickvals = list(range(0, plot_data.semester_planner.semester_length, x_tick_step))
    if (plot_data.semester_planner.semester_length - 1) not in x_tickvals:
        x_tickvals.append(plot_data.semester_planner.semester_length - 1)
    x_ticktext = [
        str(val + 1) for val in x_tickvals
    ]  # Night indices (1-indexed for display, matching birdseye)

    # Create calendar date labels for secondary x-axis (top axis)
    # Format dates as "Feb<br>01" (month and day on separate lines)
    x_ticktext_dates = []
    for day_idx in x_tickvals:
        if day_idx < len(plot_data.semester_planner.access_obj.all_dates_array):
            date_str = plot_data.semester_planner.access_obj.all_dates_array[day_idx]
            # Parse date and format as "Feb<br>01" using HTML break tag
            date_obj = datetime.strptime(date_str, "%Y-%m-%d")
            month = date_obj.strftime("%b")
            day = date_obj.strftime("%d")
            x_ticktext_dates.append(f"{month}<br>{day}")
        else:
            x_ticktext_dates.append("")

    # Calculate legend height based on number of traces
    num_traces = len(all_stars) + 2  # +2 for "Even Burn Rate" and "Total"
    legend_height = min(
        300, max(150, num_traces * 25)
    )  # Between 150-300px, 25px per trace

    yaxis_title = (
        "Time charged (% of awarded hours)" if use_time else "Visit % Complete"
    )
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
            bgcolor="rgba(255,255,255,0.7)",
            bordercolor="black",
            borderwidth=1,
            font=dict(size=labelsize - 18),
            # Standardize legend size
            itemsizing="constant",  # All legend items same size
            itemwidth=30,  # Fixed width for legend items
            # Make legend more compact
            groupclick="toggleitem",  # Click group to toggle all items
            # Standardize legend dimensions
            tracegroupgap=5,  # Gap between trace groups
            traceorder="normal",  # Keep order as traces were added
        ),
        xaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=x_tickvals,
            ticktext=x_ticktext,
            tickmode="array",
            showgrid=False,
            zeroline=False,
            anchor="y",
            side="bottom",
            range=[0, plot_data.semester_planner.semester_length - 1],  # Explicitly set range
        ),
        xaxis2=dict(
            title="",
            tickvals=x_tickvals,
            ticktext=x_ticktext_dates,
            tickmode="array",
            showgrid=False,
            side="top",
            overlaying="x",
            tickfont=dict(size=labelsize - 6),
            showticklabels=True,
            range=[
                0,
                plot_data.semester_planner.semester_length - 1,
            ],  # Match primary x-axis range
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            showgrid=False,
            zeroline=False,
        ),
        margin=dict(
            b=200, t=100
        ),  # Bottom margin for legend below, top margin for date labels
    )

    # Add an invisible trace AFTER layout to force the secondary x-axis to appear
    # This trace must be associated with xaxis='x2' to make the secondary axis visible
    fig.add_trace(
        go.Scatter(
            x=[0, len(plot_data.semester_planner.access_obj.all_dates_array) - 1],
            y=[100, 100],  # Position at top of y-axis range
            mode="markers",
            marker=dict(size=0.01, opacity=0),
            showlegend=False,
            hoverinfo="skip",
            xaxis="x2",
            name="",  # Empty name to prevent legend entry
        )
    )

    # Explicitly hide any trace with xaxis='x2' or empty name from the legend
    for trace in fig.data:
        if hasattr(trace, "xaxis") and str(trace.xaxis) == "x2":
            trace.update(showlegend=False)
        if hasattr(trace, "name") and (trace.name == "" or trace.name is None):
            trace.update(showlegend=False)

    return fig


def get_birdseye(plot_data, selection):
    """
    Produce the plotly figure showing the day/slot matrix intersection for a selection of stars

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        availability (array): a 2D array of N_slots by N_nights, binary 1/0, it is the intersection of is_allocated and is_night
        all_stars (array): a array of StarPlotter objects

    Returns:
        fig (plotly figure): a plotly figure showing the day/slot matrix intersection for a selection of stars
    """

    all_stars = plot_data.views(selection)
    availablity = plot_data.nulltime
    semester_planner = plot_data.semester_planner
    fig = go.Figure()
    # fig.update_layout(width=1200, height=800, plot_bgcolor=clear, paper_bgcolor=clear)
    fig.update_layout(plot_bgcolor=clear, paper_bgcolor=clear)

    # when multiple StarPlotter obects are submitted or a programmatic StarPlotter object,
    # show the grayed out slots from the intersection of is_allocated and is_night
    if len(all_stars) > 1 or all_stars[0].allow_mapview == False:
        fig.add_trace(
            go.Heatmap(
                z=availablity,
                colorscale=[[0, "rgba(0,0,0,0)"], [1, gray]],
                zmin=0,
                zmax=1,
                opacity=1.0,
                showscale=False,
                name="Not On Sky",
                showlegend=False,
            )
        )
    # when just one StarPlotter object is submitted, show the overlay of all maps
    else:
        colors = sns.color_palette("deep", len(all_stars[0].maps_names) + 1)
        rgb_strings = [
            f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})" for r, g, b in colors
        ]
        for m in range(len(all_stars[0].maps_names)):
            # Skip the is_observable_now map
            if all_stars[0].maps_names[m] == "is_observable_now":
                continue
            map_name = all_stars[0].maps_names[m]
            z_data = (
                1 - all_stars[0].maps[map_name].astype(int).T
            )  # Invert all other maps

            fig.add_trace(
                go.Heatmap(
                    z=z_data,
                    colorscale=[[0, "rgba(0,0,0,0)"], [1, gray]],
                    zmin=0,
                    zmax=1,
                    opacity=1.0,
                    showscale=False,
                    name=all_stars[0].maps_names[m],
                    showlegend=True,
                )
            )

    for i in range(len(all_stars)):
        fig.add_trace(
            go.Heatmap(
                z=all_stars[i].starmap,
                colorscale=[[0, "rgba(0,0,0,0)"], [1, all_stars[i].star_color_rgb]],
                zmin=0,
                zmax=1,
                opacity=1.0,
                showscale=False,
                name=all_stars[i].target,
                hovertemplate="<b>"
                + str(all_stars[i].target)
                + "</b><br><b>Date: %{x}</b><br><b>Slot: %{y}</b><br>Forecasted N_Obs: "
                + str(all_stars[i].total_observations_requested)
                + "<extra></extra>",
                showlegend=True,
            )
        )

        if all_stars[i].draw_lines:
            # Add connecting line for points with value 1
            points = np.argwhere(all_stars[i].starmap == 1)
            sorted_indices = np.argsort(points[:, 1])  # sort by x (column index)
            x_coords = points[sorted_indices, 1]
            y_coords = points[sorted_indices, 0]
            fig.add_trace(
                go.Scatter(
                    x=x_coords,
                    y=y_coords,
                    mode="lines+markers",
                    line=dict(color=all_stars[i].star_color_rgb, width=2),
                    marker=dict(size=6, color=all_stars[i].starcolor_rgb),
                    name="Connected Points",
                )
            )

    add_grid_lines = (
        False  # this takes a long time to plot. Might not be necessary/worth it.
    )
    if add_grid_lines:
        # Add vertical grid lines every slot (x)
        for x in np.arange(0.5, all_stars[i].starmap.shape[1], 1):
            fig.add_shape(
                type="line",
                x0=x,
                x1=x,
                y0=0,
                y1=all_stars[i].starmap.shape[0] - 1,
                line=dict(color="lightgray", width=1),
                layer="below",
            )

    # Add vertical dashed line denoting "today"
    today = plot_data.semester_planner.access_obj.current_night_index
    fig.add_vrect(
        x0=today
        - 1,  # The minus one is just for aesthetic purposes.
        x1=today - 1,
        annotation_text="Today",
        line_dash="dash",
        fillcolor=None,
        line_width=2,
        line_color="black",
        annotation_position="bottom left",
    )
    # X-axis: ticks every 23 days, plus the last day
    x_tick_step = 23
    x_tickvals = list(range(0, plot_data.semester_planner.semester_length, x_tick_step))
    if (plot_data.semester_planner.semester_length - 1) not in x_tickvals:
        x_tickvals.append(plot_data.semester_planner.semester_length - 1)
    x_ticktext = [str(val + 1) for val in x_tickvals]

    # Create calendar date labels for secondary x-axis (top axis)
    # Format dates as "Jan<br>15" or "Aug<br>12" (month and day on separate lines)
    x_ticktext_dates = []
    for day_idx in x_tickvals:
        if day_idx < len(plot_data.semester_planner.access_obj.all_dates_array):
            date_str = plot_data.semester_planner.access_obj.all_dates_array[day_idx]
            # Parse date and format as "Jan<br>15" or "Aug<br>12" using HTML break tag
            date_obj = datetime.strptime(date_str, "%Y-%m-%d")
            month = date_obj.strftime("%b")
            day = date_obj.strftime("%d")
            x_ticktext_dates.append(f"{month}<br>{day}")
        else:
            x_ticktext_dates.append("")

    # Y-axis: ticks every 2 hours, using slot_size
    n_slots = int(24 * 60 // plot_data.semester_planner.config.getint("semester", "slot_size"))
    slots_per_2hr = int(2 * 60 // plot_data.semester_planner.config.getint("semester", "slot_size"))
    y_tickvals = list(range(0, n_slots, slots_per_2hr))
    y_ticktext = []
    for slot in y_tickvals:
        total_minutes = slot * plot_data.semester_planner.config.getint("semester", "slot_size")
        hours = total_minutes // 60
        minutes = total_minutes % 60
        y_ticktext.append(f"{hours:02.0f}:{minutes:02.0f}")

    # Calculate legend height based on number of traces
    num_traces = len(all_stars) + (
        1
        if len(all_stars) > 1 or all_stars[0].allow_mapview == False
        else len([m for m in all_stars[0].maps_names if m != "is_observable_now"])
    )
    legend_height = min(
        300, max(150, num_traces * 25)
    )  # Between 150-300px, 25px per trace

    # Add an invisible trace to force the secondary x-axis to appear
    # This trace must be associated with xaxis='x2' to make the secondary axis visible
    n_slots = int(24 * 60 // plot_data.semester_planner.config.getint("semester", "slot_size"))
    fig.add_trace(
        go.Scatter(
            x=[0, len(plot_data.semester_planner.access_obj.all_dates_array) - 1],
            y=[n_slots + 1, n_slots + 1],  # Position just above visible area
            mode="markers",
            marker=dict(size=0.01, opacity=0),
            showlegend=False,
            legendgroup=None,
            hoverinfo="skip",
            xaxis="x2",
            name="",  # Empty name to prevent legend entry
        )
    )

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
            tickmode="array",
            showgrid=False,
            anchor="y",
            side="bottom",
            range=[0, plot_data.semester_planner.semester_length - 1],  # Explicitly set range
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            tickvals=y_tickvals,
            ticktext=y_ticktext,
            tickmode="array",
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
            font=dict(size=labelsize - 18),
            bgcolor="rgba(255,255,255,0.7)",
            bordercolor="black",
            borderwidth=1,
            # Standardize legend size
            itemsizing="constant",  # All legend items same size
            itemwidth=30,  # Fixed width for legend items
            # Make legend more compact
            groupclick="toggleitem",  # Click group to toggle all items
            # Standardize legend dimensions
            tracegroupgap=5,  # Gap between trace groups
            traceorder="normal",  # Keep order as traces were added
        ),
        xaxis2=dict(
            title="",
            tickvals=x_tickvals,
            ticktext=x_ticktext_dates,
            tickmode="array",
            showgrid=False,
            side="top",
            overlaying="x",
            tickfont=dict(size=labelsize - 6),
            showticklabels=True,
            range=[
                0,
                plot_data.semester_planner.semester_length - 1,
            ],  # Match primary x-axis range
        ),
        margin=dict(
            b=200, t=100
        ),  # Bottom margin for legend below, top margin for date labels
    )
    return fig


def get_tau_inter_line(plot_data, selection, use_program_colors=False):
    """
    Produce a plotly figure showing requested vs on sky inter-night cadences, grouped by star name.

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        all_stars (array): a array of StarPlotter objects
        use_program_colors (bool): If True, use program_color_rgb; if False, use star_color_rgb (default: False)

    Returns:
        fig (plotly figure): a plotly figure showing requested vs on sky inter-night cadences, grouped by star name.
    """

    all_stars = plot_data.views(selection)
    semester_planner = plot_data.semester_planner
    request_tau_inter = []
    onsky_tau_inter = []
    targets = []
    programs = []
    colors = []
    for starobj in all_stars:
        onsky_diffs = list(np.diff(np.where(np.diff(starobj.cume_observe) > 0)[0]))
        onsky_tau_inter.extend(onsky_diffs)
        request_tau_inter.extend([starobj.tau_inter] * len(onsky_diffs))
        targets.extend([starobj.target] * len(onsky_diffs))
        programs.extend([starobj.program] * len(onsky_diffs))
        # Choose color based on flag
        if use_program_colors:
            colors.extend([starobj.program_color_rgb] * len(onsky_diffs))
        else:
            colors.extend([starobj.star_color_rgb] * len(onsky_diffs))

    all_request_tau_inters = np.array(request_tau_inter)
    all_onsky_tau_inters = np.array(onsky_tau_inter)
    all_targets = np.array(targets)
    all_programs = np.array(programs)
    all_colors = np.array(colors)

    fig = go.Figure()

    # Build map from program to point indices
    program_to_indices = {}
    for i, prog in enumerate(all_programs):
        program_to_indices.setdefault(prog, []).append(i)

    # Create one trace per star (grouped by target)
    maxyvals = []
    # Build map from target to point indices
    target_to_indices = {}
    for i, target in enumerate(all_targets):
        target_to_indices.setdefault(target, []).append(i)

    for target, indices in target_to_indices.items():
        idx_array = np.array(indices)
        x_vals = all_request_tau_inters[idx_array]
        y_vals = all_onsky_tau_inters[idx_array]
        text_vals = [f"{all_targets[i]} in {all_programs[i]}" for i in indices]
        color_vals = all_colors[idx_array].tolist()  # Convert to list for Plotly
        maxyvals.append(np.max(y_vals))
        fig.add_trace(
            go.Scatter(
                x=x_vals,
                y=y_vals,
                mode="markers",
                name=target,  # Use target for legend
                marker=dict(size=10, color=color_vals),
                text=text_vals,
                hovertemplate="%{text}<br>X: %{x}<br>Y: %{y}<extra></extra>",
            )
        )

    # Add 1-to-1 line
    min_val = 0
    if maxyvals == []:
        max_val = 0
    else:
        max_val = max(maxyvals)
    fig.add_trace(
        go.Scatter(
            x=[min_val, max_val],
            y=[min_val, max_val],
            mode="lines",
            line=dict(color="black", dash="dash"),
            name="1-to-1 line",
            showlegend=True,
        )
    )

    fig.update_layout(
        width=1400,
        height=800,
        xaxis_title="Requested Minimum Inter-Night Cadence",
        yaxis_title="On Sky Inter-Night Cadence",
        template="plotly_white",
        xaxis=dict(
            type="log",
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            showgrid=True,
            gridcolor="lightgray",
            gridwidth=0.5,
            tickmode="array",
            tickvals=[1, 10, 100],
            ticktext=["1", "10", "100"],
            range=[
                np.log10(0.5),
                np.log10(180),
            ],  # Set range from 0.5 to 180 in log scale
        ),
        yaxis=dict(
            type="log",
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            showgrid=True,
            gridcolor="lightgray",
            gridwidth=0.5,
            tickmode="array",
            tickvals=[1, 10, 100],
            ticktext=["1", "10", "100"],
            range=[
                np.log10(0.5),
                np.log10(180),
            ],  # Set range from 0.5 to 180 in log scale
        ),
    )
    return fig


def get_rawobs(plot_data, selection, use_program_colors=False):
    """
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
    """

    all_stars = plot_data.views(selection)
    semester_planner = plot_data.semester_planner
    fig = go.Figure()
    fig.update_layout(plot_bgcolor=clear, paper_bgcolor=clear)

    # Prepare data for each star
    targets = []
    total_requested = []
    past_obs = []
    future_obs = []
    total_completed = []  # past + scheduled
    pct_complete = []
    star_colors = []

    for star in all_stars:
        targets.append(star.target)
        total = star.total_observations_requested

        # Sum past observations
        past_total = (
            sum(star.observations_past.values()) if star.observations_past else 0
        )

        # Sum future observations
        future_total = (
            sum(star.observations_future.values()) if star.observations_future else 0
        )

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
        fig.add_trace(
            go.Scatter(
                x=[total_requested[i]],
                y=[total_completed[i]],
                mode="markers",
                marker=dict(
                    size=10,
                    color=star_colors[i],  # Use each star's individual color
                    opacity=0.7,
                ),
                name=targets[i],  # Target for legend (allows toggling)
                text=[targets[i]],  # Target for hover
                hovertemplate="<b>%{text}</b><br>"
                + "Total Requested: %{x}<br>"
                + "Past: %{customdata[0]}<br>"
                + "Scheduled: %{customdata[1]}<br>"
                + "Total (Past + Scheduled): %{y}<br>"
                + "% Complete: %{customdata[2]:.1f}%<extra></extra>",
                customdata=[[past_obs[i], future_obs[i], pct_complete[i]]],
            )
        )

    # Add diagonal lines for reference (y = x for 100% complete, y = 0.5x for 50% complete)
    # For log scale, we need to use log values
    min_val = min(
        min(total_requested) if total_requested else 1,
        min(total_completed) if total_completed else 1,
    )
    max_val = max(
        max(total_requested) if total_requested else 1,
        max(total_completed) if total_completed else 1,
    )
    # Ensure min_val is at least 1 for log scale
    if min_val < 1:
        min_val = 1

    # Add 100% complete reference line (y = x) - solid black line
    fig.add_trace(
        go.Scatter(
            x=[min_val, max_val],
            y=[min_val, max_val],
            mode="lines",
            line=dict(color="black", width=1, dash="solid"),
            name="100% Complete",
            showlegend=False,  # Hide reference line from legend
            hovertemplate="100% Complete Reference Line<extra></extra>",
        )
    )

    # Add 50% complete reference line (y = 0.5x)
    fig.add_trace(
        go.Scatter(
            x=[min_val, max_val],
            y=[min_val * 0.5, max_val * 0.5],
            mode="lines",
            line=dict(color="gray", width=1, dash="dash"),
            name="50% Complete",
            showlegend=False,  # Hide reference line from legend
            hovertemplate="50% Complete Reference Line<extra></extra>",
        )
    )

    # Add annotation at the top explaining the reference lines
    fig.add_annotation(
        x=0.5,  # Center horizontally
        y=1.02,  # Just above the plot
        xref="paper",
        yref="paper",
        text="solid = 1:1<br>dashed = 1:2",
        showarrow=False,
        font=dict(size=labelsize - 8, color="black"),
        align="center",
    )

    fig.update_layout(
        width=1400,
        height=800,
        xaxis_title="Total Requested Observations",
        yaxis_title="Total Observations (Past + Scheduled)",
        template="plotly_white",
        showlegend=True,  # Show legend so stars can be toggled on/off
        xaxis=dict(
            type="log",  # Log scale for x-axis
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            showgrid=True,
            gridcolor="lightgray",
            minor=dict(
                showgrid=False,  # Hide minor grid lines
                ticks="",  # Hide minor tick marks
            ),
            dtick=1,  # Major ticks at powers of 10
        ),
        yaxis=dict(
            type="log",  # Log scale for y-axis
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            showgrid=True,
            gridcolor="lightgray",
            minor=dict(
                showgrid=False,  # Hide minor grid lines
                ticks="",  # Hide minor tick marks
            ),
            dtick=1,  # Major ticks at powers of 10
        ),
        margin=dict(b=100, t=50),
    )

    return fig


def get_timebar(
    plot_data,
    selection,
    use_program_colors=False,
    prevent_negative=False,
):
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
    all_stars = plot_data.views(selection)
    semester_planner = plot_data.semester_planner
    programmatics = plot_data.semester_planner.programs

    # Charged hours are the splan slot-based single source of truth.
    ps = plot_data.semester_planner.timeline

    total_requested_hours = 0
    programs_used = []
    for starobj in all_stars:
        total_requested_hours += starobj.total_requested_hours
        programs_used.append(starobj.program)

    uids = [str(s.unique_id) for s in all_stars]
    total_past_hours, total_future_hours = _charged_hours_from_ps(plot_data.semester_planner, ps, unique_ids=uids
    )
    total_incomplete_hours = (
        total_requested_hours - total_past_hours - total_future_hours
    )

    programs_used_unique = sorted(set(programs_used))
    program_rows = programmatics.loc[
        programmatics.index.isin(programs_used_unique)
    ]
    total_allocated_hours = program_rows["hours"].sum()
    total_allocated_nights = total_allocated_hours / hours_per_night
    max_schedulable_hours = (
        program_rows["hours"] * program_rows["max_fillfactor"]
    ).sum()

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
        "<b>Incomplete Time</b><br>(requested - past - future)<br>If you have incomplete time, <br>some of your requests are infeasible <br> consider changing them, <br> i.e. cadence or redistributing",
        "<b>Future Scheduled Time</b>",
        "<b>Past Completed Time</b>",
        "<b>Requested Time</b>",
    ]
    sum_hours = (
        total_past_hours + total_future_hours + total_incomplete_hours + unused_hours
    )
    values = [
        unused_hours,
        total_incomplete_hours,
        total_future_hours,
        total_past_hours,
        total_requested_hours,
    ]
    colors = [
        "#FF0000",
        "#F18F01",
        "#A23B72",
        "#2E86AB",
        "#00FF00",
    ]  # Red, Orange, Purple, Blue, Green

    # Create the horizontal bar chart
    # Calculate percentages based on total allocated hours for all bars
    text_labels = []
    for i, (label, val) in enumerate(zip(labels, values)):
        # Calculate percentage relative to total allocated hours
        pct = (val / total_allocated_hours * 100) if total_allocated_hours > 0 else 0
        text_labels.append(f"{val:.1f} hrs ({pct:.1f}%)")

    fig = go.Figure(
        data=[
            go.Bar(
                x=values,
                y=labels,
                orientation="h",
                marker=dict(color=colors),
                text=text_labels,
                textposition="auto",
                hovertemplate="<b>%{y}</b><br>%{x:.2f} hours<br><extra></extra>",
            )
        ]
    )

    # Adjust margin if there's a warning to display
    top_margin = 180 if total_requested_hours > total_allocated_hours else 130

    fig.update_layout(
        title_text=f"<b>Total Requested:</b> {total_requested_hours:.1f} hours ≈ {total_requested_hours / hours_per_night:.1f} nights<br><b>Total Allocated:</b> {total_allocated_hours:.1f} hours ≈ {total_allocated_nights:.1f} nights ----> w/ losses = {total_allocated_nights * 0.75:.1f} nights <br>Requested and allocated time are measured in hours ({hours_per_night:.0f} hours per night for night equivalents).<br>Past and future bars use splan charged hours (slot-based).",
        template="plotly_white",
        showlegend=False,
        height=710,  # Increased height for more vertical spacing between labels
        width=1400,
        margin=dict(t=top_margin, b=50, l=200, r=50),
        bargap=0.2,
        xaxis=dict(title="Hours", titlefont=dict(size=14), tickfont=dict(size=12)),
        yaxis=dict(title="", titlefont=dict(size=14), tickfont=dict(size=11)),
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
        yref="y",
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
        yref="y",
    )

    # Add gray vertical dashed line at max schedulable time (hours * max_fillfactor)
    fig.add_shape(
        type="line",
        x0=max_schedulable_hours,
        x1=max_schedulable_hours,
        y0=-0.5,
        y1=len(labels) - 0.5,
        line=dict(color="gray", width=2, dash="dash"),
        xref="x",
        yref="y",
    )

    # Add invisible scatter trace for hover text on the allocated time line
    # Use the same categorical labels as the bar chart to avoid numeric y-axis ticks
    fig.add_trace(
        go.Scatter(
            x=[total_allocated_hours] * len(labels),
            y=labels,  # Use categorical labels instead of numeric positions
            mode="markers",
            marker=dict(size=20, opacity=0),  # Invisible but hoverable markers
            hovertemplate=f"<b>Allocated Time</b><br>{total_allocated_hours:.2f} hours<br>This line represents the total allocated time for your program<extra></extra>",
            hoverlabel=dict(bgcolor="black", font_color="white"),
            showlegend=False,
        )
    )

    # Add invisible scatter trace for hover text on the weather loss factor line
    weather_loss_value = (
        total_allocated_hours - total_allocated_hours * weather_loss_factor
    )
    fig.add_trace(
        go.Scatter(
            x=[weather_loss_value] * len(labels),
            y=labels,  # Use categorical labels instead of numeric positions
            mode="markers",
            marker=dict(size=20, opacity=0),  # Invisible but hoverable markers
            hovertemplate=f"<b>Weather Loss Factor</b><br>{weather_loss_value:.2f} hours<br>Allocated time minus {weather_loss_factor * 100:.0f}% weather loss<br>This is only a first order estimate based on historical losses.<extra></extra>",
            hoverlabel=dict(bgcolor="gray", font_color="white"),
            showlegend=False,
        )
    )

    # Add invisible scatter trace for hover text on the max schedulable line
    fig.add_trace(
        go.Scatter(
            x=[max_schedulable_hours] * len(labels),
            y=labels,  # Use categorical labels instead of numeric positions
            mode="markers",
            marker=dict(size=20, opacity=0),  # Invisible but hoverable markers
            hovertemplate=(
                f"<b>Maximum Schedulable Time</b><br>{max_schedulable_hours:.2f} hours<br>"
                "Sum of awarded hours times per-program max_fillfactor "
                "(default 1.25). Algorithmically, you are forbidden from "
                "getting more time than this.<extra></extra>"
            ),
            hoverlabel=dict(bgcolor="gray", font_color="white"),
            showlegend=False,
        )
    )

    # Add warning annotation if requested time exceeds allocated time
    if total_requested_hours > total_allocated_hours * 1.1:
        fig.add_annotation(
            text="<b>You have requested more time than you are allocated.</b>",
            xref="paper",
            yref="paper",
            x=0.5,
            y=1.35,
            showarrow=False,
            font=dict(size=18, color="red"),
            xanchor="center",
            yanchor="middle",
        )

    return fig


def get_timebar_by_program(plot_data, selection=None, prevent_negative=False):
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
    semester_planner = plot_data.semester_planner
    programs_dict = plot_data.program_dict
    ledger = programs_ledger_for_plot(plot_data.semester_planner)

    all_programs_in_csv = set(ledger.index)
    programs_with_requests = set(programs_dict.keys())

    # Find programs in CSV that don't have any requests
    programs_without_requests = all_programs_in_csv - programs_with_requests

    # Combine all programs: those with requests and those without
    all_program_codes = sorted(
        list(programs_with_requests) + list(programs_without_requests)
    )

    # Store data for each program
    program_data = {}
    max_x_value = 0  # Track maximum x value for consistent scaling

    # Process programs with requests
    for program_code in sorted(programs_with_requests):
        program_stars = programs_dict[program_code]

        total_requested_hours = sum(
            starobj.total_requested_hours for starobj in program_stars
        )
        if program_code in ledger.index:
            row = ledger.loc[program_code]
            total_past_hours = float(row["past_hours"])
            total_future_hours = float(row["sched_hours"])
            total_allocated_hours = float(row["hours"])
        else:
            total_past_hours = 0.0
            total_future_hours = 0.0
            total_allocated_hours = 0.0
        total_incomplete_hours = (
            total_requested_hours - total_past_hours - total_future_hours
        )

        # Calculate unused hours
        unused_hours = total_allocated_hours - total_future_hours - total_past_hours

        # Apply negative value prevention if enabled
        if prevent_negative:
            total_incomplete_hours = max(0, total_incomplete_hours)
            unused_hours = max(0, unused_hours)

        program_data[program_code] = {
            "unused": unused_hours,
            "incomplete": total_incomplete_hours,
            "future": total_future_hours,
            "past": total_past_hours,
            "requested": total_requested_hours,
            "allocated": total_allocated_hours,
        }

        # Update max value for scaling
        max_x_value = max(
            max_x_value,
            total_requested_hours,
            total_allocated_hours,
            unused_hours,
            total_incomplete_hours,
            total_future_hours,
            total_past_hours,
        )

    # Process programs without requests (all bars = 0, but show allocated time)
    for program_code in sorted(programs_without_requests):
        if program_code in ledger.index:
            total_allocated_hours = float(ledger.loc[program_code, "hours"])
        else:
            total_allocated_hours = 0.0

        # All values are zero for programs with no requests
        program_data[program_code] = {
            "unused": total_allocated_hours,  # All allocated time is unused
            "incomplete": 0,
            "future": 0,
            "past": 0,
            "requested": 0,
            "allocated": total_allocated_hours,
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
        vertical_spacing=0.12,
    )

    # Colors in display order: Red, Orange, Purple, Blue, Green
    display_colors = ["#FF0000", "#F18F01", "#A23B72", "#2E86AB", "#00FF00"]
    category_names = [
        "Unused",
        "Incomplete",
        "Future Scheduled",
        "Past Completed",
        "Requested",
    ]

    # Add bars for each program in its own subplot
    for idx, program_code in enumerate(all_program_codes):
        data = program_data[program_code]

        # Calculate row and column position (1-indexed)
        row = (idx // num_cols) + 1
        col = (idx % num_cols) + 1

        # Prepare bar data for this program
        program_values = [
            data["unused"],
            data["incomplete"],
            data["future"],
            data["past"],
            data["requested"],
        ]

        # Add bars to this subplot
        fig.add_trace(
            go.Bar(
                x=program_values,
                y=category_names,
                orientation="h",
                marker=dict(color=display_colors),
                text=[f"{v:.1f}" if v > 0 else "" for v in program_values],
                textposition="auto",
                hovertemplate=f"<b>{program_code}</b><br>%{{y}}<br>%{{x:.2f}} hours<extra></extra>",
                showlegend=False,
            ),
            row=row,
            col=col,
        )

        # Add vertical dashed line for allocated time
        allocated = data["allocated"]
        # For subplots, determine the correct axis reference
        # In make_subplots, axes are numbered: x, x2, x3, ... and y, y2, y3, ...
        if idx == 0:
            xref, yref = "x", "y"
        else:
            xref, yref = f"x{idx + 1}", f"y{idx + 1}"

        fig.add_shape(
            type="line",
            x0=allocated,
            x1=allocated,
            y0=-0.5,
            y1=4.5,
            line=dict(color="black", width=2, dash="dash"),
            xref=xref,
            yref=yref,
        )

        # Add gray vertical dashed line for weather loss estimate
        weather_loss_factor = 0.2
        fig.add_shape(
            type="line",
            x0=allocated - allocated * weather_loss_factor,
            x1=allocated - allocated * weather_loss_factor,
            y0=-0.5,
            y1=4.5,
            line=dict(color="gray", width=2, dash="dash"),
            xref=xref,
            yref=yref,
        )

        # Add gray vertical dashed line at allocated * max_fillfactor
        if program_code in ledger.index:
            max_ff = float(ledger.loc[program_code, "max_fillfactor"])
        else:
            max_ff = 1.25
        max_schedulable = allocated * max_ff
        fig.add_shape(
            type="line",
            x0=max_schedulable,
            x1=max_schedulable,
            y0=-0.5,
            y1=4.5,
            line=dict(color="gray", width=2, dash="dash"),
            xref=xref,
            yref=yref,
        )

        # Add invisible scatter for hover on allocated line
        fig.add_trace(
            go.Scatter(
                x=[allocated],
                y=[category_names[2]],  # Middle bar (Future Scheduled)
                mode="markers",
                marker=dict(size=15, opacity=0),
                hovertemplate=f"<b>{program_code} Allocated Time</b><br>{allocated:.2f} hours<br>Total allocated time for this program<extra></extra>",
                hoverlabel=dict(bgcolor="black", font_color="white"),
                showlegend=False,
            ),
            row=row,
            col=col,
        )

        # Add invisible scatter for hover on weather loss line
        weather_loss_value = allocated - allocated * weather_loss_factor
        fig.add_trace(
            go.Scatter(
                x=[weather_loss_value],
                y=[category_names[2]],  # Middle bar (Future Scheduled)
                mode="markers",
                marker=dict(size=15, opacity=0),
                hovertemplate=f"<b>{program_code} Weather Loss Factor</b><br>{weather_loss_value:.2f} hours<br>Allocated time minus {weather_loss_factor * 100:.0f}% weather loss<extra></extra>",
                hoverlabel=dict(bgcolor="gray", font_color="white"),
                showlegend=False,
            ),
            row=row,
            col=col,
        )

        # Add invisible scatter for hover on max schedulable line
        fig.add_trace(
            go.Scatter(
                x=[max_schedulable],
                y=[category_names[2]],  # Middle bar (Future Scheduled)
                mode="markers",
                marker=dict(size=15, opacity=0),
                hovertemplate=(
                    f"<b>{program_code} Maximum Schedulable</b><br>"
                    f"{max_schedulable:.2f} hours<br>"
                    f"Awarded hours times max_fillfactor ({max_ff:.2f})"
                    "<extra></extra>"
                ),
                hoverlabel=dict(bgcolor="gray", font_color="white"),
                showlegend=False,
            ),
            row=row,
            col=col,
        )

        # Update x-axis for this subplot (scaled to this program's data)
        # Include max schedulable and weather loss so the gray lines are visible
        weather_loss_value = allocated - allocated * weather_loss_factor
        program_max = max(
            data["unused"],
            data["incomplete"],
            data["future"],
            data["past"],
            data["requested"],
            data["allocated"],
            max_schedulable,
            weather_loss_value,
        )
        program_max = max(program_max, 1.0)  # Ensure at least 1.0 to avoid empty scale

        fig.update_xaxes(title="Hours", range=[0, program_max * 1.1], row=row, col=col)

        # Update y-axis for this subplot (no labels)
        fig.update_yaxes(title="", showticklabels=False, row=row, col=col)

    # Update overall layout
    fig.update_layout(
        title_text="<b>Time Breakdown by Program</b><br>Each program shows 5 bars (top to bottom): Requested (green), Past Completed (blue), Future Scheduled (purple), Incomplete (orange), Unused (red)<br>Dashed vertical line represents total allocated time. Note each grid is on its own scaling.",
        template="plotly_white",
        showlegend=False,
        height=max(600, num_rows * 250),
        width=1400,
        margin=dict(t=150, b=50, l=50, r=50),
    )

    return fig


def get_football(plot_data, selection, use_program_colors=False):
    """
    Mollweide sky map: per-target scatter points layered on top of a
    semester-wide observability heatmap (nights per (ra, dec) for which at
    least one slot is dark, above the horizon, and outside the moon-avoidance
    zone). The heatmap is computed on a coarse RA/Dec grid via a fresh
    `Access` instance and cached per semester to disk.

    Parameters:
        semester_planner: the semester planner object
        all_stars (list): array of StarPlotter objects
        use_program_colors (bool): If True, use program_color_rgb; if False, use star_color_rgb (default: False)

    Returns:
        fig (plotly figure): the assembled Mollweide sky map.
    """

    all_stars = plot_data.views(selection)
    semester_planner = plot_data.semester_planner
    star_ras = [s.ra for s in all_stars]
    star_decs = [s.dec for s in all_stars]
    targets = [s.target for s in all_stars]
    programs = [s.program for s in all_stars]
    if use_program_colors:
        colors = [s.program_color_rgb for s in all_stars]
    else:
        colors = [s.star_color_rgb for s in all_stars]
    program_frame = pd.DataFrame(
        {
            "target": targets,
            "program_code": programs,
            "color": colors,
            "ra": star_ras,
            "dec": star_decs,
        }
    )

    # Equal-area-in-dec sky grid (uniform in sin(dec)).
    n_ra = 90
    n_dec = 90
    grid_ra_axis = np.linspace(0, 360, n_ra)
    grid_dec_axis = np.degrees(np.arcsin(np.linspace(-1, 1, n_dec)))
    RA_grid, DEC_grid = np.meshgrid(grid_ra_axis, grid_dec_axis)

    n_points = n_dec * n_ra
    grid_frame = pd.DataFrame(
        {
            "target": [f"noname_{i}" for i in range(n_points)],
            "ra": RA_grid.flatten(),
            "dec": DEC_grid.flatten(),
        }
    )

    semester = (
        plot_data.semester_planner.config.get("global", "semester_start_day")[:4] + plot_data.semester_planner.config.get("global", "semester")[-1]
    )
    cache_dir = _football_cache_dir(plot_data.semester_planner)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_grids_file = str(cache_dir / f"{semester}_sky_grids.npz")
    cache_image_file = str(cache_dir / f"{semester}_sky_availability_image.txt")
    semester_length = plot_data.semester_planner.semester_length

    if os.path.exists(cache_grids_file):
        cached_data = np.load(cache_grids_file)
        RA_grid = cached_data["RA_grid"]
        DEC_grid = cached_data["DEC_grid"]
        NIGHTS_grid = cached_data["NIGHTS_grid"]
    else:
        # Seasonality for the sky-grid points. A fresh Access scoped to the
        # grid keeps the planner's real access_obj pristine. Bare defaults
        # give all-True cubes for future / custom / inter / allocated /
        # clear; only altaz, moon, and night actually gate.
        seasonality_frame = grid_frame[["target", "ra", "dec"]].copy()
        seasonality_frame["unique_id"] = seasonality_frame["target"]
        grid_access = ac.Access(
            queue=plot_data.semester_planner.queue,
            request_frame=seasonality_frame,
            semester_start_date=plot_data.semester_planner.config.get("global", "semester_start_day"),
            semester_length=semester_length,
            slot_size=plot_data.semester_planner.config.getint("semester", "slot_size"),
        )
        record = grid_access.build_access()
        is_observable_now = np.logical_and.reduce(
            [record.is_altaz, record.is_moon, grid_access.compute_night()]
        )
        grid_frame["nights_observable"] = (
            is_observable_now.any(axis=2).sum(axis=1).astype(int)
        )

        NIGHTS_grid = griddata(
            points=(grid_frame.ra, grid_frame.dec),
            values=grid_frame.nights_observable,
            xi=(RA_grid, DEC_grid),
            method="linear",
        )

        np.savez(
            cache_grids_file,
            RA_grid=RA_grid,
            DEC_grid=DEC_grid,
            NIGHTS_grid=NIGHTS_grid,
        )

    if os.path.exists(cache_image_file):
        with open(cache_image_file, "r") as f:
            img_base64 = f.read()
    else:
        RA_shifted = np.radians(RA_grid - 180)
        DEC_rad = np.radians(DEC_grid)

        fig_mpl, ax = plt.subplots(
            subplot_kw={"projection": "mollweide"}, figsize=(10, 5)
        )
        im = ax.pcolormesh(
            RA_shifted,
            DEC_rad,
            NIGHTS_grid,
            cmap="gray",
            shading="nearest",
            vmin=70,
            vmax=semester_length,
        )
        ax.axis("off")

        buf = BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", pad_inches=0, dpi=150)
        plt.close()
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode()

        with open(cache_image_file, "w") as f:
            f.write(img_base64)

    fig = go.Figure()

    fig.add_layout_image(
        dict(
            source=f"data:image/png;base64,{img_base64}",
            xref="paper",
            yref="paper",
            x=0,
            y=1,
            sizex=1,
            sizey=1,
            xanchor="left",
            yanchor="top",
            sizing="stretch",
            layer="below",
            opacity=1,
        )
    )

    # Invisible Contour trace whose sole purpose is to carry the colorbar.
    # plotly has no first-class colorbar-only object; opacity=0 keeps the
    # contour itself hidden while still rendering the legend strip.
    fig.add_trace(
        go.Contour(
            z=NIGHTS_grid,
            x=RA_grid[0] - 180,
            y=DEC_grid[:, 0],
            showscale=True,
            colorscale="gray",
            contours=dict(start=70, end=semester_length, size=10),
            opacity=0,
            colorbar=dict(
                title="Observable<br>Nights",
                titleside="top",
                x=-0.15,
                len=0.75,
                thickness=15,
            ),
        )
    )

    if not program_frame.empty:
        marker = "star"
        size = 20 if len(all_stars) == 1 else 10
        grouped = program_frame.groupby("program_code")
        for program, group in grouped:
            group.reset_index(inplace=True, drop=True)
            hover = [f"{name} in {program}" for name in group["target"]]
            color = group["color"].tolist()

            fig.add_trace(
                go.Scattergeo(
                    lon=group["ra"] - 180,
                    lat=group["dec"],
                    mode="markers",
                    name=program,
                    marker=dict(symbol=marker, size=size, color=color, opacity=1),
                    text=hover,
                    hovertemplate="%{text}<br>RA: %{lon:.2f}°, Dec: %{lat:.2f}°<extra></extra>",
                )
            )

    fig.update_layout(
        shapes=[
            dict(
                type="circle",
                xref="paper",
                yref="paper",
                x0=0.0,
                y0=0.0,
                x1=1.0,
                y1=1.0,
                line=dict(color="black", width=2),
            )
        ]
    )

    # Step 5: Layout
    fig.update_layout(
        geo=dict(
            projection_type="mollweide",
            showland=False,
            showcoastlines=False,
            showframe=False,
            bgcolor="rgba(0,0,0,0)",
            lonaxis=dict(showgrid=False),
            lataxis=dict(showgrid=False),
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        template="none",
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
                font=dict(size=14),
            ),
            dict(
                text="Dec (deg)",  # Y-axis label
                x=-0.07,
                y=0.5,
                xref="paper",
                yref="paper",
                showarrow=False,
                textangle=-90,
                font=dict(size=14),
            ),
        ],
    )
    return fig


def _completion_splan_weight(req_row):
    """Numeric ``splan_weight`` for completion plots (matches semester optimizer)."""
    if req_row is None:
        return np.nan
    for col in ("splan_weight", "weight"):
        if col in req_row.index and pd.notna(req_row.get(col)):
            return float(pd.to_numeric(req_row[col], errors="coerce"))
    return np.nan


def _splan_weight_legend_label(weight):
    if pd.isna(weight):
        return "splan_weight: (missing)"
    w = float(weight)
    if w == int(w):
        return f"splan_weight: {int(w)}"
    return f"splan_weight: {w}"


def _completion_by_request_frame(plot_data, selection):
    """Per-request semester completion % joined with request.csv ``splan_weight``."""
    all_stars = plot_data.views(selection)
    semester_planner = plot_data.semester_planner
    req = plot_data.semester_planner.requests.set_index("unique_id")

    rows = []
    for star in all_stars:
        pct = (
            float(star.cume_observe_pct[-1])
            if len(star.cume_observe_pct) > 0
            else 0.0
        )
        req_row = req.loc[star.unique_id] if star.unique_id in req.index else None
        rows.append(
            {
                "unique_id": star.unique_id,
                "target": star.target,
                "program": star.program,
                "completion_pct": pct,
                "splan_weight": _completion_splan_weight(req_row),
            }
        )
    return pd.DataFrame(rows)


def _weight_legend_label(weight):
    if pd.isna(weight):
        return "weight: (missing)"
    w = float(weight)
    if w == int(w):
        return f"weight: {int(w)}"
    return f"weight: {w}"


def _sorted_weight_values(weights):
    def sort_key(w):
        if pd.isna(w):
            return (2, 0.0)
        try:
            return (0, float(w))
        except (TypeError, ValueError):
            return (1, str(w))

    return sorted(weights, key=sort_key)


def _weight_series_matches(series, weight):
    """Match weight values with numeric coercion (fixes int 1 vs float 1.0)."""
    numeric = pd.to_numeric(series, errors="coerce")
    if pd.isna(weight):
        return numeric.isna()
    return numeric == float(weight)


_COMPLETION_BIN_LABELS = [f"[{lo}, {lo + 10})" for lo in range(0, 100, 10)] + ["100"]


def _completion_bin_label(pct):
    if pct >= 100:
        return "100"
    lo = int(pct // 10) * 10
    return f"[{lo}, {lo + 10})"


def get_completion_histogram_by_weight(plot_data, selection):
    """
    Histogram of request completion rate (%), one curve per ``splan_weight``.
    """
    all_stars = plot_data.views(selection)
    semester_planner = plot_data.semester_planner
    df = _completion_by_request_frame(plot_data, selection)
    fig = go.Figure()
    weight_values = _sorted_weight_values(df["splan_weight"].unique())
    colors = sns.color_palette("deep", max(len(weight_values), 1))
    rgb_strings = [
        f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})"
        for r, g, b in colors
    ]

    for i, weight in enumerate(weight_values):
        subset = df.loc[
            _weight_series_matches(df["splan_weight"], weight), "completion_pct"
        ]
        bin_labels = subset.map(_completion_bin_label)
        counts = (
            bin_labels.value_counts()
            .reindex(_COMPLETION_BIN_LABELS, fill_value=0)
            .astype(int)
        )
        fig.add_trace(
            go.Bar(
                x=_COMPLETION_BIN_LABELS,
                y=counts,
                name=_splan_weight_legend_label(weight),
                marker_color=rgb_strings[i % len(rgb_strings)],
            )
        )

    fig.update_layout(
        width=1400,
        height=600,
        title="Completion Rate by splan_weight",
        xaxis_title="Completion Rate (%)",
        yaxis_title="Number of Requests",
        barmode="stack",
        plot_bgcolor=clear,
        paper_bgcolor=clear,
        xaxis=dict(categoryorder="array", categoryarray=_COMPLETION_BIN_LABELS),
        legend=dict(title="splan_weight"),
    )
    return fig


def get_completion_vs_target_name(plot_data, selection):
    """
    Scatter of completion rate (%) vs target name, sorted alphabetically by target.
    """
    all_stars = plot_data.views(selection)
    semester_planner = plot_data.semester_planner
    df = _completion_by_request_frame(plot_data, selection)
    df = df.sort_values("target", kind="mergesort").reset_index(drop=True)
    target_order = df["target"].tolist()

    program_colors = {star.program: star.program_color_rgb for star in all_stars}
    fig = go.Figure()
    for program in sorted(df["program"].unique()):
        sub = df[df["program"] == program]
        fig.add_trace(
            go.Scatter(
                x=sub["target"],
                y=sub["completion_pct"],
                mode="markers",
                name=program,
                marker=dict(size=8, color=program_colors.get(program, "steelblue")),
                customdata=np.stack(
                    [
                        sub["program"].to_numpy(),
                        sub["splan_weight"].map(_splan_weight_legend_label).to_numpy(),
                    ],
                    axis=-1,
                ),
                hovertemplate=(
                    "<b>%{x}</b><br>Program: %{customdata[0]}<br>"
                    "Completion: %{y:.1f}%<br>%{customdata[1]}<extra></extra>"
                ),
            )
        )

    fig.update_layout(
        width=1400,
        height=700,
        title="Completion Rate by Target",
        xaxis_title="Target",
        yaxis_title="Completion Rate (%)",
        plot_bgcolor=clear,
        paper_bgcolor=clear,
        xaxis=dict(
            categoryorder="array",
            categoryarray=target_order,
            tickangle=-45,
        ),
        yaxis=dict(range=[0, 100]),
        margin=dict(b=150),
        showlegend=True,
    )
    return fig


def get_request_frame(plot_data, selection):
    """
    Get a filtered request frame containing only the stars in all_stars.

    Args:
        semester_planner: the semester planner object
        all_stars (list): array of StarPlotter objects

    Returns:
        filtered_frame (pd.DataFrame): filtered request frame with only the specified stars
    """
    all_stars = plot_data.views(selection)
    semester_planner = plot_data.semester_planner
    # Extract targets from the StarPlotter objects
    starids = [star.unique_id for star in all_stars]

    # Filter the request frame to only include the specified stars
    filtered_frame = plot_data.semester_planner.requests[
        plot_data.semester_planner.requests["unique_id"].isin(starids)
    ].copy()

    return filtered_frame
