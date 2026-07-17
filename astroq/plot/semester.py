"""Semester-scope Plotly figures."""

from astroq.plot._common import (  # noqa: F401
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
    _charged_hours_from_ps,
    _football_cache_dir,
    programs_ledger_for_plot,
    _render_datatable,
    _TEMPLATE_ENV,
    _Path,
)
from astroq.plot.context import MAP_NAMES, PlotData, PlotSelection  # noqa: F401

def get_cof(plot_data, selection, use_time=False):
    """Cumulative Observability Function (COF) for a selection of requests/programs.

    Args:
        plot_data (PlotData): shared plot context from build_plot_data.
        selection (PlotSelection): which requests/programs to include.
        use_time (bool): normalize by awarded program hours instead of visit count.

    Returns:
        plotly.graph_objects.Figure: the COF figure.
    """
    sel = plot_data.select(selection)
    table = sel.table
    fig = go.Figure()
    fig.update_layout(plot_bgcolor=gray, paper_bgcolor=clear)

    dates = plot_data.all_dates_array
    n_nights = plot_data.n_nights
    night_indices = np.arange(n_nights)
    burn_line = np.round(np.linspace(0, 100, n_nights), 2)

    fig.add_shape(
        type="line",
        x0=night_indices[0],
        y0=burn_line[0],
        x1=night_indices[-1],
        y1=burn_line[-1],
        line=dict(color="black", width=2, dash="dash"),
        layer="below",
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            line=dict(color="black", width=2, dash="dash"),
            name="Even Burn Rate",
            showlegend=True,
            hoverinfo="skip",
        )
    )

    total_color = table["program_color"].iloc[0]
    slots_per_hour = 60 / plot_data.slot_size
    prog_hours = plot_data.programs.set_index("program")["hours"]

    if not use_time:
        total_denom = int(table["requested_visits"].sum())
        total_cume = sel.cume_visits.sum(axis=1).to_numpy(dtype=float)
        if total_denom > 0:
            total_pct = np.round(total_cume / total_denom * 100, 2)
        else:
            total_past = int(total_cume[-1]) if len(total_cume) else 0
            total_pct = (
                np.round(total_cume / total_past * 100, 2)
                if total_past > 0
                else np.zeros(n_nights)
            )
        fig.add_trace(
            go.Scatter(
                x=night_indices,
                y=total_pct,
                mode="lines",
                line=dict(color=total_color, width=2),
                name="Total",
                hovertemplate="Night: %{x}"
                + "<br>Date: "
                + "%{customdata}"
                + "<br>% Complete: %{y}"
                + "<br># Visits Requested: "
                + str(int(total_denom))
                + "<br>",
                customdata=dates,
            )
        )
    else:
        programs_in = set(table["program_code"])
        total_program_hours = prog_hours[prog_hours.index.isin(programs_in)].sum()
        cume_hours = sel.cume_slots.sum(axis=1).to_numpy(dtype=float) / slots_per_hour
        if total_program_hours > 0:
            total_pct = np.round(cume_hours / total_program_hours * 100, 2)
        else:
            total_pct = np.zeros(n_nights)
        if len(programs_in) == 1:
            total_trace_label = "<b>" + list(programs_in)[0] + "</b> (Total)<br>"
        else:
            total_trace_label = "<b>All programs (Total)</b><br>"
        fig.add_trace(
            go.Scatter(
                x=night_indices,
                y=total_pct,
                mode="lines",
                line=dict(color=total_color, width=2),
                name="Total",
                hovertemplate=total_trace_label
                + "Night: %{x}"
                + "<br>Date: "
                + "%{customdata}"
                + "<br>Time charged (% of awarded hours): %{y}"
                + "<br>Total program time: "
                + f"{total_program_hours:.1f} hours<br>"
                + "<extra></extra>",
                customdata=dates,
            )
        )

    def add_request_trace(uid):
        row = table.loc[uid]
        if use_time:
            tp = prog_hours[prog_hours.index == row["program_code"]].sum()
            cume_h = sel.cume_slots[uid].to_numpy(dtype=float) / slots_per_hour
            y_vals = np.round(cume_h / tp * 100, 2) if tp > 0 else np.zeros(n_nights)
            hovertemplate = (
                "<b>"
                + str(row["program_code"])
                + "</b><br>Night: %{x}"
                + "<br>Date: "
                + "%{customdata}"
                + "<br>Time charged (% of awarded hours): %{y}<br>Total program time: "
                + f"{tp:.1f} hours<br>"
                + "<extra></extra>"
            )
        else:
            y_vals = sel.cume_visits_pct[uid].to_numpy()
            hovertemplate = (
                "Night: %{x}"
                + "<br>Date: "
                + "%{customdata}"
                + "<br>% Complete: %{y}"
                + "<br># Visits Requested: "
                + str(int(row["requested_visits"]))
                + "<br>"
            )
        fig.add_trace(
            go.Scatter(
                x=night_indices,
                y=y_vals,
                mode="lines",
                line=dict(color=row["star_color"], width=2),
                name=row["target"],
                hovertemplate=hovertemplate,
                customdata=dates,
            )
        )
        return None

    pd.Series(sel.ids).apply(add_request_trace)

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
    """Day/slot forecast matrix for a selection of requests/programs.

    Args:
        plot_data (PlotData): shared plot context from build_plot_data.
        selection (PlotSelection): which requests/programs to include.

    Returns:
        plotly.graph_objects.Figure: the birdseye figure.
    """
    sel = plot_data.select(selection)
    table = sel.table
    availablity = plot_data.nulltime
    fig = go.Figure()
    fig.update_layout(plot_bgcolor=clear, paper_bgcolor=clear)

    # Multiple requests or a program aggregate: show the grayed-out unavailable
    # slots. A single request: overlay its per-map availability cubes.
    single_map = len(sel.ids) == 1 and bool(table["allow_mapview"].iloc[0])
    if not single_map:
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
    else:
        maps = plot_data.maps_for(sel.ids[0])
        for map_name in MAP_NAMES:
            if map_name == "is_observable_now":
                continue
            fig.add_trace(
                go.Heatmap(
                    z=1 - maps[map_name].astype(int).T,
                    colorscale=[[0, "rgba(0,0,0,0)"], [1, gray]],
                    zmin=0,
                    zmax=1,
                    opacity=1.0,
                    showscale=False,
                    name=map_name,
                    showlegend=True,
                )
            )

    def add_starmap(uid):
        row = table.loc[uid]
        fig.add_trace(
            go.Heatmap(
                z=plot_data.starmap_for(uid, is_program=sel.is_program),
                colorscale=[[0, "rgba(0,0,0,0)"], [1, row["star_color"]]],
                zmin=0,
                zmax=1,
                opacity=1.0,
                showscale=False,
                name=row["target"],
                hovertemplate="<b>"
                + str(row["target"])
                + "</b><br><b>Date: %{x}</b><br><b>Slot: %{y}</b><br>Forecasted N_Obs: "
                + str(row["total_observations_requested"])
                + "<extra></extra>",
                showlegend=True,
            )
        )
        return None

    pd.Series(sel.ids).apply(add_starmap)

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
    """Requested vs on-sky inter-night cadence, one marker series per target.

    Args:
        plot_data (PlotData): shared plot context from build_plot_data.
        selection (PlotSelection): which requests/programs to include.
        use_program_colors (bool): color by program instead of per-request.

    Returns:
        plotly.graph_objects.Figure: requested vs on-sky cadence scatter.
    """

    sel = plot_data.select(selection)
    table = sel.table
    color_col = "program_color" if use_program_colors else "star_color"

    meta = table[["target", "program_code", "tau_inter", color_col]].rename(
        columns={color_col: "color", "program_code": "program"}
    )
    order = {uid: i for i, uid in enumerate(sel.ids)}
    points = (
        plot_data.onsky_cadence[plot_data.onsky_cadence["unique_id"].isin(sel.ids)]
        .merge(meta, left_on="unique_id", right_index=True, how="inner")
        .assign(_order=lambda d: d["unique_id"].map(order))
        .sort_values("_order", kind="mergesort")
    )

    fig = go.Figure()

    def add_target_trace(group):
        fig.add_trace(
            go.Scatter(
                x=group["tau_inter"],
                y=group["onsky_tau_inter"],
                mode="markers",
                name=group.name,
                marker=dict(size=10, color=group["color"].tolist()),
                text=[
                    f"{t} in {p}"
                    for t, p in zip(group["target"], group["program"])
                ],
                hovertemplate="%{text}<br>X: %{x}<br>Y: %{y}<extra></extra>",
            )
        )
        return None

    if not points.empty:
        points.groupby("target", sort=False).apply(add_target_trace)

    # Add 1-to-1 line
    min_val = 0
    max_val = int(points["onsky_tau_inter"].max()) if not points.empty else 0
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
    """Scatter of total requested vs completed (past + scheduled) observations.

    Args:
        plot_data (PlotData): shared plot context from build_plot_data.
        selection (PlotSelection): which requests/programs to include.
        use_program_colors (bool): color by program instead of per-request.

    Returns:
        plotly.graph_objects.Figure: one marker per request.
    """

    sel = plot_data.select(selection)
    color_col = "program_color" if use_program_colors else "star_color"
    fig = go.Figure()
    fig.update_layout(plot_bgcolor=clear, paper_bgcolor=clear)

    t = sel.table.assign(
        total_completed=lambda d: d["past_visits"] + d["future_visits"]
    )
    t["pct_complete"] = np.where(
        t["total_observations_requested"] > 0,
        t["total_completed"] / t["total_observations_requested"] * 100,
        0,
    )

    def add_point(row):
        fig.add_trace(
            go.Scatter(
                x=[row["total_observations_requested"]],
                y=[row["total_completed"]],
                mode="markers",
                marker=dict(size=10, color=row[color_col], opacity=0.7),
                name=row["target"],
                text=[row["target"]],
                hovertemplate="<b>%{text}</b><br>"
                + "Total Requested: %{x}<br>"
                + "Past: %{customdata[0]}<br>"
                + "Scheduled: %{customdata[1]}<br>"
                + "Total (Past + Scheduled): %{y}<br>"
                + "% Complete: %{customdata[2]:.1f}%<extra></extra>",
                customdata=[
                    [row["past_visits"], row["future_visits"], row["pct_complete"]]
                ],
            )
        )
        return None

    t.apply(add_point, axis=1)

    # Add diagonal lines for reference (y = x for 100% complete, y = 0.5x for 50% complete)
    # For log scale, we need to use log values
    min_val = min(
        int(t["total_observations_requested"].min()) if len(t) else 1,
        int(t["total_completed"].min()) if len(t) else 1,
    )
    max_val = max(
        int(t["total_observations_requested"].max()) if len(t) else 1,
        int(t["total_completed"].max()) if len(t) else 1,
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
    """Horizontal bar chart of requested vs past vs scheduled vs allocated hours.

    Args:
        plot_data (PlotData): shared plot context from build_plot_data.
        selection (PlotSelection): which requests/programs to include.
        use_program_colors (bool): retained for API symmetry (bars use fixed colors).
        prevent_negative (bool): clamp Incomplete/Unused categories at zero.

    Returns:
        plotly.graph_objects.Figure: the time-budget bar chart.
    """
    sel = plot_data.select(selection)
    table = sel.table
    semester_planner = plot_data.semester_planner
    programmatics = plot_data.semester_planner.programs

    # Charged hours are the splan slot-based single source of truth.
    ps = plot_data.semester_planner.timeline

    total_requested_hours = float(table["total_requested_hours"].sum())

    if sel.is_program:
        total_past_hours, total_future_hours = _charged_hours_from_ps(
            semester_planner, ps, program_codes=sel.ids
        )
        programs_used_unique = sorted(set(sel.ids))
    else:
        total_past_hours, total_future_hours = _charged_hours_from_ps(
            semester_planner, ps, unique_ids=sel.ids
        )
        programs_used_unique = sorted(set(table["program_code"]))
    total_incomplete_hours = (
        total_requested_hours - total_past_hours - total_future_hours
    )

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
    A dashed vertical line represents their total allocated time. Programs are
    arranged in a 3-column grid, each subplot on its own x-scale.

    Args:
        plot_data (PlotData): shared plot context from build_plot_data.
        selection (PlotSelection): unused; every program is shown.
        prevent_negative (bool): clamp Incomplete/Unused categories at zero.

    Returns:
        plotly.graph_objects.Figure: grid of per-program time-budget bars.
    """
    semester_planner = plot_data.semester_planner
    ledger = programs_ledger_for_plot(plot_data.semester_planner)
    prog_table = plot_data.program_table

    all_program_codes = sorted(set(ledger.index) | set(prog_table.index))

    # Per-program hour breakdown. Requested comes from the aggregated request
    # table; past/scheduled/allocated from the ledger. Programs without requests
    # fall out with requested=0 and unused=allocated via the same arithmetic.
    df = pd.DataFrame(index=all_program_codes)
    df["requested"] = (
        prog_table["total_requested_hours"].reindex(all_program_codes).fillna(0.0)
    )
    df["past"] = ledger["past_hours"].reindex(all_program_codes).fillna(0.0)
    df["future"] = ledger["sched_hours"].reindex(all_program_codes).fillna(0.0)
    df["allocated"] = ledger["hours"].reindex(all_program_codes).fillna(0.0)
    df["incomplete"] = df["requested"] - df["past"] - df["future"]
    df["unused"] = df["allocated"] - df["future"] - df["past"]
    if prevent_negative:
        df["incomplete"] = df["incomplete"].clip(lower=0)
        df["unused"] = df["unused"].clip(lower=0)
    program_data = df.to_dict("index")

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

    Args:
        plot_data (PlotData): shared plot context from build_plot_data.
        selection (PlotSelection): which requests/programs to include.
        use_program_colors (bool): color by program instead of per-request.

    Returns:
        plotly.graph_objects.Figure: the assembled Mollweide sky map.
    """

    sel = plot_data.select(selection)
    table = sel.table
    color_col = "program_color" if use_program_colors else "star_color"
    program_frame = pd.DataFrame(
        {
            "target": table["target"].to_numpy(),
            "program_code": table["program_code"].to_numpy(),
            "color": table[color_col].to_numpy(),
            "ra": table["ra"].to_numpy(),
            "dec": table["dec"].to_numpy(),
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
        size = 20 if len(sel.ids) == 1 else 10

        def add_program_trace(group):
            program = group.name
            fig.add_trace(
                go.Scattergeo(
                    lon=group["ra"] - 180,
                    lat=group["dec"],
                    mode="markers",
                    name=program,
                    marker=dict(
                        symbol=marker,
                        size=size,
                        color=group["color"].tolist(),
                        opacity=1,
                    ),
                    text=[f"{name} in {program}" for name in group["target"]],
                    hovertemplate="%{text}<br>RA: %{lon:.2f}°, Dec: %{lat:.2f}°<extra></extra>",
                )
            )
            return None

        program_frame.groupby("program_code").apply(add_program_trace)

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
    t = plot_data.select(selection).table
    return pd.DataFrame(
        {
            "unique_id": list(t.index),
            "target": t["target"].to_numpy(),
            "program": t["program_code"].to_numpy(),
            "completion_pct": t["completion_pct"].to_numpy(),
            "splan_weight": t["splan_weight"].to_numpy(),
        }
    )


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
    df = _completion_by_request_frame(plot_data, selection)
    df = df.sort_values("target", kind="mergesort").reset_index(drop=True)
    target_order = df["target"].tolist()

    program_colors = plot_data.program_colors
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
    """Filtered request frame containing only the selected requests.

    Args:
        plot_data (PlotData): shared plot context from build_plot_data.
        selection (PlotSelection): which requests/programs to include.

    Returns:
        pd.DataFrame: request.csv rows for the selection (source order preserved).
    """
    sel = plot_data.select(selection)
    requests = plot_data.semester_planner.requests
    return requests[requests["unique_id"].astype(str).isin(sel.ids)].copy()
