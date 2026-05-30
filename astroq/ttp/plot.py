"""TTP-specific Plotly plotting utilities.

This module owns the plots that visualize a single-night TTP solution
(``astroq.ttp.model.TTPModel``):

- :func:`plot_path_2D_interactive` -- az/alt vs time for the chosen tour.
- :func:`get_slew_animation_plotly` -- animated polar plot of the slew.
- :func:`createTelSlewPath` -- helper that resamples the schedule onto
  animation frames.

Lifted verbatim from ``astroq.plot`` during the Stage 1 file reorg;
Stage 4 retargets the reads to the new ``TTPModel`` attributes
(``model.night_start``, ``model.observer``, ``model.wrap_limit``,
``model.requests_frame``, ``model.inaccessible_zones``) and collapses the
three hardcoded obstruction sections into one loop over
``inaccessible_zones``.
"""

# Standard library imports
import numpy as np
import pandas as pd

# Third-party imports
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
import astropy.units as u
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def _as_model(data):
    """Accept ``TTPModel`` or legacy ``[TTPModel]`` wrapper."""
    return data[0] if isinstance(data, (list, tuple)) else data


def schedule_to_ladder_frame(model):
    """Build a ladder-plot frame from ``model.schedule`` (scheduled + extras)."""
    sched = model.schedule
    on_sky = sched[~sched["is_anchor"]]
    scheduled = on_sky[on_sky["scheduled"]].sort_values("order")
    extras = on_sky[~on_sky["scheduled"]].sort_values("t_early")

    def _pack(df, *, scheduled_rows):
        starname = df.get("starname", df["unique_id"])
        return pd.DataFrame(
            {
                "Starname": df["unique_id"],
                "human_starname": starname,
                "First Available": df["t_early"],
                "Last Available": df["t_late"],
                "Start Exposure": df["t_start"] if scheduled_rows else 0.0,
                "Stop Exposure": df["t_end"] if scheduled_rows else df["t_visit"],
                "Total Exp Time (min)": df["t_visit"],
                "Exposure Time (min)": df["exptime"],
                "N_shots": df["n_exp"],
                "Priority": df["priority"],
                "Slew to Next (min)": df["t_slew"].fillna(0.0),
                "Minutes the from Start of the Night": (
                    (df["t_start"] + df["t_end"]) / 2 if scheduled_rows else 0.0
                ),
            }
        )

    parts = []
    if len(scheduled):
        parts.append(_pack(scheduled, scheduled_rows=True))
    if len(extras):
        parts.append(_pack(extras, scheduled_rows=False))
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def createTelSlewPath(stamps, changes, pointings, animationStep=120):
    """
    Correctly assign each frame of the animation to the telescope pointing at that time

    stamps (list of zeros) - the list where each element represents a frame of the animation. We manipulate and return this at the end.
    changes (list) - the times at which the telescope pointing changes (in order of the slew path)
    poitings (list) - the astropy target objects of for the stars to be observed, in order of the slew path
    animationStep (int) - the time, in seconds, between frames

    return
        stamps - now a list where element holds the pointing of the telescope (aka the star object) at that frame

    """
    minPerStep = int(animationStep / 60)
    mins = int(60 / minPerStep)

    changes = (changes - changes[0]) * 24 * mins
    for c in range(len(changes)):
        changes[c] = int(changes[c])

    for i in range(len(changes) - 1):
        for j in range(len(stamps)):
            if j >= changes[i] and j < changes[i + 1]:
                stamps[j] = pointings[i]

    if len(stamps) > 0:
        k = 0
        while k < len(stamps) and stamps[k] == 0:
            stamps[k] = pointings[0]
            k += 1
        l = len(stamps) - 1
        while l >= 0 and stamps[l] == 0:
            stamps[l] = pointings[-1]
            l -= 1

    return stamps


def _inaccessible_zone_traces(inaccessible_zones):
    """Build per-zone ``Scatterpolar`` traces from ``model.inaccessible_zones``.

    Each zone is a ``(az_min, az_max, alt_min, alt_max)`` rectangle in degrees;
    we draw it as a closed polar polygon at zenith-distance r=90-alt. Only the
    first trace carries ``showlegend=True`` so the legend is not cluttered with
    one entry per zone.
    """
    traces = []
    for idx, (az_min, az_max, alt_min, alt_max) in enumerate(inaccessible_zones or []):
        theta = np.linspace(az_min, az_max, 100)
        r_inner = np.full(100, 90 - alt_max)  # nearer to zenith
        r_outer = np.full(100, 90 - alt_min)  # nearer to horizon
        traces.append(
            go.Scatterpolar(
                r=np.concatenate([r_inner, r_outer[::-1], [r_inner[0]]]),
                theta=np.concatenate([theta, theta[::-1], [theta[0]]]),
                fill="toself",
                fillcolor="rgba(255, 0, 0, 0.7)",
                line=dict(color="rgba(255, 0, 0, 0)"),
                showlegend=(idx == 0),
                name="Excluded zone",
                hoverinfo="skip",
            )
        )
    return traces


def get_slew_animation_plotly(
    data, request_selected_path, animationStep=120, inaccessible_zones=None
):
    """Create a Plotly animated polar plot showing telescope slew path during observations.

    Args:
        data: ``TTPModel`` or ``[TTPModel]`` solution.
        request_selected_path: Path to request_selected.csv (used only to map
            ``unique_id`` -> human-readable ``starname`` for the hover text).
        animationStep (int): the time, in seconds, between animation frames. Default 120s.
        inaccessible_zones: optional list of obstruction boxes from ``Queue``.

    Returns:
        fig (plotly figure): an interactive animated figure with play/pause controls
    """

    model = _as_model(data)

    request_selected_df = pd.read_csv(request_selected_path)

    t = np.arange(
        model.night_start.jd,
        model.night_end.jd,
        TimeDelta(animationStep, format="sec").jd,
    )
    t = Time(t, format="jd")

    on_sky = model.schedule[~model.schedule["is_anchor"]]
    scheduled = on_sky[on_sky["scheduled"]].sort_values("order")
    list_targets = SkyCoord(
        scheduled.ra.values * u.deg,
        scheduled.dec.values * u.deg,
        frame="icrs",
    )
    names = scheduled["unique_id"].tolist()

    AZ = model.observer.altaz(t, list_targets, grid_times_targets=True)
    alt = np.round(AZ.az.rad, 2)
    az = 90 - np.round(AZ.alt.deg, 2)

    schedule_times = model.night_start.jd + scheduled["t_start"].to_numpy() / (24 * 60)

    stamps = [0] * len(t)
    slewPath = createTelSlewPath(stamps, schedule_times, list_targets)
    AZ1 = model.observer.altaz(t, slewPath, grid_times_targets=False)
    tel_az = np.round(AZ1.az.rad, 2)
    tel_zen = 90 - np.round(AZ1.alt.deg, 2)

    names_array = np.array(names)

    unique_id_to_starname = dict(
        zip(
            request_selected_df["unique_id"].astype(str),
            request_selected_df["starname"],
        )
    )
    human_starname_array = np.array(
        [unique_id_to_starname.get(str(uid), str(uid)) for uid in names_array]
    )

    zone_traces = _inaccessible_zone_traces(inaccessible_zones)
    n_zones = len(zone_traces)

    frames = []
    for i in range(len(t)):
        is_observed = schedule_times <= float(t[i].jd)

        # Per-frame: rebuild zone traces so the (first-frame-only) legend flag
        # is on for frame 0 and off for subsequent frames.
        if i == 0:
            zones_this_frame = zone_traces
        else:
            zones_this_frame = []
            for ztrace in zone_traces:
                ztrace_copy = go.Scatterpolar(ztrace.to_plotly_json())
                ztrace_copy.update(showlegend=False)
                zones_this_frame.append(ztrace_copy)

        frame_data = list(zones_this_frame) + [
            go.Scatterpolar(
                r=az[:, i][is_observed],
                theta=np.degrees(alt[:, i][is_observed]),
                mode="markers",
                marker=dict(size=10, color="orange", symbol="star"),
                name="Observed",
                showlegend=(i == 0),
                text=human_starname_array[is_observed],
                hovertemplate="<b>%{text}</b><br>Az: %{theta:.1f}°<br>ZD: %{r:.1f}°<extra></extra>",
            ),
            go.Scatterpolar(
                r=az[:, i][~is_observed],
                theta=np.degrees(alt[:, i][~is_observed]),
                mode="markers",
                marker=dict(size=10, color="white", symbol="star"),
                name="Scheduled",
                showlegend=(i == 0),
                text=human_starname_array[~is_observed],
                hovertemplate="<b>%{text}</b><br>Az: %{theta:.1f}°<br>ZD: %{r:.1f}°<extra></extra>",
            ),
            go.Scatterpolar(
                r=tel_zen[: i + 1] if i > 0 else tel_zen[:1],
                theta=np.degrees(tel_az[: i + 1] if i > 0 else tel_az[:1]),
                mode="lines",
                line=dict(color="orange", width=2),
                name="Telescope Path",
                showlegend=(i == 0),
            ),
        ]

        frames.append(go.Frame(data=frame_data, name=str(i)))

    fig = go.Figure(data=frames[0].data if frames else [], frames=frames)

    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                range=[0, 90],
                showticklabels=False,
                ticks="",
                showline=False,
                gridcolor="rgba(255, 255, 255, 0.2)",
                gridwidth=1,
            ),
            angularaxis=dict(
                direction="counterclockwise",
                rotation=90,
                gridcolor="rgba(255, 255, 255, 0.2)",
                gridwidth=1,
                tickfont=dict(size=18, color="black"),
                showticklabels=True,
            ),
            bgcolor="black",
        ),
        annotations=[
            dict(
                text="<b>N</b>",
                x=0.495,
                y=1.1,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=22, color="black"),
            ),
            dict(
                text="<b>W</b>",
                x=1.0,
                y=0.5,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=22, color="black"),
            ),
            dict(
                text="<b>S</b>",
                x=0.495,
                y=-0.1,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=22, color="black"),
            ),
            dict(
                text="<b>E</b>",
                x=-0.0,
                y=0.5,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=22, color="black"),
            ),
        ],
        transition={"duration": 0},
        updatemenus=[
            {
                "type": "buttons",
                "showactive": False,
                "direction": "left",
                "x": 0.35,
                "y": -0.2,
                "xanchor": "left",
                "yanchor": "bottom",
                "buttons": [
                    {
                        "label": "  \u25b6 Play  ",
                        "method": "animate",
                        "args": [
                            None,
                            {
                                "frame": {"duration": 100, "redraw": True},
                                "fromcurrent": True,
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                    {
                        "label": "  \u23f8 Pause  ",
                        "method": "animate",
                        "args": [
                            [None],
                            {
                                "frame": {"duration": 0, "redraw": False},
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                ],
                "bgcolor": "white",
                "bordercolor": "black",
                "borderwidth": 2,
                "font": {"size": 16, "color": "black", "family": "Arial"},
            }
        ],
        sliders=[
            {
                "active": 0,
                "yanchor": "top",
                "y": -0.15,
                "xanchor": "left",
                "currentvalue": {
                    "prefix": "Time: ",
                    "visible": True,
                    "xanchor": "right",
                    "font": {"size": 14, "color": "black"},
                },
                "pad": {"b": 10, "t": 50},
                "len": 0.9,
                "x": 0.1,
                "font": {"size": 12, "color": "black"},
                "steps": [
                    {
                        "args": [
                            [f.name],
                            {
                                "frame": {"duration": 100, "redraw": True},
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                        "label": t[k].datetime.strftime("%H:%M"),
                        "method": "animate",
                    }
                    for k, f in enumerate(frames)
                ],
                "transition": {"duration": 100},
            }
        ],
        width=800,
        height=800,
        title=dict(text="Telescope Slew Animation", font=dict(color="black", size=20)),
        template="plotly_white",
        paper_bgcolor="white",
        plot_bgcolor="white",
        font=dict(color="black"),
        hovermode="closest",
    )

    return fig


def plot_path_2D_interactive(data, night_start_time=None):
    """Create an interactive Plotly plot showing telescope azimuth and altitude paths with UTC times and white background.

    Args:
        data: ``TTPModel`` or ``[TTPModel]`` solution
        night_start_time: Astropy Time object representing the start of night (Minute 0) from allocation file

    Returns:
        fig (plotly figure): an interactive plot showing telescope azimuth and altitude paths with UTC times and white background.
    """

    model = _as_model(data)
    wrap = model.wrap_limit

    if night_start_time is None:
        night_start_time = model.night_start
    night_start_jd = night_start_time.jd

    on_sky = model.schedule[~model.schedule["is_anchor"]]
    scheduled = on_sky[on_sky["scheduled"]].sort_values("order")
    if scheduled.empty:
        fig = make_subplots(
            rows=2,
            cols=1,
            shared_xaxes=True,
            subplot_titles=("Azimuth Path", "Elevation Path"),
            vertical_spacing=0.1,
        )
        fig.update_layout(height=600, width=1000, template="plotly_white")
        return fig

    starname = scheduled.get("starname", scheduled["unique_id"])
    t_start = scheduled["t_start"].to_numpy()
    t_end = scheduled["t_end"].to_numpy()
    t_start_time = model.night_start + TimeDelta(t_start * 60, format="sec")
    t_end_time = model.night_start + TimeDelta(t_end * 60, format="sec")
    coords = SkyCoord(
        scheduled.ra.values * u.deg,
        scheduled.dec.values * u.deg,
        frame="icrs",
    )

    aa_start = model.observer.altaz(t_start_time, coords)
    aa_end = model.observer.altaz(t_end_time, coords)
    az_start = np.atleast_1d(aa_start.az.deg)
    alt_start = np.atleast_1d(aa_start.alt.deg)
    az_end = np.atleast_1d(aa_end.az.deg)
    alt_end = np.atleast_1d(aa_end.alt.deg)

    obs_time = np.empty(2 * len(scheduled))
    az_path = np.empty(2 * len(scheduled))
    alt_path = np.empty(2 * len(scheduled))
    names = []
    for i in range(len(scheduled)):
        obs_time[2 * i] = t_start_time[i].jd
        obs_time[2 * i + 1] = t_end_time[i].jd
        az_path[2 * i], az_path[2 * i + 1] = az_start[i], az_end[i]
        alt_path[2 * i], alt_path[2 * i + 1] = alt_start[i], alt_end[i]
        names.extend([starname.iloc[i], starname.iloc[i]])

    if len(obs_time) == 2 * len(names):
        expanded_names = []
        for name in names:
            expanded_names.append(name)
            expanded_names.append(name)
        names = expanded_names
    elif len(obs_time) != len(names):
        names = names * (len(obs_time) // len(names) + 1)
        names = names[: len(obs_time)]

    min_len = min(len(obs_time), len(az_path), len(alt_path), len(names))
    obs_time = obs_time[:min_len]
    az_path = np.array(az_path[:min_len])
    alt_path = np.array(alt_path[:min_len])
    names = names[:min_len]

    az_path = np.mod(az_path, 360)
    az_path_original = az_path.copy()

    # Values above 270° displayed as negative (subtract 360) so e.g. 350° → -10°.
    az_path_display = az_path.copy()
    az_path_display[az_path_display > 270] -= 360

    time_labels = [Time(t, format="jd").isot[11:16] for t in obs_time]

    hover_text_az = [
        f"Time: {time_labels[i]}<br>Target: {names[i]}<br>Az: {az_path_original[i]:.1f}°"
        for i in range(len(obs_time))
    ]
    hover_text_alt = [
        f"Time: {time_labels[i]}<br>Target: {names[i]}<br>Alt: {alt_path[i]:.1f}°"
        for i in range(len(obs_time))
    ]

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        subplot_titles=("Azimuth Path", "Elevation Path"),
        vertical_spacing=0.1,
    )

    fig.add_trace(
        go.Scatter(
            x=obs_time,
            y=az_path_display,
            mode="lines+markers",
            marker=dict(color="indigo"),
            name="Azimuth",
            text=hover_text_az,
            hovertemplate="%{text}<extra></extra>",
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=obs_time,
            y=alt_path,
            mode="lines+markers",
            marker=dict(color="seagreen"),
            name="Elevation",
            text=hover_text_alt,
            hovertemplate="%{text}<extra></extra>",
        ),
        row=2,
        col=1,
    )

    if wrap is not None:
        wrap_normalized = wrap % 360
        wrap_display = wrap_normalized
        if wrap_display > 270:
            wrap_display -= 360

        fig.add_shape(
            type="line",
            x0=obs_time[0],
            x1=obs_time[-1],
            y0=wrap_display,
            y1=wrap_display,
            line=dict(color="red", dash="dash", width=2),
            row=1,
            col=1,
        )
        fig.add_annotation(
            x=obs_time[-1],
            y=wrap_display,
            text=f"Wrap = {wrap_normalized}\u00b0",
            showarrow=False,
            font=dict(color="red", size=10),
            row=1,
            col=1,
        )

    # Shade Start-Exposure → Stop-Exposure intervals (minutes-from-night-start).
    if len(scheduled):
        start_exposures = t_start
        stop_exposures = t_end

        for i, (start_min, stop_min) in enumerate(zip(start_exposures, stop_exposures)):
            start_jd = night_start_jd + (start_min / 1440.0)
            stop_jd = night_start_jd + (stop_min / 1440.0)

            fig.add_vrect(
                x0=start_jd,
                x1=stop_jd,
                fillcolor="yellow",
                opacity=0.3,
                layer="below",
                line_width=0,
                row=1,
                col=1,
            )
            fig.add_vrect(
                x0=start_jd,
                x1=stop_jd,
                fillcolor="yellow",
                opacity=0.3,
                layer="below",
                line_width=0,
                row=2,
                col=1,
            )

    time_span = obs_time[-1] - obs_time[0]
    if time_span < 0.1:  # < ~2.4 h
        interval_hours = 0.5
    elif time_span < 0.3:  # < ~7 h
        interval_hours = 1.0
    else:
        interval_hours = 2.0

    interval_jd = interval_hours / 24

    start_time = obs_time[0]
    end_time = obs_time[-1]
    num_ticks = int((end_time - start_time) / interval_jd) + 2
    tick_positions = np.linspace(start_time, end_time, num_ticks)

    tick_labels = [Time(t, format="jd").isot[11:16] for t in tick_positions]

    fig.update_xaxes(
        tickmode="array",
        tickvals=tick_positions,
        ticktext=tick_labels,
        title_text="Time (UTC)",
        row=2,
        col=1,
    )

    # Range -95° to 275° with 5° buffer; values >270° displayed as negative.
    az_y_min = -95
    az_y_max = 275

    tick_interval = 45
    az_tick_positions = np.arange(-90, 271, tick_interval)

    az_tick_labels = []
    for pos in az_tick_positions:
        if pos < 0:
            label = int(pos + 360)
        else:
            label = int(pos)
        az_tick_labels.append(f"{label}\u00b0")

    fig.update_yaxes(
        tickmode="array",
        tickvals=az_tick_positions,
        ticktext=az_tick_labels,
        range=[az_y_min, az_y_max],
        title_text="Azimuth (deg)",
        row=1,
        col=1,
    )

    fig.update_yaxes(range=[0, 90], title_text="Altitude (deg)", row=2, col=1)

    fig.update_layout(height=600, width=1000, template="plotly_white")
    return fig
