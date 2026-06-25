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


def _encoder_az_display(az_deg, enc_min, enc_max):
    """Sky az -> encoder az in ``[enc_min, enc_max]`` (else nearest, for display).

    Mirrors ``astroq.queue.base.Queue._encoder_az`` but never returns NaN so a
    scheduled point always plots somewhere sensible.
    """
    az = np.asarray(az_deg, dtype=float)
    out = az.copy()
    for k in (-360.0, 0.0, 360.0):
        cand = az + k
        in_range = (cand >= enc_min) & (cand <= enc_max)
        out = np.where(in_range, cand, out)
    return out


def _as_model(data):
    """Accept ``TTPModel`` or legacy ``[TTPModel]`` wrapper."""
    return data[0] if isinstance(data, (list, tuple)) else data


def schedule_to_ladder_frame(model):
    """Build a ladder-plot frame from ``model.schedule`` (scheduled + extras)."""
    sched = model.schedule
    on_sky = sched[~sched["is_anchor"]]
    scheduled = on_sky[on_sky["scheduled"]].sort_values("order")
    extras = on_sky[~on_sky["scheduled"]].sort_values("t_earliest_start")

    def _pack(df, *, scheduled_rows):
        target = df.get("target", df["unique_id"])
        return pd.DataFrame(
            {
                "unique_id": df["unique_id"],
                "Target": target,
                "Earliest Start": df["t_earliest_start"],
                "Latest Finish": df["t_latest_finish"],
                "Start Exposure": df["t_start"] if scheduled_rows else 0.0,
                "Stop Exposure": df["t_end"] if scheduled_rows else df["t_visit"],
                "Visit Length (min)": df["t_visit"],
                "Exposure Time (min)": df["exptime"],
                "N_shots": df["n_exp"],
                "Weight": df["weight"],
                "Slew to Next (min)": df["t_slew"].fillna(0.0),
                "is_scheduled": df["scheduled"].astype(bool),
                "Scheduled (min. from start)": (
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


def _telescope_track(model, scheduled, sample_s=15.0):
    """Fine telescope trajectory ``(jd, sky_az_deg, zen_deg)`` for the tour.

    During each visit the telescope tracks the target (sidereal motion); between
    visits it slews, interpolated linearly in the telescope's *wrap frame* and
    sampled every ``sample_s`` seconds so the drawn line follows the actual
    cable-wrap route -- short south-wrap moves for the two-state model, long
    unwinds for the legacy single-cut model. Falls back to shortest-arc azimuth
    interpolation if no wrap information is available on ``model``.
    """
    n = len(scheduled)
    if n == 0:
        return np.array([]), np.array([]), np.array([])

    ns_jd = model.night_start.jd
    t_start = scheduled["t_start"].to_numpy(dtype=float)
    t_end = (scheduled["t_end"].to_numpy(dtype=float)
             if "t_end" in scheduled.columns else t_start.copy())
    t_end = np.where(np.isfinite(t_end) & (t_end > t_start), t_end, t_start)
    coords = SkyCoord(scheduled.ra.values * u.deg, scheduled.dec.values * u.deg, frame="icrs")

    wrap_states = getattr(model, "wrap_states", None)
    wrap_limit = getattr(model, "wrap_limit", None)
    has_state = (
        wrap_states is not None
        and "wrap_state" in scheduled.columns
        and scheduled["wrap_state"].notna().any()
    )
    st = scheduled["wrap_state"].to_numpy() if has_state else None

    def to_enc(az, i):
        """Sky az (deg) -> continuous encoder az for node ``i``'s wrap frame."""
        if has_state and np.isfinite(st[i]):
            lo, hi = wrap_states[int(st[i])][1], wrap_states[int(st[i])][2]
            return float(_encoder_az_display(np.array([az]), lo, hi)[0])
        if wrap_limit:
            return float(np.mod(az + (360.0 - wrap_limit), 360.0))
        return None

    def to_sky(enc, i):
        if has_state and np.isfinite(st[i]):
            return float(np.mod(enc, 360.0))
        if wrap_limit:
            return float(np.mod(enc - (360.0 - wrap_limit), 360.0))
        return float(np.mod(enc, 360.0))

    sample_jd = TimeDelta(sample_s, format="sec").jd
    seg_t, seg_az, seg_zen = [], [], []

    def track_visit(i):
        a = ns_jd + t_start[i] / (24 * 60)
        b = ns_jd + t_end[i] / (24 * 60)
        m = max(int((b - a) / sample_jd) + 1, 1) if b > a else 1
        tt = Time(np.linspace(a, b, m), format="jd")
        aa = model.observer.altaz(tt, coords[i])
        seg_t.append(np.atleast_1d(tt.jd))
        seg_az.append(np.atleast_1d(aa.az.deg))
        seg_zen.append(90.0 - np.atleast_1d(aa.alt.deg))

    track_visit(0)
    for i in range(1, n):
        a = ns_jd + t_end[i - 1] / (24 * 60)
        b = ns_jd + t_start[i] / (24 * 60)
        if b > a:
            aa0 = model.observer.altaz(Time(a, format="jd"), coords[i - 1])
            aa1 = model.observer.altaz(Time(b, format="jd"), coords[i])
            az0, alt0 = float(aa0.az.deg), float(aa0.alt.deg)
            az1, alt1 = float(aa1.az.deg), float(aa1.alt.deg)
            m = max(int((b - a) / sample_jd) + 1, 2)
            fr = np.linspace(0.0, 1.0, m)
            alt = alt0 + fr * (alt1 - alt0)
            e0, e1 = to_enc(az0, i - 1), to_enc(az1, i)
            if e0 is not None and e1 is not None:
                enc = e0 + fr * (e1 - e0)
                azs = np.array([to_sky(e, i) for e in enc])
            else:  # shortest-arc fallback
                daz = ((az1 - az0 + 180.0) % 360.0) - 180.0
                azs = np.mod(az0 + fr * daz, 360.0)
            seg_t.append(a + fr * (b - a))
            seg_az.append(azs)
            seg_zen.append(90.0 - alt)
        track_visit(i)

    T = np.concatenate(seg_t)
    A = np.concatenate(seg_az)
    Z = np.concatenate(seg_zen)
    o = np.argsort(T)
    return T[o], A[o], Z[o]


def get_slew_animation_plotly(
    data, request_selected_path, animationStep=120, inaccessible_zones=None,
    slew_sample_s=15.0,
):
    """Create a Plotly animated polar plot showing telescope slew path during observations.

    Args:
        data: ``TTPModel`` or ``[TTPModel]`` solution.
        request_selected_path: Path to request_selected.csv (used only to map
            ``unique_id`` -> human-readable ``target`` for the hover text).
        animationStep (int): the time, in seconds, between animation frames. Default 120s.
        inaccessible_zones: optional list of obstruction boxes from ``Queue``.
        slew_sample_s (float): cadence, in seconds, at which the telescope slew
            path is sampled so the drawn line traces the actual motion. Default 15s.

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

    # Actual telescope trajectory, finely sampled (wrap-aware) so the path line
    # traces the real slew motion instead of jumping between targets.
    track_jd, track_az, track_zen = _telescope_track(
        model, scheduled, sample_s=slew_sample_s
    )

    # Plot every attempted target (scheduled + considered-but-not-hit). Targets
    # that are never hit stay gray for the whole animation; hit targets turn
    # orange once their observation time passes.
    attempted = on_sky
    all_targets = SkyCoord(
        attempted.ra.values * u.deg,
        attempted.dec.values * u.deg,
        frame="icrs",
    )
    AZ = model.observer.altaz(t, all_targets, grid_times_targets=True)
    alt = np.round(AZ.az.rad, 2)
    az = 90 - np.round(AZ.alt.deg, 2)

    # Observation time per attempted target; inf (never observed) when unscheduled.
    obs_time = np.where(
        attempted["scheduled"].to_numpy(),
        model.night_start.jd + attempted["t_start"].to_numpy() / (24 * 60),
        np.inf,
    )

    names_array = np.array(attempted["unique_id"].tolist())

    unique_id_to_target = dict(
        zip(
            request_selected_df["unique_id"].astype(str),
            request_selected_df["target"],
        )
    )
    human_target_array = np.array(
        [unique_id_to_target.get(str(uid), str(uid)) for uid in names_array]
    )

    zone_traces = _inaccessible_zone_traces(inaccessible_zones)
    n_zones = len(zone_traces)

    frames = []
    for i in range(len(t)):
        is_observed = obs_time <= float(t[i].jd)

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
                r=az[:, i][~is_observed],
                theta=np.degrees(alt[:, i][~is_observed]),
                mode="markers",
                marker=dict(size=10, color="gray", symbol="star"),
                name="Attempted",
                showlegend=(i == 0),
                text=human_target_array[~is_observed],
                hovertemplate="<b>%{text}</b><br>Az: %{theta:.1f}°<br>ZD: %{r:.1f}°<extra></extra>",
            ),
            go.Scatterpolar(
                r=az[:, i][is_observed],
                theta=np.degrees(alt[:, i][is_observed]),
                mode="markers",
                marker=dict(size=10, color="orange", symbol="star"),
                name="Observed",
                showlegend=(i == 0),
                text=human_target_array[is_observed],
                hovertemplate="<b>%{text}</b><br>Az: %{theta:.1f}°<br>ZD: %{r:.1f}°<extra></extra>",
            ),
            go.Scatterpolar(
                r=track_zen[track_jd <= float(t[i].jd)],
                theta=track_az[track_jd <= float(t[i].jd)],
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


def save_slew_animation(fig, html_path, gif_path=None, **gif_kw):
    """Write interactive HTML and a GIF copy of a slew-animation figure."""
    fig.write_html(html_path)
    if gif_path is None:
        gif_path = html_path.rsplit(".", 1)[0] + ".gif"
    write_slew_animation_gif(fig, gif_path, **gif_kw)


def write_slew_animation_gif(
    fig,
    path,
    *,
    max_frames=90,
    fps=3,
    width=640,
    height=640,
):
    """Export a Plotly slew-animation figure to an animated GIF via Kaleido.

    Long nights produce hundreds of Plotly frames; this subsamples evenly to
    ``max_frames`` so GIF size stays reasonable.
    """
    import io

    from PIL import Image

    if not fig.frames:
        raise ValueError("figure has no animation frames")

    n = len(fig.frames)
    if n > max_frames:
        indices = np.unique(np.round(np.linspace(0, n - 1, max_frames)).astype(int))
    else:
        indices = np.arange(n)

    layout = fig.layout.to_plotly_json()
    layout.pop("updatemenus", None)
    layout.pop("sliders", None)

    images = []
    for i in indices:
        frame = fig.frames[int(i)]
        frame_fig = go.Figure(data=frame.data, layout=layout)
        png = frame_fig.to_image(
            format="png", width=width, height=height, engine="kaleido"
        )
        images.append(Image.open(io.BytesIO(png)))

    images[0].save(
        path,
        save_all=True,
        append_images=images[1:],
        duration=int(1000 / fps),
        loop=0,
        optimize=True,
    )


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

    target = scheduled.get("target", scheduled["unique_id"])
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

    # State-aware (cable-wrap) plotting: when the schedule carries per-node
    # wrap_state and the model knows its wrap_states, draw the azimuth path in
    # continuous encoder coordinates of the chosen winding instead of sky az.
    wrap_states = getattr(model, "wrap_states", None)
    state_aware = wrap_states is not None and "wrap_state" in scheduled.columns
    node_state = (
        scheduled["wrap_state"].to_numpy() if state_aware else None
    )

    obs_time = np.empty(2 * len(scheduled))
    az_path = np.empty(2 * len(scheduled))
    alt_path = np.empty(2 * len(scheduled))
    state_path = np.empty(2 * len(scheduled))
    names = []
    for i in range(len(scheduled)):
        obs_time[2 * i] = t_start_time[i].jd
        obs_time[2 * i + 1] = t_end_time[i].jd
        az_path[2 * i], az_path[2 * i + 1] = az_start[i], az_end[i]
        alt_path[2 * i], alt_path[2 * i + 1] = alt_start[i], alt_end[i]
        if state_aware:
            state_path[2 * i] = node_state[i]
            state_path[2 * i + 1] = node_state[i]
        names.extend([target.iloc[i], target.iloc[i]])

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

    if state_aware:
        # Encoder azimuth of each node's chosen winding (continuous, physical).
        state_path = state_path[:min_len].astype(int)
        az_path_display = az_path.copy()
        for s, (_, lo, hi) in enumerate(wrap_states):
            sel = state_path == s
            az_path_display[sel] = _encoder_az_display(az_path[sel], lo, hi)
    else:
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

    if state_aware:
        # Draw each used winding's encoder-azimuth bounds as reference lines.
        for s, (name, lo, hi) in enumerate(wrap_states):
            if not (state_path == s).any():
                continue
            for edge in (lo, hi):
                fig.add_shape(
                    type="line",
                    x0=obs_time[0], x1=obs_time[-1], y0=edge, y1=edge,
                    line=dict(color="red", dash="dash", width=1),
                    row=1, col=1,
                )
            fig.add_annotation(
                x=obs_time[-1], y=hi,
                text=f"{name}-wrap [{lo:g}, {hi:g}]\u00b0",
                showarrow=False, font=dict(color="red", size=10),
                row=1, col=1,
            )
    elif wrap is not None:
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
