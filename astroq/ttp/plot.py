"""TTP-specific Plotly plotting utilities.

This module owns the plots that visualize a single-night TTP solution
(``astroq.ttp.model.TTPModel``):

- :func:`get_ladder` -- Gantt-style night plan ladder (targets vs time).
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
from datetime import timedelta

import numpy as np
import pandas as pd

# Third-party imports
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
import astropy.units as u
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from astroq.plot._layout import (
    add_off_canvas_legend_rects,
    force_secondary_xaxis,
    interval_axis_ticks,
)


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


# --- Ladder plot (night plan Gantt) ---

_LADDER_WIDTH = 800
_VISIT_COLOR = "lightgreen"
_SLEW_COLOR = "#FFE4B5"
_IDLE_COLOR = "#FFB3B3"
_ACCESS_LINE_COLOR = "lightgreen"
_ACCESS_LINE_WIDTH = 2
_ACCESS_FILL_COLOR = "white"
_LADDER_BG = "#f4f4f4"
_SEGMENT_COLORS = {"visit": _VISIT_COLOR, "slew": _SLEW_COLOR, "idle": _IDLE_COLOR}

_LADDER_LEGEND_ITEMS = [
    (_VISIT_COLOR, "Visit"),
    (_SLEW_COLOR, "Slew"),
    (_IDLE_COLOR, "Idle"),
    (_ACCESS_FILL_COLOR, "Accessible", dict(color=_ACCESS_LINE_COLOR, width=_ACCESS_LINE_WIDTH)),
]

_LADDER_HOVER_TEMPLATE = (
    "<b>%{customdata[0]}</b><br>"
    "Scheduled (min. from start): %{x:.1f}<br>"
    "Scheduled (UTC): %{customdata[1]}<br>"
    "Earliest start (min. from start): %{customdata[2]:.1f}<br>"
    "Earliest start (UTC): %{customdata[3]}<br>"
    "Latest finish (min. from start): %{customdata[4]:.1f}<br>"
    "Latest finish (UTC): %{customdata[5]}<br>"
    "Visit Length (min): %{customdata[6]:.1f}<br>"
    "Slew to Next (min): %{customdata[7]:.1f}"
    "<extra></extra>"
)


def _min_to_utc_hhmm(night_start, minutes):
    if night_start is None or pd.isna(minutes):
        return ""
    return (night_start + TimeDelta(float(minutes) * 60, format="sec")).isot[11:16]


def _floor_utc_hour(dt):
    return dt.replace(minute=0, second=0, microsecond=0)


def _ladder_utc_ticks(night_start_time, x_min, x_max):
    """Whole-hour UTC tick positions within ``[x_min, x_max]`` (minutes from night start)."""
    start_dt = night_start_time.to_datetime()
    utc_tickvals, utc_ticktext = [], []
    tick = _floor_utc_hour(start_dt)
    while True:
        offset_min = (tick - start_dt).total_seconds() / 60.0
        if offset_min > x_max + 1e-9:
            break
        if offset_min >= x_min - 1e-9:
            utc_tickvals.append(offset_min)
            utc_ticktext.append(tick.strftime("%H:%M"))
        tick += timedelta(hours=1)
    return utc_tickvals, utc_ticktext


def _add_ladder_night_boundary(fig, x, utc_hhmm, *, side):
    """Vertical marker at night start or end with UTC time label."""
    fig.add_shape(
        type="line",
        x0=x,
        x1=x,
        y0=0,
        y1=1,
        xref="x",
        yref="paper",
        line=dict(color="black", width=1.5),
        layer="above",
    )
    if side == "start":
        label = f"Night start {utc_hhmm} (UT)"
        xanchor = "left"
    else:
        label = f"Night end {utc_hhmm} (UT)"
        xanchor = "right"
    fig.add_annotation(
        x=x,
        y=0.99,
        xref="x",
        yref="paper",
        text=f"<b>{label}</b>",
        showarrow=False,
        yanchor="top",
        xanchor=xanchor,
        font=dict(size=11, color="black"),
    )


def _night_duration_min(model):
    """Night length in minutes (works on live and HDF5-reloaded models)."""
    dur = getattr(model, "dur_min", None)
    if dur is not None:
        return float(dur)
    stats = getattr(model, "stats", None) or {}
    if "dur_min" in stats:
        return float(stats["dur_min"])
    night_start = getattr(model, "night_start", None)
    night_end = getattr(model, "night_end", None)
    if night_start is not None and night_end is not None:
        return (night_end.jd - night_start.jd) * 24 * 60
    return 600.0


def _night_timeline_segments(model):
    """Chronological (idle, visit, slew) segments from 0 .. dur_min."""
    scheduled = model.schedule[~model.schedule["is_anchor"]]
    scheduled = scheduled[scheduled["scheduled"]].sort_values("order")
    dur = _night_duration_min(model)
    cursor = 0.0
    segments = []
    for _, row in scheduled.iterrows():
        t0, t1 = float(row["t_start"]), float(row["t_end"])
        slew = float(row["t_slew"]) if pd.notna(row["t_slew"]) else 0.0
        if t0 > cursor + 1e-9:
            segments.append(("idle", cursor, t0))
        segments.append(("visit", t0, t1))
        if slew > 0:
            segments.append(("slew", t1, t1 + slew))
        cursor = t1 + slew
    if cursor < dur - 1e-9:
        segments.append(("idle", cursor, dur))
    return segments


def _night_aggregate_segments(model):
    """Contiguous visit, slew, idle totals (left-aligned) for the aggregate summary row."""
    totals = {"visit": 0.0, "slew": 0.0, "idle": 0.0}
    for kind, x0, x1 in _night_timeline_segments(model):
        totals[kind] += x1 - x0
    cursor = 0.0
    segments = []
    for kind in ("visit", "slew", "idle"):
        width = totals[kind]
        if width <= 1e-9:
            continue
        segments.append((kind, cursor, cursor + width))
        cursor += width
    return segments


def _add_ladder_required_visit_bar(fig, *, y, x0, visit_len, showlegend=False):
    """Cross-hatched bar for an unscheduled request's required visit duration."""
    x1 = x0 + visit_len
    fig.add_trace(
        go.Scatter(
            x=[x0, x1, x1, x0, x0],
            y=[y, y, y, y, y],
            mode="lines",
            fill="toself",
            fillcolor="white",
            fillpattern=dict(
                shape="/",
                bgcolor="white",
                fgcolor="#666666",
                fgopacity=0.75,
                size=8,
                solidity=0.4,
            ),
            line=dict(color="#888888", width=1),
            hoverinfo="skip",
            showlegend=showlegend,
            name="Required visit",
            legendgroup="required_visit",
        )
    )


def _synthetic_ladder_row(columns):
    """Blank ladder row with NaN numerics and empty UTC strings."""
    row = {}
    for col in columns:
        if col == "_row_kind":
            continue
        if col.endswith("(UTC)") or col == "Scheduled (UTC)":
            row[col] = ""
        elif col == "is_scheduled":
            row[col] = False
        elif col == "Target":
            row[col] = " "
        else:
            row[col] = np.nan
    return row


def _ladder_x_max(model, order_data):
    """Upper x-axis limit in minutes since night start."""
    try:
        return _night_duration_min(model)
    except (AttributeError, KeyError, TypeError, ValueError):
        pass
    if len(order_data) == 0:
        return 600.0
    target_rows = order_data[order_data["_row_kind"] == "target"]
    if target_rows.empty:
        return 600.0
    end_times = (
        target_rows["Start Exposure"]
        + target_rows["Visit Length (min)"]
        + target_rows["Slew to Next (min)"]
    )
    return float(end_times.max())


def get_ladder(data, tonight_start_time):  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    """Produce a plotly figure which illustrates the night plan solution.

    Args:
        data: a TTP data object containing the schedule information
        tonight_start_time: Astropy Time for the night start (UTC axis labels)

    Returns:
        plotly.graph_objects.Figure: the night plan ladder figure.
    """
    model = _as_model(data)
    order_data = schedule_to_ladder_frame(model)
    if order_data.empty:
        order_data = pd.DataFrame(
            columns=[
                "unique_id",
                "Target",
                "Earliest Start",
                "Latest Finish",
                "Start Exposure",
                "Stop Exposure",
                "Visit Length (min)",
                "Slew to Next (min)",
                "Scheduled (min. from start)",
            ]
        )
    if "Slew to Next (min)" not in order_data.columns:
        order_data["Slew to Next (min)"] = 0.0

    night_start = getattr(model, "night_start", None)
    if night_start is not None and len(order_data):
        order_data["Scheduled (UTC)"] = [
            _min_to_utc_hhmm(night_start, se) if se > 0 else ""
            for se in order_data["Start Exposure"]
        ]
        order_data["Earliest start (UTC)"] = [
            _min_to_utc_hhmm(night_start, t) for t in order_data["Earliest Start"]
        ]
        order_data["Latest finish (UTC)"] = [
            _min_to_utc_hhmm(night_start, t) for t in order_data["Latest Finish"]
        ]
    elif len(order_data):
        order_data["Scheduled (UTC)"] = ""
        order_data["Earliest start (UTC)"] = ""
        order_data["Latest finish (UTC)"] = ""

    on_sky = model.schedule[~model.schedule["is_anchor"]]
    n_unscheduled = int((~on_sky["scheduled"]).sum())

    order_data = order_data.iloc[::-1].reset_index(drop=True)
    order_data["_row_kind"] = "target"

    n_before_insert = len(order_data)
    summary_y = None
    aggregate_y = None
    if 0 < n_unscheduled < n_before_insert:
        unsched_header = _synthetic_ladder_row(order_data.columns)
        unsched_header["unique_id"] = "__unsched_header__"
        unsched_header["Target"] = "Unscheduled targets"
        unsched_header["_row_kind"] = "section_header"

        aggregate = _synthetic_ladder_row(order_data.columns)
        aggregate["unique_id"] = "__aggregate__"
        aggregate["Target"] = " "
        aggregate["_row_kind"] = "aggregate"

        summary = _synthetic_ladder_row(order_data.columns)
        summary["unique_id"] = "__summary__"
        summary["Target"] = "All scheduled targets"
        summary["_row_kind"] = "summary"

        sched_header = _synthetic_ladder_row(order_data.columns)
        sched_header["unique_id"] = "__sched_header__"
        sched_header["Target"] = "Scheduled targets"
        sched_header["_row_kind"] = "section_header"

        order_data = pd.concat(
            [
                order_data.iloc[:n_unscheduled],
                pd.DataFrame([unsched_header]),
                pd.DataFrame([aggregate]),
                pd.DataFrame([summary]),
                order_data.iloc[n_unscheduled:],
                pd.DataFrame([sched_header]),
            ],
            ignore_index=True,
        )
        aggregate_y = n_unscheduled + 1
        summary_y = n_unscheduled + 2

    mask = order_data["_row_kind"] != "target"
    order_data.loc[mask, "Scheduled (min. from start)"] = np.nan
    order_data["_ladder_y"] = order_data.index.astype(str)

    plot_height = max(400, 40 * len(order_data) + 200)
    categories = order_data["_ladder_y"].tolist()
    y_ticktext = [
        "" if kind in ("section_header", "summary", "aggregate") else target
        for target, kind in zip(order_data["Target"], order_data["_row_kind"])
    ]
    fig = px.scatter(
        order_data,
        x="Scheduled (min. from start)",
        y="_ladder_y",
        title="Night Plan",
        width=_LADDER_WIDTH,
        height=plot_height,
    )
    fig.update_traces(
        customdata=np.column_stack(
            [
                order_data["Target"],
                order_data["Scheduled (UTC)"].fillna(""),
                order_data["Earliest Start"],
                order_data["Earliest start (UTC)"].fillna(""),
                order_data["Latest Finish"],
                order_data["Latest finish (UTC)"].fillna(""),
                order_data["Visit Length (min)"],
                order_data["Slew to Next (min)"],
            ]
        ),
        hovertemplate=_LADDER_HOVER_TEMPLATE,
        marker=dict(size=0, opacity=0),
    )
    fig.update_layout(
        margin=dict(l=160),
        plot_bgcolor=_LADDER_BG,
        paper_bgcolor="white",
        yaxis_title="",
        yaxis=dict(
            categoryorder="array",
            categoryarray=categories,
            tickmode="array",
            tickvals=categories,
            ticktext=y_ticktext,
        ),
    )

    add_off_canvas_legend_rects(fig, _LADDER_LEGEND_ITEMS)

    new_already_processed = []
    ifixer = 0
    required_visit_bars = []
    for i in range(len(order_data)):
        if order_data["_row_kind"].iloc[i] != "target":
            continue
        if order_data["unique_id"][i] not in new_already_processed:
            indices = [
                k
                for k in range(len(order_data))
                if order_data["unique_id"][k] == order_data["unique_id"][i]
            ]
            for j, idx in enumerate(indices):
                if j == 0:
                    fig.add_shape(
                        type="rect",
                        x0=order_data["Earliest Start"][idx],
                        x1=order_data["Latest Finish"][idx],
                        y0=i + ifixer - 0.5,
                        y1=i + ifixer + 0.5,
                        fillcolor=_ACCESS_FILL_COLOR,
                        line=dict(color=_ACCESS_LINE_COLOR, width=_ACCESS_LINE_WIDTH),
                        showlegend=False,
                    )
                start_exp = float(order_data["Start Exposure"][idx])
                visit_len = float(order_data["Visit Length (min)"][idx])
                is_scheduled = bool(order_data["is_scheduled"].iloc[idx])
                if is_scheduled and start_exp > 0:
                    fig.add_shape(
                        type="rect",
                        x0=start_exp,
                        x1=start_exp + visit_len,
                        y0=i + ifixer - 0.5,
                        y1=i + ifixer + 0.5,
                        fillcolor=_VISIT_COLOR,
                        line=dict(width=0),
                    )
                elif not is_scheduled and visit_len > 0:
                    earliest = float(order_data["Earliest Start"][idx])
                    if not np.isnan(earliest):
                        required_visit_bars.append(
                            (
                                order_data["_ladder_y"].iloc[idx],
                                earliest,
                                visit_len,
                            )
                        )
                slew = float(order_data["Slew to Next (min)"][idx])
                if is_scheduled and slew > 0:
                    fig.add_shape(
                        type="rect",
                        x0=order_data["Stop Exposure"][idx],
                        x1=order_data["Stop Exposure"][idx] + slew,
                        y0=i + ifixer - 0.5,
                        y1=i + ifixer + 0.5,
                        fillcolor=_SLEW_COLOR,
                        line=dict(width=0),
                    )
            new_already_processed.append(order_data["unique_id"][i])
        else:
            ifixer -= 1

    for bar_idx, (target, x0, visit_len) in enumerate(required_visit_bars):
        _add_ladder_required_visit_bar(
            fig,
            y=target,
            x0=x0,
            visit_len=visit_len,
            showlegend=(bar_idx == 0),
        )

    if aggregate_y is not None:
        for kind, x0, x1 in _night_aggregate_segments(model):
            fig.add_shape(
                type="rect",
                x0=x0,
                x1=x1,
                y0=aggregate_y - 0.5,
                y1=aggregate_y + 0.5,
                fillcolor=_SEGMENT_COLORS[kind],
                line=dict(width=0),
                showlegend=False,
            )

    if summary_y is not None:
        for kind, x0, x1 in _night_timeline_segments(model):
            fig.add_shape(
                type="rect",
                x0=x0,
                x1=x1,
                y0=summary_y - 0.5,
                y1=summary_y + 0.5,
                fillcolor=_SEGMENT_COLORS[kind],
                line=dict(width=0),
                showlegend=False,
            )

    for idx in range(len(order_data)):
        row = order_data.iloc[idx]
        if row["_row_kind"] not in ("section_header", "summary"):
            continue
        fig.add_annotation(
            y=idx,
            xref="paper",
            x=0,
            xanchor="right",
            text=f"<b>{row['Target']}</b>",
            showarrow=False,
            font=dict(size=13, color="black"),
        )

    x_min = 0.0
    x_max = _ladder_x_max(model, order_data)
    if tonight_start_time is not None:
        utc_tickvals, utc_ticktext = _ladder_utc_ticks(
            tonight_start_time, x_min, x_max
        )
    else:
        utc_tickvals, utc_ticktext = [], []

    min_tickvals, min_ticktext = interval_axis_ticks(x_min, x_max)

    for x_line in utc_tickvals:
        fig.add_shape(
            type="line",
            x0=x_line,
            x1=x_line,
            y0=0,
            y1=1,
            xref="x",
            yref="paper",
            line=dict(color="white", width=1),
            layer="below",
        )

    if tonight_start_time is not None:
        night_end = getattr(model, "night_end", None)
        start_utc = tonight_start_time.isot[11:16]
        if night_end is not None:
            end_utc = night_end.isot[11:16]
        else:
            end_utc = _min_to_utc_hhmm(tonight_start_time, x_max)
        _add_ladder_night_boundary(fig, 0.0, start_utc, side="start")
        _add_ladder_night_boundary(fig, x_max, end_utc, side="end")

    y_ref = order_data["_ladder_y"].iloc[-1] if len(order_data) else ""
    force_secondary_xaxis(
        fig, x_min, x_max, y_ref, marker_size=0.001
    )
    fig.update_layout(
        xaxis=dict(
            title="time since start [min]",
            range=[x_min, x_max],
            tickmode="array",
            tickvals=min_tickvals,
            ticktext=min_ticktext,
            showgrid=False,
        ),
        xaxis2=dict(
            title=dict(text="time [UTC]", standoff=0),
            overlaying="x",
            side="top",
            range=[x_min, x_max],
            tickmode="array",
            tickvals=utc_tickvals,
            ticktext=utc_ticktext,
            showgrid=False,
            showline=True,
            mirror=True,
        ),
    )

    return fig


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
