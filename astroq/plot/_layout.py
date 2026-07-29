"""Shared Plotly layout presets and helpers for semester figures."""

from datetime import datetime

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

from astroq.plot._common import clear, labelsize

FIG_WIDTH = 1400
TIMELINE_HEIGHT = 1000
LOG_SCATTER_HEIGHT = 800

TIMEBAR_COLORS = ["#FF0000", "#F18F01", "#A23B72", "#2E86AB", "#00FF00"]
TIMEBAR_CATEGORY_NAMES = [
    "Unused",
    "Incomplete",
    "Future Scheduled",
    "Past Completed",
    "Requested",
]
TIMEBAR_WEATHER_LOSS_FACTOR = 0.2

_HORIZONTAL_LEGEND = dict(
    orientation="h",
    x=0.5,
    y=-0.15,
    xanchor="center",
    yanchor="top",
    bgcolor="rgba(255,255,255,0.7)",
    bordercolor="black",
    borderwidth=1,
    font=dict(size=labelsize - 18),
    itemsizing="constant",
    itemwidth=30,
    groupclick="toggleitem",
    tracegroupgap=5,
    traceorder="normal",
)

pio.templates["astroq_semester"] = go.layout.Template(
    layout=go.Layout(
        template="plotly_white",
        width=FIG_WIDTH,
        plot_bgcolor=clear,
        paper_bgcolor=clear,
        font=dict(size=labelsize),
        xaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            showgrid=False,
        ),
        yaxis=dict(
            title_font=dict(size=labelsize),
            tickfont=dict(size=labelsize - 4),
            showgrid=False,
        ),
    )
)


def semester_night_ticks(semester_planner, step=23):
    """Night-index tick positions plus calendar labels for the date axis."""
    semester_length = semester_planner.semester_length
    tickvals = list(range(0, semester_length, step))
    if (semester_length - 1) not in tickvals:
        tickvals.append(semester_length - 1)
    ticktext = [str(val + 1) for val in tickvals]

    ticktext_dates = []
    dates = semester_planner.access_obj.all_dates_array
    for day_idx in tickvals:
        if day_idx < len(dates):
            date_obj = datetime.strptime(dates[day_idx], "%Y-%m-%d")
            ticktext_dates.append(
                f"{date_obj.strftime('%b')}<br>{date_obj.strftime('%d')}"
            )
        else:
            ticktext_dates.append("")
    return tickvals, ticktext, ticktext_dates


def semester_night_xaxis(semester_planner, tickvals, ticktext, *, zeroline=False):
    """Primary bottom x-axis: night index in semester."""
    axis = dict(
        title_font=dict(size=labelsize),
        tickfont=dict(size=labelsize - 4),
        tickvals=tickvals,
        ticktext=ticktext,
        tickmode="array",
        showgrid=False,
        anchor="y",
        side="bottom",
        range=[0, semester_planner.semester_length - 1],
    )
    if zeroline:
        axis["zeroline"] = False
    return axis


def semester_date_xaxis2(semester_planner, tickvals, ticktext_dates):
    """Secondary top x-axis: calendar dates aligned to night ticks."""
    return dict(
        title="",
        tickvals=tickvals,
        ticktext=ticktext_dates,
        tickmode="array",
        showgrid=False,
        side="top",
        overlaying="x",
        tickfont=dict(size=labelsize - 6),
        showticklabels=True,
        range=[0, semester_planner.semester_length - 1],
    )


def horizontal_legend():
    """Standard horizontal legend below semester timeline figures."""
    return dict(_HORIZONTAL_LEGEND)


def timeline_layout(
    semester_planner,
    *,
    yaxis_title,
    yaxis,
    xaxis_title="Night in Semester",
    tickvals=None,
    ticktext=None,
    ticktext_dates=None,
    **overrides,
):
    """Layout preset for COF / birdseye semester-timeline figures."""
    if tickvals is None:
        tickvals, ticktext, ticktext_dates = semester_night_ticks(semester_planner)
    layout = {
        "width": FIG_WIDTH,
        "height": TIMELINE_HEIGHT,
        "template": "astroq_semester",
        "plot_bgcolor": clear,
        "paper_bgcolor": clear,
        "xaxis_title": xaxis_title,
        "yaxis_title": yaxis_title,
        "showlegend": True,
        "legend": horizontal_legend(),
        "xaxis": semester_night_xaxis(semester_planner, tickvals, ticktext),
        "xaxis2": semester_date_xaxis2(
            semester_planner, tickvals, ticktext_dates
        ),
        "yaxis": yaxis,
        "margin": dict(b=200, t=100),
    }
    layout.update(overrides)
    return layout


def today_vrect(fig, night_index):
    """Vertical dashed line marking the current semester night."""
    fig.add_vrect(
        x0=night_index,
        x1=night_index,
        annotation_text="Today",
        line_dash="dash",
        fillcolor=None,
        line_width=2,
        line_color="black",
        annotation_position="bottom left",
    )


def force_secondary_xaxis(fig, x_min, x_max, y_dummy, *, marker_size=0.01):
    """Invisible scatter trace that forces an overlay x-axis to render."""
    fig.add_trace(
        go.Scatter(
            x=[x_min, x_max],
            y=[y_dummy, y_dummy],
            mode="markers",
            marker=dict(size=marker_size, opacity=0),
            showlegend=False,
            hoverinfo="skip",
            xaxis="x2",
            name="",
        )
    )


def interval_axis_ticks(x_min, x_max, interval=60):
    """Tick positions at regular intervals (e.g. 0, 60, 120 minutes)."""
    first = 0 if x_min <= 0 else int(np.ceil(x_min / interval)) * interval
    tickvals = []
    value = float(first)
    while value <= x_max + 1e-9:
        if value >= x_min - 1e-9:
            tickvals.append(value)
        value += interval
    return tickvals, [str(int(round(v))) for v in tickvals]


def add_off_canvas_legend_rects(
    fig,
    items,
    *,
    x0=-100,
    x1=-80,
    y0=-0.5,
    y1=0.5,
):
    """Legend entries as off-canvas rectangles (used by the night ladder plot).

    Each item is ``(fillcolor, name)`` or ``(fillcolor, name, line_dict)``.
    """
    for item in items:
        if len(item) == 3:
            fillcolor, name, line = item
        else:
            fillcolor, name = item
            line = dict(width=0)
        fig.add_shape(
            type="rect",
            x0=x0,
            x1=x1,
            y0=y0,
            y1=y1,
            fillcolor=fillcolor,
            line=line,
            showlegend=True,
            name=name,
        )


def hide_x2_legend_entries(fig):
    """Keep invisible secondary-axis traces out of the legend."""
    for trace in fig.data:
        if hasattr(trace, "xaxis") and str(trace.xaxis) == "x2":
            trace.update(showlegend=False)
        if hasattr(trace, "name") and (trace.name == "" or trace.name is None):
            trace.update(showlegend=False)


def birdseye_slot_yaxis(semester_planner):
    """Y-axis tick labels every two hours of slot index."""
    slot_size = semester_planner.config.getint("semester", "slot_size")
    n_slots = int(24 * 60 // slot_size)
    slots_per_2hr = int(2 * 60 // slot_size)
    y_tickvals = list(range(0, n_slots, slots_per_2hr))
    y_ticktext = []
    for slot in y_tickvals:
        total_minutes = slot * slot_size
        y_ticktext.append(f"{total_minutes // 60:02.0f}:{total_minutes % 60:02.0f}")
    yaxis = dict(
        title_font=dict(size=labelsize),
        tickfont=dict(size=labelsize - 4),
        tickvals=y_tickvals,
        ticktext=y_ticktext,
        tickmode="array",
        showgrid=False,
    )
    return yaxis, n_slots


def log_axis_cadence():
    """Log axis preset for tau_inter cadence scatter (0.5–180 days)."""
    return dict(
        type="log",
        title_font=dict(size=labelsize),
        tickfont=dict(size=labelsize - 4),
        showgrid=True,
        gridcolor="lightgray",
        gridwidth=0.5,
        tickmode="array",
        tickvals=[1, 10, 100],
        ticktext=["1", "10", "100"],
        range=[np.log10(0.5), np.log10(180)],
    )


def log_axis_counts():
    """Log axis preset for rawobs observation-count scatter."""
    return dict(
        type="log",
        title_font=dict(size=labelsize),
        tickfont=dict(size=labelsize - 4),
        showgrid=True,
        gridcolor="lightgray",
        minor=dict(showgrid=False, ticks=""),
        dtick=1,
    )


def log_scatter_layout(*, xaxis_title, yaxis_title, xaxis, yaxis, **overrides):
    """Layout preset for log-scale scatter semester figures."""
    layout = {
        "width": FIG_WIDTH,
        "height": LOG_SCATTER_HEIGHT,
        "template": "astroq_semester",
        "plot_bgcolor": clear,
        "paper_bgcolor": clear,
        "xaxis_title": xaxis_title,
        "yaxis_title": yaxis_title,
        "xaxis": xaxis,
        "yaxis": yaxis,
    }
    layout.update(overrides)
    return layout


def timebar_category_labels():
    """Rich HTML y-axis labels for the single-program timebar chart."""
    return [
        (
            "<b>Unused Time</b><br>(allocation - past - future)<br>"
            "If you have positive unused time, <br>"
            "consider adding or changing requests"
        ),
        (
            "<b>Incomplete Time</b><br>(requested - past - future)<br>"
            "If you have incomplete time, <br>"
            "some of your requests are infeasible <br> consider changing them, "
            "<br> i.e. cadence or redistributing"
        ),
        "<b>Future Scheduled Time</b>",
        "<b>Past Completed Time</b>",
        "<b>Requested Time</b>",
    ]


def timebar_clip_hours(unused, incomplete, prevent_negative):
    """Optionally clamp incomplete/unused categories at zero."""
    if prevent_negative:
        incomplete = max(0, incomplete)
        unused = max(0, unused)
    return unused, incomplete


def timebar_text_labels(values, total_allocated_hours):
    """Bar annotation text: hours and percent of allocation."""
    labels = []
    for val in values:
        pct = (val / total_allocated_hours * 100) if total_allocated_hours > 0 else 0
        labels.append(f"{val:.1f} hrs ({pct:.1f}%)")
    return labels


def timebar_single_title(
    total_requested_hours,
    total_allocated_hours,
    total_allocated_nights,
    hours_per_night,
):
    """HTML title block for the aggregate timebar chart."""
    return (
        f"<b>Total Requested:</b> {total_requested_hours:.1f} hours ≈ "
        f"{total_requested_hours / hours_per_night:.1f} nights<br>"
        f"<b>Total Allocated:</b> {total_allocated_hours:.1f} hours ≈ "
        f"{total_allocated_nights:.1f} nights ----> w/ losses = "
        f"{total_allocated_nights * 0.75:.1f} nights <br>"
        "Requested and allocated time are measured in hours "
        f"({hours_per_night:.0f} hours per night for night equivalents).<br>"
        "Past and future bars use splan charged hours (slot-based)."
    )


def timebar_layout(title_text, top_margin, **overrides):
    """Layout preset for the single-program timebar chart."""
    layout = {
        "title_text": title_text,
        "template": "astroq_semester",
        "showlegend": False,
        "height": 710,
        "width": FIG_WIDTH,
        "margin": dict(t=top_margin, b=50, l=200, r=50),
        "bargap": 0.2,
        "xaxis": dict(
            title="Hours",
            title_font=dict(size=14),
            tickfont=dict(size=12),
        ),
        "yaxis": dict(
            title="",
            title_font=dict(size=14),
            tickfont=dict(size=11),
        ),
    }
    layout.update(overrides)
    return layout


def add_timebar_guide_lines(
    fig,
    labels,
    *,
    total_allocated_hours,
    max_schedulable_hours,
    weather_loss_factor=TIMEBAR_WEATHER_LOSS_FACTOR,
):
    """Vertical reference lines and invisible hover scatters for timebar charts."""
    y0, y1 = -0.5, len(labels) - 0.5

    fig.add_shape(
        type="line",
        x0=total_allocated_hours,
        x1=total_allocated_hours,
        y0=y0,
        y1=y1,
        line=dict(color="black", width=2, dash="dash"),
        xref="x",
        yref="y",
    )

    weather_loss_value = total_allocated_hours * (1 - weather_loss_factor)
    fig.add_shape(
        type="line",
        x0=weather_loss_value,
        x1=weather_loss_value,
        y0=y0,
        y1=y1,
        line=dict(color="gray", width=2, dash="dash"),
        xref="x",
        yref="y",
    )

    fig.add_shape(
        type="line",
        x0=max_schedulable_hours,
        x1=max_schedulable_hours,
        y0=y0,
        y1=y1,
        line=dict(color="gray", width=2, dash="dash"),
        xref="x",
        yref="y",
    )

    fig.add_trace(
        go.Scatter(
            x=[total_allocated_hours] * len(labels),
            y=labels,
            mode="markers",
            marker=dict(size=20, opacity=0),
            hovertemplate=(
                f"<b>Allocated Time</b><br>{total_allocated_hours:.2f} hours<br>"
                "This line represents the total allocated time for your program"
                "<extra></extra>"
            ),
            hoverlabel=dict(bgcolor="black", font_color="white"),
            showlegend=False,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=[weather_loss_value] * len(labels),
            y=labels,
            mode="markers",
            marker=dict(size=20, opacity=0),
            hovertemplate=(
                f"<b>Weather Loss Factor</b><br>{weather_loss_value:.2f} hours<br>"
                f"Allocated time minus {weather_loss_factor * 100:.0f}% weather loss<br>"
                "This is only a first order estimate based on historical losses."
                "<extra></extra>"
            ),
            hoverlabel=dict(bgcolor="gray", font_color="white"),
            showlegend=False,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=[max_schedulable_hours] * len(labels),
            y=labels,
            mode="markers",
            marker=dict(size=20, opacity=0),
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


def add_timebar_subplot_guides(
    fig,
    *,
    row,
    col,
    xref,
    yref,
    program_code,
    allocated,
    max_schedulable,
    category_names,
    max_fillfactor=1.25,
    weather_loss_factor=TIMEBAR_WEATHER_LOSS_FACTOR,
):
    """Reference lines and hover scatters for one timebar-by-program subplot."""
    y0, y1 = -0.5, 4.5
    hover_y = category_names[2]

    fig.add_shape(
        type="line",
        x0=allocated,
        x1=allocated,
        y0=y0,
        y1=y1,
        line=dict(color="black", width=2, dash="dash"),
        xref=xref,
        yref=yref,
    )

    weather_loss_value = allocated * (1 - weather_loss_factor)
    fig.add_shape(
        type="line",
        x0=weather_loss_value,
        x1=weather_loss_value,
        y0=y0,
        y1=y1,
        line=dict(color="gray", width=2, dash="dash"),
        xref=xref,
        yref=yref,
    )

    fig.add_shape(
        type="line",
        x0=max_schedulable,
        x1=max_schedulable,
        y0=y0,
        y1=y1,
        line=dict(color="gray", width=2, dash="dash"),
        xref=xref,
        yref=yref,
    )

    fig.add_trace(
        go.Scatter(
            x=[allocated],
            y=[hover_y],
            mode="markers",
            marker=dict(size=15, opacity=0),
            hovertemplate=(
                f"<b>{program_code} Allocated Time</b><br>{allocated:.2f} hours<br>"
                "Total allocated time for this program<extra></extra>"
            ),
            hoverlabel=dict(bgcolor="black", font_color="white"),
            showlegend=False,
        ),
        row=row,
        col=col,
    )

    fig.add_trace(
        go.Scatter(
            x=[weather_loss_value],
            y=[hover_y],
            mode="markers",
            marker=dict(size=15, opacity=0),
            hovertemplate=(
                f"<b>{program_code} Weather Loss Factor</b><br>"
                f"{weather_loss_value:.2f} hours<br>"
                f"Allocated time minus {weather_loss_factor * 100:.0f}% weather loss"
                "<extra></extra>"
            ),
            hoverlabel=dict(bgcolor="gray", font_color="white"),
            showlegend=False,
        ),
        row=row,
        col=col,
    )

    fig.add_trace(
        go.Scatter(
            x=[max_schedulable],
            y=[hover_y],
            mode="markers",
            marker=dict(size=15, opacity=0),
            hovertemplate=(
                f"<b>{program_code} Maximum Schedulable</b><br>"
                f"{max_schedulable:.2f} hours<br>"
                f"Awarded hours times max_fillfactor ({max_fillfactor:.2f})"
                "<extra></extra>"
            ),
            hoverlabel=dict(bgcolor="gray", font_color="white"),
            showlegend=False,
        ),
        row=row,
        col=col,
    )


def timebar_grid_layout(title_text, num_rows, **overrides):
    """Layout preset for the per-program timebar grid."""
    layout = {
        "title_text": title_text,
        "template": "astroq_semester",
        "showlegend": False,
        "height": max(600, num_rows * 250),
        "width": FIG_WIDTH,
        "margin": dict(t=150, b=50, l=50, r=50),
    }
    layout.update(overrides)
    return layout


def completion_layout(**overrides):
    """Minimal layout preset for completion-rate figures."""
    layout = {
        "width": FIG_WIDTH,
        "template": "astroq_semester",
        "plot_bgcolor": clear,
        "paper_bgcolor": clear,
    }
    layout.update(overrides)
    return layout
