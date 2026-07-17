"""
Module for constructing the standard AstroQ plots. All plots are returned as html strings.
From there, they can be used as is or saved as png files.
"""

# Standard library imports
from collections import defaultdict
from datetime import datetime, timedelta
from html import escape as html_escape
from urllib.parse import quote
import os
import base64
import re
from io import BytesIO

# Third-party imports
import numpy as np
import pandas as pd
import seaborn as sns
import astropy.units as u
from astropy.coordinates import SkyCoord
import astroplan as apl
import jinja2
import matplotlib
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from astropy.time import TimeDelta
from scipy.interpolate import griddata

# Local imports
import astroq.access as ac

# Jinja templates ship with the webapp subpackage. ``importlib.resources`` keeps
# the lookup correct for both source checkouts and installed wheels.
from importlib.resources import files as _resource_files
from pathlib import Path as _Path

_TEMPLATE_DIR = str(_resource_files("astroq.webapp").joinpath("templates"))


def _football_cache_dir(semester_planner):
    """Directory holding the cached sky-availability grids for ``get_football``.

    One cache per workdir keeps the on-disk artifacts adjacent to the run that
    produced them and avoids polluting the installed package. Tests monkeypatch
    this function to redirect the cache into a tmp dir.
    """
    return _Path(semester_planner.config.get("global", "workdir")) / "cache"


_TEMPLATE_ENV = jinja2.Environment(
    loader=jinja2.FileSystemLoader(_TEMPLATE_DIR),
    # We intentionally embed pandas-rendered <table> HTML and CSS/JS verbatim;
    # values that need escaping (e.g. tooltips) are escaped explicitly via html_escape.
    autoescape=False,
    trim_blocks=True,
    lstrip_blocks=True,
)

# Configure matplotlib for headless rendering
matplotlib.use("Agg")

# just used for color reproducibility
np.random.seed(24)

# Global variables from dynamic.py
gray = "rgb(210,210,210)"
clear = "rgba(255,255,255,1)"
labelsize = 38
hours_per_night = 12.0


def cumulative_by_night(ps, n_nights, *, group_col, group_val, metric="visits"):
    """Cumulative visits or slot sum aligned to night index 0..n_nights-1."""
    sub = ps.loc[ps[group_col] == group_val]
    if sub.empty:
        return np.zeros(n_nights, dtype=float)
    if metric == "visits":
        daily = sub.groupby(sub.index).size()
    else:
        daily = sub.groupby(sub.index)["t_visit_slots"].sum()
    return daily.reindex(range(n_nights), fill_value=0).cumsum().to_numpy(dtype=float)


def daily_visits_by_night(ps, n_nights, *, group_col, group_val):
    """Per-night visit counts (non-cumulative) aligned to night index."""
    sub = ps.loc[ps[group_col] == group_val]
    if sub.empty:
        return np.zeros(n_nights, dtype=int)
    return (
        sub.groupby(sub.index)
        .size()
        .reindex(range(n_nights), fill_value=0)
        .to_numpy(dtype=int)
    )


def _visit_counts_by_date(
    ps, *, group_col, group_val, past, today_idx, all_dates_array
):
    """Map calendar date -> visit count for one star or program."""
    sub = ps.loc[ps[group_col] == group_val]
    if past:
        sub = sub.loc[sub.index < today_idx]
    else:
        sub = sub.loc[sub.index >= today_idx]
    if sub.empty:
        return {}
    return {
        all_dates_array[int(d)]: int(n)
        for d, n in sub.groupby(sub.index).size().items()
    }


def _charged_hours_from_ps(semester_planner, ps, *, program_codes=None, unique_ids=None):
    """Past and scheduled charged hours from ``timeline``."""
    slot_size = semester_planner.config.getfloat("semester", "slot_size")
    slots_per_hour = 60 / slot_size
    sub = ps
    if program_codes is not None:
        sub = sub[sub["program_code"].isin(program_codes)]
    if unique_ids is not None:
        uids = {str(u) for u in unique_ids}
        sub = sub[sub["unique_id"].isin(uids)]
    today_idx = semester_planner.access_obj.current_night_index
    past_h = (
        sub.loc[sub.index < today_idx, "t_visit_slots"].sum() / slots_per_hour
    )
    sched_h = (
        sub.loc[sub.index >= today_idx, "t_visit_slots"].sum() / slots_per_hour
    )
    return float(past_h), float(sched_h)


def _cof_pct_curve(
    semester_planner,
    ps,
    n_nights,
    *,
    group_col,
    group_val,
    use_time,
    denominator,
):
    """Cumulative COF % array for one program or request."""
    metric = "slots" if use_time else "visits"
    cume = cumulative_by_night(
        ps, n_nights, group_col=group_col, group_val=group_val, metric=metric
    )
    if use_time:
        slot_size = semester_planner.config.getfloat("semester", "slot_size")
        cume = cume / (60 / slot_size)
    if denominator > 0:
        return np.round(cume / denominator * 100, 2)
    if not use_time and cume[-1] > 0:
        return np.round(cume / cume[-1] * 100, 2)
    return np.zeros(n_nights, dtype=float)


def _cof_group_for_star(star):
    """Return (group_col, group_val) for slicing ``timeline``."""
    if getattr(star, "allow_mapview", True) is False:
        return "program_code", star.program
    return "unique_id", str(star.unique_id)


def _visit_denominator(star):
    return getattr(star, "requested_visits", star.total_observations_requested)


def programs_ledger_for_plot(semester_planner):
    """Timeline charged hours per program without requiring a live Gurobi model."""
    if semester_planner.schedule is None:
        raise RuntimeError("No schedule created")
    sph = semester_planner.slots_per_hour
    idx = semester_planner.programs.index
    ps = semester_planner.timeline
    today = semester_planner.access_obj.current_night_index

    def hours(mask):
        return (
            ps.loc[mask, "t_visit_slots"]
            .groupby(ps.loc[mask, "program_code"])
            .sum()
            .reindex(idx, fill_value=0)
            / sph
        )

    return semester_planner.programs.assign(
        past_hours=hours(ps.index < today),
        sched_hours=hours(ps.index >= today),
    )


def _render_datatable(
    df,
    *,
    template_name,
    table_id,
    variant,
    column_widths,
    column_defs,
    tooltips=None,
    page_size=25,
    sort_column=0,
    sort_dir="asc",
    numeric_cols=(),
    time_cols=(),
    has_band_padding=False,
    has_column_filters=False,
    has_dt_controls_styling=False,
    has_init_complete_header_style=False,
    filter_placeholder="",
    auto_width=False,
    table_layout=None,
    scroll_x=None,
    responsive=None,
    add_tfoot=False,
):
    """Render a DataFrame to a self-contained DataTables HTML fragment.

    All public table-rendering helpers in this module funnel through this function.
    The Jinja templates live under ``templates/`` and ``templates/partials/``.

    Args:
        df (pd.DataFrame): table to render. Must already have its display column
            names and HTML cell contents finalized (the renderer doesn't transform
            cell values).
        template_name (str): page template, e.g. ``'request_table.html.j2'``.
        table_id (str): DOM id for the rendered ``<table>``.
        variant (str): ``'compact'`` (request/nightplan) or ``'card'`` (legacy generic).
        column_widths (list[str]): one width spec per column, e.g. ``'12ch'`` or
            ``'100px'``. Used by the CSS partial in ``compact`` variant.
        column_defs (list[dict]): list of ``{'target': i, 'width': w}`` dicts that
            DataTables consumes via ``columnDefs``.
        tooltips (list[str] | None): per-column tooltip strings. Any non-empty
            entry triggers a ``data-tooltip`` attribute on the matching ``<th>`` and
            enables the tooltip CSS block.
        page_size, sort_column, sort_dir: DataTables init values.
        numeric_cols, time_cols: column indices that filter as numbers / HH:MM.
        has_band_padding: enable tighter padding for the band columns
            (request_frame layout only).
        has_column_filters: render per-column filter inputs in ``<tfoot>``.
        has_dt_controls_styling: style DataTables length/buttons/paginate (card variant).
        has_init_complete_header_style: re-apply header CSS in ``initComplete``
            (legacy generic table; redundant with the CSS but matches prior output).
        filter_placeholder (str): placeholder for column filter inputs.
        auto_width, table_layout, scroll_x, responsive: passthroughs to DataTables.
        add_tfoot (bool): inject an empty ``<tfoot>`` so column-filter inputs have
            a row to mount onto.

    Returns:
        str: HTML fragment containing ``<style>`` + ``<table>`` + ``<script>``.
    """
    tooltips = list(tooltips) if tooltips else []
    has_tooltips = any(bool(t) for t in tooltips)

    table_html = df.to_html(
        classes="table table-striped table-hover",
        index=False,
        escape=False,
        table_id=table_id,
    )

    if has_tooltips:
        # pandas emits one <th>label</th> per column header; replace each in order.
        idx = [0]

        def _add_th_tooltip(m):
            i = idx[0]
            idx[0] += 1
            t = tooltips[i] if i < len(tooltips) else ""
            return (
                f'<th data-tooltip="{html_escape(t)}">{m.group(1)}</th>'
                if t
                else m.group(0)
            )

        table_html = re.sub(
            r"<th>([^<]*)</th>",
            _add_th_tooltip,
            table_html,
            count=len(df.columns),
        )

    if add_tfoot:
        tfoot_cells = "".join("<th></th>" for _ in df.columns)
        table_html = table_html.replace(
            "</tbody>", "</tbody><tfoot><tr>" + tfoot_cells + "</tr></tfoot>"
        )

    template = _TEMPLATE_ENV.get_template(template_name)
    ctx = {
        "table_id": table_id,
        "table_html": table_html,
        "variant": variant,
        "column_widths": list(column_widths),
        "column_defs": list(column_defs),
        "page_size": int(page_size),
        "sort_column": int(sort_column),
        "sort_dir": sort_dir,
        "numeric_cols": list(numeric_cols),
        "time_cols": list(time_cols),
        "has_tooltips": has_tooltips,
        "has_band_padding": has_band_padding,
        "has_column_filters": has_column_filters,
        "has_dt_controls_styling": has_dt_controls_styling,
        "has_init_complete_header_style": has_init_complete_header_style,
        "filter_placeholder": filter_placeholder,
        "auto_width": auto_width,
        "table_layout": table_layout,
    }
    if scroll_x is not None:
        ctx["scroll_x"] = scroll_x
    if responsive is not None:
        ctx["responsive"] = responsive
    return template.render(**ctx)

