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

    def get_stats(self, row, slot_size, queue):
        """
        Grab the observational stategy information for a given star from the requests.csv file.

        Args:
            row (pd.Series): A row from the requests.csv file as a DataFrame
            slot_size (int): The slot size in minutes
            queue (astroq.queue.base.Queue): provides ``slew_overhead_mean``
                and ``readout_time`` for the per-request seconds accounting.

        Returns:
            expected_nobs_per_night (int): how many exposures we expect to take
            total_observations_requested (int): sum of observational strategie values
            exposure_time (int): exposure time of single shot
            slots_per_night (int): number of slots required to complete all exposures in a night
            program (str): the program code
        """
        # Access row data directly instead of filtering the entire DataFrame (PERFORMANCE OPTIMIZATION)
        self.target = row["target"]
        self.inactive = row["inactive"]
        self.ra = float(row["ra"])
        self.dec = float(row["dec"])
        self.program = str(row["program_code"])
        self.exptime = int(row["exptime"])
        self.n_exp = int(row["n_exp"])
        self.n_intra_max = int(row["n_intra_max"])
        self.n_intra_min = int(row["n_intra_min"])
        self.tau_intra = int(row["tau_intra"])
        self.n_inter_max = int(row["n_inter_max"])
        self.tau_inter = int(row["tau_inter"])
        self.total_observations_requested = (
            self.n_exp * self.n_intra_max * self.n_inter_max
        )
        self.total_requested_seconds = (
            self.total_observations_requested * self.exptime
            + queue.readout_time * (self.n_exp - 1) * self.n_inter_max
            + queue.slew_overhead_mean * self.n_intra_max * self.n_inter_max
        )
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
        star_obs_past = past[past["target"] == str(self.unique_id)]
        # Parse date from timestamp and group by date
        star_obs_past = star_obs_past.copy()
        star_obs_past["date"] = star_obs_past["timestamp"].str[:10]
        observations_past = star_obs_past.groupby("date").size().to_dict()
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
        star_rows = forecast_df[forecast_df["unique_id"] == str(self.unique_id)]
        # Count number of slots scheduled per night (d)
        observations_future = {}
        for d, group in star_rows.groupby("d"):
            # d may be int or str; ensure it's int for indexing
            date = all_dates_array[int(d)]
            n_slots = len(
                group
            )  # this is the number of starting slots in given to this target in this night
            observations_future[date] = (
                n_slots  # no need to multiply by nexp here because we do it later in timebar; so that COF has right values.
            )
        self.observations_future = observations_future

    def get_map(self, semester_planner, forecast_df):
        """
        Build the 2D d/s matrix starmap for teh given star using semester_plan.csv.
        Only set starmap[d, s] = 1 if sched['unique_id'] == self.unique_id.

        Args:
            semester_planner: The semester planner object
            forecast_df (pd.DataFrame): Pre-loaded forecast DataFrame with columns ['unique_id', 'd', 's']
        """
        n_nights = semester_planner.semester_length
        n_slots = int((24 * 60) / semester_planner.config.getint("semester", "slot_size"))
        starmap = np.zeros((n_nights, n_slots), dtype=int)

        # Filter to only this star's rows
        star_forecast = forecast_df[forecast_df["unique_id"] == str(self.unique_id)]

        if len(star_forecast) > 0:
            # Vectorized approach: extract d,s values as numpy arrays and set all at once (PERFORMANCE OPTIMIZATION)
            d_values = star_forecast["d"].values.astype(int)
            s_values = star_forecast["s"].values.astype(int)

            # Set the primary slots
            starmap[d_values, s_values] = 1

            # Set the reserve slots
            rf = semester_planner.requests_frame
            row = rf.loc[rf["unique_id"] == str(self.unique_id)]
            reserve_slots = int(row["t_visit_slots"].iloc[0]) if len(row) else 1
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
        nulltime (array): a 2D array of N_slots by N_nights, binary 1/0, it is the intersection of is_allocated and is_night
    """

    # Create a starmap of the times when we cannot observe due to twilight and allocation constraints
    # Used in the birdseye view plot to blackout the unavailble squares

    # Use the stored access record from the semester planner instead of recomputing
    access = semester_planner.access_record
    nulltime = access["is_allocated"][0]
    nulltime = 1 - nulltime
    nulltime = np.array(nulltime).T

    forecast_df = semester_planner.schedule
    forecast_df["unique_id"] = forecast_df["unique_id"].astype(str)

    # Per-visit overhead scalars come from the queue (single source of truth).
    queue = semester_planner.queue
    slew_overhead = queue.slew_overhead_mean
    readout_overhead = queue.readout_time

    targets = semester_planner.requests_frame_all["target"].unique()
    programs = semester_planner.requests_frame_all["program_code"].unique()

    # Make colors consistent for all stars in each program
    colors = sns.color_palette("deep", len(programs))
    rgb_strings = [
        f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})" for r, g, b in colors
    ]
    program_colors_rgb_vals = dict(zip(programs, rgb_strings))

    # Per-(uid, night) past-observation aggregates derived from past_df. Keys
    # are UT calendar dates (timestamp[:10]). Empty dicts when past_df is empty.
    past_df = semester_planner.past_df
    if not past_df.empty:
        pdf = past_df.assign(
            unique_id=past_df["unique_id"].astype(str),
            night=past_df["timestamp"].astype(str).str[:10],
        )
        n_obs_by_uid = (
            pdf.groupby("unique_id").apply(
                lambda g: g.groupby("night").size().to_dict()
            ).to_dict()
        )
        n_visits_by_uid = (
            pdf.groupby("unique_id").apply(
                lambda g: g.groupby("night")["timestamp"].nunique().to_dict()
            ).to_dict()
        )
    else:
        n_obs_by_uid = {}
        n_visits_by_uid = {}

    all_stars = []
    i = 0
    for i, row in semester_planner.requests_frame_all.iterrows():
        # Create a StarPlotter object for each request, fill and compute relavant information
        newstar = StarPlotter(row["unique_id"])
        newstar.get_map(semester_planner, forecast_df)
        newstar.get_stats(row, semester_planner.config.getint("semester", "slot_size"), queue)
        uid = str(newstar.unique_id)
        newstar.observations_past = n_visits_by_uid.get(uid, {})
        newstar.observations_past_exposures = n_obs_by_uid.get(uid, {})
        newstar.get_future(forecast_df, semester_planner.all_dates_array)

        # Create COF arrays for each request
        combined_set = set(
            list(newstar.observations_past.keys())
            + list(newstar.observations_future.keys())
        )
        # For inactive stars, only include past observations; for active stars, include both past and future
        if newstar.inactive == False:
            newstar.dates_observe = [
                newstar.observations_past[date]
                if date in newstar.observations_past.keys()
                else (
                    newstar.observations_future[date] * newstar.n_exp
                    if date in combined_set
                    else 0
                )
                for date in semester_planner.all_dates_array
            ]
            newstar.dates_observe_time = [
                (
                    newstar.observations_past_exposures[date] * newstar.exptime
                    + readout_overhead * (newstar.observations_past[date] - 1)
                    + slew_overhead * (newstar.observations_past[date] - 1)
                )
                / 3600
                if date in newstar.observations_past_exposures.keys()
                else (
                    (
                        newstar.observations_future[date]
                        * newstar.n_exp
                        * newstar.exptime
                        + readout_overhead
                        * (newstar.n_exp - 1)
                        * newstar.observations_future[date]
                        + slew_overhead * newstar.observations_future[date]
                    )
                    / 3600
                    if date in combined_set
                    else 0
                )
                for date in semester_planner.all_dates_array
            ]
        else:
            # For inactive stars, only show past observations
            newstar.dates_observe = [
                newstar.observations_past[date]
                if date in newstar.observations_past.keys()
                else 0
                for date in semester_planner.all_dates_array
            ]
            newstar.dates_observe_time = [
                (
                    newstar.observations_past_exposures[date] * newstar.exptime
                    + readout_overhead * (newstar.observations_past[date] - 1)
                    + slew_overhead * (newstar.observations_past[date] - 1)
                )
                / 3600
                if date in newstar.observations_past_exposures.keys()
                else 0
                for date in semester_planner.all_dates_array
            ]

        newstar.cume_observe = np.cumsum(newstar.dates_observe)
        newstar.cume_observe_time = np.cumsum(newstar.dates_observe_time)  # in hours

        if newstar.inactive:
            newstar.total_observations_requested = np.max(newstar.cume_observe)
            newstar.total_requested_seconds = (
                newstar.total_observations_requested * newstar.exptime
                + slew_overhead * newstar.total_observations_requested
            )
            newstar.total_requested_hours = newstar.total_requested_seconds / 3600
            newstar.total_requested_nights = (
                newstar.total_requested_hours / hours_per_night
            )

        # Handle division by zero for inactive stars (total_observations_requested = 0)
        if newstar.total_observations_requested > 0:
            newstar.cume_observe_pct = np.round(
                (
                    np.cumsum(newstar.dates_observe)
                    / newstar.total_observations_requested
                )
                * 100.0,
                3,
            )
        else:
            # For inactive stars, show percentage based on total past observations if any exist
            total_past_obs = (
                sum(newstar.observations_past.values())
                if newstar.observations_past
                else 0
            )
            if total_past_obs > 0:
                newstar.cume_observe_pct = np.round(
                    (np.cumsum(newstar.dates_observe) / total_past_obs) * 100.0, 3
                )
            else:
                newstar.cume_observe_pct = np.zeros(
                    len(semester_planner.all_dates_array)
                )

        # Create consistent colors across programs, and random colors for each star within programs
        newstar.program_color_rgb = program_colors_rgb_vals[newstar.program]
        # Ensure rgb_strings has at least one element before random selection
        if len(rgb_strings) > 1:
            newstar.star_color_rgb = rgb_strings[
                np.random.randint(0, len(rgb_strings) - 1)
            ]
        else:
            newstar.star_color_rgb = rgb_strings[0]
        newstar.draw_lines = False
        newstar.maps_names = [
            "is_allocated",
            "is_custom",
            "is_altaz",
            "is_moon",
            "is_inter",
            "is_future",
            "is_clear",
            "is_observable_now",
        ]
        # Find the target index for this star in the access record
        # For inactive targets, they won't be in requests_frame, so create zero maps
        try:
            target_idx = np.where(
                semester_planner.requests_frame["unique_id"] == newstar.unique_id
            )[0][0]
            # Extract the 2D slice for this specific target from each 3D map
            newstar.maps = {
                name: access[name][target_idx] for name in newstar.maps_names
            }
            newstar.allow_mapview = True
        except (IndexError, KeyError):
            # Target is inactive (not in access record) - create zero maps with appropriate shape
            n_nights = semester_planner.semester_length
            n_slots = int((24 * 60) / semester_planner.config.getint("semester", "slot_size"))
            newstar.maps = {
                name: np.zeros((n_nights, n_slots), dtype=bool)
                for name in newstar.maps_names
            }
            newstar.allow_mapview = False

        all_stars.append(newstar)
        i += 1

    # Now create StarPlotter objects for each program, as it were one star.
    # These will not have all the attributes, but we only need these for the admin COF plot
    # These StarPlotter objects cannot be used to create a birdseye plot, they don't have all attributes
    programmatics = pd.read_csv(
        os.path.join(semester_planner.config.get("global", "workdir"), "programs.csv")
    )

    unique_programs = sorted(set(star.program for star in all_stars))
    programs_as_stars = {}
    for i in range(len(unique_programs)):
        prog_indices = [
            j for j, star in enumerate(all_stars) if star.program == unique_programs[i]
        ]
        prog_objs = [
            star
            for j, star in enumerate(all_stars)
            if star.program == unique_programs[i]
        ]

        # This is the quasi-StarPlotter object definition
        programmatic_star = StarPlotter(all_stars[prog_indices[0]].program)
        programmatic_star.target = all_stars[prog_indices[0]].program
        programmatic_star.program = all_stars[prog_indices[0]].program

        # Compute the COF data for all stars in the given program
        cume_observe = [all_stars[k].cume_observe for k in prog_indices]
        programmatic_star.cume_observe = np.sum(
            [all_stars[k].cume_observe for k in prog_indices], axis=0
        )
        stars_stacked = np.vstack(cume_observe)
        summed_cumulative = np.sum(stars_stacked, axis=0)
        max_value = np.sum(
            [all_stars[k].total_observations_requested for k in prog_indices]
        )
        programmatic_star.cume_observe_pct = np.round(
            summed_cumulative / max_value * 100, 2
        )

        # Compute the cumulative observe time for all stars in the given program
        cume_observe_time = [all_stars[k].cume_observe_time for k in prog_indices]
        stars_stacked_time = np.vstack(cume_observe_time)
        summed_cumulative_time = np.sum(stars_stacked_time, axis=0)
        total_requested_prog = np.sum(
            [all_stars[k].total_requested_hours for k in prog_indices]
        )
        allocated = programmatics[programmatics["program"] == unique_programs[i]][
            "hours"
        ].sum()
        # Use requested as divisor when requested < allocated, else allocated
        max_value_time = min(total_requested_prog, allocated)
        # summed_cumulative_time and max_value_time are both in hours
        if max_value_time > 0:
            programmatic_star.cume_observe_time_pct = np.round(
                summed_cumulative_time / max_value_time * 100, 2
            )
        else:
            programmatic_star.cume_observe_time_pct = np.zeros(
                len(semester_planner.all_dates_array)
            )
        programmatic_star.cume_observe_time = summed_cumulative_time  # in hours

        # Handle division by zero for programs with only inactive stars
        if max_value > 0:
            programmatic_star.cume_observe_pct = np.round(
                summed_cumulative / max_value * 100, 2
            )
        else:
            # For inactive-only programs, use total past observations as denominator
            total_past_obs = sum(
                sum(all_stars[k].observations_past.values())
                if all_stars[k].observations_past
                else 0
                for k in prog_indices
            )
            if total_past_obs > 0:
                programmatic_star.cume_observe_pct = (
                    summed_cumulative / total_past_obs * 100
                )
            else:
                programmatic_star.cume_observe_pct = np.zeros(
                    len(semester_planner.all_dates_array)
                )

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

        programmatic_star.total_observations_requested = np.sum(
            [all_stars[k].total_observations_requested for k in prog_indices]
        )
        programmatic_star.total_requested_hours = np.sum(
            [all_stars[k].total_requested_hours for k in prog_indices]
        )
        programmatic_star.draw_lines = False
        programmatic_star.allow_mapview = False

        # Set colors to match program color
        programmatic_star.program_color_rgb = all_stars[
            prog_indices[0]
        ].program_color_rgb
        programmatic_star.star_color_rgb = all_stars[prog_indices[0]].program_color_rgb

        # Create list of "stars" objects which are really the programmatic overview
        programs_as_stars[all_stars[prog_indices[0]].program] = programmatic_star

    # Group stars into lists by program indexed by a dictionary
    program_dict = defaultdict(list)
    for obj in all_stars:
        program_dict[obj.program].append(obj)

    return program_dict, programs_as_stars, nulltime


def get_cof(semester_planner, all_stars, use_time=False):
    """
    Produce a plotly figure showing the Cumulative Observability Function (COF) for a selection of stars

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        all_stars (array): a array of StarPlotter objects
        use_time (bool): if True, use the cumulative observe time percentage instead of the cumulative observe percentage

    Returns:
        fig (plotly figure): a plotly figure showing the COF for a selection of stars
    """

    fig = go.Figure()
    fig.update_layout(
        plot_bgcolor=gray, paper_bgcolor=clear
    )  # autosize=True,margin=dict(l=40, r=40, t=40, b=40),

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
            total_past_obs = sum(
                sum(star.observations_past.values()) if star.observations_past else 0
                for star in all_stars
            )
            if total_past_obs > 0:
                cume_observe_pct = (cume_observe / total_past_obs) * 100
            else:
                cume_observe_pct = np.zeros(len(semester_planner.all_dates_array))

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
                + "<br># Obs Requested: "
                + str(max_value)
                + "<br>",
                customdata=semester_planner.all_dates_array,
            )
        )
    else:
        # use_time=True: normalize by program hours from programs.csv
        programmatics_cof = pd.read_csv(
            os.path.join(semester_planner.config.get("global", "workdir"), "programs.csv")
        )
        programs_in_stars = set(
            getattr(s, "program", getattr(s, "target", None)) for s in all_stars
        )
        programs_in_stars = {p for p in programs_in_stars if p is not None}
        summed_cume_time = np.sum(
            [
                getattr(
                    s,
                    "cume_observe_time",
                    np.zeros(len(semester_planner.all_dates_array)),
                )
                for s in all_stars
            ],
            axis=0,
        )
        total_program_hours = programmatics_cof[
            programmatics_cof["program"].isin(programs_in_stars)
        ]["hours"].sum()
        # summed_cume_time and total_program_hours are both in hours
        if total_program_hours > 0:
            cume_time_pct = np.round(summed_cume_time / total_program_hours * 100, 2)
        else:
            cume_time_pct = np.zeros(len(semester_planner.all_dates_array))

        # Add the Total trace (time-based)
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
                + "<br>Time % Complete: %{y}"
                + "<br>Total program time: "
                + f"{total_program_hours:.1f} hours<br>"
                + "<extra></extra>",
                customdata=semester_planner.all_dates_array,
            )
        )

    # Then add individual star traces (so they appear above the Total trace)
    for i in range(len(all_stars)):
        if use_time:
            y_vals = getattr(all_stars[i], "cume_observe_time_pct", None)
            prog_for_star = getattr(all_stars[i], "program", all_stars[i].target)
            total_prog_hours = (
                programmatics_cof.loc[
                    programmatics_cof["program"] == prog_for_star, "hours"
                ].iloc[0]
                if prog_for_star in programmatics_cof["program"].values
                else 0.0
            )
            if y_vals is None:
                # Individual stars: compute from cume_observe_time (hours) / program hours
                y_vals = (
                    np.round(all_stars[i].cume_observe_time / total_prog_hours * 100, 2)
                    if total_prog_hours > 0
                    else np.zeros(len(semester_planner.all_dates_array))
                )
            hovertemplate = (
                "<b>"
                + str(prog_for_star)
                + "</b><br>Night: %{x}"
                + "<br>Date: "
                + "%{customdata}"
                + "<br>Time % Complete: %{y}<br>Total program time: "
                + f"{total_prog_hours:.1f} hours<br>"
                + "<extra></extra>"
            )
        else:
            y_vals = all_stars[i].cume_observe_pct
            hovertemplate = (
                "Night: %{x}"
                + "<br>Date: "
                + "%{customdata}"
                + "<br>% Complete: %{y}"
                + "<br># Obs Requested: "
                + str(all_stars[i].total_observations_requested)
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
                customdata=semester_planner.all_dates_array,
            )
        )
        last_pct = float(np.round(y_vals[-1], 2)) if len(y_vals) else 0
        lines.append(str(all_stars[i].target) + "," + str(last_pct))

    # Find the night index for "today" (current_day)
    try:
        today_night_index = semester_planner.all_dates_array.index(
            semester_planner.config.get("global", "current_day")
        )
    except (ValueError, AttributeError):
        # Fallback to today_starting_night if available, otherwise use 0
        today_night_index = getattr(semester_planner, "today_starting_night", 0) - 1

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
    x_tickvals = list(range(0, semester_planner.semester_length, x_tick_step))
    if (semester_planner.semester_length - 1) not in x_tickvals:
        x_tickvals.append(semester_planner.semester_length - 1)
    x_ticktext = [
        str(val + 1) for val in x_tickvals
    ]  # Night indices (1-indexed for display, matching birdseye)

    # Create calendar date labels for secondary x-axis (top axis)
    # Format dates as "Feb<br>01" (month and day on separate lines)
    x_ticktext_dates = []
    for day_idx in x_tickvals:
        if day_idx < len(semester_planner.all_dates_array):
            date_str = semester_planner.all_dates_array[day_idx]
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
        "Time % Complete (vs program hours)" if use_time else "Request % Complete"
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
            range=[0, semester_planner.semester_length - 1],  # Explicitly set range
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
                semester_planner.semester_length - 1,
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
            x=[0, len(semester_planner.all_dates_array) - 1],
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


def get_birdseye(semester_planner, availablity, all_stars):
    """
    Produce the plotly figure showing the day/slot matrix intersection for a selection of stars

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        availability (array): a 2D array of N_slots by N_nights, binary 1/0, it is the intersection of is_allocated and is_night
        all_stars (array): a array of StarPlotter objects

    Returns:
        fig (plotly figure): a plotly figure showing the day/slot matrix intersection for a selection of stars
    """

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
    fig.add_vrect(
        x0=semester_planner.today_starting_night
        - 1,  # The minus one is just for aesthetic purposes.
        x1=semester_planner.today_starting_night - 1,
        annotation_text="Today",
        line_dash="dash",
        fillcolor=None,
        line_width=2,
        line_color="black",
        annotation_position="bottom left",
    )
    # X-axis: ticks every 23 days, plus the last day
    x_tick_step = 23
    x_tickvals = list(range(0, semester_planner.semester_length, x_tick_step))
    if (semester_planner.semester_length - 1) not in x_tickvals:
        x_tickvals.append(semester_planner.semester_length - 1)
    x_ticktext = [str(val + 1) for val in x_tickvals]

    # Create calendar date labels for secondary x-axis (top axis)
    # Format dates as "Jan<br>15" or "Aug<br>12" (month and day on separate lines)
    x_ticktext_dates = []
    for day_idx in x_tickvals:
        if day_idx < len(semester_planner.all_dates_array):
            date_str = semester_planner.all_dates_array[day_idx]
            # Parse date and format as "Jan<br>15" or "Aug<br>12" using HTML break tag
            date_obj = datetime.strptime(date_str, "%Y-%m-%d")
            month = date_obj.strftime("%b")
            day = date_obj.strftime("%d")
            x_ticktext_dates.append(f"{month}<br>{day}")
        else:
            x_ticktext_dates.append("")

    # Y-axis: ticks every 2 hours, using slot_size
    n_slots = int(24 * 60 // semester_planner.config.getint("semester", "slot_size"))
    slots_per_2hr = int(2 * 60 // semester_planner.config.getint("semester", "slot_size"))
    y_tickvals = list(range(0, n_slots, slots_per_2hr))
    y_ticktext = []
    for slot in y_tickvals:
        total_minutes = slot * semester_planner.config.getint("semester", "slot_size")
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
    n_slots = int(24 * 60 // semester_planner.config.getint("semester", "slot_size"))
    fig.add_trace(
        go.Scatter(
            x=[0, len(semester_planner.all_dates_array) - 1],
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
            range=[0, semester_planner.semester_length - 1],  # Explicitly set range
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
                semester_planner.semester_length - 1,
            ],  # Match primary x-axis range
        ),
        margin=dict(
            b=200, t=100
        ),  # Bottom margin for legend below, top margin for date labels
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


def get_rawobs(semester_planner, all_stars, use_program_colors=False):
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
    semester_planner, all_stars, use_program_colors=False, prevent_negative=False
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
    programmatics = pd.read_csv(
        os.path.join(semester_planner.config.get("global", "workdir"), "programs.csv")
    )

    # Per-visit overhead scalars come from the queue (single source of truth).
    slew_overhead = semester_planner.queue.slew_overhead_mean
    readout_overhead = semester_planner.queue.readout_time

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
            total_past += visits * (
                starobj.exptime * starobj.n_exp
                + readout_overhead * (starobj.n_exp - 1)
                + slew_overhead
            )
        # Future: same day-by-day formula (a) exposures*visits, (b) readout*(n_exp-1)*visits, (c) slew*visits
        for visits in starobj.observations_future.values():
            total_future += visits * (
                starobj.exptime * starobj.n_exp
                + readout_overhead * (starobj.n_exp - 1)
                + slew_overhead
            )
        total_requested_hours += starobj.total_requested_hours
        programs_used.append(starobj.program)

    # Convert to hours for better readability
    total_past_hours = total_past / 3600
    total_future_hours = total_future / 3600
    total_incomplete_hours = (
        total_requested_hours - total_past_hours - total_future_hours
    )

    if len(programs_used) > 1:
        program_rows = programmatics[programmatics["program"].isin(programs_used)]
        total_allocated_hours = program_rows["hours"].sum()
        total_allocated_nights = program_rows["nights"].sum()
    else:
        program_rows = programmatics[programmatics["program"] == programs_used[0]]
        total_allocated_hours = program_rows["hours"].sum()
        total_allocated_nights = program_rows["nights"].sum()

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
        title_text=f"<b>Total Requested:</b> {total_requested_hours:.1f} hours ≈ {total_requested_hours / hours_per_night:.1f} nights<br><b>Total Allocated:</b> {total_allocated_hours:.1f} hours = {total_allocated_nights:.1f} nights ----> w/ losses = {total_allocated_nights * 0.75:.1f} nights <br>Requested time is measured in hours. Allocated time is measured in nights. Conversion is 12 hours per night.<br>All bars include exposure times and standard overheads.",
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

    # Add gray vertical dashed line at total_allocated_hours * throttle_grace
    grace_factor = semester_planner.config.getfloat("semester", "throttle_grace")
    fig.add_shape(
        type="line",
        x0=total_allocated_hours * grace_factor,
        x1=total_allocated_hours * grace_factor,
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

    # Add invisible scatter trace for hover text on the throttle grace line
    grace_value = total_allocated_hours * grace_factor
    fig.add_trace(
        go.Scatter(
            x=[grace_value] * len(labels),
            y=labels,  # Use categorical labels instead of numeric positions
            mode="markers",
            marker=dict(size=20, opacity=0),  # Invisible but hoverable markers
            hovertemplate=f"<b>Maximum Schedulable Time</b><br>{grace_value:.2f} hours<br>We allow for over-filled requests by a factor of up to {grace_factor:.2f} your allocation<br>Algorithmically, you are forbidden from getting more time than this.<extra></extra>",
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
    programmatics = pd.read_csv(
        os.path.join(semester_planner.config.get("global", "workdir"), "programs.csv")
    )

    # Per-visit overhead scalars come from the queue (single source of truth).
    slew_overhead = semester_planner.queue.slew_overhead_mean
    readout_overhead = semester_planner.queue.readout_time

    # Get all programs from programs.csv
    all_programs_in_csv = set(programmatics["program"].unique())
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

        # Calculate times for this program (same logic as get_timebar)
        total_past = 0
        total_future = 0
        total_requested_hours = 0

        for starobj in program_stars:
            # Past: day-by-day sum of (exposure) + (readout) + (slew) per visit; per date: visits * (exptime*n_exp + readout*(n_exp-1) + slew)
            for visits in starobj.observations_past.values():
                total_past += visits * (
                    starobj.exptime * starobj.n_exp
                    + readout_overhead * (starobj.n_exp - 1)
                    + slew_overhead
                )
            # Future: same day-by-day formula
            for visits in starobj.observations_future.values():
                total_future += visits * (
                    starobj.exptime * starobj.n_exp
                    + readout_overhead * (starobj.n_exp - 1)
                    + slew_overhead
                )
            total_requested_hours += starobj.total_requested_hours

        # Convert to hours
        total_past_hours = total_past / 3600
        total_future_hours = total_future / 3600
        total_incomplete_hours = (
            total_requested_hours - total_past_hours - total_future_hours
        )

        # Get allocated hours for this program
        program_row = programmatics[programmatics["program"] == program_code]
        if len(program_row) > 0:
            total_allocated_hours = program_row["hours"].sum()
        else:
            total_allocated_hours = 0

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
        # Get allocated hours for this program from programs.csv
        program_row = programmatics[programmatics["program"] == program_code]
        if len(program_row) > 0:
            total_allocated_hours = program_row["hours"].sum()
        else:
            total_allocated_hours = 0

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
            yref=yref,
        )

        # Add gray vertical dashed line at allocated * throttle_grace
        grace_factor = semester_planner.config.getfloat("semester", "throttle_grace")
        fig.add_shape(
            type="line",
            x0=allocated * grace_factor,
            x1=allocated * grace_factor,
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

        # Add invisible scatter for hover on throttle grace line
        grace_value = allocated * grace_factor
        fig.add_trace(
            go.Scatter(
                x=[grace_value],
                y=[category_names[2]],  # Middle bar (Future Scheduled)
                mode="markers",
                marker=dict(size=15, opacity=0),
                hovertemplate=f"<b>{program_code} Throttle Grace</b><br>{grace_value:.2f} hours<br>Allocated time times throttle grace factor ({grace_factor:.2f})<extra></extra>",
                hoverlabel=dict(bgcolor="gray", font_color="white"),
                showlegend=False,
            ),
            row=row,
            col=col,
        )

        # Update x-axis for this subplot (scaled to this program's data)
        # Include allocated*grace and weather loss so the gray lines are visible when they exceed the bars
        weather_loss_value = allocated - allocated * weather_loss_factor
        program_max = max(
            data["unused"],
            data["incomplete"],
            data["future"],
            data["past"],
            data["requested"],
            data["allocated"],
            allocated * grace_factor,
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


def get_football(semester_planner, all_stars, use_program_colors=False):
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
        semester_planner.config.get("global", "semester_start_day")[:4] + semester_planner.config.get("global", "semester")[-1]
    )
    cache_dir = _football_cache_dir(semester_planner)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_grids_file = str(cache_dir / f"{semester}_sky_grids.npz")
    cache_image_file = str(cache_dir / f"{semester}_sky_availability_image.txt")
    semester_length = semester_planner.semester_length

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
            queue=semester_planner.queue,
            request_frame=seasonality_frame,
            semester_start_date=semester_planner.config.get("global", "semester_start_day"),
            semester_length=semester_length,
            slot_size=semester_planner.config.getint("semester", "slot_size"),
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


def get_request_frame(semester_planner, all_stars):
    """
    Get a filtered request frame containing only the stars in all_stars.

    Args:
        semester_planner: the semester planner object
        all_stars (list): array of StarPlotter objects

    Returns:
        filtered_frame (pd.DataFrame): filtered request frame with only the specified stars
    """
    # Extract targets from the StarPlotter objects
    starids = [star.unique_id for star in all_stars]

    # Filter the request frame to only include the specified stars
    filtered_frame = semester_planner.requests_frame_all[
        semester_planner.requests_frame_all["unique_id"].isin(starids)
    ].copy()

    return filtered_frame


def get_ladder(data, tonight_start_time):
    """Produce a plotly figure which illustrates the night plan solution.

    Args:
        data (obj): a TTP data object containing the schedule information

    Returns:
        fig (plotly figure): a plotly figure illustrating the night plan solution.
    """

    from astroq.ttp.plot import schedule_to_ladder_frame, _as_model

    model = _as_model(data)
    orderData = schedule_to_ladder_frame(model)
    if orderData.empty:
        orderData = pd.DataFrame(
            columns=[
                "unique_id",
                "human_target",
                "First Available",
                "Last Available",
                "Start Exposure",
                "Stop Exposure",
                "Total Exp Time (min)",
                "Slew to Next (min)",
                "Minutes the from Start of the Night",
            ]
        )
    if "Slew to Next (min)" not in orderData.columns:
        orderData["Slew to Next (min)"] = 0.0

    if model.night_start is not None and len(orderData):
        orderData["UTC Start Time"] = [
            (model.night_start + TimeDelta(se * 60, format="sec")).isot[11:16]
            if se > 0
            else ""
            for se in orderData["Start Exposure"]
        ]

    on_sky = model.schedule[~model.schedule["is_anchor"]]
    n_unscheduled = int((~on_sky["scheduled"]).sum())

    # reverse so the plot flows top -> bottom with time; after reversal,
    # the lowest indices (bottom of plot) hold the unscheduled block.
    orderData = orderData.iloc[::-1].reset_index(drop=True)

    # Each priority gets a different color. Keys are integer tiers (1–10); values may be float.
    colordict = {
        "10": "red",
        "9": "tomato",
        "8": "darkorange",
        "9": "sandybrown",
        "7": "gold",
        "6": "olive",
        "5": "green",
        "4": "cyan",
        "3": "darkviolet",
        "2": "magenta",
        "1": "blue",
    }

    def _priority_color(priority):
        key = str(int(round(float(priority))))
        return colordict.get(key, "gray")

    hover_cols = [
        "First Available",
        "Last Available",
        "Exposure Time (min)",
        "N_shots",
        "Total Exp Time (min)",
        "Slew to Next (min)",
        "UTC Start Time",
    ]
    fig = px.scatter(
        orderData,
        x="Minutes the from Start of the Night",
        y="human_target",
        hover_data=hover_cols,
        title="Night Plan",
        width=800,
        height=1000,
    )  # color='Program'
    fig.update_layout(yaxis_title="")
    fig.add_shape(
        type="rect",
        x0=-100,
        x1=-80,
        y0=-0.5,
        y1=0.5,
        fillcolor="red",
        showlegend=True,
        name="Exposure",
    )
    fig.add_shape(
        type="rect",
        x0=-100,
        x1=-80,
        y0=-0.5,
        y1=0.5,
        fillcolor="dimgray",
        showlegend=True,
        name="Slew",
    )
    fig.add_shape(
        type="rect",
        x0=-100,
        x1=-80,
        y0=-0.5,
        y1=0.5,
        fillcolor="lime",
        opacity=0.3,
        showlegend=True,
        name="Accessible",
    )

    new_already_processed = []
    ifixer = 0  # for multi-visit targets, it throws off the one row per target plotting...this fixes it
    for i in range(len(orderData["unique_id"])):
        if orderData["unique_id"][i] not in new_already_processed:
            indices = [
                k
                for k in range(len(orderData["unique_id"]))
                if orderData["unique_id"][k] == orderData["unique_id"][i]
            ]
            for j in range(len(indices)):
                if j == 0:
                    # only do this once, otherwise the green bar gets discolored compared to other rows
                    fig.add_shape(
                        type="rect",
                        x0=orderData["First Available"][indices[j]],
                        x1=orderData["Last Available"][indices[j]],
                        y0=i + ifixer - 0.5,
                        y1=i + ifixer + 0.5,
                        fillcolor="lime",
                        opacity=0.3,
                        showlegend=False,
                    )
                fig.add_shape(
                    type="rect",
                    x0=orderData["Start Exposure"][indices[j]],
                    x1=orderData["Start Exposure"][indices[j]]
                    + orderData["Total Exp Time (min)"][indices[j]],
                    y0=i + ifixer - 0.5,
                    y1=i + ifixer + 0.5,
                    fillcolor=_priority_color(orderData["Priority"][indices[j]]),
                )
                slew = float(orderData["Slew to Next (min)"][indices[j]])
                if slew > 0:
                    fig.add_shape(
                        type="rect",
                        x0=orderData["Stop Exposure"][indices[j]],
                        x1=orderData["Stop Exposure"][indices[j]] + slew,
                        y0=i + ifixer - 0.5,
                        y1=i + ifixer + 0.5,
                        fillcolor="dimgray",
                        line=dict(width=0),
                    )
            new_already_processed.append(orderData["unique_id"][i])
        else:
            # if we already did this star, it is a multi-visit star and we need to adjust the row counter for plotting purposes
            ifixer -= 1

    if n_unscheduled and n_unscheduled < len(orderData):
        sep_y = n_unscheduled - 0.5
        fig.add_hline(y=sep_y, line_color="black", line_width=1, line_dash="solid")

    x_min = 0
    night_start = getattr(model, "night_start", None)
    night_end = getattr(model, "night_end", None)
    if night_start is not None and night_end is not None:
        x_max = (night_end.jd - night_start.jd) * 24 * 60
    elif len(orderData) > 0:
        end_times = (
            orderData["Start Exposure"]
            + orderData["Total Exp Time (min)"]
            + orderData["Slew to Next (min)"]
        )
        x_max = end_times.max()
    else:
        x_max = 600
    fig.update_layout(xaxis_range=[x_min, x_max])
    for x_line, label in [(x_min, "start"), (x_max, "end")]:
        fig.add_vline(
            x=x_line,
            line_color="black",
            line_width=1,
            annotation_text=label,
            annotation_position="top",
        )
    # Add secondary x-axis with UTC time
    start_time = tonight_start_time.to_datetime()
    # Create tick positions (every 60 minutes or so, adjust as needed)
    tick_interval = 60  # minutes
    tick_positions = list(range(0, int(x_max) + tick_interval, tick_interval))
    tick_labels = [
        (start_time + timedelta(minutes=pos)).strftime("%H:%M")
        for pos in tick_positions
    ]
    # Add secondary x-axis
    # Add an invisible trace to force the secondary axis to appear
    fig.add_trace(
        go.Scatter(
            x=[x_min, x_max],
            y=["unique_id", "unique_id"],  # Place just below the visible range
            mode="markers",
            marker=dict(size=0.1, opacity=0),
            showlegend=False,
            hoverinfo="skip",
            xaxis="x2",
        )
    )
    # Create the secondary x-axis configuration
    fig.update_layout(
        xaxis2=dict(
            title=dict(text="UTC Time", standoff=0),
            overlaying="x",
            side="top",
            range=[x_min, x_max],
            tickmode="array",
            tickvals=tick_positions,
            ticktext=tick_labels,
            showgrid=False,
            showline=True,
            mirror=True,
        )
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
    request_selected_path = os.path.join(
        night_planner.output_directory, "request_selected.csv"
    )

    if not os.path.exists(request_selected_path):
        raise FileNotFoundError(
            f"request_selected.csv not found at {request_selected_path}"
        )

    # Read the request_selected.csv file
    request_selected_df = pd.read_csv(request_selected_path)
    solution = night_planner.solution
    on_sky = solution.schedule[~solution.schedule["is_anchor"]]
    scheduled = on_sky[on_sky["scheduled"]].sort_values("order")

    merged_df = request_selected_df.merge(
        scheduled[["unique_id", "t_start", "t_early", "t_late"]],
        on="unique_id",
        how="inner",
    )
    merged_df = merged_df.rename(
        columns={
            "t_start": "Start Exposure",
            "t_early": "First Available",
            "t_late": "Last Available",
        }
    )

    # Select and reorder only the specific columns requested
    # desired_columns = [
    #     'Start Exposure', 'unique_id', 'target', 'program_code', 'ra', 'dec',
    #     'exptime', 'n_exp', 'n_intra_max', 'tau_intra', 'weather_band_1', 'weather_band_2', 'weather_band_3', 'teff',
    #     'jmag', 'Vmag', 'epoch', 'gaia_id', 'First Available', 'Last Available'
    # ]
    desired_columns = [
        "First Available",
        "Start Exposure",
        "Last Available",
        "unique_id",
        "target",
        "program_code",
        "ra",
        "dec",
        "exptime",
        "n_exp",
        "n_intra_max",
        "tau_intra",
        "jmag",
        "Vmag",
    ]

    # Keep only the columns that exist in the merged dataframe
    available_columns = [col for col in desired_columns if col in merged_df.columns]

    # Reorder columns to match the desired structure
    final_df = merged_df[available_columns].copy()

    # Round numeric fields to appropriate decimal places
    if "ra" in final_df.columns:
        # Ensure ra is numeric before rounding, handle 'None' strings
        final_df["ra"] = final_df["ra"].replace("None", pd.NA)
        final_df["ra"] = pd.to_numeric(final_df["ra"], errors="coerce").round(1)

    if "dec" in final_df.columns:
        # Ensure dec is numeric before rounding, handle 'None' strings
        final_df["dec"] = final_df["dec"].replace("None", pd.NA)
        final_df["dec"] = pd.to_numeric(final_df["dec"], errors="coerce").round(1)

    if "jmag" in final_df.columns:
        # Ensure jmag is numeric before rounding, handle 'None' strings
        final_df["jmag"] = final_df["jmag"].replace("None", pd.NA)
        final_df["jmag"] = pd.to_numeric(final_df["jmag"], errors="coerce").round(1)

    if "Vmag" in final_df.columns:
        final_df["Vmag"] = final_df["Vmag"].replace("None", pd.NA)
        final_df["Vmag"] = pd.to_numeric(final_df["Vmag"], errors="coerce").round(1)

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
            night_planner.allocation_file, night_planner.current_day
        )

        # Convert the time columns to HST timestamps
        if "Start Exposure" in final_df.columns:
            final_df["Start Exposure"] = final_df["Start Exposure"].apply(
                lambda x: (
                    str(TimeDelta(x * 60, format="sec") + night_start_time)[11:16]
                    if pd.notna(x)
                    else ""
                )
            )

        if "First Available" in final_df.columns:
            final_df["First Available"] = final_df["First Available"].apply(
                lambda x: (
                    str(TimeDelta(x * 60, format="sec") + night_start_time)[11:16]
                    if pd.notna(x)
                    else ""
                )
            )

        if "Last Available" in final_df.columns:
            final_df["Last Available"] = final_df["Last Available"].apply(
                lambda x: (
                    str(TimeDelta(x * 60, format="sec") + night_start_time)[11:16]
                    if pd.notna(x)
                    else ""
                )
            )

    except Exception as e:
        print(f"Warning: Could not convert time fields to HST timestamps: {e}")
        print("Time fields will remain as minutes from start of night")

    # Handle missing values and 'None' strings
    final_df = final_df.replace(["", "NoGaiaName", "None"], pd.NA)

    # Ensure DataFrame is clean and properly structured for DataTables
    final_df = final_df.reset_index(drop=True)
    # Remove duplicate column names if any exist
    final_df = final_df.loc[:, ~final_df.columns.duplicated(keep="first")]
    # Fill NaN values with empty strings to ensure consistent structure
    final_df = final_df.fillna("")
    # Ensure all columns have consistent data types (convert objects to strings)
    for col in final_df.columns:
        if final_df[col].dtype == "object":
            final_df[col] = (
                final_df[col]
                .astype(str)
                .replace("nan", "")
                .replace("None", "")
                .replace("", "")
            )

    return final_df


REQUEST_FRAME_COLUMNS = [
    "target",
    "unique_id",
    "program_code",
    "ra",
    "dec",
    "exptime",
    "n_exp",
    "n_inter_max",
    "tau_inter",
    "n_intra_max",
    "n_intra_min",
    "tau_intra",
    "weather_band_1",
    "weather_band_2",
    "weather_band_3",
    "inactive",
    "comments",
]
BOOLEAN_COLUMNS = {
    "weather_band_1": "Band1",
    "weather_band_2": "Band2",
    "weather_band_3": "Band3",
    "inactive": "Inactive",
}
REQUEST_FRAME_DISPLAY_NAMES = {
    "target": "Target",
    "unique_id": "ID",
    "program_code": "Program",
    "ra": "RA",
    "dec": "Dec",
    "exptime": "ExpTime",
    "comments": "Comments",
}
# Tooltips shown when hovering over column headers.
REQUEST_FRAME_COLUMN_TOOLTIPS = {
    "Star": "Name of the star",
    "ID": "Keck OB database unique ID",
    "Program": "Program Code",
    "RA": "RA in decimal degrees",
    "Dec": "Declination in decimal degrees",
    "ExpTime": "Exposure time in seconds",
    "n_exp": "Number of Exposures per Visit",
    "n_inter_max": "Maximum number of unique nights to observe the star",
    "tau_inter": "The minimum inter-night cadence between unique night observations",
    "n_intra_max": "The desired number of visits to the star in each night it is observed",
    "n_intra_min": "The accepted minimum number of visits to the star in each night it is observed",
    "tau_intra": "The minimum intra-night cadence between visits within a night in hours",
    "Band1": "Allowed to observe in Band1?",
    "Band2": "Allowed to observe in Band2?",
    "Band3": "Allowed to observe in Band3?",
    "Inactive": "Is the star set to inactive?",
    "Comments": "Observer notes (e.g. from Keck star list)",
}


_REQUEST_BAND_COLS = ("Band1", "Band2", "Band3", "Inactive")
_REQUEST_NO_PAD_COLS = (
    "n_inter_max",
    "tau_inter",
    "n_intra_max",
    "n_intra_min",
    "tau_intra",
)
# Numeric columns (post-rename): RA(3), Dec(4), ExpTime(5), n_exp(6),
# n_inter_max(7), tau_inter(8), n_intra_max(9), n_intra_min(10), tau_intra(11).
_REQUEST_NUMERIC_COLS = [3, 4, 5, 6, 7, 8, 9, 10, 11]


def _is_true(val):
    """Coerce a CSV-ish value to bool with permissive parsing."""
    if pd.isna(val) or val == "":
        return False
    s = str(val).lower()
    if s in ("true", "1", "yes"):
        return True
    if s in ("false", "0", "no"):
        return False
    try:
        return bool(float(val))
    except (ValueError, TypeError):
        return False


def _yn_cell(val, *, inactive_semantics):
    """Render a boolean cell as Y/N with green=good / red=bad background."""
    green_bg = "rgba(34, 139, 34, 0.25)"
    red_bg = "rgba(220, 53, 69, 0.25)"
    truth = _is_true(val)
    if inactive_semantics:
        y_n, bg = ("Y", red_bg) if truth else ("N", green_bg)
    else:
        y_n, bg = ("Y", green_bg) if truth else ("N", red_bg)
    return (
        f'<span style="background:{bg};padding:2px 6px;border-radius:4px;">{y_n}</span>'
    )


def _visible_len(s):
    """Length of ``s`` with HTML tags stripped (used for content-fit column widths)."""
    return len(re.sub(r"<[^>]+>", "", str(s)).strip())


def request_frame_to_html(
    request_df,
    semester_code=None,
    date=None,
    band=None,
    table_id="request-table",
    page_size=25,
):
    """
    Convert a request frame (from request.csv) to HTML for admin/program/star pages.

    Displays only: target, unique_id, program_code, ra, dec, exptime, n_exp,
    n_inter_max, tau_inter, n_intra_max, n_intra_min, tau_intra, Band1, Band2, Band3, Inactive, Comments.
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
    df = request_df.copy().reset_index(drop=True)
    cols = [c for c in REQUEST_FRAME_COLUMNS if c in df.columns]
    df = df[cols].copy().fillna("")

    for coord in ("ra", "dec"):
        if coord in df.columns:
            df[coord] = pd.to_numeric(df[coord], errors="coerce")
            df[coord] = df[coord].apply(lambda x: f"{x:.2f}" if pd.notna(x) else "")

    if (
        semester_code
        and date
        and band
        and "program_code" in df.columns
        and "target" in df.columns
    ):
        df["target"] = df.apply(
            lambda row: (
                f'<a href="/{semester_code}/{date}/{band}/'
                f'{quote(str(row["program_code"]))}/{quote(str(row["target"]))}">'
                f'{row["target"]}</a>'
            ),
            axis=1,
        )

    for orig in BOOLEAN_COLUMNS:
        if orig not in df.columns:
            continue
        is_inactive = orig == "inactive"
        df[orig] = df[orig].apply(
            lambda v, _i=is_inactive: _yn_cell(v, inactive_semantics=_i)
        )

    df = df.rename(columns={**BOOLEAN_COLUMNS, **REQUEST_FRAME_DISPLAY_NAMES})
    for col in df.columns:
        if df[col].dtype == "object" and col not in BOOLEAN_COLUMNS.values():
            df[col] = df[col].astype(str).replace("nan", "").replace("None", "")

    widths = []
    for col in df.columns:
        content_max = max((_visible_len(c) for c in df[col]), default=0)
        header_len = len(str(col))
        pad = 3 if col == "n_exp" else (0 if col in _REQUEST_NO_PAD_COLS else 2)
        ch_width = max(content_max, header_len, 1) + pad
        if col in _REQUEST_BAND_COLS:
            ch_width = max(ch_width, 3)
        widths.append(f"{ch_width}ch")

    column_defs = [{"target": i, "width": w} for i, w in enumerate(widths)]
    tooltips = [REQUEST_FRAME_COLUMN_TOOLTIPS.get(col, "") for col in df.columns]

    return _render_datatable(
        df,
        template_name="request_table.html.j2",
        table_id=table_id,
        variant="compact",
        column_widths=widths,
        column_defs=column_defs,
        tooltips=tooltips,
        page_size=page_size,
        sort_column=0,
        numeric_cols=_REQUEST_NUMERIC_COLS,
        has_band_padding=True,
        has_column_filters=True,
        filter_placeholder="Filter... (use > < >= <= for numbers)",
        add_tfoot=True,
    )


NIGHTPLAN_COLUMNS = [
    "First Available",
    "Start Exposure",
    "Last Available",
    "unique_id",
    "target",
    "program_code",
    "ra",
    "dec",
    "exptime",
    "n_exp",
    "n_intra_max",
    "tau_intra",
    "jmag",
    "Vmag",
]
NIGHTPLAN_COLUMN_TOOLTIPS = {
    "First Available": "First available time to observe (HH:MM). Use > < >= <= with HH:MM to filter.",
    "Start Exposure": "Scheduled start time (HH:MM). Use > < >= <= with HH:MM to filter.",
    "Last Available": "Last available time to observe (HH:MM). Use > < >= <= with HH:MM to filter.",
    "unique_id": "Keck OB database unique ID",
    "target": "Name of the target",
    "program_code": "Program Code",
    "ra": "Right ascension in decimal degrees",
    "dec": "Declination in decimal degrees",
    "exptime": "Exposure time in seconds",
    "n_exp": "Number of exposures per visit",
    "n_intra_max": "Maximum intra-night visits",
    "tau_intra": "Minimum intra-night cadence in hours",
    "jmag": "J-band magnitude",
    "Vmag": "V-band magnitude",
}


# Numeric (post-select) cols: ra(6), dec(7), exptime(8), n_exp(9),
# n_intra_max(10), tau_intra(11), jmag(12), Vmag(13). Time cols: 0,1,2.
_NIGHTPLAN_NUMERIC_COLS = [6, 7, 8, 9, 10, 11, 12, 13]
_NIGHTPLAN_TIME_COLS = [0, 1, 2]


def nightplan_table_to_html(script_df, table_id="script-table", page_size=100):
    """
    Convert nightplan script DataFrame to HTML with same styling as request_frame_to_html.

    Same colors, fonts, fontsize, filtering (partial match, numeric > < >= <=), hover tooltips.
    Displays: First Available, Start Exposure, Last Available, unique_id, target, program_code,
    ra, dec, exptime, n_exp, n_intra_max, tau_intra, jmag, Vmag.
    """
    df = script_df.copy().reset_index(drop=True)
    cols = [c for c in NIGHTPLAN_COLUMNS if c in df.columns]
    df = df[cols].copy().fillna("")
    for col in df.columns:
        if df[col].dtype == "object":
            df[col] = df[col].astype(str).replace("nan", "").replace("None", "")

    widths = []
    for col in df.columns:
        content_max = max((_visible_len(c) for c in df[col]), default=0)
        header_len = len(str(col))
        widths.append(f"{max(content_max, header_len, 1) + 2}ch")

    column_defs = [{"target": i, "width": w} for i, w in enumerate(widths)]
    tooltips = [NIGHTPLAN_COLUMN_TOOLTIPS.get(col, "") for col in df.columns]

    return _render_datatable(
        df,
        template_name="nightplan_table.html.j2",
        table_id=table_id,
        variant="compact",
        column_widths=widths,
        column_defs=column_defs,
        tooltips=tooltips,
        page_size=page_size,
        sort_column=1,
        numeric_cols=_NIGHTPLAN_NUMERIC_COLS,
        time_cols=_NIGHTPLAN_TIME_COLS,
        has_column_filters=True,
        filter_placeholder="Filter... (> < for HH:MM or numbers)",
        add_tfoot=True,
    )


# Default per-column widths for `dataframe_to_html` (legacy generic table).
_GENERIC_WIDTH_MAP = {
    "First Available": "80px",
    "Start Exposure": "80px",
    "Last Available": "80px",
    "unique_id": "200px",
    "target": "200px",
    "program_code": "120px",
    "ra": "100px",
    "dec": "100px",
    "exptime": "80px",
    "n_exp": "60px",
    "n_intra_max": "80px",
    "tau_intra": "80px",
    "jmag": "60px",
    "Vmag": "60px",
}


def dataframe_to_html(dataframe, sort_column=2, page_size=10, table_id="request-table"):
    """
    Convert a pandas dataframe into an HTML string for rendering
    on the webapp pages.

    Args:
        dataframe (pd.DataFrame): The dataframe to convert
        sort_column (int): Column index to sort by (default: 2 for target)
        page_size (int): Default number of rows per page (default: 25)
        table_id (str): Unique ID for the table (default: 'request-table')

    Returns:
        table_html (str): HTML string with table and DataTables initialization
    """
    df = dataframe.reset_index(drop=True)
    df = df.loc[:, ~df.columns.duplicated(keep="first")]
    df = df.fillna("")
    for col in df.columns:
        if df[col].dtype == "object":
            df[col] = df[col].astype(str).replace("nan", "").replace("None", "")

    if sort_column >= len(df.columns):
        sort_column = 0

    widths = [_GENERIC_WIDTH_MAP.get(col, "100px") for col in df.columns]
    column_defs = [{"target": i, "width": w} for i, w in enumerate(widths)]

    return _render_datatable(
        df,
        template_name="generic_table.html.j2",
        table_id=table_id,
        variant="card",
        column_widths=widths,
        column_defs=column_defs,
        page_size=page_size,
        sort_column=sort_column,
        has_dt_controls_styling=True,
        has_init_complete_header_style=True,
        table_layout="auto",
        scroll_x=False,
        responsive=False,
    )
