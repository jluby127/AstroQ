"""HTML table renderers for webapp and CLI."""

import re
from urllib.parse import quote

import pandas as pd

from astroq.plot._common import _render_datatable

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
    n_inter_max, tau_inter, n_intra_max, n_intra_min, tau_intra, Band1, Band2,
    Band3, Inactive, Comments.
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
    "Earliest Start",
    "Start Exposure",
    "Latest Finish",
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
    "Earliest Start": "Earliest allowed start time (HH:MM). Use > < >= <= with HH:MM to filter.",
    "Start Exposure": "Scheduled start time (HH:MM). Sorted from local noon→next noon. Use > < >= <= with HH:MM to filter.",
    "Latest Finish": "Latest allowed finish time (HH:MM). Use > < >= <= with HH:MM to filter.",
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
    Displays: Earliest Start, Start Exposure, Latest Finish, unique_id, target, program_code,
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
    "Earliest Start": "80px",
    "Start Exposure": "80px",
    "Latest Finish": "80px",
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
