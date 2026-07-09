"""HIRES-CPS data ingestion.

Pulls request tabs from per-program Google Sheets, the Keck observing schedule,
and JUMP past-history. Consumed by :func:`astroq.driver.hirescps_prep` to
produce the on-disk ``request.csv``, ``custom.csv``, ``allocation.csv``, and
``past.csv`` used by the planner.
"""

# Standard library imports
import csv
import hashlib
import io
import logging
import os
import math
import re
import urllib.parse

# Third-party imports
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u

logs = logging.getLogger(__name__)


# =============================================================================
# Google Sheets — request.csv and custom.csv
# =============================================================================
# Per-program CPS request tabs (export?format=csv). Produces on-disk request
# rows plus custom-window rows parsed from bracketed start/stop columns.

# Shared request fields through ``priority`` (exptime/maxtime are seconds; Keck / MAGIQ convention).
REQUEST_COLS_CORE = [
    "program_code",
    "target",
    "unique_id",
    "ra",
    "dec",
    "exptime",
    "maxtime",
    "n_exp",
    "n_inter_max",
    "tau_inter",
    "n_intra_max",
    "n_intra_min",
    "tau_intra",
    "minimum_elevation",
    "minimum_moon_separation",
    "weather_band_1",
    "weather_band_2",
    "weather_band_3",
    "gaia_id",
    "teff",
    "jmag",
    "Vmag",
    "pmra",
    "pmdec",
    "epoch",
    "exp_meter_threshold",
    "inactive",
    "decker",
    "cell in/out?",
    "priority",
]

# On-disk ``request.csv`` from prep (no start/stop columns): comments last.
REQUEST_COLS = REQUEST_COLS_CORE + ["comments"]

# Full Google Sheet ``requests`` tab (CPS template): … priority, start, stop, comments.
REQUEST_COLS_READ = REQUEST_COLS_CORE + ["start", "stop", "comments"]

# Column definitions for custom dataframe (built from start/stop on requests)
CUSTOM_COLS = ["unique_id", "target", "start", "stop"]


def _customs_from_requests_df(req_df):
    """Build ``custom.csv`` rows from bracketed ``start`` / ``stop`` sheet columns.

    Each cell holds a bracketed, comma-separated list, e.g.
    ``[2026-02-01 12:00, 2026-03-01 12:00]``. Start/stop lists are paired by
    index; empty or invalid cells yield no rows for that field.
    """
    empty = pd.DataFrame(columns=CUSTOM_COLS)
    if req_df is None or req_df.empty:
        return empty
    for col in ("start", "stop"):
        if col not in req_df.columns:
            raise ValueError(f"requests DataFrame missing required column: {col}")

    def _bracket_list(value):
        if not isinstance(value, str):
            return []
        text = value.strip()
        if len(text) < 2 or text[0] != "[" or text[-1] != "]":
            return []
        inner = text[1:-1].strip()
        return [part.strip() for part in inner.split(",") if part.strip()] if inner else []

    rows = []
    for _, row in req_df.iterrows():
        uid = row.get("unique_id", "")
        target = row.get("target", "")
        for start, stop in zip(
            _bracket_list(row.get("start")),
            _bracket_list(row.get("stop")),
        ):
            rows.append(
                {"unique_id": uid, "target": target, "start": start, "stop": stop}
            )
    return pd.DataFrame(rows, columns=CUSTOM_COLS) if rows else empty


_SHEET_ID_RE = re.compile(r"/spreadsheets/d/([^/?#]+)")
_GID_RE = re.compile(r"[?#&]gid=(\d+)")


def _fetch_sheet_dataframe(url, skip_rows=3):
    """
    Fetch one tab of a HIRES-CPS Google Sheet as a DataFrame.

    Given the address-bar URL of the relevant tab (must include both the
    workbook ID and the per-tab ``gid``, e.g.
    ``https://docs.google.com/spreadsheets/d/<SHEET_ID>/edit?gid=<GID>#gid=<GID>``),
    download via Google's public ``export?format=csv&gid=<GID>`` endpoint --
    the same as File > Download > CSV in the UI -- and parse using the fixed
    HIRES-CPS layout: rows 1-3 are template labels, row 4 is the canonical
    header, row 5+ is data.

    Returns a DataFrame with exactly ``REQUEST_COLS_READ`` columns.

    Raises ``ValueError`` for malformed URLs, HTML responses (sharing not set
    to "Anyone with the link"), missing ``program_code`` header, or any
    missing required columns.
    """
    url = (url or "").strip()
    sheet_match = _SHEET_ID_RE.search(url)
    if not sheet_match:
        raise ValueError(f"No Google Sheet ID in URL: {url[:80]}...")
    gid_match = _GID_RE.search(url)
    if not gid_match:
        raise ValueError(
            "URL must include a tab-specific gid (e.g. .../edit?gid=12345#gid=12345). "
            f"Got: {url[:120]}..."
        )
    sheet_id, gid = sheet_match.group(1), gid_match.group(1)

    csv_url = (
        f"https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=csv&gid={gid}"
    )
    print(f"downloading requests from {url}")
    resp = requests.get(csv_url, timeout=15)
    resp.raise_for_status()

    # Log label from Content-Disposition (``<workbook> - <tab>.csv``), else the sheet URL.
    cd = resp.headers.get("Content-Disposition", "")
    m = re.search(r"filename\*=UTF-8''([^;]+)", cd, re.I) or re.search(
        r'filename="([^"]+)"', cd, re.I
    )
    title = url
    if m:
        name = urllib.parse.unquote(m.group(1).strip())
        if name.endswith(".csv"):
            name = name[:-4]
        if " - " in name:
            name = name.rsplit(" - ", 1)[0]
        title = name

    text = resp.text
    stripped = text.lstrip()
    if not stripped or stripped.startswith("<!") or "<html" in stripped[:200].lower():
        raise ValueError(
            f"Sheet {sheet_id} (gid={gid}): export endpoint returned HTML, "
            "not CSV. Verify sharing is set to 'Anyone with the link'."
        )

    df = pd.read_csv(io.StringIO(text), skiprows=skip_rows, dtype=str)
    df = df.dropna(how="all")
    df.columns = [str(c).strip() for c in df.columns]
    # Google Sheets still use legacy column name; normalize to canonical schema.
    if "starname" in df.columns and "target" not in df.columns:
        df = df.rename(columns={"starname": "target"})
    elif "starname" in df.columns and "target" in df.columns:
        raise ValueError(
            "CSV has both 'starname' and 'target' columns; remove the duplicate."
        )
    if "comments" not in df.columns:
        df["comments"] = ""
    if "program_code" not in df.columns:
        raise ValueError(
            "CSV header row does not contain 'program_code' -- check that "
            "row 4 of the requests_new tab matches the canonical column "
            f"template. Found columns: {list(df.columns)[:8]}..."
        )
    df = df[df["program_code"].astype(str).str.strip() != ""]
    df = df[df["program_code"].astype(str).str.lower() != "nan"]
    missing = set(REQUEST_COLS_READ) - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {sorted(missing)}")
    print(f"read {len(df)} records from {title}\n")
    return df[REQUEST_COLS_READ].copy()


def _dedup_requests_by_hash(requests_df, custom_df):
    """
    Deduplicate ``requests_df`` rows that share the same ``unique_id`` across programs.

    For each duplicate group the winner is the row with the lowest SHA-256 digest of
    ``f"{program_code}__{unique_id}"``; ties (theoretically impossible for distinct
    canonical strings) break lexicographically on ``program_code``. Losers are dropped
    entirely and matching rows in ``custom_df`` are filtered out.

    A ``logs.warning`` is emitted for every duplicate group, naming the kept program and
    urging PIs to remove duplicates upstream.

    Args:
        requests_df (pd.DataFrame): Concatenated request rows.
        custom_df (pd.DataFrame | None): Associated custom-window rows keyed by
            ``unique_id``.

    Returns:
        tuple: ``(requests_df, custom_df)`` with duplicates removed.
    """
    if (
        requests_df is None
        or requests_df.empty
        or "unique_id" not in requests_df.columns
    ):
        return requests_df, custom_df

    df = requests_df.copy()
    canonical = df["program_code"].astype(str) + "__" + df["unique_id"].astype(str)
    scores = canonical.map(lambda s: hashlib.sha256(s.encode("utf-8")).hexdigest())
    df = df.assign(precedence_score=scores)

    keep_idx = []
    duplicate_blocks = []
    star_col = "target" if "target" in df.columns else "unique_id"
    for uid, grp in df.groupby("unique_id", sort=False):
        if len(grp) == 1:
            keep_idx.append(grp.index[0])
            continue
        ranked = grp.sort_values(
            ["precedence_score", "program_code"], ascending=[True, True]
        )
        winner_idx = ranked.index[0]
        keep_idx.append(winner_idx)
        block_lines = []
        for i, row_idx in enumerate(ranked.index):
            pc = str(ranked.loc[row_idx, "program_code"])
            sn = str(ranked.loc[row_idx, star_col])
            marker = "*" if i == 0 else " "
            block_lines.append(f"  {pc} {sn} {marker}".rstrip())
        duplicate_blocks.append("\n".join(block_lines))

    if duplicate_blocks:
        msg = (
            "Duplicate rows exist!\n"
            "\n"
            "- Requests selected based on hash scheme.\n"
            "- * indicates selected target\n"
            "- Coordinate with PIs to resolve duplicates\n"
            "\n" + "\n\n".join(duplicate_blocks)
        )
        logs.warning(msg)

    out = (
        df.loc[sorted(keep_idx)]
        .drop(columns=["precedence_score"])
        .reset_index(drop=True)
    )

    if (
        custom_df is not None
        and not custom_df.empty
        and "unique_id" in custom_df.columns
    ):
        kept_uids = set(out["unique_id"].astype(str))
        custom_df = custom_df[
            custom_df["unique_id"].astype(str).isin(kept_uids)
        ].reset_index(drop=True)

    return out, custom_df


def pull_requests(request_urls_path):
    """
    Pull HIRES-CPS request and custom-window data from the per-program Google
    Sheets listed in ``request_urls_path`` (column ``url``).

    For each URL, fetches the relevant tab via :func:`_fetch_sheet_dataframe`
    (Google's ``export?format=csv&gid=<GID>`` endpoint) and parses it using
    the fixed HIRES-CPS layout. Each tab must expose the columns in
    ``REQUEST_COLS_READ`` (canonical fields through ``priority`` plus
    ``start``, ``stop``, ``comments``); a missing ``comments`` column is
    backfilled as empty strings.

    Args:
        request_urls_path: Path to CSV with a ``url`` column listing sheet URLs.

    Returns:
        tuple: ``(requests_df, custom_df)`` where ``requests_df`` has
        ``REQUEST_COLS`` and ``custom_df`` has
        ``[unique_id, target, start, stop]``.

    Raises:
        ValueError: If ``request_urls_path`` is unset or empty.
    """
    if not request_urls_path:
        raise ValueError("request_urls_path is required.")
    sheet_urls = pd.read_csv(request_urls_path)["url"].tolist()

    request_dfs = []
    custom_dfs = []
    for url in sheet_urls:
        url = (url or "").strip()
        df = _fetch_sheet_dataframe(url)
        request_dfs.append(df[REQUEST_COLS])
        custom_dfs.append(_customs_from_requests_df(df))
    requests_df = (
        pd.concat(request_dfs, ignore_index=True)
        if request_dfs
        else pd.DataFrame(columns=REQUEST_COLS)
    )
    custom_df = (
        pd.concat(custom_dfs, ignore_index=True)
        if custom_dfs
        else pd.DataFrame(columns=CUSTOM_COLS)
    )
    # Convert ra (HH:MM:SS.ss) and dec (+/-DD:MM:SS.s) from sexagesimal to decimal degrees
    if (
        not requests_df.empty
        and "ra" in requests_df.columns
        and "dec" in requests_df.columns
    ):
        c = SkyCoord(
            ra=requests_df["ra"].astype(str),
            dec=requests_df["dec"].astype(str),
            unit=(u.hourangle, u.deg),
        )
        requests_df = requests_df.copy()
        requests_df["ra"] = c.ra.deg
        requests_df["dec"] = c.dec.deg
    requests_df, custom_df = _dedup_requests_by_hash(requests_df, custom_df)
    return requests_df, custom_df


# =============================================================================
# Keck schedule — allocation.csv
# =============================================================================
# Keck tel schedule query form (HIRESr nights). Crossmatched against
# request_urls.csv program codes to build allocation blocks for AstroQ.

KECK_SCHEDULE_QUERY_URL = (
    "https://www2.keck.hawaii.edu/observing/keckSchedule/queryForm.php"
)
KECK_SCHEDULE_INSTRUMENT = "HIRESr"


def pull_all_scheduled(start_date, end_date, output_path=None, timeout=60):
    """Query the Keck schedule form for HIRESr and return a DataFrame.

    Date bounds come from config ``semester_start_day`` / ``semester_end_day``.

    Columns: ``Date, Time, Dark, TelNr, Instrument, Account, PI, Institution, ProjCode``.
    If ``output_path`` is given, also write the same DataFrame to CSV.
    """
    payload = {
        "doQuery": "1",
        "table": "schedule",
        "Date": f"between {start_date} and {end_date}",
        "Instrument": KECK_SCHEDULE_INSTRUMENT,
        "cb_Date": "on",
        "cb_TelNr": "on",
        "cb_Instrument": "on",
        "cb_Account": "on",
        "cb_Principal": "on",
        "cb_Institution": "on",
        "cb_ProjCode": "on",
        "excel": "on",
        "sched": "Query Tel Schedule",
    }
    response = requests.post(KECK_SCHEDULE_QUERY_URL, data=payload, timeout=timeout)
    response.raise_for_status()
    text = response.text.strip()
    if not text.startswith("Date,"):
        snippet = text[:200].replace("\n", " ")
        raise RuntimeError(f"Unexpected response from schedule form: {snippet}")

    df = pd.DataFrame(csv.DictReader(io.StringIO(text)))
    if df.empty:
        return df
    df = df.sort_values(
        ["Date", "Time", "Instrument", "ProjCode"], kind="mergesort"
    ).reset_index(drop=True)
    if output_path is not None:
        df.to_csv(output_path, index=False)
    return df


def koip_block_to_utc_iso(date, start_time, end_time):
    """Build UTC ``allocation.csv`` timestamps from Keck KOIP schedule fields.

    ``date`` is the HST civil observing-night label (``YYYY-MM-DD``).
    ``start_time`` and ``end_time`` are UTC clock times (``HH:MM``).
    Keck convention: the UTC calendar date is the civil date plus one day.
    """
    civil = pd.to_datetime(date.astype(str).str.strip())
    utc_day = civil + pd.Timedelta(days=1)
    start = pd.to_datetime(
        utc_day.dt.strftime("%Y-%m-%d")
        + " "
        + start_time.astype(str).str.strip(),
        utc=True,
    )
    stop = pd.to_datetime(
        utc_day.dt.strftime("%Y-%m-%d")
        + " "
        + end_time.astype(str).str.strip(),
        utc=True,
    )
    stop = stop.where(stop >= start, stop + pd.Timedelta(days=1))
    return (
        start.dt.strftime("%Y-%m-%dT%H:%M"),
        stop.dt.strftime("%Y-%m-%dT%H:%M"),
    )


def crossmatch_allocation(scheduled_df, request_urls_path, semester, output_path=None):
    """Filter the all-scheduled Keck DataFrame to rows whose ProjCode appears in
    ``request_urls_<sem>.csv`` (column ``program_code``, e.g. ``2026A_C364``).

    Splits the schedule's ``Time`` cell (``"05:03 - 13:22 ( 75%)"``) into
    UTC ``StartTime`` / ``EndTime`` and emits AstroQ allocation columns from
    HST civil ``Date`` plus one day. If ``output_path`` is given, also writes
    the result as CSV.
    """
    req = pd.read_csv(request_urls_path)
    req["ProjCode"] = req["program_code"].str.removeprefix(f"{semester}_")
    matched = scheduled_df.merge(req[["ProjCode"]], on="ProjCode", how="inner")
    matched[["StartTime", "EndTime"]] = matched["Time"].str.extract(
        r"(\d{1,2}:\d{2})\s*-\s*(\d{1,2}:\d{2})"
    )
    matched["start"], matched["stop"] = koip_block_to_utc_iso(
        matched["Date"], matched["StartTime"], matched["EndTime"]
    )
    matched["comment"] = ""
    cols = [
        "Date",
        "StartTime",
        "EndTime",
        "start",
        "stop",
        "PI",
        "Instrument",
        "ProjCode",
        "comment",
    ]
    matched = (
        matched[cols]
        .sort_values(["Date", "StartTime", "ProjCode"])
        .reset_index(drop=True)
    )
    if output_path is not None:
        matched.to_csv(output_path, index=False)
    return matched


# =============================================================================
# JUMP — past.csv
# =============================================================================
# Authenticated explorer download of HIRES frame history; collapsed to one
# row per visit attempt for the semester planner.

JUMP_BASE_URL = "https://jump.caltech.edu"
# Parameterized JUMP explorer query returning all HIRES observations between
# ``start_date`` and ``end_date`` (UTC). Replaces the legacy named saved query
# "HIRES2026A - All Observations".
JUMP_HIRES_PAST_EXPLORER_ID = 285
JUMP_PAST_QUERY_TMP_FILENAME = "past_jump-query-tmp.csv"


def exposures_to_visits(data, request_csv_path=None):
    """Exposures to visits.
    
    Jump store HIRES exposures. However, we want to interpret them as visits since
    that's how AstroQ interprets them. 
    
    - We identify templates by searching for (B1/B3 decker, iodine out). These are
      labeled with a _t suffix.
    
    - We identify groups of successive exposures on the same target as groups. If the
      group has at least 50% of the exposures required for the target, we count it as a
      successful visit
    
    Convert raw JUMP frame rows to ``past.csv`` visit rows.
    """
    cols = ["unique_id", "target", "timestamp", "exposure_time"]
    if data.empty:
        return pd.DataFrame(columns=cols)

    n_exp_by_uid = {}
    if request_csv_path and os.path.isfile(request_csv_path):
        req = pd.read_csv(request_csv_path)
        if "unique_id" in req.columns and "n_exp" in req.columns:
            n_exp_by_uid = (
                req.set_index("unique_id")["n_exp"]
                .astype(int)
                .to_dict()
            )
        else:
            logs.warning(
                "request.csv missing unique_id/n_exp columns; assuming n_exp=1"
            )
    elif request_csv_path:
        logs.warning(
            "request.csv not found at %r; assuming n_exp=1 for all targets",
            request_csv_path,
        )

    df = data.copy()
    if {"decker", "iodine_in"}.issubset(df.columns):
        is_template = (
            df["decker"].astype(str).str.upper().isin(["B1", "B3"])
            & (df["iodine_in"] == False)
        )
        df.loc[is_template, "target"] = (
            df.loc[is_template, "target"].astype(str) + "_t"
        )

    df["target"] = df["target"].astype(str)
    df["unique_id"] = df["target"]
    df["exposure_time"] = (
        pd.to_numeric(df["exposure_time"], errors="coerce").fillna(0).astype(int)
    )
    df["_ts"] = pd.to_datetime(df["timestamp"], errors="coerce")
    df = df.loc[df["_ts"].notna()].copy()
    if df.empty:
        return pd.DataFrame(columns=cols)

    df["_night"] = df["_ts"].dt.strftime("%Y-%m-%d")
    rows = []
    rejected_groups = []
    for (uid, _night), grp in df.groupby(["unique_id", "_night"], sort=False):
        n_exp = int(n_exp_by_uid.get(uid, 1))
        min_frames = math.ceil(0.5 * n_exp)
        grp = grp.sort_values("_ts")
        n_frames = len(grp)
        if n_frames < min_frames:
            rejected_groups.append(
                f"  {uid} {_night}: {n_frames}/{n_exp} frames (need {min_frames})"
            )
            continue
        first = grp.iloc[0]
        rows.append(
            {
                "unique_id": uid,
                "target": uid,
                "timestamp": first["timestamp"],
                "exposure_time": int(grp["exposure_time"].sum()),
            }
        )

    summary = (
        f"{len(df)} raw frame(s) -> {len(rows)} accounted visit(s) "
        f"({len(rejected_groups)} group(s) below 50% threshold)"
    )
    if rejected_groups:
        summary += ":\n" + "\n".join(rejected_groups)
    print(summary)
    return pd.DataFrame(rows, columns=cols)


def get_hires_past_history(
    path_to_csv,
    semester_start_day=None,
    semester_end_day=None,
    request_csv_path=None,
):
    """Pull HIRES past history from JUMP and write processed ``path_to_csv``.

    Fetches the parameterized JUMP explorer query
    :data:`JUMP_HIRES_PAST_EXPLORER_ID` (285). Downloads the raw explorer CSV
    to ``past_jump-query-tmp.csv`` beside ``path_to_csv`` (preserved for
    inspection), then groups frames into visit attempts. A group counts toward
    ``past.csv`` only when ``len(frames) >= ceil(0.5 * n_exp)`` where ``n_exp``
    comes from ``request_csv_path``.

    Rows with ``decker`` B1/B3 and ``iodine_in`` False get ``_t`` appended to
    ``target`` so they match template request ``unique_id``s.

    Output schema: one row per accounted visit —
    ``unique_id, target, timestamp, exposure_time`` (``timestamp`` = first frame).

    Args:
        path_to_csv (str): Output CSV path.
        semester_start_day (str): ``YYYY-MM-DD`` from config
            ``[global] semester_start_day``; passed to JUMP as ``start_date``.
        semester_end_day (str): ``YYYY-MM-DD`` from config
            ``[global] semester_end_day``; passed to JUMP as ``end_date``.
        request_csv_path (str, optional): ``request.csv`` path supplying
            ``unique_id`` / ``n_exp`` for visit-collapse thresholds.

    Raises:
        ValueError: if ``semester_start_day`` or ``semester_end_day`` is
            missing (both are required to build the query window).
    """
    if not semester_start_day or not semester_end_day:
        raise ValueError(
            "get_hires_past_history requires both semester_start_day and "
            "semester_end_day to build the JUMP explorer query window."
        )

    raw_path = os.path.join(os.path.dirname(path_to_csv), JUMP_PAST_QUERY_TMP_FILENAME)

    login_url = f"{JUMP_BASE_URL}/user/login/"
    username = os.environ.get("KPFCC_JUMP_USERNAME")
    password = os.environ.get("KPFCC_JUMP_PASSWORD")
    if not username or not password:
        raise RuntimeError(
            "Missing JUMP credentials. Set KPFCC_JUMP_USERNAME and "
            "KPFCC_JUMP_PASSWORD in the environment or in the workspace .env file."
        )

    session = requests.Session()
    login_page = session.get(login_url, timeout=60)
    login_page.raise_for_status()
    login_resp = session.post(
        login_url,
        data={
            "action": "login",
            "username": username,
            "password": password,
            "csrfmiddlewaretoken": session.cookies["csrftoken"],
        },
        headers={"Referer": login_url},
        timeout=60,
    )
    login_resp.raise_for_status()
    if login_resp.url.rstrip("/") == login_url.rstrip("/"):
        raise RuntimeError(
            "JUMP login appears to have failed: still on login page after submitting credentials. "
            "Check KPFCC_JUMP_USERNAME/KPFCC_JUMP_PASSWORD."
        )

    param_str = (
        f"start_date:{semester_start_day}|end_date:{semester_end_day}"
    )
    param_slug = urllib.parse.quote(param_str, safe="")
    # Build the download URL directly. Scraping the page's download <a href>
    # corrupts the query string: HTML parsers decode ``&params`` as ``&para``
    # -> ``¶``, yielding ``¶ms=`` and a 500 from JUMP.
    download_url = (
        f"{JUMP_BASE_URL}/explorer/{JUMP_HIRES_PAST_EXPLORER_ID}/download"
        f"?format=csv&params={param_slug}"
    )
    print(f"JUMP download URL: {download_url}")

    csv_response = session.get(download_url, timeout=60)
    csv_response.raise_for_status()
    with open(raw_path, "wb") as f:
        f.write(csv_response.content)
    print(f"Raw JUMP pull saved to {raw_path}")

    data = pd.read_csv(raw_path)
    if "starname" in data.columns and "target" not in data.columns:
        data = data.rename(columns={"starname": "target"})
    elif "starname" in data.columns:
        data = data.rename(columns={"starname": "target"})

    visits = exposures_to_visits(data, request_csv_path)
    visits.to_csv(path_to_csv, index=False)
    print(f"Processed past history saved to {path_to_csv}")
