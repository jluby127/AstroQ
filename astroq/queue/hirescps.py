"""
Module for preparing all data specific to from Keck Observatory's custom made Observing Block (OB) database.
This is specific to the KPF-CC program and the observatory's infrastructure as way to power the prep kpfcc command.
New observatories should write their own module to connect to a new "prep <your observatory>" command.
"""

# Standard library imports
import csv
import hashlib
import io
import logging
import os
import urllib.parse
from datetime import date

# Third-party imports
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
import astropy.units as u
import re
from bs4 import BeautifulSoup

# Local imports
from astroq.queue.base import Queue

logs = logging.getLogger(__name__)


class HIRESCPS(Queue):
    """HIRES-CPS on Keck-I.

    HIRES is permanently installed on Keck-I, so the telescope/instrument
    pairing is unique and this single class fully describes the queue.

    Pointing geometry references the Keck-I limits page:
    https://www2.keck.hawaii.edu/inst/common/TelLimits.html

    The upper elevation limit (85 deg) matches AstroQ's accessibility policy
    rather than TTP's historical Keck1 zenith limit (84 deg); the physical
    zenith limit is 88.9 deg.
    """

    slew_rate = 0.6  # deg/s; matches TTP Keck1 (6./10.)
    wrap_limit = 235.0  # deg azimuth
    nSlots = 4  # TTP slew-slot granularity
    readout_time = 45.0  # seconds; per-shot detector readout
    slew_overhead_mean = (
        60.0  # seconds; mean per-visit slew + acquisition (splan-only estimate)
    )

    # Inaccessible (alt, az) boxes, degrees. (az_min, az_max, alt_min, alt_max).
    # A sky point is excluded iff it lies inside ANY box. See Queue.is_accessible.
    # Duplicated independently in HIRESCPS and KPFCC; the two queues may
    # legitimately diverge on elevation policy.
    inaccessible_zones = [
        (5.3, 146.2, 0.0, 33.3),  # Nasmyth deck obstruction
        (0.0, 360.0, 0.0, 18.0),  # below 18 deg elevation clamp
        (0.0, 360.0, 85.0, 90.0),  # above 85 deg elevation clamp
    ]

    # Constraints `Access` should compute for HIRES-CPS. Omits ``"clear"``;
    # weather loss is handled separately and the cube defaults to all-True.
    access_constraints = (
        "altaz", "future", "moon", "custom", "inter", "allocated",
    )

    def __init__(self):
        import astroplan as apl

        self.observatory = apl.Observer.at_site(
            "Keck Observatory", name="Keck", timezone="US/Hawaii"
        )

    def write_starlist(self, *args, **kwargs):
        return write_starlist(*args, **kwargs)


# Shared request fields through ``priority`` (exptime/maxtime are seconds; Keck / MAGIQ convention).
REQUEST_COLS_CORE = [
    "program_code",
    "starname",
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
CUSTOM_COLS = ["unique_id", "starname", "start", "stop"]


def _parse_bracket_array(s):
    """
    Parse a string like '[2026-02-01 12:00, 2026-03-01 12:00]' into a list of strings.
    Handles '[]', single element '[2026-02-01 12:00]', and multiple. Returns [] if invalid.
    """
    if s is None or not isinstance(s, str):
        return []
    s = s.strip()
    if not s or len(s) < 2 or s[0] != "[" or s[-1] != "]":
        return []
    inner = s[1:-1].strip()
    if not inner:
        return []
    return [part.strip() for part in inner.split(",") if part.strip()]


def _customs_from_requests_df(req_df):
    """
    Build customs DataFrame from requests DataFrame that has start/stop columns.
    start/stop are strings like '[2026-02-01 12:00, 2026-03-01 12:00]'. Pairs by index.
    """
    if req_df is None or req_df.empty:
        return pd.DataFrame(columns=CUSTOM_COLS)
    for col in ("start", "stop"):
        if col not in req_df.columns:
            raise ValueError(f"requests DataFrame missing required column: {col}")
    rows = []
    for _, r in req_df.iterrows():
        uid = r.get("unique_id", "")
        star = r.get("starname", "")
        starts = _parse_bracket_array(r.get("start"))
        stops = _parse_bracket_array(r.get("stop"))
        n = min(len(starts), len(stops))
        for i in range(n):
            rows.append(
                {
                    "unique_id": uid,
                    "starname": star,
                    "start": starts[i],
                    "stop": stops[i],
                }
            )
    return (
        pd.DataFrame(rows, columns=CUSTOM_COLS)
        if rows
        else pd.DataFrame(columns=CUSTOM_COLS)
    )


_SHEET_ID_RE = re.compile(r"/spreadsheets/d/([^/?#]+)")
_GID_RE = re.compile(r"[?#&]gid=(\d+)")
_FILENAME_STAR_RE = re.compile(r"filename\*=UTF-8''([^;]+)", re.IGNORECASE)
_FILENAME_RE = re.compile(r'filename="([^"]+)"', re.IGNORECASE)


def _workbook_title_from_response(resp):
    """
    Extract the Google Sheets workbook title from a CSV-export response.

    Google's ``export?format=csv`` endpoint returns ``Content-Disposition``
    of the form ``attachment; filename="<workbook> - <tab>.csv"; filename*=UTF-8''...``.
    Prefer the RFC 5987 ``filename*=UTF-8''...`` form (which preserves the
    space separator after URL-decoding) and split off the trailing ``- <tab>``.
    Returns ``None`` if no filename is present in the header.
    """
    cd = resp.headers.get("Content-Disposition", "")
    m = _FILENAME_STAR_RE.search(cd)
    if m:
        name = urllib.parse.unquote(m.group(1).strip())
    else:
        m = _FILENAME_RE.search(cd)
        if not m:
            return None
        name = m.group(1)
    if name.endswith(".csv"):
        name = name[:-4]
    if " - " in name:
        name = name.rsplit(" - ", 1)[0]
    return name


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
    title = _workbook_title_from_response(resp) or url
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
    star_col = "starname" if "starname" in df.columns else "unique_id"
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


def pull_requests():
    """
    Pull HIRES-CPS request and custom-window data from the per-program Google
    Sheets listed in the CSV at ``$HIRES_PROGRAM_SHEET_URLS_CSV`` (column
    ``url``).

    For each URL, fetches the relevant tab via :func:`_fetch_sheet_dataframe`
    (Google's ``export?format=csv&gid=<GID>`` endpoint) and parses it using
    the fixed HIRES-CPS layout. Each tab must expose the columns in
    ``REQUEST_COLS_READ`` (canonical fields through ``priority`` plus
    ``start``, ``stop``, ``comments``); a missing ``comments`` column is
    backfilled as empty strings.

    Returns:
        tuple: ``(requests_df, custom_df)`` where ``requests_df`` has
        ``REQUEST_COLS`` and ``custom_df`` has
        ``[unique_id, starname, start, stop]``.

    Raises:
        ValueError: If ``HIRES_PROGRAM_SHEET_URLS_CSV`` is unset or empty.
    """
    path = os.environ.get("HIRES_PROGRAM_SHEET_URLS_CSV")
    if not path:
        raise ValueError(
            "HIRES_PROGRAM_SHEET_URLS_CSV is not set. "
            "Point it at a CSV with a 'url' column, e.g. request_urls_2026A.csv."
        )
    sheet_urls = pd.read_csv(path)["url"].tolist()

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


KECK_SCHEDULE_QUERY_URL = (
    "https://www2.keck.hawaii.edu/observing/keckSchedule/queryForm.php"
)
DEFAULT_INSTRUMENTS = ("KPF", "KPF-CC", "HIRES")


def canonicalize_instrument(name):
    """Map Keck schedule instrument variants (HIRESr/HIRESb/...) to a canonical name."""
    value = (name or "").strip()
    lower = value.lower()
    if lower.startswith("hires"):
        return "HIRES"
    if lower.startswith("kpf-cc"):
        return "KPF-CC"
    if lower == "kpf":
        return "KPF"
    return value


def infer_current_semester(today=None):
    """Return the Keck semester (e.g. ``"2026A"``) for ``today``.

    Convention: months 02-07 are A, 08+ are B, 01 is the previous year's B.
    """
    today = today or date.today()
    if 2 <= today.month <= 7:
        return f"{today.year}A"
    if today.month >= 8:
        return f"{today.year}B"
    return f"{today.year - 1}B"


def semester_date_range(semester):
    """Return ``(start_date, end_date)`` ISO strings for a semester like ``"2026A"``."""
    semester = semester.strip().upper()
    if (
        len(semester) != 5
        or semester[-1] not in {"A", "B"}
        or not semester[:4].isdigit()
    ):
        raise ValueError(
            f"Invalid semester format: {semester!r}. Expected forms like 2026A or 2026B."
        )
    year = int(semester[:4])
    if semester[-1] == "A":
        return f"{year}-02-01", f"{year}-07-31"
    return f"{year}-08-01", f"{year + 1}-01-31"


def fetch_schedule_csv(semester, instrument, timeout=60):
    """Submit the Keck schedule query form for one instrument and return matching rows."""
    start_date, end_date = semester_date_range(semester)
    payload = {
        "doQuery": "1",
        "table": "schedule",
        "Date": f"between {start_date} and {end_date}",
        "Instrument": instrument,
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

    reader = csv.DictReader(io.StringIO(text))
    rows = [dict(row) for row in reader]

    canonical_requested = canonicalize_instrument(instrument).lower()
    exact_rows = []
    for row in rows:
        normalized = canonicalize_instrument(row.get("Instrument") or "")
        if normalized.lower() == canonical_requested:
            new_row = dict(row)
            new_row["Instrument"] = normalized
            exact_rows.append(new_row)
    return exact_rows


def pull_all_scheduled(semester, instruments=DEFAULT_INSTRUMENTS, output_path=None):
    """Query the Keck schedule form for each instrument and return one combined DataFrame.

    Columns: ``Date, Time, Dark, TelNr, Instrument, Account, PI, Institution, ProjCode``.
    If ``output_path`` is given, also write the same DataFrame to CSV.
    """
    rows = []
    for inst in instruments:
        rows.extend(fetch_schedule_csv(semester, inst))
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["Instrument"] = df["Instrument"].map(canonicalize_instrument)
    df = df.sort_values(
        ["Date", "Time", "Instrument", "ProjCode"], kind="mergesort"
    ).reset_index(drop=True)
    if output_path is not None:
        df.to_csv(output_path, index=False)
    return df


def crossmatch_allocation(scheduled_df, request_urls_path, semester, output_path=None):
    """Filter the all-scheduled Keck DataFrame to rows whose ProjCode appears in
    ``request_urls_<sem>.csv`` (column ``program_code``, e.g. ``2026A_C364``).

    Splits the schedule's ``Time`` cell (``"05:03 - 13:22 ( 75%)"``) into
    ``StartTime`` / ``EndTime`` and emits the AstroQ allocation columns. If
    ``output_path`` is given, also writes the result as CSV.
    """
    req = pd.read_csv(request_urls_path)
    req["ProjCode"] = req["program_code"].str.removeprefix(f"{semester}_")
    matched = scheduled_df.merge(req[["ProjCode"]], on="ProjCode", how="inner")
    matched[["StartTime", "EndTime"]] = matched["Time"].str.extract(
        r"(\d{1,2}:\d{2})\s*-\s*(\d{1,2}:\d{2})"
    )
    matched["start"] = matched["Date"] + "T" + matched["StartTime"]
    matched["stop"] = matched["Date"] + "T" + matched["EndTime"]
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


def login_JUMP():
    login_url = "https://jump.caltech.edu/user/login/"
    s = requests.session()
    login_page = s.get(login_url)
    login_page.raise_for_status()
    csrftoken = s.cookies["csrftoken"]
    username = os.environ.get("KPFCC_JUMP_USERNAME")
    password = os.environ.get("KPFCC_JUMP_PASSWORD")
    if not username or not password:
        raise RuntimeError(
            "Missing JUMP credentials. Set KPFCC_JUMP_USERNAME and "
            "KPFCC_JUMP_PASSWORD in the environment or in the workspace .env file."
        )
    payload = {
        "action": "login",
        "username": username,
        "password": password,
        "csrfmiddlewaretoken": csrftoken,
    }
    new_login = s.post(login_url, data=payload, headers=dict(Referer=login_url))
    new_login.raise_for_status()
    if new_login.url.rstrip("/") == login_url.rstrip("/"):
        raise RuntimeError(
            "JUMP login appears to have failed: still on login page after submitting credentials. "
            "Check KPFCC_JUMP_USERNAME/KPFCC_JUMP_PASSWORD."
        )
    return s


def get_database_explorer(
    name, path_for_csv, url="https://jump.caltech.edu/explorer/", links=None
):
    # log into JUMP and go to DataBase page
    if links is None:
        links = []
    session = login_JUMP()
    response = session.get(url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")
    # find the table of queries
    table = soup.find("tbody", attrs={"class": "list"})
    if table is None:
        title = (
            soup.title.string.strip()
            if soup.title and soup.title.string
            else "unknown title"
        )
        preview = " ".join(soup.get_text(" ", strip=True).split()[:40])
        raise RuntimeError(
            "JUMP explorer page did not contain the expected query table. "
            f"Fetched URL: {response.url!r}. Page title: {title!r}. "
            f"Response preview: {preview!r}"
        )
    # find the correct row by the query name
    for row in table.find_all("tr"):
        tab = row.find("td", attrs={"class": "name"})
        try:
            x = tab.a.string
        except:
            pass
        else:
            # once it finds the match, save the download link
            if x == name:
                for l in row.find_all("a", href=re.compile("download")):
                    links.append("/".join(url.split("/")[:3]) + l.get("href"))
    # make sure it finds it before saving
    if links != []:
        response = session.get(links[0])
        response.raise_for_status()
        # pretty sure you can just read this directly into pandas if that's what you want
        open(path_for_csv, "wb").write(response.content)
    else:
        raise RuntimeError(
            f"JUMP explorer did not contain a download link for query {name!r}. "
            "The query may have been renamed or the explorer markup may have changed."
        )
    session.keep_alive = False
    return


def get_hires_past_history(path_to_csv, semester_start_day=None):
    """
    Pull HIRES past history from Jump and write ``path_to_csv``.

    Args:
        path_to_csv (str): Output CSV path.
        semester_start_day (str, optional): ``YYYY-MM-DD``; rows with ``timestamp`` strictly
            before this calendar instant are dropped so ``past.csv`` matches the current
            semester (avoids KeyError in internight logic when last obs is outside planned nights).
    """
    name = "HIRES2026A - All Observations"
    # comment this line out when playing with synthetic schedules
    get_database_explorer(name, path_to_csv)
    print("All KPF observations pulled from Jump. Saved to csv: " + path_to_csv)

    data = pd.read_csv(path_to_csv)

    data = data.rename(columns={"starname": "target", "program_name": "semid"})

    if semester_start_day and "timestamp" in data.columns and len(data) > 0:
        ts = pd.to_datetime(data["timestamp"], errors="coerce")
        cutoff = pd.Timestamp(semester_start_day)
        n_before = len(data)
        valid_ts = ts.notna()
        keep = valid_ts & (ts >= cutoff)
        data = data.loc[keep].copy()
        dropped = n_before - len(data)
        if dropped:
            print(
                f"Dropped {dropped} past-history row(s) with timestamp before semester start "
                f"{semester_start_day}."
            )

    # Add dummy columns
    data["exposure_start_time"] = data["timestamp"]
    data["observer"] = ""

    # Create unique_id column
    data["id"] = data["target"]

    # make ints
    data["exposure_time"] = data["exposure_time"].astype(int)

    data.to_csv(path_to_csv, index=False)
    print("Data cleaned. Done.")


def write_starlist(
    frame,
    schedule,
    night_start_time,
    filler_stars,
    current_day,
    outputdir,
    version="nominal",
    all_active_requests=None,
    past_history=None,
):
    """
    Generate the nightly script in the correct format.

    Args:
        frame (dataframe): the request.csv in dataframe format for just the targets that were selected to be observed tonight
        schedule (dataframe): TTP ``schedule`` DataFrame (parallel to ``nodes``)
        night_start_time (astropy time object): Beginning of observing interval
        filler_stars (array): star names of the stars added in the bonus round
        current_day (str): today's date in format YYYY-MM-DD
        outputdir (str): the directory to save the script file
        version (str): a tag for thescript (e.g. nominal, slowdown, backups, etc)
        all_active_requests (pd.DataFrame | None): full active request frame for
            the semester. When provided together with ``past_history``, a third
            ``BACKUPS`` section is appended listing every active request along
            with an ``obs=N/M`` token (past nights observed / requested).
        past_history (dict | None): mapping ``unique_id -> StarHistory``. Used to
            compute the ``obs=N/M`` token for the BACKUPS section.

    Returns:
        list[str]: the lines written to the script file.
    """
    frame["starname"] = frame["starname"].astype(str)

    on_sky = schedule[~schedule["is_anchor"]]
    scheduled_df = on_sky[on_sky["scheduled"]].sort_values("order")
    extras_df = on_sky[~on_sky["scheduled"]]

    total_exptime = 0
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    script_file = os.path.join(
        outputdir, "script_{}_{}.txt".format(current_day, version)
    )

    lines = []
    for _, srow in scheduled_df.iterrows():
        uid = str(srow["unique_id"])
        filler_flag = uid in filler_stars
        row = frame.loc[frame["unique_id"] == uid]
        row.reset_index(inplace=True)
        total_exptime += float(row["exptime"].iloc[0])

        start_exposure_hst = str(
            TimeDelta(srow["t_start"] * 60, format="sec") + night_start_time
        )[11:16]
        first_available_hst = str(
            TimeDelta(srow["t_early"] * 60, format="sec") + night_start_time
        )[11:16]
        last_available_hst = str(
            TimeDelta(srow["t_late"] * 60, format="sec") + night_start_time
        )[11:16]
        lines.append(
            format_hires_row(
                row,
                start_exposure_hst,
                first_available_hst,
                last_available_hst,
                current_day,
                filler_flag=filler_flag,
            )
        )

    lines.append("")
    lines.append("X" * 45 + "EXTRAS" + "X" * 45)
    lines.append("")

    for _, erow in extras_df.iterrows():
        uid = str(erow["unique_id"])
        filler_flag = uid in filler_stars
        row = frame.loc[frame["unique_id"] == uid]
        row.reset_index(inplace=True)
        first_available_hst = str(
            TimeDelta(erow["t_early"] * 60, format="sec") + night_start_time
        )[11:16]
        last_available_hst = str(
            TimeDelta(erow["t_late"] * 60, format="sec") + night_start_time
        )[11:16]
        lines.append(
            format_hires_row(
                row,
                "24:00",
                first_available_hst,
                last_available_hst,
                current_day,
                filler_flag,
                True,
            )
        )

    if all_active_requests is not None and past_history is not None:
        backup_df = all_active_requests.copy()
        backup_df["_ra_float"] = pd.to_numeric(backup_df["ra"], errors="coerce")
        backup_df["_vmag_float"] = pd.to_numeric(
            backup_df.get("Vmag", pd.Series(dtype=float)), errors="coerce"
        )
        backup_df = backup_df.sort_values("_ra_float", kind="mergesort")

        def emit_block(header, sub_df):
            lines.append("")
            lines.append(header)
            lines.append("")
            for _, req_row in sub_df.iterrows():
                uid = req_row["unique_id"]
                n_done = (
                    past_history[uid].total_n_unique_nights
                    if uid in past_history
                    else 0
                )
                n_req_raw = req_row.get("n_inter_max", 0)
                n_req = int(n_req_raw) if pd.notna(n_req_raw) else 0
                obs_token = f"obs={n_done}/{n_req}"

                row = sub_df.loc[sub_df["unique_id"] == uid].head(1)
                row = row.reset_index()
                lines.append(
                    format_hires_row(
                        row,
                        None,
                        None,
                        None,
                        current_day,
                        filler_flag=False,
                        extra=False,
                        omit_timing=True,
                        obs_token=obs_token,
                    )
                )

        emit_block(
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 2026A - Requests - All XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
            backup_df,
        )
        emit_block(
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 2026A - Requests - V < 8 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
            backup_df[backup_df["_vmag_float"] < 8],
        )

    # add buffer lines to end of file
    lines.append("")
    lines.append("")

    with open(script_file, "w") as f:
        f.write("\n".join(lines))
    # logs.info("Total Open Shutter Time Scheduled: " + str(np.round((total_exptime/3600),2)) + " hours")
    return lines


def format_hires_row(
    row,
    obs_time,
    first_available,
    last_available,
    current_day,
    filler_flag=False,
    extra=False,
    omit_timing=False,
    obs_token=None,
):
    """
    Format request data in the specific way needed for the script (relates to the Keck "Magiq"
    software's data ingestion requirements).

    Args:
        row (dataframe): a single row from the requests sheet dataframe
        obs_time (str): the timestamp of the night to begin the exposure according to the TTP.
                        In format HH:MM in HST timezone
        first_available (str): the timestamp of the night where the star is first accessible.
                                In format HH:MM in HST timezone.
        last_available (str): the timestamp of the night where the star is last accessible.
                                In format HH:MM in HST timezone.
        filler_flag (boolean): True of the target was added in the bonus round
        extra (boolean): is this an "extra" target
        omit_timing (boolean): if True, drop the trailing
            ``<obs_time> <first_available> <last_available>`` triplet (used for
            the BACKUPS section, where targets are not tied to a specific time).
            ``obs_time`` / ``first_available`` / ``last_available`` may be passed
            as ``None`` in this mode.
        obs_token (str | None): optional trailing token (e.g. ``"obs=3/10"``)
            appended at the very end of the line. Used by the BACKUPS section to
            show past-vs-requested observation counts.

    Returns:
        line (str): the properly formatted string to be included in the script file
    """

    equinox = "2000"
    # Treat missing/NaN proper motions as zero so apply_space_motion stays finite.
    pmra_raw = row.get("pmra", pd.Series([0.0])).iloc[0] if "pmra" in row else 0.0
    pmdec_raw = row.get("pmdec", pd.Series([0.0])).iloc[0] if "pmdec" in row else 0.0
    pmra = 0.0 if pd.isna(pmra_raw) else float(pmra_raw)
    pmdec = 0.0 if pd.isna(pmdec_raw) else float(pmdec_raw)

    current_time = Time(current_day)
    coord = SkyCoord(
        ra=row["ra"].iloc[0] * u.deg,
        dec=row["dec"].iloc[0] * u.deg,
        pm_ra_cosdec=pmra * u.mas / u.yr,
        pm_dec=pmdec * u.mas / u.yr,
        obstime=Time(f"J{equinox}"),
    )
    new_coord = coord.apply_space_motion(new_obstime=current_time)
    updated_ra = new_coord.ra.to_string(
        unit=u.hourangle, sep=" ", pad=True, precision=1
    )
    updated_dec = new_coord.dec.to_string(
        unit=u.deg, sep=" ", pad=True, precision=0
    )
    if updated_dec[0] != "-":
        updated_dec = "+" + updated_dec

    # Annotate the script line with the observation epoch when proper motion was actually
    # applied, so the observer can tell the coords are propagated (frame/equinox is still J2000).
    epoch_token = (
        f"epoch={current_time.jyear:.1f}"
        if (pmra != 0.0 or pmdec != 0.0)
        else None
    )

    starname_str = str(row["starname"].iloc[0])
    namestring = " " * (16 - len(starname_str[:16])) + starname_str[:16]

    # Handle missing columns with default values
    vmag_val = row.get("Vmag", [15.0])[0] if "Vmag" in row else 15.0

    try:
        vmag_val = float(vmag_val) if vmag_val is not None else 15.0
    except (ValueError, TypeError):
        vmag_val = 25.0

    exposurestring = (
        " " * (4 - len(str(int(row["exptime"].iloc[0]))))
        + str(int(row["exptime"].iloc[0]))
        + "/"
        + str(int(row["maxtime"].iloc[0]))
        + " " * (4 - len(str(int(row["maxtime"].iloc[0]))))
    )

    ofstring = "1of" + str(int(row["n_intra_max"].iloc[0]))

    numstring = str(int(row["n_exp"].iloc[0])) + "x"
    vmagstring = (
        "vmag="
        + str(np.round(float(vmag_val), 1))
        + " " * (4 - len(str(np.round(float(vmag_val), 1))))
    )

    programstring = row["program_code"].iloc[0]
    priostring = row["priority"].iloc[0]
    deckerstring = row["decker"].iloc[0]
    cellstring = row["cell in/out?"].iloc[0]
    exp_meter_thresholdstring = row["exp_meter_threshold"].iloc[0]

    if extra == False:
        timestring2 = str(obs_time)
    else:
        # designate a nonsense time
        timestring2 = "24:00"

    line = (
        namestring
        + " "
        + updated_ra
        + " "
        + updated_dec
        + " "
        + str(equinox)
        + " "
        + vmagstring
        + " "
        + exposurestring
        + " "
        + exp_meter_thresholdstring
        + " "
        + deckerstring
        + " "
        + numstring
        + " "
        + cellstring
        + " "
        + priostring
        + " CC "
        + programstring
    )

    if not omit_timing:
        line += " " + timestring2 + " " + first_available + " " + last_available

    if epoch_token is not None:
        line += " " + epoch_token

    # Handle missing Observing Notes column
    observing_notes = (
        row.get("Observing Notes", [""])[0] if "Observing Notes" in row else ""
    )
    if observing_notes and not pd.isnull(observing_notes):
        line += " " + str(observing_notes)

    if obs_token is not None:
        line += " " + str(obs_token)

    return line

