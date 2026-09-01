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

from astroq.queue.prep_common import standardize_request_strategy
from astroq.queue.hirescps import jump as jump_client

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
    "jump_name",
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

# On-disk ``request.csv`` from prep (no start/stop columns).
REQUEST_COLS = REQUEST_COLS_CORE + ["splan_weight", "comments"]

# Sheet rows before ``attach_splan_weight`` (no ``splan_weight`` column).
SHEET_REQUEST_COLS = REQUEST_COLS_CORE + ["comments"]

PRIORITY_TO_SPLAN_WEIGHT = {
    "p1": 1,
    "p2": 2,
    "p3": 3,
    "p4": 3,
    "p5": 3,
}
_DEFAULT_SPLAN_WEIGHT_NO_PRIORITIES = 2
_UNKNOWN_PRIORITY_SPLAN_WEIGHT = 3

# Full Google Sheet ``requests`` tab (CPS template): … priority, start, stop, comments.
# ``jump_name`` is optional on sheets (defaults to ``unique_id`` at load time).
REQUEST_COLS_READ = REQUEST_COLS_CORE + ["start", "stop", "comments"]
_REQUEST_COLS_READ_REQUIRED = set(REQUEST_COLS_READ) - {"jump_name"}

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

# Row-1 human labels (CPS template) → canonical column names when row-4 machine
# headers are blank (common for the trailing ``comments`` column).
_HUMAN_LABEL_TO_CANONICAL = {
    "comments": "comments",
}


def _human_header_labels(raw_text: str) -> list[str]:
    """Return stripped human labels from CPS template row 1."""
    row = next(csv.reader(io.StringIO(raw_text)), [])
    return [str(c).strip() for c in row]


def _is_blank_sheet_header(name: str) -> bool:
    text = str(name).strip()
    return not text or text.startswith("Unnamed:")


def _canonicalize_sheet_columns(
    df: pd.DataFrame, human_labels: list[str]
) -> pd.DataFrame:
    """Resolve blank/Unnamed machine headers using row-1 human labels."""
    cols = list(df.columns)
    for i, col in enumerate(cols):
        if not _is_blank_sheet_header(col):
            continue
        if i >= len(human_labels):
            continue
        canonical = _HUMAN_LABEL_TO_CANONICAL.get(human_labels[i].strip().lower())
        if canonical:
            cols[i] = canonical
    df.columns = cols

    for alias, canonical in (("Comments", "comments"), ("Observing Notes", "comments")):
        if alias not in df.columns:
            continue
        if canonical not in df.columns:
            df = df.rename(columns={alias: canonical})
            continue
        empty = df[canonical].isna() | (df[canonical].astype(str).str.strip() == "")
        df.loc[empty, canonical] = df.loc[empty, alias]
        df = df.drop(columns=[alias])
    return df


def backfill_jump_name(requests_df):
    """Ensure ``jump_name`` is populated (defaults to ``unique_id``)."""
    if requests_df is None or requests_df.empty:
        return requests_df
    df = requests_df.copy()
    if "jump_name" not in df.columns:
        df["jump_name"] = df["unique_id"].astype(str)
        return df
    empty = df["jump_name"].isna() | (df["jump_name"].astype(str).str.strip() == "")
    df.loc[empty, "jump_name"] = df.loc[empty, "unique_id"].astype(str)
    return df


def _parse_cps_sheet_text(text: str, title: str, skip_rows: int = 3) -> pd.DataFrame:
    """Parse CPS template CSV text into ``REQUEST_COLS_READ`` columns."""
    human_labels = _human_header_labels(text)
    df = pd.read_csv(io.StringIO(text), skiprows=skip_rows, dtype=str)
    df = df.dropna(how="all")
    df.columns = [str(c).strip() for c in df.columns]
    df = _canonicalize_sheet_columns(df, human_labels)
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
    missing = _REQUEST_COLS_READ_REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {sorted(missing)}")
    if "jump_name" not in df.columns:
        df["jump_name"] = ""
    print(f"read {len(df)} records from {title}\n")
    return df[REQUEST_COLS_READ].copy()


def _load_request_tab(source: str) -> pd.DataFrame:
    """Load one CPS request tab from a Google Sheet URL or local path."""
    source = (source or "").strip()
    if source.startswith("http"):
        url = source
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
        resp = requests.get(csv_url, timeout=60)
        resp.raise_for_status()
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
    else:
        path = os.path.abspath(source)
        print(f"reading requests from local file {path}")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        title = os.path.basename(path)

    stripped = text.lstrip()
    if stripped.startswith("<!") or "<html" in stripped[:200].lower():
        raise ValueError(f"Expected CSV for {title}, got HTML.")
    return _parse_cps_sheet_text(text, title)


def _fetch_sheet_dataframe(url, skip_rows=3):
    """Fetch one tab of a HIRES-CPS Google Sheet as a DataFrame."""
    del skip_rows  # fixed CPS layout: skip 3 label rows in :func:`_parse_cps_sheet_text`
    return _load_request_tab(url)


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


def _normalize_priority_token(value):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip().lower()
    if not text or text == "nan":
        return None
    return text


def priority_token_to_splan_weight(token):
    """Map observer priority token (``p1`` … ``p5``) to semester MILP weight."""
    key = _normalize_priority_token(token)
    if key is None:
        return None
    weight = PRIORITY_TO_SPLAN_WEIGHT.get(key)
    if weight is None:
        logs.warning(
            "Unknown priority token %r; using splan_weight=%s",
            token,
            _UNKNOWN_PRIORITY_SPLAN_WEIGHT,
        )
        return _UNKNOWN_PRIORITY_SPLAN_WEIGHT
    return weight


def attach_splan_weight(requests_df):
    """Add ``splan_weight`` from ``priority`` with empty-row fill rules."""
    if requests_df is None or requests_df.empty:
        out = (
            requests_df.copy()
            if requests_df is not None
            else pd.DataFrame(columns=REQUEST_COLS)
        )
        if "splan_weight" not in out.columns:
            out["splan_weight"] = pd.Series(dtype=int)
        return out

    df = requests_df.copy()
    if "priority" not in df.columns:
        df["priority"] = ""

    tokens = df["priority"].map(_normalize_priority_token)
    splan = tokens.map(priority_token_to_splan_weight)

    specified = splan[tokens.notna()]
    if specified.empty:
        splan = splan.fillna(_DEFAULT_SPLAN_WEIGHT_NO_PRIORITIES)
    else:
        empty_fallback = int(specified.max())
        splan = splan.fillna(empty_fallback)

    df["splan_weight"] = splan.astype(int)
    return df


def pull_requests(request_urls_path):
    """
    Pull HIRES-CPS request and custom-window data from Google Sheets or local
    CPS template CSVs listed in ``request_urls_path``.

    Each row must have ``url`` (Google Sheet) and/or ``local_path`` (file).
    """
    ru = pd.read_csv(request_urls_path)
    urls_dir = os.path.dirname(os.path.abspath(request_urls_path))
    if "url" not in ru.columns and "local_path" not in ru.columns:
        raise ValueError("request_urls.csv must include 'url' and/or 'local_path'.")

    request_dfs = []
    custom_dfs = []
    for _, row in ru.iterrows():
        source = ""
        if "local_path" in ru.columns and pd.notna(row.get("local_path")):
            source = str(row["local_path"]).strip()
            if source and not os.path.isabs(source):
                source = os.path.join(urls_dir, source)
        if not source and "url" in ru.columns and pd.notna(row.get("url")):
            source = str(row["url"]).strip()
        if not source:
            continue
        df = _load_request_tab(source)
        if df.empty:
            continue

        c = SkyCoord(
            ra=df["ra"].astype(str),
            dec=df["dec"].astype(str),
            unit=(u.hourangle, u.deg),
        )
        df["ra"] = c.ra.deg
        df["dec"] = c.dec.deg

        request_dfs.append(df[SHEET_REQUEST_COLS])
        custom_dfs.append(_customs_from_requests_df(df))

    if not request_dfs:
        raise ValueError("No requests found.")

    requests_df = pd.concat(request_dfs, ignore_index=True)
    custom_df = pd.concat(custom_dfs, ignore_index=True)
    requests_df = backfill_jump_name(requests_df)
    requests_df, custom_df = _dedup_requests_by_hash(requests_df, custom_df)
    requests_df = attach_splan_weight(requests_df)
    requests_df = standardize_request_strategy(requests_df)
    return requests_df[REQUEST_COLS], custom_df


# =============================================================================
# Keck schedule — allocation.csv
# =============================================================================
# Keck tel schedule query form (HIRESr + KPF-CC nights). Crossmatched against
# request_urls.csv program codes to build allocation blocks for AstroQ.

KECK_SCHEDULE_QUERY_URL = (
    "https://www2.keck.hawaii.edu/observing/keckSchedule/queryForm.php"
)
KECK_SCHEDULE_INSTRUMENTS = ("HIRESr", "KPF-CC")


def _query_keck_schedule_form(instrument, start_date, end_date, timeout=60):
    """Query the Keck tel schedule form for a single instrument."""
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
        raise RuntimeError(
            f"Unexpected response from schedule form ({instrument}): {snippet}"
        )
    return pd.DataFrame(csv.DictReader(io.StringIO(text)))


def pull_all_scheduled(start_date, end_date, output_path=None, timeout=60):
    """Query the Keck schedule form for HIRESr and KPF-CC nights.

    Date bounds come from config ``semester_start_day`` / ``semester_end_day``.

    Columns: ``Date, Time, Dark, TelNr, Instrument, Account, PI, Institution, ProjCode``.
    If ``output_path`` is given, also write the same DataFrame to CSV.

    A failed query for one instrument (empty or malformed form response) is
    logged and skipped so other instruments can still contribute rows.
    """
    frames = []
    for inst in KECK_SCHEDULE_INSTRUMENTS:
        try:
            frame = _query_keck_schedule_form(
                inst, start_date, end_date, timeout=timeout
            )
        except RuntimeError as exc:
            logs.warning(
                "Keck schedule query for %s failed (%s); skipping.",
                inst,
                exc,
            )
            continue
        if not frame.empty:
            frames.append(frame)
    if not frames:
        df = pd.DataFrame()
    else:
        df = pd.concat(frames, ignore_index=True)
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
    scheduled_df = scheduled_df.copy()

    # 2026B_U258 (Ryan Foley): Nov 16 half-night swapped for the Dec 28 first half.
    # Out: 2026-11-16 04:32-10:07 (50%) HIRESr hires1 Ryan Foley UCSC U258
    # In:  first half of the Keck MOSFIRE slot 2026-12-28 04:43-16:04 (100%)
    #      mosfire6 Hennawi U238, split at the 10:23 midpoint
    #      → credited as HIRESr hires1 Ryan Foley UCSC U258 for CPS.
    # Half for half: the swap must leave the CPS total allocation unchanged.
    if semester == "2026B":
        drop_nov16 = (scheduled_df["Date"] == "2026-11-16") & (
            scheduled_df["ProjCode"].str.strip(",") == "U258"
        )
        if drop_nov16.any():
            logs.info(
                "Dropping %d Keck U258 Nov 16 half-night row(s) (swapped to Dec 28)",
                int(drop_nov16.sum()),
            )
            scheduled_df = scheduled_df.loc[~drop_nov16].reset_index(drop=True)

        proj = scheduled_df["ProjCode"].str.strip(",")
        dec28_u258 = (scheduled_df["Date"] == "2026-12-28") & (proj == "U258")
        if not dec28_u258.any():
            injected = pd.DataFrame(
                [
                    {
                        "Date": "2026-12-28",
                        "Time": "04:43 - 10:23 (50%)",
                        "Dark": "37",
                        "TelNr": "1",
                        "Instrument": "HIRESr",
                        "Account": "hires1",
                        "PI": "Ryan Foley",
                        "Institution": "UCSC",
                        "ProjCode": "U258",
                    }
                ]
            )
            scheduled_df = pd.concat([scheduled_df, injected], ignore_index=True)
            logs.info(
                "Injected U258 Dec 28 first-half row (MOSFIRE U238 slot → HIRES CPS)"
            )

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

    # 2026B_C362 (Yapeng Zhang): CPS queue is Aug 3 half-night only (PI confirm 2026-07-27).
    # That half night is charged to 2026B_C275 (Heather Knutson).
    if semester == "2026B":
        _C362_DROP_DATES = ("2026-08-02", "2026-08-30")
        drop = (matched["ProjCode"] == "C362") & matched["Date"].isin(_C362_DROP_DATES)
        if drop.any():
            logs.info(
                "Dropping %d Keck C362 block(s) not in CPS queue: %s",
                int(drop.sum()),
                sorted(matched.loc[drop, "Date"].unique()),
            )
            matched = matched.loc[~drop].reset_index(drop=True)

        recode = (matched["ProjCode"] == "C362") & (matched["Date"] == "2026-08-03")
        if recode.any():
            n = int(recode.sum())
            logs.info(
                "Recoding %d C362 Aug 3 half-night block(s) to C275 (2026B queue assignment)",
                n,
            )
            matched.loc[recode, "ProjCode"] = "C275"

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


def jump_query_to_past(data, request_csv_path, semester_start_day=None):
    """Convert raw JUMP frame rows to ``past.csv`` visit rows.

    JUMP ``starname`` values are matched to ``request.csv`` via ``jump_name``
    (defaulting to ``unique_id`` when blank). One JUMP visit fans out to every
    request row sharing that ``jump_name``. Template frames (B1/B3 decker,
    iodine out) get a ``_t`` suffix before matching.
    """
    cols = ["unique_id", "target", "timestamp", "exposure_time"]
    if data.empty:
        return pd.DataFrame(columns=cols)

    df = data.rename(columns={"starname": "target"}).copy()
    df["jump_name"] = df["target"].astype(str)
    is_template = (
        df["decker"].astype(str).str.upper().isin(["B1", "B3"])
        & (df["iodine_in"] == False)
    )
    df.loc[is_template, "jump_name"] += "_t"

    req = pd.read_csv(request_csv_path)
    req["unique_id"] = req["unique_id"].astype(str)
    req = backfill_jump_name(req)
    n_exp_by_uid = req.groupby("unique_id")["n_exp"].max().astype(int).to_dict()

    jump_canon = {}
    jump_to_uids = {}
    for _, row in req.iterrows():
        jn = str(row["jump_name"]).strip()
        uid = str(row["unique_id"]).strip()
        jump_canon[jn.lower()] = jn
        jump_to_uids.setdefault(jn.lower(), [])
        if uid not in jump_to_uids[jn.lower()]:
            jump_to_uids[jn.lower()].append(uid)

    def _resolve_jump(name):
        key = str(name).lower()
        return jump_canon.get(key, str(name))

    df["jump_key"] = df["jump_name"].map(_resolve_jump)
    renames = {
        jn: jk
        for jn, jk in zip(df["jump_name"], df["jump_key"])
        if str(jn) != str(jk)
    }
    if renames:
        print("JUMP rename (jump_name -> request jump_name):")
        for jn, jk in sorted(renames.items()):
            print(f"  {jn} -> {jk}")
        n_lines = int((df["jump_name"] != df["jump_key"]).sum())
        print(f"JUMP rename: {n_lines} frame line(s) where jump_name != jump_key")
    else:
        print("JUMP rename: no jump_name required renaming")

    df["exposure_time"] = (
        pd.to_numeric(df["exposure_time"], errors="coerce").fillna(0).astype(int)
    )
    df["_ts"] = pd.to_datetime(df["timestamp"], errors="coerce")
    df = df.loc[df["_ts"].notna()].copy()
    if df.empty:
        return pd.DataFrame(columns=cols)
    df["_night"] = df["_ts"].dt.strftime("%Y-%m-%d")

    if semester_start_day:
        n_before = len(df)
        df = df.loc[df["_night"] >= semester_start_day].copy()
        n_dropped = n_before - len(df)
        if n_dropped:
            print(
                f"JUMP past filter: dropped {n_dropped} frame(s) before "
                f"semester_start_day={semester_start_day}"
            )
        if df.empty:
            return pd.DataFrame(columns=cols)

    rows, rejected = [], []
    for (jump_key, night), grp in df.groupby(["jump_key", "_night"], sort=False):
        uids = jump_to_uids.get(str(jump_key).lower(), [str(jump_key)])
        first = grp.sort_values("_ts").iloc[0]
        total_exptime = int(grp["exposure_time"].sum())
        for uid in uids:
            n_exp = n_exp_by_uid.get(uid, 1)
            min_frames = math.ceil(0.5 * n_exp)
            if len(grp) < min_frames:
                rejected.append(
                    f"  {uid} {night}: {len(grp)}/{n_exp} frames (need {min_frames})"
                )
                continue
            rows.append(
                {
                    "unique_id": uid,
                    "target": uid,
                    "timestamp": first["timestamp"],
                    "exposure_time": total_exptime,
                }
            )

    summary = (
        f"{len(df)} raw frame(s) -> {len(rows)} accounted visit(s) "
        f"({len(rejected)} group(s) below 50% threshold)"
    )
    if rejected:
        summary += ":\n" + "\n".join(rejected)
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

    session = jump_client.jump_session()
    param_str = (
        f"start_date:{semester_start_day}|end_date:{semester_end_day}"
    )
    param_slug = urllib.parse.quote(param_str, safe="")
    download_url = (
        f"{jump_client.JUMP_BASE_URL}/explorer/{jump_client.JUMP_HIRES_PAST_EXPLORER_ID}/download"
        f"?format=csv&params={param_slug}"
    )
    print(f"JUMP download URL: {download_url}")

    csv_response = session.get(download_url, timeout=60)
    csv_response.raise_for_status()
    with open(raw_path, "wb") as f:
        f.write(csv_response.content)
    print(f"Raw JUMP pull saved to {raw_path}")

    data = pd.read_csv(raw_path)

    visits = jump_query_to_past(
        data, request_csv_path, semester_start_day=semester_start_day
    )
    visits.to_csv(path_to_csv, index=False)
    print(f"Processed past history saved to {path_to_csv}")
