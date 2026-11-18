"""Authenticated access to JUMP.

jump2 (cannon.caltech.edu, via the optional ``jump2-client`` package) supplies
the HIRES frame history behind ``past.csv``. The jump1 helpers
(django-sql-explorer saved queries and ad-hoc SQL) remain for one-off scripts
such as ``ops/HIRES/2026B/select_recommissioning.py``.
"""

from __future__ import annotations

import io
import os
import re
import urllib.parse

import pandas as pd
import requests

JUMP_BASE_URL = "https://jump.caltech.edu"
JUMP_PAST_QUERY_TMP_FILENAME = "past_jump2-query-tmp.csv"

JUMP2_BASE_URL = "https://cannon.caltech.edu"
JUMP2_INSTALL_HINT = (
    "pip install -U --extra-index-url https://cannon.caltech.edu/client/simple/ "
    "jump2-client  (or run `bash ops/HIRES/setup_env.sh`)"
)
JUMP2_AUTH_HINT = (
    "jump2 API key missing or expired. Create one under Profile -> API Keys on "
    "https://cannon.caltech.edu and set JUMP_API_KEY, or run `jump2 scripts login`."
)
# core_observation.instrument value for HIRES frames (the j chip carries the
# science and template exposures used by the planner).
JUMP2_HIRES_INSTRUMENT = "hires_j"
# Row cap of the SQL endpoint used by jump2-client < 0.17, which truncates
# without raising. Frame queries are batched by night to stay under it.
JUMP2_SQL_ROW_CAP = 1000


def jump_session() -> requests.Session:
    """Return an authenticated ``requests.Session`` for JUMP."""
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
            "JUMP login appears to have failed: still on login page after submitting "
            "credentials. Check KPFCC_JUMP_USERNAME/KPFCC_JUMP_PASSWORD."
        )
    return session


def run_saved_query(
    session: requests.Session,
    query_id: int,
    params: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Download a saved explorer query as a CSV ``DataFrame``."""
    url = f"{JUMP_BASE_URL}/explorer/{query_id}/download?format=csv"
    if params:
        param_str = "|".join(f"{k}:{v}" for k, v in params.items())
        url += f"&params={urllib.parse.quote(param_str, safe='')}"
    resp = session.get(url, timeout=120)
    resp.raise_for_status()
    return pd.read_csv(io.StringIO(resp.text))


def run_sql(session: requests.Session, sql: str) -> pd.DataFrame:
    """Execute ad-hoc read-only SQL via ``/explorer/play/`` and parse the HTML table."""
    play_url = f"{JUMP_BASE_URL}/explorer/play/"
    session.get(play_url, timeout=60)
    resp = session.post(
        play_url,
        data={
            "sql": sql,
            "csrfmiddlewaretoken": session.cookies["csrftoken"],
            "title": "",
            "connection": "",
        },
        headers={"Referer": play_url},
        timeout=300,
    )
    resp.raise_for_status()
    headers = [
        re.sub(r"<[^>]+>", "", cell).strip()
        for cell in re.findall(r"<th[^>]*>(.*?)</th>", resp.text, re.S)
    ]
    headers = [h for h in headers if h and not h.startswith("Sum:")]
    rows = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", resp.text, re.S):
        cells = [
            re.sub(r"<[^>]+>", "", cell).strip()
            for cell in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)
        ]
        if cells:
            rows.append(cells[1:] if len(cells) > 1 and cells[0].isdigit() else cells)
    if not rows:
        return pd.DataFrame()
    ncols = max(len(r) for r in rows)
    if headers and len(headers) >= ncols:
        col_names = headers[:ncols]
    else:
        col_names = [f"col{i}" for i in range(ncols)]
    padded = [r + [""] * (ncols - len(r)) for r in rows]
    return pd.DataFrame(padded, columns=col_names)


# =============================================================================
# jump2 — HIRES frame history for past.csv
# =============================================================================


def jump2_client():
    """Return a ``jump2_client.Client`` for cannon.caltech.edu.

    Credentials come from the client's own lookup: ``JUMP_API_KEY``, else
    ``~/.config/jump2/credentials.toml`` (written by ``jump2 scripts login``).
    ``JUMP_BASE_URL`` overrides the server. ``jump2-client`` is only needed by
    the HIRES queue, so it is imported here rather than at module load.
    """
    try:
        from jump2_client import Client
    except ImportError as exc:
        raise ImportError(
            f"jump2-client is required to pull HIRES history. Install it with: "
            f"{JUMP2_INSTALL_HINT}"
        ) from exc
    return Client(base_url=os.environ.get("JUMP_BASE_URL") or JUMP2_BASE_URL)


def _jump2_query(client, sql: str) -> pd.DataFrame:
    """Run ``sql`` through ``client.data.query``, turning auth failures into a
    renew-your-key message."""
    try:
        from jump2_client import AuthError
    except ImportError:  # a test double stands in for the client
        AuthError = ()
    try:
        return client.data.query(sql)
    except AuthError as exc:
        raise RuntimeError(f"{JUMP2_AUTH_HINT} ({exc})") from exc


def hires_frames(client, start_day: str, stop_day: str) -> pd.DataFrame:
    """HIRES science frames with ``start_day <= utctime < stop_day`` (UTC).

    Junk frames and calibration frames (anything not tagged ``stars``) are
    dropped in the query. Frames without exposure-meter counts are kept.

    Returns the columns ``jump_query_to_past`` expects: ``starname, timestamp,
    exposure_time, decker, iodine_in, counts``. ``timestamp`` is jump2's
    mid-exposure UTC time as ``YYYY-MM-DD HH:MM``.

    A per-night count comes first; frames are then fetched in batches of whole
    nights small enough for the 1,000-row cap of older clients. A batch that
    still reaches the cap raises rather than returning a truncated history.
    """
    def where(lo, hi):
        return f"""
WHERE core_observation.instrument ILIKE '{JUMP2_HIRES_INSTRUMENT}'
  AND core_observation.utctime >= TIMESTAMPTZ '{lo} 00:00+00'
  AND core_observation.utctime < TIMESTAMPTZ '{hi} 00:00+00'
  AND core_observation.isjunk IS NOT TRUE
  -- Science frames only: drops bias/flat/ThAr/iodine/focus frames, including
  -- calibrations taken under the last target's name.
  AND 'stars' = ANY(core_observation.exposure_types)"""

    from_ = """
FROM core_star
INNER JOIN core_observation
  ON core_star.id = core_observation.star_id"""

    counts = _jump2_query(
        client,
        "SELECT (core_observation.utctime AT TIME ZONE 'UTC')::date AS night,"
        f"\n  COUNT(DISTINCT core_observation.id) AS n{from_}{where(start_day, stop_day)}"
        "\nGROUP BY 1\nORDER BY 1",
    )

    columns = ["starname", "timestamp", "exposure_time", "decker", "iodine_in", "counts"]
    if counts.empty:
        return pd.DataFrame(columns=columns)

    nights = [str(d)[:10] for d in counts["night"]]
    sizes = pd.to_numeric(counts["n"]).astype(int).tolist()
    batches, batch, total = [], [], 0
    for d, n in zip(nights, sizes):
        if batch and total + n >= JUMP2_SQL_ROW_CAP:
            batches.append(batch)
            batch, total = [], 0
        batch.append(d)
        total += n
    batches.append(batch)

    frames = []
    for batch in batches:
        lo = batch[0]
        hi = (pd.Timestamp(batch[-1]) + pd.Timedelta(days=1)).strftime("%Y-%m-%d")
        sql = f"""SELECT DISTINCT ON (core_observation.id)
  core_star.name AS starname,
  TO_CHAR(core_observation.utctime AT TIME ZONE 'UTC', 'YYYY-MM-DD HH24:MI') AS timestamp,
  core_observation.hdr_exposure_time AS exposure_time,
  core_observation.hdr_decker AS decker,
  core_observation.iodine_in AS iodine_in,
  core_observation.hdr_counts AS counts{from_}{where(lo, hi)}
ORDER BY core_observation.id"""
        df = _jump2_query(client, sql)
        if len(df) >= JUMP2_SQL_ROW_CAP:
            raise RuntimeError(
                f"jump2 returned {len(df)} HIRES frames for {lo}..{batch[-1]}, at the "
                f"{JUMP2_SQL_ROW_CAP}-row cap; the history may be truncated."
            )
        frames.append(df)
    return pd.concat(frames, ignore_index=True)[columns]
