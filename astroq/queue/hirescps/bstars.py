"""B-star calibration block for HIRES nightly scripts.

Resolves a name-only reference list via SIMBAD, caches catalog fields on disk,
and formats propagated coordinates for the MAGIQ script footer.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

from astroq.queue.hirescps.script_columns import (
    format_exposure_token,
    format_section_header,
    format_vmag_token,
)

logs = logging.getLogger(__name__)

_BSTARS_FILE = Path(__file__).with_name("bstars.txt")
_CACHE_FILENAME = "bstars_simbad.csv"
_BSTARS_HEADER = format_section_header("B-Stars-Coordinates-Advanced")
_BSTARS_NAME_WIDTH = 17
_CACHE_COLS = [
    "name",
    "ra_deg",
    "dec_deg",
    "pmra",
    "pmdec",
    "vmag",
    "sp_type",
    "rot_vel",
    "simbad_id",
]
_SIMBAD_BATCH = 25


def _load_bstar_names(path: Path | None = None) -> pd.DataFrame:
    """Read ``bstars.txt``; return DataFrame with ``name`` and optional ``comment``."""
    path = path or _BSTARS_FILE
    names, comments = [], []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split("\t", 1)
        name = parts[0].strip().lower()
        comment = parts[1].strip() if len(parts) > 1 else ""
        names.append(name)
        comments.append(comment)
    return pd.DataFrame({"name": names, "comment": comments})


def _simbad_query_id(name: str) -> str:
    """Map script name ``hr9098`` to SIMBAD identifier ``HR 9098``."""
    m = re.fullmatch(r"hr(\d+)", name.lower())
    if not m:
        return name
    return f"HR {int(m.group(1))}"


def _name_from_simbad_id(simbad_id: str) -> str:
    m = re.search(r"HR\s*(\d+)", str(simbad_id), flags=re.IGNORECASE)
    if not m:
        return str(simbad_id).strip().lower()
    return f"hr{int(m.group(1))}"


def _float_or_zero(val) -> float:
    try:
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return 0.0
        return float(val)
    except (TypeError, ValueError):
        return 0.0


def _float_or_nan(val) -> float:
    try:
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return float("nan")
        return float(val)
    except (TypeError, ValueError):
        return float("nan")


def _simbad_table_to_rows(table, requested_names: list[str]) -> list[dict]:
    """Collapse SIMBAD mesRot join duplicates to one row per requested HR name."""
    if table is None or len(table) == 0:
        return []

    df = table.to_pandas()
    id_col = "user_specified_id" if "user_specified_id" in df.columns else "matched_id"
    if "mesrot.mespos" in df.columns:
        df = df.sort_values("mesrot.mespos", kind="mergesort")
    df = df.groupby(id_col, as_index=False).first()

    rows_by_name = {}
    for _, row in df.iterrows():
        name = _name_from_simbad_id(row[id_col])
        vsini_col = "mesrot.vsini" if "mesrot.vsini" in row.index else None
        rows_by_name[name] = {
            "name": name,
            "ra_deg": float(row["ra"]),
            "dec_deg": float(row["dec"]),
            "pmra": _float_or_zero(row.get("pmra")),
            "pmdec": _float_or_zero(row.get("pmdec")),
            "vmag": _float_or_nan(row.get("V")),
            "sp_type": str(row.get("sp_type", "") or "").strip(),
            "rot_vel": _float_or_nan(row[vsini_col]) if vsini_col else float("nan"),
            "simbad_id": str(row.get("main_id", "") or "").strip(),
        }

    out = []
    for name in requested_names:
        out.append(rows_by_name.get(name, {"name": name}))
    return out


def _query_simbad(names: list[str]) -> pd.DataFrame:
    """Batch SIMBAD lookup for HR star names."""
    from astroquery.simbad import Simbad

    simbad = Simbad()
    simbad.reset_votable_fields()
    simbad.add_votable_fields("ra", "dec", "pmra", "pmdec", "V", "sp_type", "mesRot", "main_id")

    rows = []
    ids = [_simbad_query_id(n) for n in names]
    for start in range(0, len(ids), _SIMBAD_BATCH):
        batch_ids = ids[start : start + _SIMBAD_BATCH]
        batch_names = names[start : start + _SIMBAD_BATCH]
        logs.info("Querying SIMBAD for %d B-star(s).", len(batch_ids))
        table = simbad.query_objects(batch_ids)
        rows.extend(_simbad_table_to_rows(table, batch_names))

    return pd.DataFrame(rows)


def _read_cache(cache_path: Path) -> pd.DataFrame:
    if not cache_path.is_file():
        return pd.DataFrame(columns=_CACHE_COLS)
    df = pd.read_csv(cache_path)
    for col in _CACHE_COLS:
        if col not in df.columns:
            df[col] = np.nan if col in ("vmag", "rot_vel", "ra_deg", "dec_deg") else ""
    return df[_CACHE_COLS]


def _write_cache(cache_path: Path, df: pd.DataFrame) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    df[_CACHE_COLS].to_csv(cache_path, index=False)


def _resolve_bstars(names_df: pd.DataFrame, cache_dir: str | Path) -> pd.DataFrame:
    """Return cached + freshly queried SIMBAD rows for all names."""
    cache_path = Path(cache_dir) / _CACHE_FILENAME
    cache = _read_cache(cache_path)
    cached_names = set(cache["name"].astype(str)) if not cache.empty else set()
    missing = [n for n in names_df["name"].astype(str) if n not in cached_names]
    if missing:
        fresh = _query_simbad(missing)
        cache = pd.concat([cache, fresh], ignore_index=True)
        _write_cache(cache_path, cache)
        logs.info("Wrote %d B-star row(s) to %s", len(fresh), cache_path)
    else:
        logs.info("Using cached B-star catalog (%s).", cache_path)

    merged = names_df.merge(cache, on="name", how="left")
    missing_rows = merged[merged["ra_deg"].isna()]
    if not missing_rows.empty:
        logs.warning(
            "B-stars missing from cache/SIMBAD: %s",
            ", ".join(missing_rows["name"].astype(str)),
        )
    return merged


def _format_coord_strings(coord: SkyCoord) -> tuple[str, str]:
    ra = coord.ra.to_string(unit=u.hourangle, sep=" ", pad=True, precision=1)
    dec = coord.dec.to_string(unit=u.deg, sep=" ", pad=True, precision=0)
    if dec[0] != "-":
        dec = "+" + dec
    return ra, dec


def propagate_bstar_coord(
    ra_deg: float,
    dec_deg: float,
    pmra: float,
    pmdec: float,
    current_day: str,
) -> SkyCoord:
    """Propagate J2000 coordinates to ``current_day`` using proper motion."""
    coord = SkyCoord(
        ra=ra_deg * u.deg,
        dec=dec_deg * u.deg,
        pm_ra_cosdec=pmra * u.mas / u.yr,
        pm_dec=pmdec * u.mas / u.yr,
        obstime=Time("J2000"),
    )
    return coord.apply_space_motion(new_obstime=Time(current_day))


def _display_name(name: str) -> str:
    """Map ``hr9098`` to script name ``HR9098``."""
    m = re.fullmatch(r"hr(\d+)", name.lower())
    if m:
        return f"HR{m.group(1)}"
    return name.upper()


def _format_vmag_token(vmag: float) -> str:
    """Match ``format_hires_row`` vmag padding."""
    return format_vmag_token(vmag)


def _format_bstar_exposure() -> str:
    """Fixed B-star exposure token matching ``format_hires_row`` spacing."""
    return format_exposure_token(5, 500)


def format_bstar_row(
    name: str,
    coord: SkyCoord,
    vmag: float,
    pmra: float,
    pmdec: float,
    current_day: str,
    *,
    comment: str = "",
) -> str:
    """Format one B-star line like a request row with ``XX`` instead of program."""
    ra_str, dec_str = _format_coord_strings(coord)
    display_name = _display_name(name)
    name_str = " " * (_BSTARS_NAME_WIDTH - len(display_name[:_BSTARS_NAME_WIDTH])) + display_name[:_BSTARS_NAME_WIDTH]
    line = (
        f"{name_str} {ra_str} {dec_str} 2000 {_format_vmag_token(vmag)}"
        f"{_format_bstar_exposure()} 250k B5 1x  in p3 XX"
    )
    line += f"  epoch={Time(current_day).jyear:.1f}"
    if comment:
        line += f", {comment}"
    return line


def build_bstars_section(
    current_day: str,
    cache_dir: str | Path,
    *,
    names_path: Path | None = None,
    query_fn=None,
) -> list[str]:
    """Build header + B-star rows for append to the nightly script."""
    names_df = _load_bstar_names(names_path)
    resolve = query_fn or _resolve_bstars
    table = resolve(names_df, cache_dir)

    lines = [_BSTARS_HEADER]
    for _, row in table.iterrows():
        if not np.isfinite(row.get("ra_deg", np.nan)):
            continue
        pmra = float(row.get("pmra", 0.0) or 0.0)
        pmdec = float(row.get("pmdec", 0.0) or 0.0)
        coord = propagate_bstar_coord(
            float(row["ra_deg"]),
            float(row["dec_deg"]),
            pmra,
            pmdec,
            current_day,
        )
        lines.append(
            format_bstar_row(
                str(row["name"]),
                coord,
                float(row.get("vmag", np.nan)),
                pmra,
                pmdec,
                current_day,
                comment=str(row.get("comment", "") or ""),
            )
        )
    return lines
