"""CSV input contracts and validated reads for AstroQ semester inputs."""

import logging
import os

import pandas as pd
from astropy.time import Time

logs = logging.getLogger(__name__)

REQUEST_SCHEMA = {
    "unique_id": str,
    "target": str,
    "program_code": str,
    "ra": float,            # deg
    "dec": float,           # deg
    "exptime": float,       # seconds
    "n_exp": int,
    "n_inter_max": int,
    "tau_inter": int,       # days
    "n_intra_min": int,
    "n_intra_max": int,
    "tau_intra": float,     # hours
    "inactive": bool,
    "splan_weight": float,
}

PAST_SCHEMA = {
    "unique_id": str,
    "target": str,
    "timestamp": str,        # UT ISO
    "exposure_time": float,  # seconds
}

ALLOCATION_SCHEMA = {"start": Time, "stop": Time}

CUSTOM_SCHEMA = {"unique_id": str, "target": str, "start": Time, "stop": Time}

PROGRAMS_SCHEMA = {"program": str, "hours": float}

DEFAULT_MIN_FILLFACTOR = 0.0
DEFAULT_MAX_FILLFACTOR = 1.25

REQUEST_COLS = list(REQUEST_SCHEMA)
PAST_COLS = list(PAST_SCHEMA)



def read_csv(path, type):
    """Read ``path`` and validate/coerce against the schema for ``type``.

    Args:
        path (str): CSV location.
        type (str): one of ``"request"``, ``"past"``, ``"programs"``,
            ``"allocation"``, ``"custom"``.

    Returns:
        pandas.DataFrame with every schema column coerced; extra columns pass through.

    Raises:
        ValueError: unknown ``type``, missing columns, nulls, bad dtypes, duplicates.
        FileNotFoundError: missing file when the type requires one.
    """
    if type == "past":
        df = _load_frame(path, PAST_SCHEMA, "past.csv", empty_ok=True)

    elif type == "request":
        df = _load_frame(path, REQUEST_SCHEMA, "request.csv")
        active = df[~df["inactive"]]
        dup_mask = active["unique_id"].duplicated(keep=False)
        if dup_mask.any():
            dup_ids = sorted(active.loc[dup_mask, "unique_id"].unique())
            raise ValueError(
                f"Duplicate unique_id among active requests: {dup_ids}. "
                f"Remove or merge duplicate rows so each active request has one row."
            )

    elif type == "programs":
        df = _load_frame(path, PROGRAMS_SCHEMA, "programs.csv", key="program")
        if "min_fillfactor" not in df.columns:
            df["min_fillfactor"] = DEFAULT_MIN_FILLFACTOR
        else:
            df["min_fillfactor"] = (
                pd.to_numeric(df["min_fillfactor"], errors="coerce")
                .fillna(DEFAULT_MIN_FILLFACTOR)
                .astype(float)
            )
        if "max_fillfactor" not in df.columns:
            df["max_fillfactor"] = DEFAULT_MAX_FILLFACTOR
        else:
            df["max_fillfactor"] = (
                pd.to_numeric(df["max_fillfactor"], errors="coerce")
                .fillna(DEFAULT_MAX_FILLFACTOR)
                .astype(float)
            )
        df = df.set_index("program")

    elif type == "allocation":
        df = _load_frame(path, ALLOCATION_SCHEMA, "allocation.csv")

    elif type == "custom":
        df = _load_frame(path, CUSTOM_SCHEMA, "custom.csv", empty_ok=True)

    else:
        raise ValueError(f"read_csv type {type!r} not implemented")

    return df


def _load_frame(path, schema, name, *, empty_ok=False, key=None):
    """Read ``path`` and validate/coerce against ``schema``.

    Args:
        path (str): CSV location.
        schema (dict): column -> target dtype.
        name (str): label used in error messages (e.g. ``"request.csv"``).
        empty_ok (bool): missing/zero-byte/header-only files return an empty
            frame with the schema's columns instead of raising.
        key (str, optional): column whose values must be unique.

    Returns:
        pandas.DataFrame with every schema column coerced; extra columns pass through.

    Raises:
        FileNotFoundError: missing file when ``empty_ok=False``.
        ValueError: missing columns, null values, bad dtypes, or duplicate ``key``.
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        if empty_ok:
            return pd.DataFrame({c: pd.Series(dtype=object) for c in schema})
        raise FileNotFoundError(f"{name} not found: {path}")
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        if empty_ok:
            return pd.DataFrame({c: pd.Series(dtype=object) for c in schema})
        raise

    missing = [c for c in schema if c not in df.columns]
    if missing:
        raise ValueError(f"{name} missing required column(s): {missing}")

    nulls = [c for c in schema if df[c].isna().any()]
    if nulls:
        raise ValueError(
            f"{name} has null values in required column(s): {nulls}. "
            f"Inputs must arrive clean from the prep stage."
        )

    for col, dtype in schema.items():
        if dtype is bool:
            if df[col].dtype != bool:
                raise ValueError(
                    f"{name} column {col!r} must be boolean (True/False); "
                    f"got dtype {df[col].dtype}."
                )
        elif dtype is int:
            vals = pd.to_numeric(df[col])
            if (vals % 1 != 0).any():
                raise ValueError(f"{name} column {col!r} has non-integer values.")
            df[col] = vals.astype(int)
        elif dtype is float:
            df[col] = pd.to_numeric(df[col]).astype(float)
        elif dtype is Time:
            df[col] = df[col].apply(Time)
        else:
            df[col] = df[col].astype(str)

    if key is not None:
        dup = df[key].duplicated(keep=False)
        if dup.any():
            dup_vals = sorted(df.loc[dup, key].unique())
            raise ValueError(f"{name} has duplicate {key!r} values: {dup_vals}")

    return df.reset_index(drop=True)
