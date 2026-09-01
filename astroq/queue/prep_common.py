"""Shared request.csv standardization used by every ``<queue>.prep`` module.

The semester planner (:mod:`astroq.splan`) expects a clean ``request.csv``: it
validates the schema and does not repair values. Producers are responsible for
supplying valid strategy fields, so both the HIRES-CPS and KPF-CC prep paths
funnel their finished request frame through :func:`standardize_request_strategy`
before writing it to disk.
"""

import os

import numpy as np
import pandas as pd

import astroq.io

# Intra-night strategy columns and the value to use when the source is blank
# or the legacy early-HIRES-CPS webform emitted the literal string ``"None"``.
_INTRA_DEFAULTS = (
    ("n_intra_max", 1),
    ("n_intra_min", 1),
    ("tau_intra", 0),
)

PROGRAMS_COLS = ("program", "hours", "min_fill", "max_fill", "max_feasible_fill")


def write_programs_csv(programs_df, savepath):
    """Write ``programs.csv`` sorted by program with hours rounded to 0.01.

    Emits the full column set so operators can hand-edit the fill bounds:
    ``min_fill`` and ``max_fill`` default to the :mod:`astroq.io` values, and
    ``max_feasible_fill`` is left empty (NaN) for ``astroq compute-max-fill``
    to populate.
    """
    out = programs_df.copy()
    if "hours" in out.columns:
        out["hours"] = out["hours"].round(2)
    for col, default in (
        ("min_fill", astroq.io.DEFAULT_MIN_FILL),
        ("max_fill", astroq.io.DEFAULT_MAX_FILL),
        ("max_feasible_fill", np.nan),
    ):
        if col not in out.columns:
            out[col] = default
    out = out.sort_values("program", kind="mergesort").reset_index(drop=True)
    extra = [c for c in out.columns if c not in PROGRAMS_COLS]
    out = out[list(PROGRAMS_COLS) + extra]
    out.to_csv(os.path.join(savepath, "programs.csv"), index=False)
    return out


def merge_program_overrides_from_request_urls(programmatics, request_urls_path):
    """Apply optional ``hours``, ``min_fill``, ``max_fill`` from ``request_urls.csv``.

    Programs listed in the URLs file but absent from Keck allocation are added
    when ``hours`` is set; otherwise they receive 0 hours (legacy behavior).
    """
    ru = pd.read_csv(request_urls_path)
    if "program_code" not in ru.columns:
        return programmatics

    out = programmatics.copy()
    fill_cols = ("hours", "min_fill", "max_fill")
    for _, row in ru.iterrows():
        prog = row.get("program_code")
        if pd.isna(prog) or str(prog).strip() == "":
            continue
        prog = str(prog).strip()
        mask = out["program"].astype(str) == prog
        if mask.any():
            for col in fill_cols:
                if col not in ru.columns or col not in out.columns:
                    continue
                val = row[col] if col in row.index else row.get(col)
                if pd.notna(val):
                    out.loc[mask, col] = float(val)
            continue

        hours = (
            float(row["hours"])
            if "hours" in ru.columns and pd.notna(row.get("hours"))
            else 0.0
        )
        new_row = {"program": prog, "hours": hours}
        for col in ("min_fill", "max_fill"):
            if col in ru.columns and pd.notna(row.get(col)):
                new_row[col] = float(row[col])
        out = pd.concat([out, pd.DataFrame([new_row])], ignore_index=True)
    return out.copy()


def standardize_request_strategy(df):
    """Fill intra-night strategy defaults on a request frame in place.

    Tolerates legacy ``"None"`` strings and missing values, applying
    ``n_intra_max -> 1``, ``n_intra_min -> 1``, ``tau_intra -> 0``. No-op on
    already-clean inputs. Mutates and returns ``df``.
    """
    for col, default in _INTRA_DEFAULTS:
        df[col] = df[col].replace("None", np.nan).fillna(default)
    return df
