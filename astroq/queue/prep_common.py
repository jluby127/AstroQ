"""Shared request.csv standardization used by every ``<queue>.prep`` module.

The semester planner (:mod:`astroq.splan`) expects a clean ``request.csv``: it
validates the schema and does not repair values. Producers are responsible for
supplying valid strategy fields, so both the HIRES-CPS and KPF-CC prep paths
funnel their finished request frame through :func:`standardize_request_strategy`
before writing it to disk.
"""

import os

import numpy as np

# Intra-night strategy columns and the value to use when the source is blank
# or the legacy early-HIRES-CPS webform emitted the literal string ``"None"``.
_INTRA_DEFAULTS = (
    ("n_intra_max", 1),
    ("n_intra_min", 1),
    ("tau_intra", 0),
)


def write_programs_csv(programs_df, savepath):
    """Write ``programs.csv`` sorted by program with hours rounded to 0.01."""
    out = programs_df.copy()
    if "hours" in out.columns:
        out["hours"] = out["hours"].round(2)
    out = out.sort_values("program", kind="mergesort").reset_index(drop=True)
    out.to_csv(os.path.join(savepath, "programs.csv"), index=False)
    return out


def standardize_request_strategy(df):
    """Fill intra-night strategy defaults on a request frame in place.

    Tolerates legacy ``"None"`` strings and missing values, applying
    ``n_intra_max -> 1``, ``n_intra_min -> 1``, ``tau_intra -> 0``. No-op on
    already-clean inputs. Mutates and returns ``df``.
    """
    for col, default in _INTRA_DEFAULTS:
        df[col] = df[col].replace("None", np.nan).fillna(default)
    return df
