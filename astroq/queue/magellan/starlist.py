"""Minimal Magellan nightly starlist writer (placeholder format)."""

import os

import pandas as pd


def write_starlist(
    frame,
    schedule,
    night_start_time,
    filler_stars,
    current_day,
    outputdir,
    version="nominal",
    all_active_requests=None,
    **kwargs,
):
    """Write a simple CSV starlist from the TTP schedule.

    This is a placeholder until a Magellan-specific observer format is defined.
    """
    del filler_stars, all_active_requests, kwargs

    os.makedirs(outputdir, exist_ok=True)
    out_path = os.path.join(outputdir, f"starlist_{current_day}_{version}.csv")

    rows = []
    for i, row in frame.iterrows():
        start_time = ""
        if schedule is not None and i < len(schedule):
            start_time = schedule.iloc[i].get("start_time", "")
        rows.append(
            {
                "unique_id": row.get("unique_id", ""),
                "target": row.get("target", ""),
                "start_time": start_time,
            }
        )

    pd.DataFrame(rows).to_csv(out_path, index=False)
    return out_path
