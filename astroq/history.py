"""
Module for processing the past.csv observation history into per-target
``StarHistory`` records consumed by the semester planner. The instrument-
specific code that builds ``past.csv`` lives in ``astroq.queue.*``.
"""

# Standard library imports
import os
from collections import namedtuple

# Third-party imports
import pandas as pd
from astropy.time import Time, TimeDelta

# Global timezone offset: hours difference between local time and UT
# Positive values mean local time is ahead of UT (e.g., +8 for UTC+8)
# Negative values mean local time is behind UT (e.g., -10 for UTC-10)
TIMEZONE_OFFSET_HOURS = -10  # Adjust this value for your timezone

# Define StarHistory namedtuple at module level for proper pickling
StarHistory = namedtuple(
    "StarHistory",
    [
        "name",
        "date_last_observed",
        "total_n_exposures",
        "total_n_visits",
        "total_n_unique_nights",
        "total_open_shutter_time",
        "n_obs_on_nights",
        "n_visits_on_nights",
    ],
)


def process_star_history(filename):
    """
    Process the past.csv file to return a dict of star histories.

    Args:
        filename (str): the path to the past.csv file

    Returns:
        result (dict): a dictionary of star histories where keys are 'id', values are objects with attributes: date_last_observed, total_n_exposures, total_n_visits, total_n_unique_nights, total_open_shutter_time
    """

    # Check if file is empty or doesn't exist
    if not os.path.exists(filename) or os.path.getsize(filename) == 0:
        # Create empty DataFrame with expected columns
        df = pd.DataFrame(
            columns=[
                "id",
                "target",
                "semid",
                "timestamp",
                "exposure_start_time",
                "exposure_time",
                "observer",
            ]
        )
    else:
        try:
            df = pd.read_csv(filename)
        except pd.errors.EmptyDataError:
            # Create empty DataFrame with expected columns
            df = pd.DataFrame(
                columns=[
                    "id",
                    "target",
                    "semid",
                    "timestamp",
                    "exposure_start_time",
                    "exposure_time",
                    "observer",
                ]
            )
    result = {}

    # Group by target (star)
    for unique_id, star_df in df.groupby("id"):
        starname = df[df["id"] == unique_id]["target"].values[0]
        # Group by visit (timestamp)
        visits = []
        for ts, visit_df in star_df.groupby("timestamp"):
            # Each visit is a DataFrame of exposures
            exposure_start_time = list(visit_df["exposure_start_time"])
            exposure_time = list(visit_df["exposure_time"])
            junk = list(visit_df["junk"]) if "junk" in visit_df.columns else []
            # Convert junk to bool if needed
            junk = [bool(j) for j in junk]
            # Apply junk logic: skip if >= half are True
            if junk and sum(junk) >= len(junk) / 2:
                continue
            visits.append(
                {
                    "id": visit_df["id"].iloc[0],
                    "exposure_start_time": exposure_start_time,
                    "exposure_time": exposure_time,
                    "timestamp": ts,
                }
            )
        if not visits:
            continue
        star_id = visits[0]["id"]
        all_times = [Time(t) for v in visits for t in v["exposure_start_time"]]
        # Adjust UT time to local timezone
        ut_time = max(all_times) if all_times else None
        if ut_time:
            # Subtract timezone offset (negative offset means local time is behind UT)
            local_time = ut_time - TimeDelta(
                abs(TIMEZONE_OFFSET_HOURS) / 24, format="jd"
            )
            date_last_observed = local_time.isot[:10]
        else:
            date_last_observed = ""
        total_n_exposures = sum(len(v["exposure_start_time"]) for v in visits)
        total_n_visits = len(visits)
        unique_nights = set(
            Time(t).isot[:10] for v in visits for t in v["exposure_start_time"]
        )
        total_n_unique_nights = len(unique_nights)
        total_open_shutter_time = int(
            round(sum(sum(map(float, v["exposure_time"])) for v in visits))
        )
        # Build n_obs_on_nights and n_visits_on_nights
        n_obs_on_nights = {}
        n_visits_on_nights = {}
        for v in visits:
            visit_date = v["timestamp"][:10]
            n_obs_on_nights[visit_date] = n_obs_on_nights.get(visit_date, 0) + len(
                v["exposure_start_time"]
            )
            n_visits_on_nights[visit_date] = n_visits_on_nights.get(visit_date, 0) + 1
        result[unique_id] = StarHistory(
            name=starname,
            date_last_observed=date_last_observed,
            total_n_exposures=total_n_exposures,
            total_n_visits=total_n_visits,
            total_n_unique_nights=total_n_unique_nights,
            total_open_shutter_time=total_open_shutter_time,
            n_obs_on_nights=n_obs_on_nights,
            n_visits_on_nights=n_visits_on_nights,
        )
    return result
