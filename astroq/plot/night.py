"""Night-plan table prep (script plan DataFrame for HTML tables)."""

import os

import pandas as pd
from astropy.time import TimeDelta

from astroq.nplan import get_nightly_times_from_allocation
from astroq.plot.tables import NIGHTPLAN_COLUMNS


def get_script_plan(night_planner):
    """Generate script plan DataFrame from semester planner and night planner objects.

    This function reads the request_selected.csv file from the semester planner's output directory,
    merges it with the night planner's solution data, and returns a properly formatted DataFrame
    with the same column structure as the original get_script_plan function.

    Args:
        night_planner: NightPlanner object containing solution attribute

    Returns:
        final_df (pd.DataFrame): a formatted observing plan DataFrame
    """

    # Read the request_selected.csv file from the semester planner's output directory
    request_selected_path = os.path.join(
        night_planner.output_directory, "request_selected.csv"
    )

    if not os.path.exists(request_selected_path):
        raise FileNotFoundError(
            f"request_selected.csv not found at {request_selected_path}"
        )

    # Read the request_selected.csv file
    request_selected_df = pd.read_csv(request_selected_path)
    solution = night_planner.solution
    on_sky = solution.schedule[~solution.schedule["is_anchor"]]
    scheduled = on_sky[on_sky["scheduled"]].sort_values("order")

    merged_df = request_selected_df.merge(
        scheduled[["unique_id", "t_start", "t_earliest_start", "t_latest_finish"]],
        on="unique_id",
        how="inner",
    )
    merged_df = merged_df.rename(
        columns={
            "t_start": "Start Exposure",
            "t_earliest_start": "Earliest Start",
            "t_latest_finish": "Latest Finish",
        }
    )

    desired_columns = NIGHTPLAN_COLUMNS

    # Keep only the columns that exist in the merged dataframe
    available_columns = [col for col in desired_columns if col in merged_df.columns]

    # Reorder columns to match the desired structure
    final_df = merged_df[available_columns].copy()

    # Round numeric fields to appropriate decimal places
    if "ra" in final_df.columns:
        # Ensure ra is numeric before rounding, handle 'None' strings
        final_df["ra"] = final_df["ra"].replace("None", pd.NA)
        final_df["ra"] = pd.to_numeric(final_df["ra"], errors="coerce").round(1)

    if "dec" in final_df.columns:
        # Ensure dec is numeric before rounding, handle 'None' strings
        final_df["dec"] = final_df["dec"].replace("None", pd.NA)
        final_df["dec"] = pd.to_numeric(final_df["dec"], errors="coerce").round(1)

    if "jmag" in final_df.columns:
        # Ensure jmag is numeric before rounding, handle 'None' strings
        final_df["jmag"] = final_df["jmag"].replace("None", pd.NA)
        final_df["jmag"] = pd.to_numeric(final_df["jmag"], errors="coerce").round(1)

    if "Vmag" in final_df.columns:
        final_df["Vmag"] = final_df["Vmag"].replace("None", pd.NA)
        final_df["Vmag"] = pd.to_numeric(final_df["Vmag"], errors="coerce").round(1)

    # Convert time fields from "minutes from start of night" to HST timestamps
    try:
        night_start_time, _ = get_nightly_times_from_allocation(
            night_planner.allocation_file, night_planner.current_day
        )

        def _minutes_to_hhmm(minutes):
            if pd.notna(minutes):
                return str(TimeDelta(minutes * 60, format="sec") + night_start_time)[
                    11:16
                ]
            return ""

        for col in ("Start Exposure", "Earliest Start", "Latest Finish"):
            if col in final_df.columns:
                final_df[col] = final_df[col].apply(_minutes_to_hhmm)

    except (AttributeError, OSError, TypeError, ValueError) as e:
        print(f"Warning: Could not convert time fields to HST timestamps: {e}")
        print("Time fields will remain as minutes from start of night")

    # Handle missing values and 'None' strings
    final_df = final_df.replace(["", "NoGaiaName", "None"], pd.NA)

    # Ensure DataFrame is clean and properly structured for DataTables
    final_df = final_df.reset_index(drop=True)
    # Remove duplicate column names if any exist
    final_df = final_df.loc[:, ~final_df.columns.duplicated(keep="first")]
    # Fill NaN values with empty strings to ensure consistent structure
    final_df = final_df.fillna("")
    # Ensure all columns have consistent data types (convert objects to strings)
    for col in final_df.columns:
        if final_df[col].dtype == "object":
            final_df[col] = (
                final_df[col]
                .astype(str)
                .replace("nan", "")
                .replace("None", "")
                .replace("", "")
            )

    return final_df
