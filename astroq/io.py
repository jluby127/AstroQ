"""
Module for the building and writing of reports for human readability and debugging of schedule solution.
"""

# Standard library imports
import os

# Third-party imports
import numpy as np
import pandas as pd
from jinja2 import Template


def validate_active_request_unique_ids(requests_frame, request_path):
    """
    Raise ValueError if active requests contain duplicate ``unique_id`` values.

    Duplicate IDs produce duplicate (unique_id, day, slot) keys for Gurobi variables
    and fail with a cryptic ``Duplicate keys in Model.addVars()`` error.

    Args:
        requests_frame: DataFrame of active requests (inactive rows already removed).
        request_path: Path to the request CSV, for error messages.
    """
    uid = requests_frame["unique_id"]
    dup_mask = uid.duplicated(keep=False)
    if not dup_mask.any():
        return
    dup_ids = sorted(uid[dup_mask].unique())
    raise ValueError(
        f"Duplicate unique_id among active requests in {request_path!r}: {dup_ids}. "
        "Remove or merge duplicate rows so each active request has one row."
    )


def serialize_schedule(Yrds, semester_planner):
    """
    Turns the ragged matrix of the Gurobi solution Yrds into a human readable
    solution by filling in the slots where a star's exposure is started.

    Args:
        Yrds (array): the Gurobi solution with keys of (unique_id, day, slot) and values 1 or 0.
        semester_planner (obj): a SemesterPlanner object from splan.py

    Returns:
        None
    """
    df = pd.DataFrame(Yrds.keys(), columns=["r", "d", "s"])
    df["value"] = [Yrds[k].x for k in Yrds.keys()]
    sparse = df.query("value>0").copy()
    sparse.drop(columns=["value"], inplace=True)
    sparse["name"] = (
        sparse["r"]
        .map(
            dict(
                zip(
                    semester_planner.requests_frame["unique_id"],
                    semester_planner.requests_frame["starname"],
                )
            )
        )
        .fillna("NO MATCHING NAME")
    )
    sparse.to_csv(
        os.path.join(semester_planner.output_directory, "semester_plan.csv"),
        index=False,
        na_rep="",
    )
    semester_planner.future_forecast = "semester_plan.csv"

    day, slot = np.mgrid[
        : semester_planner.semester_length, : semester_planner.n_slots_in_night
    ]
    dense1 = pd.DataFrame(dict(d=day.flatten(), s=slot.flatten()))
    dense1 = pd.merge(
        dense1, sparse, left_on=["d", "s"], right_on=["d", "s"], how="left"
    )
    dense1["r"] = dense1["r"].fillna("")
    # dense1 has keys for all days and slots, where no star was scheduled to start its observation, the r column is blank
    dense1.to_csv(
        os.path.join(
            semester_planner.output_directory, "serialized_outputs_dense_v1.csv"
        ),
        index=False,
        na_rep="",
    )

    dense2 = dense1.copy()
    # Use the stored access record from the first run (no need to recompute)
    access = semester_planner.access_record
    isAlloc = access["is_allocated"].flatten()
    isClear = access["is_clear"].flatten()
    # have to go backwards otherwise you're adding stars into slots and then testing if the star is in the next slot
    for slot in range(semester_planner.n_slots_in_semester - 1, -1, -1):
        name_string = ""
        if isAlloc[slot] == 0:
            name_string += "X"
        if isClear[slot] == 0:
            name_string += "W"
        dense2.loc[slot, "r"] = name_string + str(dense2.loc[slot, "r"])
        if dense2.loc[slot, "r"] in list(semester_planner.requests_frame["unique_id"]):
            slots_needed = semester_planner.slots_needed_for_exposure_dict[
                dense2.loc[slot, "r"]
            ]
            if slots_needed > 1:
                for t in range(1, slots_needed):
                    dense2.loc[slot + t, "r"] = str(dense2.loc[slot + t, "r"]) + str(
                        dense2.loc[slot, "r"]
                    )
    # dense2 has keys for all days and slots, manually fill in the reserved slots for each observation and fill in Past/Twilight/Weather info
    dense2.to_csv(
        os.path.join(
            semester_planner.output_directory, "serialized_outputs_dense_v2.csv"
        ),
        index=False,
        na_rep="",
    )

    # Generate the fullness report
    build_fullness_report(semester_planner, "Round1")
    return sparse


def compute_program_statistics(semester_planner, schedule_df):
    """
    Compute awarded, requested, and scheduled statistics per program.

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        schedule_df (DataFrame): DataFrame containing scheduled observations

    Returns:
        dict: Dictionary keyed by program name, with values containing:
            - awarded_nights, awarded_hours, awarded_slots
            - past_nights, past_hours, past_slots
            - requested_nights, requested_hours, requested_slots
            - scheduled_nights, scheduled_hours, scheduled_slots
            - fullness1 ((scheduled+past)/requested %)
            - fullness2 ((scheduled+past)/awarded %)
    """
    # Get conversion factors
    slot_size = semester_planner.slot_size  # in minutes
    hours_per_night = semester_planner.hours_per_night  # in hours
    slots_per_hour = 60 / slot_size
    slots_per_night = hours_per_night * slots_per_hour

    program_frame = pd.read_csv(semester_planner.programs_file)
    program_stats = {}
    for _, prog_row in program_frame.iterrows():
        program = prog_row["program"]
        awarded_nights = prog_row["nights"]
        awarded_hours = awarded_nights * hours_per_night
        awarded_slots = awarded_nights * slots_per_night

        # Calculate requested slots for this program
        program_requests = semester_planner.requests_frame[
            semester_planner.requests_frame["program_code"] == program
        ]

        # Map slots_needed from dict and calculate requested slots using vectorized operations
        program_requests = program_requests.copy()
        program_requests["slots_needed"] = (
            program_requests["unique_id"]
            .map(semester_planner.slots_needed_for_exposure_dict)
            .fillna(0)
        )
        requested_slots = (
            program_requests["slots_needed"]
            * program_requests["n_intra_max"]
            * program_requests["n_inter_max"]
        ).sum()
        requested_hours = requested_slots / slots_per_hour
        requested_nights = requested_hours / hours_per_night

        # Calculate past slots for this program (already observed)
        past_slots = 0
        for _, req_row in program_requests.iterrows():
            star_id = req_row["unique_id"]
            if star_id in semester_planner.past_history:
                past_hist = semester_planner.past_history[star_id]
                # Get slots per exposure for this star
                slots_per_exposure = (
                    semester_planner.slots_needed_for_exposure_dict.get(star_id, 1)
                )
                # Total past slots = number of exposures * slots per exposure
                past_slots += past_hist.total_n_exposures * slots_per_exposure
        past_hours = past_slots / slots_per_hour
        past_nights = past_hours / hours_per_night

        # Calculate scheduled slots for this program
        # Merge schedule_df with requests_frame to get program_code for each scheduled star
        schedule_with_program = schedule_df.merge(
            semester_planner.requests_frame[["unique_id", "program_code"]],
            left_on="r",
            right_on="unique_id",
            how="inner",
        )
        # Filter for current program and map slots_needed, then sum
        program_scheduled = schedule_with_program[
            schedule_with_program["program_code"] == program
        ]
        scheduled_slots = (
            program_scheduled["r"]
            .map(semester_planner.slots_needed_for_exposure_dict)
            .fillna(1)
            .sum()
        )
        scheduled_hours = scheduled_slots / slots_per_hour
        scheduled_nights = scheduled_hours / hours_per_night

        # Calculate fullness percentages based on (scheduled + past)
        total_scheduled_and_past_slots = scheduled_slots + past_slots
        fullnessA = (
            (total_scheduled_and_past_slots * 100.0 / requested_slots)
            if requested_slots > 0
            else 0.0
        )
        fullnessB = (
            (total_scheduled_and_past_slots * 100.0 / awarded_slots)
            if awarded_slots > 0
            else 0.0
        )
        fullnessC = (
            (requested_slots * 100.0 / awarded_slots) if awarded_slots > 0 else 0.0
        )

        program_stats[program] = {
            "awarded_nights": awarded_nights,
            "awarded_hours": awarded_hours,
            "awarded_slots": awarded_slots,
            "past_nights": past_nights,
            "past_hours": past_hours,
            "past_slots": past_slots,
            "requested_nights": requested_nights,
            "requested_hours": requested_hours,
            "requested_slots": requested_slots,
            "scheduled_nights": scheduled_nights,
            "scheduled_hours": scheduled_hours,
            "scheduled_slots": scheduled_slots,
            "fullnessA": fullnessA,
            "fullnessB": fullnessB,
            "fullnessC": fullnessC,
        }

    return program_stats


def build_fullness_report(semester_planner, round_info):
    """
    Determine how full the schedule is: slots available, slots scheduled, and slots required

    Args:
        semester_planner (obj): a SemesterPlanner object from splan.py
        round_info (str): information about the optimization round (e.g. Round1, Round2)

    Returns:
        None
    """
    file_path = os.path.join(semester_planner.output_directory, "runReport.txt")
    # logs.info(f"Writing runReport.txt to: {file_path}")

    print(semester_planner.output_directory)
    print(file_path)

    # Read the semester plan CSV file
    semester_plan_path = os.path.join(
        semester_planner.output_directory, "semester_plan.csv"
    )
    if not os.path.exists(semester_plan_path):
        # logs.info(f"Warning: semester_plan.csv not found at {semester_plan_path}")
        return

    schedule_df = pd.read_csv(semester_plan_path)

    # Get access record data
    access = semester_planner.access_record
    is_allocated = access["is_allocated"]

    # Calculate statistics
    # is_allocated is 3D array, but all 2D arrays are the same, so use the first one
    is_alloc_2d = is_allocated[0]  # Take the first 2D array
    total_slots_in_semester = (
        is_alloc_2d.shape[0] * is_alloc_2d.shape[1]
    )  # Total slots = n_nights * n_slots_per_night
    allocated_slots = np.sum(is_alloc_2d)  # Slots that are allocated (not X)

    # Calculate scheduled slots (starting slots only)
    scheduled_starting_slots = len(
        schedule_df
    )  # Number of starting slots with stars scheduled

    # Calculate reserved slots (slots beyond starting slots for multi-slot exposures)
    reserved_slots = 0
    for _, row in schedule_df.iterrows():
        star_name = row["r"]
        if star_name in semester_planner.slots_needed_for_exposure_dict:
            slots_needed = semester_planner.slots_needed_for_exposure_dict[star_name]
            reserved_slots += (
                slots_needed - 1
            )  # Subtract 1 because the starting slot is already counted

    total_scheduled_slots = scheduled_starting_slots + reserved_slots
    empty_slots = allocated_slots - total_scheduled_slots

    # Calculate total slots requested using slots_needed_for_exposure_dict and n_intra_max
    total_slots_requested = 0
    for i in range(len(semester_planner.requests_frame)):
        star_id = semester_planner.requests_frame["unique_id"][i]
        if star_id in semester_planner.slots_needed_for_exposure_dict:
            slots_needed = semester_planner.slots_needed_for_exposure_dict[star_id]
            total_slots_requested += (
                slots_needed
                * semester_planner.requests_frame["n_intra_max"][i]
                * semester_planner.requests_frame["n_inter_max"][i]
            )

    # Calculate percentages
    percentage_of_available = (
        np.round((total_scheduled_slots * 100) / allocated_slots, 3)
        if allocated_slots > 0
        else 0
    )
    percentage_of_requested = (
        np.round((total_scheduled_slots * 100) / total_slots_requested, 3)
        if total_slots_requested > 0
        else 0
    )

    # Compute program statistics
    program_stats = compute_program_statistics(semester_planner, schedule_df)

    with open(file_path, "w") as file:
        file.write("Stats for " + str(round_info) + "\n")
        file.write("------------------------------------------------------" + "\n")
        file.write("N slots in semester:" + str(total_slots_in_semester) + "\n")
        file.write("N available slots:" + str(allocated_slots) + "\n")
        file.write(
            "N starting slots scheduled: " + str(scheduled_starting_slots) + "\n"
        )
        file.write("N reserved slots: " + str(reserved_slots) + "\n")
        file.write("N total slots scheduled: " + str(total_scheduled_slots) + "\n")
        file.write("N slots left empty: " + str(empty_slots) + "\n")
        file.write("N slots requested (total): " + str(total_slots_requested) + "\n")
        file.write(
            "Utilization (% of available slots): "
            + str(percentage_of_available)
            + "%"
            + "\n"
        )
        file.write(
            "Utilization (% of requested slots): "
            + str(percentage_of_requested)
            + "%"
            + "\n"
        )
        file.write("\n")
        file.write("Program Statistics:" + "\n")
        file.write("------------------------------------------------------" + "\n")
        if program_stats:
            # Template for program statistics
            program_template = Template("""{{ program }}
 -- Awarded {{ "%.2f"|format(stats['awarded_nights']) }} nights = {{ "%.2f"|format(stats['awarded_hours']) }} hours = {{ "%.1f"|format(stats['awarded_slots']) }} slots.
 -- Requested {{ "%.2f"|format(stats['requested_nights']) }} nights = {{ "%.2f"|format(stats['requested_hours']) }} hours = {{ "%.1f"|format(stats['requested_slots']) }} slots
 ------ Fullness of requested to awarded: {{ "%.2f"|format(stats['fullnessC']) }}%
 -- Past {{ "%.2f"|format(stats['past_nights']) }} nights = {{ "%.2f"|format(stats['past_hours']) }} hours = {{ "%.1f"|format(stats['past_slots']) }} slots.
 -- Scheduled {{ "%.2f"|format(stats['scheduled_nights']) }} nights = {{ "%.2f"|format(stats['scheduled_hours']) }} hours = {{ "%.1f"|format(stats['scheduled_slots']) }} slots
 ------ Fullness of past/scheduled to requested: {{ "%.2f"|format(stats['fullnessA']) }}%
 ------ Fullness of past/scheduled to awarded: {{ "%.2f"|format(stats['fullnessB']) }}%


""")
            # Sort by program name for consistent output
            for program in sorted(program_stats.keys()):
                stats = program_stats[program]
                rendered = program_template.render(program=program, stats=stats)
                file.write(rendered)

        else:
            file.write("  No program information available" + "\n")
        file.close()
