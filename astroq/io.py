"""
Module for the building and writing of reports for human readability and debugging of schedule solution. 
"""

# Standard library imports
import math
import os
import re
import time

# Third-party imports
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time, TimeDelta
from astropy import units as u
import numpy as np
import pandas as pd

# Local imports
import astroq.history as hs
import astroq.access as ac

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
    df = pd.DataFrame(Yrds.keys(),columns=['r','d','s'])
    df['value'] = [Yrds[k].x for k in Yrds.keys()]
    sparse = df.query('value>0').copy()
    sparse.drop(columns=['value'], inplace=True)
    sparse['name'] = sparse['r'].map(
        dict(zip(semester_planner.requests_frame['unique_id'], semester_planner.requests_frame['starname']))
    ).fillna("NO MATCHING NAME")
    sparse.to_csv(os.path.join(semester_planner.output_directory, "semester_plan.csv"), index=False, na_rep="")
    semester_planner.future_forecast = "semester_plan.csv"

    day, slot = np.mgrid[:semester_planner.semester_length,:semester_planner.n_slots_in_night]
    dense1 = pd.DataFrame(dict(d=day.flatten(), s=slot.flatten()))
    dense1 = pd.merge(dense1, sparse, left_on=['d','s'],right_on=['d','s'],how='left')
    dense1['r'] = dense1['r'].fillna('')
    # dense1 has keys for all days and slots, where no star was scheduled to start its observation, the r column is blank
    dense1.to_csv(os.path.join(semester_planner.output_directory, "serialized_outputs_dense_v1.csv"), index=False, na_rep="")

    dense2 = dense1.copy()
    # Use the stored access record from the first run (no need to recompute)
    access = semester_planner.access_record
    isAlloc = access['is_alloc'].flatten()
    isClear = access['is_clear'].flatten()
    # have to go backwards otherwise you're adding stars into slots and then testing if the star is in the next slot
    for slot in range(semester_planner.n_slots_in_semester-1, -1, -1):
        name_string = ""
        if isAlloc[slot] == 0:
            name_string += "X"
        if isClear[slot] == 0:
            name_string += "W"
        dense2.loc[slot, 'r'] = name_string + str(dense2.loc[slot, 'r'])
        if dense2.loc[slot, 'r'] in list(semester_planner.requests_frame['unique_id']):
            slots_needed = semester_planner.slots_needed_for_exposure_dict[dense2.loc[slot, 'r']]
            if slots_needed > 1:
                for t in range(1, slots_needed):
                    dense2.loc[slot + t, 'r'] = str(dense2.loc[slot + t, 'r']) + str(dense2.loc[slot, 'r'])
    # dense2 has keys for all days and slots, manually fill in the reserved slots for each observation and fill in Past/Twilight/Weather info
    dense2.to_csv(os.path.join(semester_planner.output_directory, "serialized_outputs_dense_v2.csv"), index=False, na_rep="")
    
    # Generate the fullness report
    build_fullness_report(semester_planner, "Round1")
    return sparse

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
    print(f"Writing runReport.txt to: {file_path}")
    
    # Read the semester plan CSV file
    semester_plan_path = os.path.join(semester_planner.output_directory, "semester_plan.csv")
    if not os.path.exists(semester_plan_path):
        print(f"Warning: semester_plan.csv not found at {semester_plan_path}")
        return
    
    schedule_df = pd.read_csv(semester_plan_path)
    
    # Get access record data
    access = semester_planner.access_record
    is_alloc = access['is_alloc']
    
    # Calculate statistics
    # is_alloc is 3D array, but all 2D arrays are the same, so use the first one
    is_alloc_2d = is_alloc[0]  # Take the first 2D array
    total_slots_in_semester = is_alloc_2d.shape[0] * is_alloc_2d.shape[1]  # Total slots = n_nights * n_slots_per_night
    allocated_slots = np.sum(is_alloc_2d)  # Slots that are allocated (not X)
    
    # Calculate scheduled slots (starting slots only)
    scheduled_starting_slots = len(schedule_df)  # Number of starting slots with stars scheduled
    
    # Calculate reserved slots (slots beyond starting slots for multi-slot exposures)
    reserved_slots = 0
    for _, row in schedule_df.iterrows():
        star_name = row['r']
        if star_name in semester_planner.slots_needed_for_exposure_dict:
            slots_needed = semester_planner.slots_needed_for_exposure_dict[star_name]
            reserved_slots += slots_needed - 1  # Subtract 1 because the starting slot is already counted
    
    total_scheduled_slots = scheduled_starting_slots + reserved_slots
    empty_slots = allocated_slots - total_scheduled_slots
    
    # Calculate total slots requested using slots_needed_for_exposure_dict and n_intra_max
    total_slots_requested = 0
    for i in range(len(semester_planner.requests_frame)):
        star_id = semester_planner.requests_frame['unique_id'][i]
        if star_id in semester_planner.slots_needed_for_exposure_dict:
            slots_needed = semester_planner.slots_needed_for_exposure_dict[star_id]
            total_slots_requested += slots_needed * semester_planner.requests_frame['n_intra_max'][i] * semester_planner.requests_frame['n_inter_max'][i]
    
    # Calculate percentages
    percentage_of_available = np.round((total_scheduled_slots * 100) / allocated_slots, 3) if allocated_slots > 0 else 0
    percentage_of_requested = np.round((total_scheduled_slots * 100) / total_slots_requested, 3) if total_slots_requested > 0 else 0
    
    with open(os.path.join(semester_planner.output_directory, "runReport.txt"), "w") as file:
        file.write("Stats for " + str(round_info) + "\n")
        file.write("------------------------------------------------------" + "\n")
        file.write("N slots in semester:" + str(total_slots_in_semester) + "\n")
        file.write("N available slots:" + str(allocated_slots) + "\n")
        file.write("N starting slots scheduled: " + str(scheduled_starting_slots) + "\n")
        file.write("N reserved slots: " + str(reserved_slots) + "\n")
        file.write("N total slots scheduled: " + str(total_scheduled_slots) + "\n")
        file.write("N slots left empty: " + str(empty_slots) + "\n")
        file.write("N slots requested (total): " + str(total_slots_requested) + "\n")
        file.write("Utilization (% of available slots): " + str(percentage_of_available) + "%" + "\n")
        file.write("Utilization (% of requested slots): " + str(percentage_of_requested) + "%" + "\n")
        file.close()
