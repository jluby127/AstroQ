"""
Module defining the admin object, which stores, manipulates, and passes information across
all requests easily.

Example usage:
    import admin_functions as af
"""
import os
import math

import numpy as np
from astropy.time import Time
from astropy.time import TimeDelta

class data_admin(object):
    """A Data Admin object, from which we can easily pass information."""

    def __init__(self, output_directory, current_day):

        self.output_directory = output_directory
        # Suggest your output directory be something so that it doesn't autosave
        # to the same directory as the run files and crowds up the GitHub repo.
        check = os.path.isdir(self.output_directory)
        if not check:
            os.makedirs(self.output_directory)

        self.current_day = current_day






    def run_admin(self, n_quarters_in_night=4, n_hours_in_night=14):
        """
        Given today's date, collate all important information about the semester.

        Args:
            current_day (string): today's date in format "YYYY-MM-DD"
            n_quarters_in_night (int): the number of quarters in a night, at Keck Observatory this is 4
            n_hours_in_night (int): the number of hours in a night, this should be safely larger than
                                    the longest night of the year, so 14 hours is a good number.

        Returns:
            semester_start_date (string): first day of the semester (civil HST) in format "YYYY-MM-DD"
            semester_length (int): total number of days in the full semester
            all_dates_dict (dictionary): keys are the calendar dates and values are the day of the semester
            all_dates_array (array): the keys of the all_dates_dict dictionary
            n_slots_in_night (int): the number of slots in a night
            n_nights_in_semester (int): the number of nights remaining in the semester
            today_starting_slot (int): the slot number which corresponds to the first slot of today's date
            today_starting_night (int): the day number in semester associated with today's date
        """

        # Get semester parameters and define important quantities
        semester_start_date, semester_end_date, semester_length, semester_year, semester_letter = \
            get_semester_info(current_day)
        all_dates_dict = build_date_dictionary(semester_start_date, semester_length)
        all_dates_array = list(all_dates_dict.keys())
        n_nights_in_semester = current_day_tracker(current_day, all_dates_dict)
        print("Total semester length: ", semester_length)
        print("There are " + str(n_nights_in_semester) + " calendar nights remaining in the semester.")

        # Compute slot grid sizes
        n_slots_in_quarter = int(((n_hours_in_night*60)/n_quarters_in_night)/slot_size)
        n_slots_in_night = n_slots_in_quarter*n_quarters_in_night
        n_slots_in_semester = n_slots_in_night*n_nights_in_semester

        # Define the slot and night represents the today's date
        today_starting_slot = all_dates_dict[current_day]*n_slots_in_night
        today_starting_night =  all_dates_dict[current_day]
        print("There are " + str(n_slots_in_semester) + " slots remaining in the semester.")

        return semester_start_date, semester_length, all_dates_dict, all_dates_array, n_slots_in_night \
                n_nights_in_semester, today_starting_slot, today_starting_night
