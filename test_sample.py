import astroq.driver as dr
import astroq.benchmarking as bn
import argparse
from configparser import ConfigParser
import os
import astroq.splan as splan
import astroq.plot as pl
import astroq.nplan as nplan
from astroq.ttp.model import TTPModel
import unittest
import pandas as pd
import numpy as np
from astropy.time import Time
from pathlib import Path
from io import BytesIO
import imageio.v3 as iio

import multiprocessing
import time
import requests
import os
from astroq.webapp import launch_app, app
import threading

class TestClass(unittest.TestCase):

    def test01_helloworld(self):
        dr.plan_semester(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini', run_band3=False))

    def test02_round2_weather(self):
        dr.plan_semester(argparse.Namespace(config_file='examples/hello_world/config_hello_world_bonus_weather.ini', run_band3=False))

    def test03_plan_night(self):
        dr.plan_night(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    def test04_bench(self):
        dr.bench(argparse.Namespace(config_file='examples/bench/config_benchmark.ini', number_slots=12, thin=10))
        dr.plan_night(argparse.Namespace(config_file='examples/bench/config_benchmark.ini'))
        dr.plot(argparse.Namespace(config_file='examples/bench/config_benchmark.ini'))

    def test05_generic_prep(self):
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world_prep.ini', allo_source='examples/hello_world/prepped/observatory_schedule.csv', past_source='examples/hello_world/prepped/jump_past_history.csv', request_source='examples/hello_world/prepped/request.csv', filler_programs='2025B_E473', band_number=1, is_full_band=False))
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world_prep.ini', allo_source='examples/hello_world/prepped/observatory_schedule.csv', past_source='examples/hello_world/prepped/jump_past_history.csv', request_source='examples/hello_world/prepped/request.csv', filler_programs='2025B_E473', band_number=3, is_full_band=True))

    def test06_kpfcc_prep(self):
       dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world_prep.ini', allo_source='db', past_source='db', request_source='db', filler_programs='2025B_E473', band_number=1, is_full_band=True))

    def test07_plot(self):
        dr.plot(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    def test08_requests_vs_schedule(self):
        sch = 'examples/hello_world/2018B/2018-08-05/band1/outputs/semester_plan.csv'
        dr.requests_vs_schedule(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini', schedule_file=sch))

    # # this is not working right now.
    # def test_simulate_history(self):
    #     dr.make_simulated_history(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    # # we don't care about Optimal Instrument Allocation yet
    # def test_oia(self):
    #     dr.kpfcc_prep(argparse.Namespace(config_file='examples/recreate_paper/oia1/config_oia1.ini'))
    #     dr.kpfcc_build(argparse.Namespace(config_file='examples/recreate_paper/oia1/config_oia1.ini'))
    #     dr.schedule(argparse.Namespace(request_file="examples/recreate_paper/oia1/outputs/2024-08-02/request_set.json", config_file='examples/recreate_paper/oia1/config_oia1.ini'))

    def test09_hdf5_validation(self):
        """Test loading and validating SemesterPlanner and NightPlanner from HDF5 files."""
        import os
        import pandas as pd
        import numpy as np
        from astropy.time import Time
        
        outputs_dir = 'examples/hello_world/2018B/2018-08-05/band1/outputs'
        semester_planner_h5 = os.path.join(outputs_dir, 'semester_planner.h5')
        night_planner_h5 = os.path.join(outputs_dir, 'night_planner.h5')
        
        # Check files exist
        self.assertTrue(os.path.exists(semester_planner_h5), f"Semester planner HDF5 file not found: {semester_planner_h5}")
        self.assertTrue(os.path.exists(night_planner_h5), f"Night planner HDF5 file not found: {night_planner_h5}")
        
        # Load SemesterPlanner
        semester_planner = splan.SemesterPlanner.from_hdf5(semester_planner_h5)
        
        # Validate SemesterPlanner scalar/string attributes
        scalar_attrs = [
            'current_day', 'semester_start_date', 'semester_length', 'semester_letter',
            'slot_size', 'n_slots_in_night', 'n_nights_in_semester', 'n_slots_in_semester',
            'today_starting_slot', 'today_starting_night', 'run_band3', 'observatory',
            'output_directory', 'run_weather_loss', 'solve_time_limit', 'gurobi_output',
            'solve_max_gap', 'max_bonus', 'run_bonus_round', 'semester_directory',
            'custom_file', 'allocation_file'
        ]
        
        for attr in scalar_attrs:
            self.assertTrue(hasattr(semester_planner, attr), f"SemesterPlanner missing attribute: {attr}")
            value = getattr(semester_planner, attr)
            self.assertIsNotNone(value, f"SemesterPlanner attribute {attr} is None")
        
        # Validate SemesterPlanner DataFrames
        self.assertIsNotNone(semester_planner.requests_frame, "requests_frame is None")
        self.assertIsInstance(semester_planner.requests_frame, pd.DataFrame, "requests_frame is not a DataFrame")
        self.assertFalse(semester_planner.requests_frame.empty, "requests_frame is empty")
        
        self.assertIsNotNone(semester_planner.serialized_schedule, "serialized_schedule is None")
        self.assertIsInstance(semester_planner.serialized_schedule, pd.DataFrame, "serialized_schedule is not a DataFrame")
        
        # Validate SemesterPlanner dictionaries
        dict_attrs = ['all_dates_dict', 'slots_needed_for_exposure_dict', 'past_nights_observed_dict']
        for attr in dict_attrs:
            self.assertTrue(hasattr(semester_planner, attr), f"SemesterPlanner missing attribute: {attr}")
            value = getattr(semester_planner, attr)
            self.assertIsNotNone(value, f"SemesterPlanner attribute {attr} is None")
            self.assertIsInstance(value, dict, f"SemesterPlanner attribute {attr} is not a dict")
        
        # Validate past_history
        self.assertTrue(hasattr(semester_planner, 'past_history'), "SemesterPlanner missing attribute: past_history")
        self.assertIsNotNone(semester_planner.past_history, "past_history is None")
        self.assertIsInstance(semester_planner.past_history, dict, "past_history is not a dict")
        
        # Validate all_dates_array
        self.assertTrue(hasattr(semester_planner, 'all_dates_array'), "SemesterPlanner missing attribute: all_dates_array")
        self.assertIsNotNone(semester_planner.all_dates_array, "all_dates_array is None")
        self.assertIsInstance(semester_planner.all_dates_array, list, "all_dates_array is not a list")
        self.assertGreater(len(semester_planner.all_dates_array), 0, "all_dates_array is empty")
        
        # Validate access_record
        self.assertTrue(hasattr(semester_planner, 'access_record'), "SemesterPlanner missing attribute: access_record")
        self.assertIsNotNone(semester_planner.access_record, "access_record is None")
        self.assertIsInstance(semester_planner.access_record, np.recarray, "access_record is not a recarray")
        
        # Validate access_obj
        self.assertTrue(hasattr(semester_planner, 'access_obj'), "SemesterPlanner missing attribute: access_obj")
        self.assertIsNotNone(semester_planner.access_obj, "access_obj is None")
        
        # Load NightPlanner
        night_planner = nplan.NightPlanner.from_hdf5(night_planner_h5)
        
        # Validate NightPlanner scalar/string attributes
        nightplanner_attrs = [
            'upstream_path', 'semester_directory', 'current_day', 'output_directory',
            'reports_directory', 'max_solve_gap', 'max_solve_time', 'show_gurobi_output',
            'allocation_file', 'filler_file', 'custom_file'
        ]
        
        for attr in nightplanner_attrs:
            self.assertTrue(hasattr(night_planner, attr), f"NightPlanner missing attribute: {attr}")
            value = getattr(night_planner, attr)
            self.assertIsNotNone(value, f"NightPlanner attribute {attr} is None")
        
        # Validate solution exists (schema v3: bare TTPModel, not [TTPModel])
        self.assertTrue(hasattr(night_planner, 'solution'), "NightPlanner missing attribute: solution")
        self.assertIsNotNone(night_planner.solution, "solution is None")
        solution = night_planner.solution
        self.assertIsInstance(solution, TTPModel, "solution is not a TTPModel")

        self.assertIsInstance(solution.night_start, Time)
        self.assertIsInstance(solution.night_end, Time)

        self.assertIsNotNone(solution.schedule)
        self.assertIsInstance(solution.schedule, pd.DataFrame)
        self.assertFalse(solution.schedule.empty)
        for col in ('scheduled', 't_start', 't_end', 'unique_id', 'is_anchor'):
            self.assertIn(col, solution.schedule.columns, f"schedule missing column: {col}")

        self.assertIsNotNone(solution.stats)
        self.assertIsInstance(solution.stats, dict)
        for key in ('dur', 'n_requested', 'n_scheduled', 't_first_start', 't_last_end',
                    't_visit_sum', 't_slew_sum'):
            self.assertIn(key, solution.stats, f"stats missing key: {key}")

        on_sky = solution.schedule[~solution.schedule['is_anchor']]
        scheduled = on_sky[on_sky['scheduled']]
        if len(scheduled):
            self.assertTrue(scheduled['order'].notna().all())
            self.assertTrue((scheduled['t_end'] >= scheduled['t_start']).all())
        
        # Validate semester_planner reference
        self.assertTrue(hasattr(night_planner, 'semester_planner'), "NightPlanner missing attribute: semester_planner")
        self.assertIsNotNone(night_planner.semester_planner, "semester_planner is None")
        self.assertIsInstance(night_planner.semester_planner, splan.SemesterPlanner, "semester_planner is not a SemesterPlanner instance")
        
        # Validate data consistency between semester_planner and night_planner
        self.assertEqual(semester_planner.current_day, night_planner.current_day, 
                        "current_day mismatch between SemesterPlanner and NightPlanner")
        self.assertEqual(semester_planner.semester_directory, night_planner.semester_directory,
                        "semester_directory mismatch between SemesterPlanner and NightPlanner")

if __name__=="__main__":
    unittest.main()
