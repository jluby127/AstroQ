import astroq.driver as dr
import astroq.benchmarking as bn
import argparse
from configparser import ConfigParser
import os
import astroq.splan as splan
import astroq.plot as pl
import astroq.nplan as nplan
import unittest
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
        
        # Validate solution exists
        self.assertTrue(hasattr(night_planner, 'solution'), "NightPlanner missing attribute: solution")
        self.assertIsNotNone(night_planner.solution, "solution is None")
        self.assertIsInstance(night_planner.solution, list, "solution is not a list")
        self.assertGreater(len(night_planner.solution), 0, "solution list is empty")
        
        solution = night_planner.solution[0]
        
        # Validate solution attributes
        self.assertTrue(hasattr(solution, 'plotly'), "solution missing attribute: plotly")
        self.assertIsNotNone(solution.plotly, "solution.plotly is None")
        self.assertIsInstance(solution.plotly, dict, "solution.plotly is not a dict")
        self.assertGreater(len(solution.plotly), 0, "solution.plotly is empty")
        # Check plotly has expected keys
        expected_plotly_keys = ['Starname', 'Start Exposure', 'Minutes the from Start of the Night']
        for key in expected_plotly_keys:
            if key in solution.plotly:
                self.assertIsNotNone(solution.plotly[key], f"solution.plotly['{key}'] is None")
        
        self.assertTrue(hasattr(solution, 'times'), "solution missing attribute: times")
        self.assertIsNotNone(solution.times, "solution.times is None")
        self.assertIsInstance(solution.times, list, "solution.times is not a list")
        if len(solution.times) > 0:
            self.assertIsInstance(solution.times[0], Time, "solution.times elements are not Time objects")
        
        self.assertTrue(hasattr(solution, 'nightstarts'), "solution missing attribute: nightstarts")
        self.assertIsNotNone(solution.nightstarts, "solution.nightstarts is None")
        self.assertIsInstance(solution.nightstarts, Time, "solution.nightstarts is not a Time object")
        
        self.assertTrue(hasattr(solution, 'nightends'), "solution missing attribute: nightends")
        self.assertIsNotNone(solution.nightends, "solution.nightends is None")
        self.assertIsInstance(solution.nightends, Time, "solution.nightends is not a Time object")
        
        self.assertTrue(hasattr(solution, 'schedule'), "solution missing attribute: schedule")
        self.assertIsNotNone(solution.schedule, "solution.schedule is None")
        self.assertIsInstance(solution.schedule, dict, "solution.schedule is not a dict")
        
        self.assertTrue(hasattr(solution, 'stars'), "solution missing attribute: stars")
        self.assertIsNotNone(solution.stars, "solution.stars is None")
        self.assertIsInstance(solution.stars, list, "solution.stars is not a list")
        if len(solution.stars) > 0:
            star = solution.stars[0]
            self.assertTrue(hasattr(star, 'name'), "star missing attribute: name")
            self.assertIsNotNone(star.name, "star.name is None")
            self.assertTrue(hasattr(star, 'target'), "star missing attribute: target")
            self.assertIsNotNone(star.target, "star.target is None")
        
        self.assertTrue(hasattr(solution, 'az_path'), "solution missing attribute: az_path")
        self.assertIsNotNone(solution.az_path, "solution.az_path is None")
        self.assertIsInstance(solution.az_path, np.ndarray, "solution.az_path is not a numpy array")
        self.assertGreater(solution.az_path.size, 0, "solution.az_path is empty")
        
        self.assertTrue(hasattr(solution, 'alt_path'), "solution missing attribute: alt_path")
        self.assertIsNotNone(solution.alt_path, "solution.alt_path is None")
        self.assertIsInstance(solution.alt_path, np.ndarray, "solution.alt_path is not a numpy array")
        self.assertGreater(solution.alt_path.size, 0, "solution.alt_path is empty")
        
        self.assertTrue(hasattr(solution, 'extras'), "solution missing attribute: extras")
        self.assertIsNotNone(solution.extras, "solution.extras is None")
        # extras can be DataFrame or dict, both are valid
        self.assertTrue(isinstance(solution.extras, (pd.DataFrame, dict)), 
                       f"solution.extras is not DataFrame or dict, got {type(solution.extras)}")
        
        self.assertTrue(hasattr(solution, 'observatory'), "solution missing attribute: observatory")
        self.assertIsNotNone(solution.observatory, "solution.observatory is None")
        
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
