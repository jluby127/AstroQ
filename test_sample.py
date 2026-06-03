import argparse
import os
import unittest

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time

import astroq.driver as dr
import astroq.nplan as nplan
import astroq.splan as splan
from astroq.ttp.model import TTPModel


class TestClass(unittest.TestCase):

    def test01_helloworld(self):
        dr.plan_semester(
            argparse.Namespace(
                config_file="examples/hello_world/config_hello_world.ini",
                run_band3=False,
            )
        )

    def test02_round2_weather(self):
        dr.plan_semester(
            argparse.Namespace(
                config_file="examples/hello_world/config_hello_world_bonus_weather.ini",
                run_band3=False,
            )
        )

    def test03_plan_night(self):
        dr.plan_night(
            argparse.Namespace(
                config_file="examples/hello_world/config_hello_world.ini",
            )
        )

    def test04_bench(self):
        dr.bench(
            argparse.Namespace(
                config_file="examples/bench/config_benchmark.ini",
                number_slots=12,
                thin=10,
            )
        )
        dr.plan_night(
            argparse.Namespace(
                config_file="examples/bench/config_benchmark.ini",
            )
        )
        dr.plot(
            argparse.Namespace(
                config_file="examples/bench/config_benchmark.ini",
            )
        )

    def test05_generic_prep(self):
        dr.kpfcc_prep(
            argparse.Namespace(
                config_file="examples/hello_world/config_hello_world_prep.ini",
                allo_source="examples/hello_world/prepped/observatory_schedule.csv",
                past_source="examples/hello_world/prepped/jump_past_history.csv",
                request_source="examples/hello_world/prepped/request.csv",
                filler_programs="2025B_E473",
                band_number=1,
                is_full_band=False,
            )
        )
        dr.kpfcc_prep(
            argparse.Namespace(
                config_file="examples/hello_world/config_hello_world_prep.ini",
                allo_source="examples/hello_world/prepped/observatory_schedule.csv",
                past_source="examples/hello_world/prepped/jump_past_history.csv",
                request_source="examples/hello_world/prepped/request.csv",
                filler_programs="2025B_E473",
                band_number=3,
                is_full_band=True,
            )
        )

    def test06_kpfcc_prep(self):
        dr.kpfcc_prep(
            argparse.Namespace(
                config_file="examples/hello_world/config_hello_world_prep.ini",
                allo_source="db",
                past_source="db",
                request_source="db",
                filler_programs="2025B_E473",
                band_number=1,
                is_full_band=True,
            )
        )

    def test07_plot(self):
        dr.plot(
            argparse.Namespace(
                config_file="examples/hello_world/config_hello_world.ini",
            )
        )

    def test08_requests_vs_schedule(self):
        sch = "examples/hello_world/2018B/2018-08-05/band1/outputs/semester_plan.csv"
        dr.requests_vs_schedule(
            argparse.Namespace(
                config_file="examples/hello_world/config_hello_world.ini",
                schedule_file=sch,
            )
        )

    def test09_hdf5_validation(self):
        """NightPlanner schema v5: config_ini_text + TTP solution round-trip."""
        outputs_dir = "examples/hello_world/2018B/2018-08-05/band1/outputs"
        semester_planner_h5 = os.path.join(outputs_dir, "semester_planner.h5")
        night_planner_h5 = os.path.join(outputs_dir, "night_planner.h5")

        self.assertTrue(os.path.exists(semester_planner_h5), semester_planner_h5)
        self.assertTrue(os.path.exists(night_planner_h5), night_planner_h5)

        semester_planner = splan.SemesterPlanner.from_hdf5(semester_planner_h5)
        night_planner = nplan.NightPlanner.from_hdf5(night_planner_h5)

        self.assertTrue(hasattr(night_planner, "config"))
        self.assertTrue(hasattr(night_planner, "_config_ini_text"))
        self.assertEqual(
            night_planner.current_day,
            night_planner.config.get("global", "current_day"),
        )
        self.assertEqual(
            semester_planner.current_day,
            night_planner.current_day,
        )
        self.assertEqual(
            semester_planner.semester_directory,
            night_planner.semester_directory,
        )

        solution = night_planner.solution
        self.assertIsInstance(solution, TTPModel)
        self.assertIsInstance(solution.night_start, Time)
        self.assertIsInstance(solution.night_end, Time)
        self.assertFalse(solution.schedule.empty)
        self.assertIsInstance(solution.requests_frame["coord"].iloc[0], SkyCoord)
        self.assertIsInstance(solution.requests_frame["first_available"].iloc[0], Time)

    def test10_nightly_availability_windows(self):
        """Access helper returns Time scalars from slotmidpoints."""
        outputs_dir = "examples/hello_world/2018B/2018-08-05/band1/outputs"
        sp = splan.SemesterPlanner.from_hdf5(
            os.path.join(outputs_dir, "semester_planner.h5"),
        )
        night_d = sp.all_dates_dict[sp.current_day]
        uids = sp.requests_frame["unique_id"].iloc[:3]
        windows = sp.access_obj.nightly_availability_windows(
            night_d,
            uids,
            sp.access_record,
        )
        self.assertEqual(len(windows), 3)
        self.assertIsInstance(windows["first_available"].iloc[0], Time)
        self.assertIsInstance(windows["last_available"].iloc[0], Time)

        obs = sp.access_obj.observability(
            sp.requests_frame,
            access=sp.access_record,
        )
        night = obs.loc[obs["d"] == night_d]
        uid = uids.iloc[0]
        s_min = night.loc[night["unique_id"] == uid, "s"].min()
        s_max = night.loc[night["unique_id"] == uid, "s"].max()
        row = windows.loc[windows["unique_id"] == uid].iloc[0]
        self.assertEqual(row["first_available"], sp.access_obj.slotmidpoints[night_d, s_min])
        self.assertEqual(row["last_available"], sp.access_obj.slotmidpoints[night_d, s_max])

    def test11_get_nightly_times_missing_day(self):
        allo = "examples/hello_world/2018B/2018-08-05/band1/allocation.csv"
        with self.assertRaises(ValueError):
            nplan.get_nightly_times_from_allocation(allo, "1900-01-01")


if __name__ == "__main__":
    unittest.main()
