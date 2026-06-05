import argparse
import os
import unittest
import warnings

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from tables.exceptions import DataTypeWarning

import astroq.driver as dr
import astroq.nplan as nplan
import astroq.splan as splan
from astroq.ttp.model import TTPModel


class TestClass(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # python -m unittest installs warnings.simplefilter("default") at
        # runner startup, which wipes the filter installed by astroq/__init__.py.
        # Re-install it inside unittest's catch_warnings() context so it sticks
        # for the duration of the test class.
        warnings.filterwarnings("ignore", category=DataTypeWarning)

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
        self.assertIsInstance(solution.requests["coord"], SkyCoord)
        self.assertIsInstance(solution.requests["first_available"], Time)

    def test10_nightly_availability_windows(self):
        """Access exposes first/last_available as Time arrays from slotmidpoints."""
        outputs_dir = "examples/hello_world/2018B/2018-08-05/band1/outputs"
        sp = splan.SemesterPlanner.from_hdf5(
            os.path.join(outputs_dir, "semester_planner.h5"),
        )
        # build_access populates first_available, last_available,
        # has_observable. Not done by from_hdf5; trigger explicitly so
        # the test does not depend on h5 contents.
        access_record = sp.access_obj.build_access()

        night_d = sp.all_dates_dict[sp.current_day]
        uids = sp.requests_frame["unique_id"].iloc[:3]

        # Cross-check first/last_available against the long-form observability
        # table for each observable uid tonight.
        req_index = sp.access_obj.request_frame.set_index("unique_id").index
        row_idx = req_index.get_indexer(uids)
        first_available = sp.access_obj.first_available[row_idx, night_d]
        last_available = sp.access_obj.last_available[row_idx, night_d]
        has_obs = sp.access_obj.has_observable[row_idx, night_d]

        self.assertEqual(len(first_available), 3)
        self.assertIsInstance(first_available, Time)
        self.assertIsInstance(last_available, Time)

        obs = sp.access_obj.observability(access_record.is_observable)
        night = obs.loc[obs["d"] == night_d]
        for k, uid in enumerate(uids):
            if not has_obs[k]:
                continue
            slots = night.loc[night["unique_id"] == uid, "s"]
            s_min, s_max = int(slots.min()), int(slots.max())
            self.assertEqual(
                first_available[k], sp.access_obj.slotmidpoints[night_d, s_min]
            )
            self.assertEqual(
                last_available[k], sp.access_obj.slotmidpoints[night_d, s_max]
            )

    def test11_get_nightly_times_missing_day(self):
        allo = "examples/hello_world/2018B/2018-08-05/band1/allocation.csv"
        with self.assertRaises(ValueError):
            nplan.get_nightly_times_from_allocation(allo, "1900-01-01")

    def test12_webapp(self):
        """Drive every webapp route via Flask's test client against the hello_world
        fixture. The plot cache dir is redirected to a tmp dir so the cache-miss
        branch of get_football's seasonality grid actually runs. Rendered HTML
        is dumped to the same tmp dir for visual inspection.
        """
        import tempfile
        from unittest.mock import patch

        import astroq.plot as pl
        import astroq.webapp.app as wa

        tmp = tempfile.mkdtemp(prefix="astroq_webapp_")
        print(f"webapp HTML dumps -> {tmp}")

        # Pick a real (program_code, starname) pair from the fixture so the URL
        # actually resolves inside data_astroq[0].
        request_selected = pd.read_csv(
            "examples/hello_world/2018B/2018-08-05/band1/"
            "outputs/request_selected.csv"
        )
        row = request_selected.iloc[0]
        program_code = str(row["program_code"])
        starname = str(row["starname"])

        # Redirect get_football's on-disk cache so the test does not write into
        # the committed data/ directory and so the cache-miss branch executes.
        from pathlib import Path
        with patch.object(pl, "_football_cache_dir", lambda sp: Path(tmp)):
            wa.uptree_path = "examples/hello_world"
            client = wa.app.test_client()

            # (path, expected_codes, dump_name)
            # expected_codes is a set; download_nightplan may legitimately 500
            # when the stored output_directory does not resolve from this cwd.
            routes = [
                ("/", {200}, "homepage.html"),
                ("/2018B/2018-08-05/band1/admin", {200}, "admin.html"),
                (
                    f"/2018B/2018-08-05/band1/{program_code}",
                    {200},
                    f"program_{program_code}.html",
                ),
                (
                    f"/2018B/2018-08-05/band1/{program_code}/{starname}",
                    {200},
                    f"star_{starname}.html",
                ),
                ("/2018B/2018-08-05/band1/nightplan", {200}, "nightplan.html"),
                (
                    "/2018B/2018-08-05/band1/download_nightplan",
                    {200, 500},
                    "download_nightplan.txt",
                ),
                # Invalid band -> 400 from abort()
                ("/2018B/2018-08-05/bogus/admin", {400}, None),
                # Valid band but missing date -> 404 from load_data_for_path
                ("/2018B/1900-01-01/band1/admin", {404}, None),
            ]
            for path, expected_codes, dump_name in routes:
                resp = client.get(path)
                msg = f"{path} -> {resp.status_code}: {resp.data[:500]!r}"
                self.assertIn(resp.status_code, expected_codes, msg=msg)
                if dump_name is not None and resp.status_code == 200:
                    with open(os.path.join(tmp, dump_name), "wb") as f:
                        f.write(resp.data)


if __name__ == "__main__":
    unittest.main()
