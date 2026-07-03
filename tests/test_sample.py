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
            )
        )

    def test02_full_mode(self):
        """mode=full exercises Rounds 1-5 on the symmetric toy model."""
        dr.plan_semester(
            argparse.Namespace(
                config_file="examples/priorities/symmetric_toy_model/config_benchmark.ini",
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

    def test08_hdf5_validation(self):
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
            semester_planner.config.get("global", "current_day"),
            night_planner.current_day,
        )
        self.assertEqual(
            semester_planner.config.get("global", "workdir"),
            night_planner.semester_directory,
        )

        solution = night_planner.solution
        self.assertIsInstance(solution, TTPModel)
        self.assertIsInstance(solution.night_start, Time)
        self.assertIsInstance(solution.night_end, Time)
        self.assertFalse(solution.schedule.empty)
        self.assertIsInstance(solution.requests["coord"], SkyCoord)
        self.assertIsInstance(solution.requests["time_earliest_start"], Time)

    def test10_nightly_availability_windows(self):
        """Night windows from is_observable_now first/last clear slot."""
        outputs_dir = "examples/hello_world/2018B/2018-08-05/band1/outputs"
        sp = splan.SemesterPlanner.from_hdf5(
            os.path.join(outputs_dir, "semester_planner.h5"),
        )
        access_record = sp.access_obj.build_access()

        night_d = sp.all_dates_dict[sp.config.get("global", "current_day")]
        uids = sp.requests_frame["unique_id"].iloc[:3]

        req_index = sp.access_obj.request_frame.set_index("unique_id").index
        row_idx = req_index.get_indexer(uids)
        now = access_record.is_observable_now[row_idx, night_d, :]
        slotmid = sp.access_obj.slotmidpoints[night_d]
        time_earliest_start = slotmid[now.argmax(1)]
        time_latest_finish = slotmid[now.shape[1] - 1 - now[:, ::-1].argmax(1)]

        self.assertEqual(len(time_earliest_start), 3)
        self.assertIsInstance(time_earliest_start, Time)
        self.assertIsInstance(time_latest_finish, Time)

        obs = sp.access_obj.observability(access_record.is_observable_now)
        night = obs.loc[obs["d"] == night_d]
        for k, uid in enumerate(uids):
            slots = night.loc[night["unique_id"] == uid, "s"]
            if slots.empty:
                continue
            s_min, s_max = int(slots.min()), int(slots.max())
            self.assertEqual(
                time_earliest_start[k], sp.access_obj.slotmidpoints[night_d, s_min]
            )
            self.assertEqual(
                time_latest_finish[k], sp.access_obj.slotmidpoints[night_d, s_max]
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

        # Pick a real (program_code, target) pair from the fixture so the URL
        # actually resolves inside data_astroq[0].
        request_selected = pd.read_csv(
            "examples/hello_world/2018B/2018-08-05/band1/"
            "outputs/request_selected.csv"
        )
        row = request_selected.iloc[0]
        program_code = str(row["program_code"])
        target = str(row["target"])

        # Redirect get_football's on-disk cache so the test does not write into
        # the committed data/ directory and so the cache-miss branch executes.
        from pathlib import Path
        with patch.object(pl, "_football_cache_dir", lambda sp: Path(tmp)):
            wa._uptree_path = "examples/hello_world"
            client = wa.app.test_client()

            # (path, expected_codes, dump_name)
            routes = [
                ("/", {200}, "homepage.html"),
                ("/2018B/2018-08-05/band1/admin", {200}, "admin.html"),
                (
                    f"/2018B/2018-08-05/band1/{program_code}",
                    {200},
                    f"program_{program_code}.html",
                ),
                (
                    f"/2018B/2018-08-05/band1/{program_code}/{target}",
                    {200},
                    f"star_{target}.html",
                ),
                ("/2018B/2018-08-05/band1/nightplan", {200}, "nightplan.html"),
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

        # Flat routes via -rp (single run directory)
        workdir = "examples/hello_world/2018B/2018-08-05/band1"
        with patch.object(pl, "_football_cache_dir", lambda sp: Path(tmp)):
            wa._run_path = workdir
            wa._uptree_path = None
            client = wa.app.test_client()
            flat_routes = [
                ("/admin", {200, 404}, "flat_admin.html"),
                ("/2018B/2018-08-05/band1/admin", {404}, None),
            ]
            for path, expected_codes, dump_name in flat_routes:
                resp = client.get(path)
                msg = f"{path} -> {resp.status_code}: {resp.data[:500]!r}"
                self.assertIn(resp.status_code, expected_codes, msg=msg)
                if dump_name is not None and resp.status_code == 200:
                    with open(os.path.join(tmp, dump_name), "wb") as f:
                        f.write(resp.data)

        # Parent -rp: /{run_name}/admin
        parent = "examples/hello_world/2018B/2018-08-05"
        with patch.object(pl, "_football_cache_dir", lambda sp: Path(tmp)):
            wa._run_path = parent
            wa._uptree_path = None
            client = wa.app.test_client()
            parent_routes = [
                ("/admin", {404}, None),
                ("/band1/admin", {200, 404}, "parent_admin.html"),
            ]
            for path, expected_codes, dump_name in parent_routes:
                resp = client.get(path)
                msg = f"{path} -> {resp.status_code}: {resp.data[:500]!r}"
                self.assertIn(resp.status_code, expected_codes, msg=msg)
                if dump_name is not None and resp.status_code == 200:
                    with open(os.path.join(tmp, dump_name), "wb") as f:
                        f.write(resp.data)

    def test13_archive(self):
        """Export admin, nightplan, and program static HTML via astroq archive."""
        import tempfile
        from unittest.mock import patch

        import astroq.plot as pl

        cf = "examples/hello_world/config_hello_world.ini"
        workdir = "examples/hello_world/2018B/2018-08-05/band1"
        archive_dir = os.path.join(workdir, "outputs", "webapp_archive")
        admin_path = os.path.join(archive_dir, "admin.html")
        night_path = os.path.join(archive_dir, "nightplan.html")
        programs_dir = os.path.join(archive_dir, "programs")

        tmp = tempfile.mkdtemp(prefix="astroq_archive_")
        from pathlib import Path

        with patch.object(pl, "_football_cache_dir", lambda sp: Path(tmp)):
            dr.archive(argparse.Namespace(config_file=cf))

        self.assertTrue(os.path.isfile(admin_path), f"missing {admin_path}")
        self.assertTrue(os.path.isfile(night_path), f"missing {night_path}")
        self.assertTrue(os.path.isdir(programs_dir), f"missing {programs_dir}")
        program_html = [
            f for f in os.listdir(programs_dir) if f.endswith(".html")
        ]
        self.assertGreater(len(program_html), 0, "no program archive pages written")

        with open(admin_path, encoding="utf-8") as f:
            admin_html = f.read()
        with open(night_path, encoding="utf-8") as f:
            night_html = f.read()
        with open(
            os.path.join(programs_dir, program_html[0]), encoding="utf-8"
        ) as f:
            program_page = f.read()

        self.assertIn("Admin Dashboard", admin_html)
        self.assertIn("plotly-graph-div", admin_html)
        self.assertNotIn('href="/2018B/2018-08-05/band1/', admin_html)

        self.assertIn("Night Plan", night_html)
        self.assertIn("plotly-graph-div", night_html)
        self.assertNotIn("download_nightplan", night_html)

        self.assertIn("Semester Plan", program_page)
        self.assertIn("plotly-graph-div", program_page)
        self.assertNotIn('href="/2018B/2018-08-05/band1/', program_page)

    def test14_jump_query_to_past(self):
        """Visit groups need >=50% of n_exp frames; one row per accepted visit."""
        import tempfile

        import astroq.queue.hirescps.prep as prep

        frames = pd.DataFrame(
            {
                "starname": ["T1", "T1", "T1", "T2", "T2", "T2"],
                "timestamp": [
                    "2026-03-02 14:15",
                    "2026-03-02 14:17",
                    "2026-03-02 14:20",
                    "2026-04-04 06:20",
                    "2026-04-05 07:13",
                    "2026-04-05 07:15",
                ],
                "exposure_time": [95, 128, 169, 95, 16, 16],
                "decker": ["C2"] * 6,
                "iodine_in": [False] * 6,
            }
        )
        tmp = tempfile.mkdtemp(prefix="astroq_visits_")
        req_csv = os.path.join(tmp, "request.csv")
        pd.DataFrame({"unique_id": ["T1", "T2"], "n_exp": [3, 3]}).to_csv(
            req_csv, index=False
        )
        out = prep.jump_query_to_past(frames, req_csv)
        self.assertEqual(len(out), 2)
        t1 = out.loc[out["unique_id"] == "T1"].iloc[0]
        self.assertEqual(t1["timestamp"], "2026-03-02 14:15")
        self.assertEqual(t1["exposure_time"], 392)
        t2 = out.loc[out["unique_id"] == "T2"].iloc[0]
        self.assertEqual(t2["timestamp"], "2026-04-05 07:13")
        self.assertEqual(t2["exposure_time"], 32)
        # Lone 1/3 frame on 2026-04-04 must not appear.
        self.assertNotIn("2026-04-04 06:20", out["timestamp"].tolist())

    def test15_hires_past_history_query(self):
        """JUMP pull writes tmp file and collapses frames to visit rows."""
        import tempfile
        from unittest.mock import MagicMock, patch

        import astroq.queue.hirescps.prep as prep

        raw_csv = (
            "starname,timestamp,exposure_time,decker,iodine_in,counts\n"
            "109358,2026-03-02 14:15,95,C2,True,250000\n"
            "109358,2026-03-02 14:17,128,C2,True,250000\n"
            "109358,2026-03-02 14:20,169,C2,True,250000\n"
            "109358,2026-04-04 06:20,95,C2,True,250000\n"
            "156079,2026-05-31 13:46,85,B3,False,2000\n"
            "T004478,2026-05-31 13:46,1799,B3,False,2000\n"
        )
        captured = {}

        login_resp = MagicMock()
        login_resp.raise_for_status = lambda: None

        post_resp = MagicMock()
        post_resp.raise_for_status = lambda: None
        post_resp.url = f"{prep.JUMP_BASE_URL}/explorer/"

        download_resp = MagicMock()
        download_resp.content = raw_csv.encode("utf-8")
        download_resp.raise_for_status = lambda: None

        def fake_get(url, *args, **kwargs):
            if "/download" in url:
                captured["url"] = url
                return download_resp
            return login_resp

        fake_session = MagicMock()
        fake_session.cookies = {"csrftoken": "test-token"}
        fake_session.get.side_effect = fake_get
        fake_session.post.return_value = post_resp

        tmp = tempfile.mkdtemp(prefix="astroq_past_")
        out_csv = os.path.join(tmp, "past.csv")
        req_csv = os.path.join(tmp, "request.csv")
        pd.DataFrame(
            {
                "unique_id": ["109358", "156079_t", "T004478_t"],
                "n_exp": [3, 1, 1],
            }
        ).to_csv(req_csv, index=False)

        with patch.object(prep.requests, "Session", return_value=fake_session):
            prep.get_hires_past_history(
                out_csv,
                semester_start_day="2026-02-01",
                semester_end_day="2026-07-31",
                request_csv_path=req_csv,
            )

        tmp_csv = os.path.join(tmp, prep.JUMP_PAST_QUERY_TMP_FILENAME)
        self.assertTrue(os.path.isfile(tmp_csv))
        self.assertIn(
            f"/explorer/{prep.JUMP_HIRES_PAST_EXPLORER_ID}/download", captured["url"]
        )
        # Download URL must carry a real ``&params=`` (regression: BeautifulSoup
        # turned it into ``¶ms=``).
        self.assertIn("&params=", captured["url"])
        self.assertNotIn("\u00b6", captured["url"])
        self.assertIn("start_date%3A2026-02-01", captured["url"])
        self.assertIn("end_date%3A2026-07-31", captured["url"])

        out = pd.read_csv(out_csv)
        self.assertEqual(
            list(out.columns), ["unique_id", "target", "timestamp", "exposure_time"]
        )
        self.assertEqual(out["unique_id"].tolist(), out["target"].tolist())
        self.assertIn("109358", out["target"].tolist())
        self.assertIn("156079_t", out["target"].tolist())
        self.assertIn("T004478_t", out["target"].tolist())
        self.assertNotIn("156079", out["target"].tolist())
        self.assertEqual(len(out.loc[out["unique_id"] == "109358"]), 1)
        self.assertEqual(
            out.loc[out["unique_id"] == "109358", "timestamp"].iloc[0],
            "2026-03-02 14:15",
        )

    def test16_hires_past_history_requires_dates(self):
        """Live JUMP pull needs both semester dates."""
        import astroq.queue.hirescps.prep as prep

        with self.assertRaises(ValueError):
            prep.get_hires_past_history("unused.csv", semester_start_day="2026-02-01")


if __name__ == "__main__":
    unittest.main()
