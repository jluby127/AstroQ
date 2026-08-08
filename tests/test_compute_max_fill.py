"""Tests for compute-max-fill driver edge cases."""

import glob
import os
import shutil
import tempfile
import unittest
from argparse import Namespace
from configparser import ConfigParser
from unittest.mock import MagicMock, patch

import pandas as pd

from astroq.driver import compute_max_fill


def _requests(program_codes):
    return pd.DataFrame(
        {
            "program_code": list(program_codes),
            "unique_id": [f"R{i}" for i in range(len(program_codes))],
            "inactive": [False] * len(program_codes),
        }
    )


class ComputeMaxFillCase(unittest.TestCase):
    """Builds a run dir with one config, programs.csv and request.csv."""

    def build_run_dir(self, programs, requests):
        self.tmp = tempfile.mkdtemp(prefix="astroq_max_fill_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.workdir = os.path.join(self.tmp, "run")
        os.makedirs(self.workdir)

        cfg = ConfigParser()
        cfg.optionxform = str
        cfg["global"] = {"workdir": self.workdir}
        cfg["data"] = {"request_file": "request.csv", "programs_file": "programs.csv"}
        self.config_path = os.path.join(self.tmp, "config.ini")
        with open(self.config_path, "w", encoding="utf-8") as fh:
            cfg.write(fh)

        self.programs_path = os.path.join(self.workdir, "programs.csv")
        programs.to_csv(self.programs_path, index=False)
        requests.to_csv(os.path.join(self.workdir, "request.csv"), index=False)
        return requests

    def run_driver(self):
        compute_max_fill(Namespace(config_file=self.config_path))
        return pd.read_csv(self.programs_path)

    def stub_planner(self, mock_planner_cls, *, request_slots, fill=None, program=None):
        planner = MagicMock()
        planner.request_slots = request_slots
        if fill is None:
            planner.F = {}
        else:
            f_var = MagicMock()
            f_var.X = fill
            planner.F = {program: f_var}
        mock_planner_cls.return_value = planner
        return planner


class TestComputeMaxFill(ComputeMaxFillCase):
    @patch("astroq.driver.splan.SemesterPlanner")
    @patch("astroq.driver.astroq.io.read_csv")
    def test_no_request_rows_sets_zero(self, mock_read_csv, mock_planner_cls):
        mock_read_csv.return_value = self.build_run_dir(
            pd.DataFrame({"program": ["P1"], "hours": [10.0]}),
            _requests(["P2"]),
        )

        programs = self.run_driver()

        self.assertEqual(programs.loc[0, "max_feasible_fill"], 0.0)
        mock_planner_cls.assert_not_called()

    @patch("astroq.driver.splan.SemesterPlanner")
    @patch("astroq.driver.astroq.io.read_csv")
    def test_zero_observable_slots_sets_zero(self, mock_read_csv, mock_planner_cls):
        mock_read_csv.return_value = self.build_run_dir(
            pd.DataFrame({"program": ["P2"], "hours": [10.0]}),
            _requests(["P2"]),
        )
        planner = self.stub_planner(mock_planner_cls, request_slots=pd.DataFrame())

        programs = self.run_driver()

        self.assertEqual(programs.loc[0, "max_feasible_fill"], 0.0)
        mock_planner_cls.assert_called_once()
        planner.build_model.assert_not_called()
        planner.solve_shortfall.assert_not_called()

    @patch("astroq.driver.splan.SemesterPlanner")
    @patch("astroq.driver.astroq.io.read_csv")
    def test_no_fill_variable_sets_zero(self, mock_read_csv, mock_planner_cls):
        mock_read_csv.return_value = self.build_run_dir(
            pd.DataFrame({"program": ["P3"], "hours": [10.0]}),
            _requests(["P3"]),
        )
        planner = self.stub_planner(
            mock_planner_cls, request_slots=pd.DataFrame({"r": ["R0"]})
        )

        programs = self.run_driver()

        self.assertEqual(programs.loc[0, "max_feasible_fill"], 0.0)
        planner.build_model.assert_called_once()
        planner.solve_shortfall.assert_called_once()

    @patch("astroq.driver.splan.SemesterPlanner")
    @patch("astroq.driver.astroq.io.read_csv")
    def test_solved_fill_is_saved(self, mock_read_csv, mock_planner_cls):
        mock_read_csv.return_value = self.build_run_dir(
            pd.DataFrame({"program": ["P4"], "hours": [10.0]}),
            _requests(["P4"]),
        )
        planner = self.stub_planner(
            mock_planner_cls,
            request_slots=pd.DataFrame({"r": ["R0"]}),
            fill=0.85,
            program="P4",
        )

        programs = self.run_driver()

        self.assertEqual(programs.loc[0, "max_feasible_fill"], 0.85)
        planner.build_model.assert_called_once()
        planner.solve_shortfall.assert_called_once()
        planner.run_model_shortfall.assert_not_called()

    @patch("astroq.driver.splan.SemesterPlanner")
    @patch("astroq.driver.astroq.io.read_csv")
    def test_requests_passed_in_memory(self, mock_read_csv, mock_planner_cls):
        mock_read_csv.return_value = self.build_run_dir(
            pd.DataFrame({"program": ["P5"], "hours": [10.0]}),
            _requests(["P5", "P6", "P5"]),
        )
        self.stub_planner(
            mock_planner_cls,
            request_slots=pd.DataFrame({"r": ["R0"]}),
            fill=1.0,
            program="P5",
        )

        self.run_driver()

        _, kwargs = mock_planner_cls.call_args
        self.assertTrue(kwargs["defer_model"])
        passed = kwargs["requests"]
        self.assertIsInstance(passed, pd.DataFrame)
        self.assertEqual(set(passed["program_code"]), {"P5"})
        self.assertEqual(len(passed), 2)

    @patch("astroq.driver.splan.SemesterPlanner")
    @patch("astroq.driver.astroq.io.read_csv")
    def test_no_request_csv_files_written(self, mock_read_csv, mock_planner_cls):
        mock_read_csv.return_value = self.build_run_dir(
            pd.DataFrame(
                {
                    "program": ["P1", "P2", "P3"],
                    "hours": [10.0, 10.0, 0.0],
                }
            ),
            _requests(["P2", "P3"]),
        )
        self.stub_planner(
            mock_planner_cls,
            request_slots=pd.DataFrame({"r": ["R0"]}),
            fill=0.5,
            program="P2",
        )

        programs = self.run_driver()

        self.assertEqual(
            dict(zip(programs["program"], programs["max_feasible_fill"])),
            {"P1": 0.0, "P2": 0.5, "P3": 0.0},
        )
        written = {
            os.path.basename(p)
            for p in glob.glob(os.path.join(self.workdir, "request*.csv"))
        }
        self.assertEqual(written, {"request.csv"})
        self.assertEqual(
            glob.glob(os.path.join(self.tmp, "request*.csv")),
            [],
        )


if __name__ == "__main__":
    unittest.main()
