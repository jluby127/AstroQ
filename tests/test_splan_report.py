"""Tests for SemesterPlanner.timeline, programs_ledger, and plot COF helpers."""

import unittest
from configparser import ConfigParser
from unittest.mock import MagicMock

import numpy as np
import pandas as pd

from astroq.plot.context import _cume_matrix
from astroq.splan import TIMELINE_VALUE_COLS, SemesterPlanner


def _stub_planner(
    *, past=None, schedule=None, requests=None, programs=None, current_night_index=2
):
    sp = SemesterPlanner.__new__(SemesterPlanner)
    sp.config = ConfigParser()
    sp.config.optionxform = str
    sp.config.read_string(
        """
        [semester]
        slot_size = 2
        hours_per_night = 12
        """
    )
    sp.schedule = schedule
    sp.requests = requests.copy()
    sp.requests["r"] = sp.requests["unique_id"]
    if "inactive" not in sp.requests.columns:
        sp.requests["inactive"] = False
    for col in ("n_intra_max", "n_inter_max"):
        if col not in sp.requests.columns:
            sp.requests[col] = 1
    sp.past = past.copy() if past is not None else pd.DataFrame(
        columns=["unique_id", "target", "timestamp", "exposure_time"]
    )
    sp.past["r"] = sp.past["unique_id"]
    sp.programs = programs.copy()
    sp._add_program_columns()
    sp.access_obj = MagicMock()
    sp.access_obj.all_dates_array = [
        "2026-02-01",
        "2026-02-02",
        "2026-02-03",
        "2026-02-04",
    ]
    sp.access_obj.current_night_index = current_night_index
    return sp


def _attach_mock_fill(sp, *, fill=None):
    fill = fill or {}
    sp.F = {p: MagicMock() for p in sp.programs.index}
    sp.model = MagicMock()
    sp.model.getAttr.side_effect = lambda attr, var: {
        "X": {p: fill.get(p, (0.0, 0.0, 1.0))[0] for p in sp.programs.index},
        "LB": {p: fill.get(p, (0.0, 0.0, 1.0))[1] for p in sp.programs.index},
        "UB": {p: fill.get(p, (0.0, 0.0, 1.0))[2] for p in sp.programs.index},
    }[attr]


class TestTimeline(unittest.TestCase):
    def test_timeline_concat(self):
        requests = pd.DataFrame(
            {
                "unique_id": ["R1", "R2"],
                "program_code": ["P1", "P1"],
                "t_visit_slots": [3, 2],
            }
        )
        requests["r"] = requests["unique_id"]
        past = pd.DataFrame(
            {
                "unique_id": ["R1", "R1"],
                "target": ["R1", "R1"],
                "timestamp": ["2026-02-01 12:00", "2026-02-02 12:00"],
                "exposure_time": [100.0, 100.0],
            }
        )
        past["r"] = past["unique_id"]
        schedule = pd.DataFrame(
            {
                "unique_id": ["R2", "R2"],
                "d": [2, 3],
                "s": [0, 5],
                "target": ["R2", "R2"],
            }
        )
        programs = pd.DataFrame({"hours": [10.0]}, index=pd.Index(["P1"], name="program"))
        sp = _stub_planner(
            past=past, schedule=schedule, requests=requests, programs=programs
        )

        ps = sp.timeline
        self.assertEqual(ps.index.name, "d")
        self.assertEqual(list(ps.columns), list(TIMELINE_VALUE_COLS))
        self.assertEqual(len(ps), 4)
        today = sp.access_obj.current_night_index
        self.assertEqual(ps.loc[ps.index < today, "t_visit_slots"].sum(), 6)
        self.assertEqual(ps.loc[ps.index >= today, "t_visit_slots"].sum(), 4)

    def test_past_slots_matches_legacy(self):
        requests = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "program_code": ["P1"],
                "t_visit_slots": [3],
            }
        )
        requests["r"] = requests["unique_id"]
        past = pd.DataFrame(
            {
                "unique_id": ["R1", "R1"],
                "target": ["R1", "R1"],
                "timestamp": ["2026-02-01 12:00", "2026-02-02 12:00"],
                "exposure_time": [100.0, 100.0],
            }
        )
        past["r"] = past["unique_id"]
        programs = pd.DataFrame({"hours": [5.0]}, index=pd.Index(["P1"], name="program"))
        sp = _stub_planner(past=past, schedule=None, requests=requests, programs=programs)
        sp._add_program_columns()
        self.assertEqual(int(sp.programs.loc["P1", "past_slots"]), 6)

    def test_programs_ledger_charged_hours(self):
        requests = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "program_code": ["P1"],
                "t_visit_slots": [6],
                "n_intra_max": [1],
                "n_inter_max": [2],
            }
        )
        requests["r"] = requests["unique_id"]
        past = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "target": ["R1"],
                "timestamp": ["2026-02-01 12:00"],
                "exposure_time": [100.0],
            }
        )
        past["r"] = past["unique_id"]
        schedule = pd.DataFrame(
            {"unique_id": ["R1"], "d": [1], "s": [0], "target": ["R1"]}
        )
        programs = pd.DataFrame({"hours": [5.0]}, index=pd.Index(["P1"], name="program"))
        sp = _stub_planner(
            past=past,
            schedule=schedule,
            requests=requests,
            programs=programs,
            current_night_index=1,
        )
        _attach_mock_fill(sp)
        ledger = sp.programs_ledger
        self.assertAlmostEqual(ledger.loc["P1", "past_hours"], 0.2)
        self.assertAlmostEqual(ledger.loc["P1", "sched_hours"], 0.2)
        self.assertAlmostEqual(ledger.loc["P1", "proj_hours"], 0.4)
        self.assertAlmostEqual(ledger.loc["P1", "requested_hours"], 0.4)

    def test_programs_ledger_requires_schedule(self):
        requests = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "program_code": ["P1"],
                "t_visit_slots": [1],
            }
        )
        requests["r"] = requests["unique_id"]
        programs = pd.DataFrame({"hours": [1.0]}, index=pd.Index(["P1"], name="program"))
        empty_past = pd.DataFrame(
            columns=["unique_id", "target", "timestamp", "exposure_time"]
        )
        sp = _stub_planner(
            past=empty_past, schedule=None, requests=requests, programs=programs
        )
        with self.assertRaises(RuntimeError):
            _ = sp.programs_ledger

    def test_timeline_requires_build_schedule(self):
        requests = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "program_code": ["P1"],
                "t_visit_slots": [1],
            }
        )
        requests["r"] = requests["unique_id"]
        programs = pd.DataFrame({"hours": [1.0]}, index=pd.Index(["P1"], name="program"))
        empty_past = pd.DataFrame(
            columns=["unique_id", "target", "timestamp", "exposure_time"]
        )
        sp = _stub_planner(
            past=empty_past, schedule=None, requests=requests, programs=programs
        )
        with self.assertRaises(RuntimeError):
            _ = sp.timeline

    def test_cache_invalidated_on_build_schedule(self):
        requests = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "program_code": ["P1"],
                "t_visit_slots": [1],
            }
        )
        requests["r"] = requests["unique_id"]
        programs = pd.DataFrame({"hours": [1.0]}, index=pd.Index(["P1"], name="program"))
        empty_past = pd.DataFrame(
            columns=["unique_id", "target", "timestamp", "exposure_time"]
        )
        schedule = pd.DataFrame(columns=["unique_id", "d", "s", "target"])
        sp = _stub_planner(
            past=empty_past, schedule=schedule, requests=requests, programs=programs
        )
        _ = sp.timeline
        self.assertIn("timeline", sp.__dict__)
        sp._invalidate_timeline()
        self.assertNotIn("timeline", sp.__dict__)


class TestProgramsLedger(unittest.TestCase):
    def test_zero_award_program_pct_columns(self):
        requests = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "program_code": ["P0"],
                "t_visit_slots": [6],
            }
        )
        past = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "target": ["R1"],
                "timestamp": ["2026-02-01 12:00"],
                "exposure_time": [100.0],
            }
        )
        schedule = pd.DataFrame(
            {"unique_id": ["R1"], "d": [1], "s": [0], "target": ["R1"]}
        )
        programs = pd.DataFrame({"hours": [0.0]}, index=pd.Index(["P0"], name="program"))
        sp = _stub_planner(
            past=past,
            schedule=schedule,
            requests=requests,
            programs=programs,
            current_night_index=1,
        )
        _attach_mock_fill(sp)
        ledger = sp.programs_ledger
        self.assertAlmostEqual(ledger.loc["P0", "past_hours"], 0.2)
        self.assertAlmostEqual(ledger.loc["P0", "proj_hours"], 0.4)

    def test_fill_columns_from_model(self):
        requests = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "program_code": ["P1"],
                "t_visit_slots": [30],
            }
        )
        schedule = pd.DataFrame(
            {"unique_id": ["R1"], "d": [0], "s": [0], "target": ["R1"]}
        )
        programs = pd.DataFrame({"hours": [10.0]}, index=pd.Index(["P1"], name="program"))
        sp = _stub_planner(
            past=pd.DataFrame(
                columns=["unique_id", "target", "timestamp", "exposure_time"]
            ),
            schedule=schedule,
            requests=requests,
            programs=programs,
        )
        _attach_mock_fill(sp, fill={"P1": (0.5, 0.1, 1.25)})
        ledger = sp.programs_ledger
        self.assertAlmostEqual(ledger.loc["P1", "fill_proj"], 0.5)
        self.assertAlmostEqual(ledger.loc["P1", "fill_min"], 0.1)
        self.assertAlmostEqual(ledger.loc["P1", "fill_max"], 1.25)

    def test_program_without_requests_has_zero_requested_hours(self):
        requests = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "program_code": ["P1"],
                "t_visit_slots": [6],
            }
        )
        schedule = pd.DataFrame(
            {"unique_id": ["R1"], "d": [0], "s": [0], "target": ["R1"]}
        )
        programs = pd.DataFrame(
            {"hours": [5.0, 3.0]},
            index=pd.Index(["P1", "P2"], name="program"),
        )
        sp = _stub_planner(
            past=pd.DataFrame(
                columns=["unique_id", "target", "timestamp", "exposure_time"]
            ),
            schedule=schedule,
            requests=requests,
            programs=programs,
        )
        _attach_mock_fill(sp)
        ledger = sp.programs_ledger
        self.assertEqual(ledger.loc["P2", "requested_hours"], 0.0)


class TestPlotHelpers(unittest.TestCase):
    def test_cume_matrix_visits(self):
        ps = pd.DataFrame(
            {
                "unique_id": ["R1", "R1", "R1"],
                "program_code": ["P1", "P1", "P1"],
                "t_visit_slots": [3, 3, 3],
            },
            index=pd.Index([0, 1, 2], name="d"),
        )
        cume = _cume_matrix(ps, 4, "unique_id", ["R1"], "visits")
        np.testing.assert_array_equal(cume["R1"].to_numpy(), [1, 2, 3, 3])

    def test_cume_matrix_slots(self):
        ps = pd.DataFrame(
            {
                "unique_id": ["R1", "R1"],
                "program_code": ["P1", "P1"],
                "t_visit_slots": [30, 30],
            },
            index=pd.Index([0, 1], name="d"),
        )
        cume = _cume_matrix(ps, 2, "program_code", ["P1"], "slots")
        np.testing.assert_array_equal(cume["P1"].to_numpy(), [30, 60])


if __name__ == "__main__":
    unittest.main()
