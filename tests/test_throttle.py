"""Tests for the program throttle's past-slot accounting.

The throttle must count past observations on ALL request rows, including
inactive ones, so a PI cannot reclaim budget by flipping a target inactive.
"""

import unittest
from configparser import ConfigParser

import pandas as pd

from astroq.splan import SemesterPlanner


def _make_planner(requests, past):
    """Build a minimal SemesterPlanner stub (no Gurobi, no I/O).

    Sets only the attributes ``_add_program_columns`` depends on and
    attaches program columns the same way the real constructor does.
    """
    sp = SemesterPlanner.__new__(SemesterPlanner)
    sp.requests = requests.copy()
    sp.requests["r"] = sp.requests["unique_id"]
    sp.past = past.copy()
    sp.past["r"] = sp.past["unique_id"]
    sp.config = ConfigParser()
    sp.config.read_string(
        """
        [semester]
        slot_size = 20
        hours_per_night = 12
        """
    )
    sp.programs = pd.DataFrame(
        {"hours": [60.0]},
        index=pd.Index(["2026A_X001"], name="program"),
    )
    sp._add_program_columns()
    return sp


def _requests():
    """One program, one active + one inactive target, identical strategy.

    ``t_visit_slots`` is supplied directly so the test does not depend on the
    queue's overhead model.
    """
    return pd.DataFrame(
        {
            "unique_id": ["ACT", "INACT"],
            "program_code": ["2026A_X001", "2026A_X001"],
            "inactive": [False, True],
            "t_visit_slots": [3, 3],
        }
    )


class TestPastSlotsByProgram(unittest.TestCase):
    def test_includes_inactive_past(self):
        rfa = _requests()
        past = pd.DataFrame(
            {
                # 2 exposures on the active target, 4 on the inactive one.
                "unique_id": ["ACT", "ACT", "INACT", "INACT", "INACT", "INACT"],
                "timestamp": ["2026-02-01T10:00"] * 6,
            }
        )
        sp = _make_planner(rfa, past)
        # (2 active + 4 inactive) exposures * 3 slots/visit = 18 slots.
        self.assertEqual(sp.programs.loc["2026A_X001", "past_slots"], 18)

    def test_inactive_only_still_counts(self):
        rfa = _requests()
        past = pd.DataFrame(
            {
                "unique_id": ["INACT", "INACT"],
                "timestamp": ["2026-02-01T10:00", "2026-02-02T10:00"],
            }
        )
        sp = _make_planner(rfa, past)
        # If inactive past were ignored this would be 0.
        self.assertEqual(sp.programs.loc["2026A_X001", "past_slots"], 6)

    def test_empty_past(self):
        rfa = _requests()
        past = pd.DataFrame(columns=["unique_id", "timestamp"])
        sp = _make_planner(rfa, past)
        self.assertEqual(sp.programs.loc["2026A_X001", "past_slots"], 0)

    def test_integer_unique_id_join(self):
        """Integer unique_ids in past.csv are stringified at load (PAST_SCHEMA),
        so they match string request ids without repair in the throttle."""
        import os
        import tempfile

        from astroq.io import read_csv

        rfa = pd.DataFrame(
            {
                "unique_id": ["101", "202"],
                "program_code": ["2026A_X001", "2026A_X001"],
                "inactive": [False, True],
                "t_visit_slots": [2, 2],
            }
        )
        raw = pd.DataFrame(
            {
                "unique_id": [101, 202, 202],
                "target": ["a", "b", "b"],
                "timestamp": ["t"] * 3,
                "exposure_time": [60, 60, 60],
            }
        )
        path = os.path.join(tempfile.mkdtemp(prefix="astroq_past_"), "past.csv")
        raw.to_csv(path, index=False)
        past = read_csv(path, "past")
        self.assertEqual(past["unique_id"].tolist(), ["101", "202", "202"])

        sp = _make_planner(rfa, past)
        # (1 active + 2 inactive) * 2 slots = 6.
        self.assertEqual(sp.programs.loc["2026A_X001", "past_slots"], 6)


if __name__ == "__main__":
    unittest.main()
