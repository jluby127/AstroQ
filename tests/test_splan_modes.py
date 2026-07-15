"""Tests for the [semester] mode switch and pipeline dispatch."""

import unittest
from configparser import ConfigParser, NoOptionError
from unittest.mock import patch

from astroq.splan import _ALLOWED_MODES, _MODE_PIPELINES, SemesterPlanner


def _planner_with_config(text):
    """Minimal planner stub: only ``config`` is needed by ``run_model``."""
    sp = SemesterPlanner.__new__(SemesterPlanner)
    sp.config = ConfigParser()
    sp.config.optionxform = str
    sp.config.read_string(text)
    return sp


class TestRunModelDispatch(unittest.TestCase):
    def test_mode_pipelines(self):
        self.assertEqual(
            _MODE_PIPELINES,
            {
                "shortfall": "run_model_shortfall",
                "shortfall,balance,prioritize,fill-empty,fill-current-day": (
                    "run_model_shortfall_balance_prioritize_fillempty_fillcurrentday"
                ),
            },
        )
        self.assertEqual(
            _ALLOWED_MODES,
            [
                "shortfall",
                "shortfall,balance,prioritize,fill-empty,fill-current-day",
            ],
        )

    def test_missing_mode_raises(self):
        sp = _planner_with_config("[semester]\n")
        with self.assertRaises(NoOptionError):
            sp.run_model()

    def test_shortfall_mode(self):
        sp = _planner_with_config("[semester]\nmode = shortfall\n")
        with patch.object(
            sp, "run_model_shortfall_balance_prioritize_fillempty_fillcurrentday"
        ) as full:
            with patch.object(sp, "run_model_shortfall") as shortfall:
                sp.run_model()
                shortfall.assert_called_once()
                full.assert_not_called()

    def test_full_pipeline_mode(self):
        sp = _planner_with_config(
            f"[semester]\nmode = {_ALLOWED_MODES[1]}\n"
        )
        with patch.object(
            sp, "run_model_shortfall_balance_prioritize_fillempty_fillcurrentday"
        ) as full:
            with patch.object(sp, "run_model_shortfall") as shortfall:
                sp.run_model()
                full.assert_called_once()
                shortfall.assert_not_called()

    def test_unknown_mode_raises(self):
        sp = _planner_with_config("[semester]\nmode = bonus\n")
        with self.assertRaisesRegex(ValueError, "mode='bonus' invalid"):
            sp.run_model()

    def test_legacy_round1_raises(self):
        sp = _planner_with_config("[semester]\nmode = round1\n")
        with self.assertRaisesRegex(ValueError, "mode='round1' invalid"):
            sp.run_model()


if __name__ == "__main__":
    unittest.main()
