"""Tests for the [semester] mode switch and pipeline dispatch."""

import unittest
from configparser import ConfigParser, NoOptionError
from unittest.mock import MagicMock, patch

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

    def test_full_pipeline_calls_all_steps(self):
        sp = _planner_with_config(
            """
            [global]
            current_day = 2018-08-05
            [semester]
            mode = shortfall,balance,prioritize,fill-empty,fill-current-day
            [semester.fill-current-day]
            global_shortfall_slack = 1.1
            """
        )
        step_order = []

        def record(step):
            step_order.append(step)

        sp.model = MagicMock()
        sp.model.SolCount = 1
        sp.model.ObjVal = 0.0
        sp.F = MagicMock()
        sp.F.sum.return_value = 0
        sp.theta = {}
        sp.access_obj = MagicMock(current_night_index=4)

        with patch.object(sp, "optimize_model", side_effect=record):
            with patch.object(sp, "build_schedule"):
                with patch.object(sp, "log_report"):
                    with patch.object(sp, "write_request_selected"):
                        with patch.object(sp, "to_hdf5"):
                            with patch.object(sp, "_constraint_fillfactor"):
                                with patch.object(
                                    sp,
                                    "_objective_prioritize_intra",
                                    return_value=0,
                                ):
                                    with patch.object(
                                        sp,
                                        "_objective_minimize_empty_slots",
                                        return_value=0,
                                    ):
                                        with patch.object(
                                            sp,
                                            "_objective_weighted_theta",
                                            return_value=0,
                                        ):
                                            with patch.object(
                                                sp,
                                                "_objective_slots_used_tonight",
                                                return_value=0,
                                            ):
                                                sp.run_model_shortfall_balance_prioritize_fillempty_fillcurrentday()

        self.assertEqual(
            step_order,
            [
                "shortfall",
                "balance",
                "prioritize",
                "fill-empty",
                "fill-current-day",
            ],
        )


class TestShortfallPipeline(unittest.TestCase):
    """``solve_shortfall`` is the side-effect-free core of the shortfall run."""

    def test_solve_shortfall_writes_nothing(self):
        sp = _planner_with_config("[semester]\n")
        sp.model = MagicMock()
        with patch.object(sp, "_constraint_fillfactor"), patch.object(
            sp, "_objective_weighted_theta", return_value=0
        ), patch.object(sp, "optimize_model") as optimize, patch.object(
            sp, "build_schedule"
        ) as build_schedule, patch.object(
            sp, "log_report"
        ) as log_report, patch.object(
            sp, "write_request_selected"
        ) as write_selected, patch.object(
            sp, "to_hdf5"
        ) as to_hdf5:
            sp.solve_shortfall()

            optimize.assert_called_once_with("shortfall")
            build_schedule.assert_not_called()
            log_report.assert_not_called()
            write_selected.assert_not_called()
            to_hdf5.assert_not_called()

    def test_run_model_shortfall_delegates_then_writes(self):
        sp = _planner_with_config("[semester]\n")
        with patch.object(sp, "solve_shortfall") as solve, patch.object(
            sp, "build_schedule"
        ) as build_schedule, patch.object(
            sp, "log_report"
        ) as log_report, patch.object(
            sp, "write_request_selected"
        ) as write_selected, patch.object(
            sp, "to_hdf5"
        ) as to_hdf5:
            sp.run_model_shortfall()

            solve.assert_called_once_with()
            build_schedule.assert_called_once_with()
            log_report.assert_called_once_with("shortfall")
            write_selected.assert_called_once_with()
            to_hdf5.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
