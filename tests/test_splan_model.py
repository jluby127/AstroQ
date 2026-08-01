"""Tests for astroq.splan: mode dispatch, fill-factor bounds, balance, throttle."""

import unittest
from configparser import ConfigParser, NoOptionError
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from gurobipy import GRB

from astroq.splan import _ALLOWED_MODES, _MODE_PIPELINES, SemesterPlanner

CONFIG = """
[global]
current_day = 2018-08-05
[semester]
mode = shortfall,balance,prioritize,fill-empty,fill-current-day
"""


def _planner_with_config(text):
    sp = SemesterPlanner.__new__(SemesterPlanner)
    sp.config = ConfigParser()
    sp.config.optionxform = str
    sp.config.read_string(text)
    return sp


def _planner_with_programs(programs_df, requests_df):
    sp = SemesterPlanner.__new__(SemesterPlanner)
    sp.config = ConfigParser()
    sp.config.optionxform = str
    sp.config.read_string(
        """
        [semester]
        slot_size = 20
        hours_per_night = 12
        """
    )
    sp.requests = requests_df.copy()
    sp.requests["r"] = sp.requests["unique_id"]
    sp.past = pd.DataFrame(columns=["unique_id", "target", "timestamp", "exposure_time"])
    sp.past["r"] = sp.past["unique_id"]
    sp.programs = programs_df.copy()
    sp._add_program_columns()
    return sp


def _stub_balance_planner(max_feasible_fill):
    sp = SemesterPlanner.__new__(SemesterPlanner)
    sp.config = ConfigParser()
    sp.config.optionxform = str
    sp.config.read_string(CONFIG)
    sp.programs = pd.DataFrame(
        {"max_feasible_fill": list(max_feasible_fill.values())},
        index=pd.Index(list(max_feasible_fill), name="program"),
    )
    sp.program_keys = sp.programs.index
    sp.F = {p: MagicMock() for p in max_feasible_fill}
    sp.theta = {}
    sp.model = MagicMock()
    sp.model.SolCount = 1
    sp.model.ObjVal = 0.0
    sp.model.getAttr.return_value = {p: 0.0 for p in max_feasible_fill}
    fsf_var = MagicMock()
    fsf_var.__ge__ = MagicMock(return_value=MagicMock())
    sp.model.addVar.return_value = fsf_var
    sp.access_obj = MagicMock(current_night_index=4)
    return sp


def _run_balance_pipeline(sp):
    with patch.object(sp, "optimize_model"), patch.object(
        sp, "build_schedule"
    ), patch.object(sp, "log_report"), patch.object(
        sp, "write_request_selected"
    ), patch.object(
        sp, "to_hdf5"
    ), patch.object(
        sp, "_constraint_fillfactor"
    ), patch.object(
        sp, "_objective_prioritize_intra", return_value=0
    ), patch.object(
        sp, "_objective_minimize_empty_slots", return_value=0
    ), patch.object(
        sp, "_objective_weighted_theta", return_value=0
    ), patch.object(
        sp, "_objective_slots_used_tonight", return_value=0
    ):
        sp.run_model_shortfall_balance_prioritize_fillempty_fillcurrentday()


def _make_throttle_planner(requests, past):
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


def _throttle_requests():
    return pd.DataFrame(
        {
            "unique_id": ["ACT", "INACT"],
            "program_code": ["2026A_X001", "2026A_X001"],
            "inactive": [False, True],
            "t_visit_slots": [3, 3],
        }
    )


class TestRunModelDispatch(unittest.TestCase):
    def test_mode_pipelines(self):
        self.assertEqual(
            _MODE_PIPELINES,
            {
                "shortfall": "run_model_shortfall",
                "shortfall,balance,prioritize,fill-empty,fill-current-day": (
                    "run_model_shortfall_balance_prioritize_fillempty_fillcurrentday"
                ),
                "shortfall,prioritize,fill-empty,fill-current-day": (
                    "run_model_shortfall_prioritize_fillempty_fillcurrentday"
                ),
            },
        )
        self.assertEqual(
            _ALLOWED_MODES,
            [
                "shortfall",
                "shortfall,balance,prioritize,fill-empty,fill-current-day",
                "shortfall,prioritize,fill-empty,fill-current-day",
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
        sp = _planner_with_config(f"[semester]\nmode = {_ALLOWED_MODES[1]}\n")
        with patch.object(
            sp, "run_model_shortfall_balance_prioritize_fillempty_fillcurrentday"
        ) as full:
            with patch.object(sp, "run_model_shortfall") as shortfall:
                sp.run_model()
                full.assert_called_once()
                shortfall.assert_not_called()

    def test_no_balance_mode(self):
        sp = _planner_with_config(f"[semester]\nmode = {_ALLOWED_MODES[2]}\n")
        with patch.object(sp, "_run_pipeline") as pipeline:
            sp.run_model()
            pipeline.assert_called_once_with(balance=False)

    def test_full_pipeline_mode_enables_balance(self):
        sp = _planner_with_config(f"[semester]\nmode = {_ALLOWED_MODES[1]}\n")
        with patch.object(sp, "_run_pipeline") as pipeline:
            sp.run_model()
            pipeline.assert_called_once_with(balance=True)

    def test_unknown_mode_raises(self):
        sp = _planner_with_config("[semester]\nmode = bonus\n")
        with self.assertRaisesRegex(ValueError, "mode='bonus' invalid"):
            sp.run_model()

    def test_legacy_round1_raises(self):
        sp = _planner_with_config("[semester]\nmode = round1\n")
        with self.assertRaisesRegex(ValueError, "mode='round1' invalid"):
            sp.run_model()


class TestShortfallPipeline(unittest.TestCase):
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


class TestValidateProgramCoverage(unittest.TestCase):
    def test_missing_program_raises(self):
        sp = _planner_with_programs(
            pd.DataFrame({"hours": [10.0]}, index=pd.Index(["P1"], name="program")),
            pd.DataFrame(
                {
                    "unique_id": ["R1"],
                    "program_code": ["P2"],
                    "inactive": [False],
                    "t_visit_slots": [1],
                }
            ),
        )
        with self.assertRaisesRegex(ValueError, "missing from programs.csv"):
            sp._validate_program_coverage()


class TestConstraintFillfactor(unittest.TestCase):
    def _stub_planner(self, *, past_slots, awarded_slots, min_ff=0.0, max_ff=1.25):
        sp = SemesterPlanner.__new__(SemesterPlanner)
        sp.programs = pd.DataFrame(
            {
                "min_fill": [min_ff],
                "max_fill": [max_ff],
                "awarded_slots": [awarded_slots],
                "past_slots": [past_slots],
            },
            index=pd.Index(["P1"], name="program"),
        )
        f_var = MagicMock()
        f_var.LB = 0.0
        f_var.UB = 1.25
        sp.F = {"P1": f_var}
        return sp, f_var

    def test_csv_baseline_bounds(self):
        sp, f_var = self._stub_planner(past_slots=0, awarded_slots=100)
        sp._constraint_fillfactor()
        self.assertEqual(f_var.LB, 0.0)
        self.assertEqual(f_var.UB, 1.25)

    def test_max_fill_override(self):
        sp, f_var = self._stub_planner(past_slots=0, awarded_slots=100)
        sp._constraint_fillfactor(max_fill=1.0)
        self.assertEqual(f_var.UB, 1.0)

    def test_min_fill_override(self):
        sp, f_var = self._stub_planner(past_slots=0, awarded_slots=100)
        sp._constraint_fillfactor(min_fill=0.5)
        self.assertEqual(f_var.LB, 0.5)

    def test_overpast_clamp_pins_ub_and_logs_warning(self):
        sp, f_var = self._stub_planner(
            past_slots=150, awarded_slots=100, max_ff=1.0
        )
        with self.assertLogs("astroq.splan", level="WARNING") as logs:
            sp._constraint_fillfactor(max_fill=1.0)
        self.assertEqual(f_var.UB, 1.5)
        self.assertTrue(any("over ceiling from past alone" in m for m in logs.output))

    def test_lb_capped_to_ub_when_floor_exceeds_clamped_ceiling(self):
        sp, f_var = self._stub_planner(
            past_slots=150, awarded_slots=100, min_ff=2.0, max_ff=1.0
        )
        sp._constraint_fillfactor(max_fill=1.0)
        self.assertEqual(f_var.UB, 1.5)
        self.assertEqual(f_var.LB, 1.5)


class TestZeroAwardExcluded(unittest.TestCase):
    def test_zero_award_not_in_program_keys(self):
        requests = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "program_code": ["P1"],
                "inactive": [False],
                "t_visit_slots": [1],
                "n_intra_max": [1],
                "n_intra_min": [1],
                "n_inter_max": [1],
                "tau_inter": [1],
                "tau_intra_slots": [0],
            }
        )
        programs = pd.DataFrame(
            {"hours": [0.0, 10.0]},
            index=pd.Index(["P0", "P1"], name="program"),
        )
        sp = _planner_with_programs(programs, requests)
        self.assertNotIn("P0", sp.programs.index[sp.programs["awarded_slots"] > 0])


class TestFillFactorIntegration(unittest.TestCase):
    def test_shortfall_caps_ub_at_one(self):
        sp = SemesterPlanner("examples/hello_world/config_hello_world.ini")
        sp._constraint_fillfactor(max_fill=1.0)
        ub = pd.Series(sp.model.getAttr("UB", sp.F))
        self.assertTrue((ub <= 1.0 + 1e-9).all())

    def test_fill_factor_linking_after_shortfall(self):
        sp = SemesterPlanner("examples/hello_world/config_hello_world.ini")
        sp._constraint_fillfactor(max_fill=1.0)
        sp.model.setObjective(sp._objective_weighted_theta(), GRB.MINIMIZE)
        sp.optimize_model("shortfall")
        self.assertGreater(sp.model.SolCount, 0)
        f_proj = pd.Series(sp.model.getAttr("X", sp.F))
        awarded = sp.programs["awarded_slots"]
        past = sp.programs["past_slots"]
        sched = pd.Series(
            {
                p: expr.getValue()
                for p, expr in sp.sched_slots_by_program.items()
            }
        ).reindex(f_proj.index).fillna(0.0)
        expected = (past.reindex(f_proj.index) + sched) / awarded.reindex(
            f_proj.index
        )
        pd.testing.assert_series_equal(f_proj, expected, atol=1e-6)


@pytest.mark.parametrize(
    "max_feasible_fill,program_keys,expect_error,msg_parts",
    [
        ({"P1": np.nan, "P2": np.nan}, None, ValueError, ["['P1', 'P2']"]),
        ({"P1": 0.8, "P2": np.nan}, None, ValueError, ["['P2']", "compute-max-fill"]),
        ({"P1": 0.0}, None, None, []),
        ({"P1": 0.8, "P_noaward": np.nan}, ["P1"], None, []),
    ],
)
def test_balance_guard_max_feasible_fill(
    max_feasible_fill, program_keys, expect_error, msg_parts
):
    sp = _stub_balance_planner(max_feasible_fill)
    if program_keys is not None:
        sp.program_keys = pd.Index(program_keys, name="program")
        sp.F = {p: sp.F[p] for p in program_keys}
    if expect_error is not None:
        with pytest.raises(expect_error) as exc_info:
            _run_balance_pipeline(sp)
        for part in msg_parts:
            assert part in str(exc_info.value)
    else:
        _run_balance_pipeline(sp)


class TestBalanceFormulation(unittest.TestCase):
    def test_objective_is_minimize_fsf_max(self):
        sp = _stub_balance_planner({"P1": 0.8, "P2": 0.5})
        _run_balance_pipeline(sp)

        sense = [c.args[1] for c in sp.model.setObjective.call_args_list]
        self.assertEqual(sense[1], GRB.MINIMIZE)
        self.assertIs(sp.model.setObjective.call_args_list[1].args[0], sp.fsf_max)

    def test_one_fsf_constraint_per_program(self):
        sp = _stub_balance_planner({"P1": 0.8, "P2": 0.5})
        _run_balance_pipeline(sp)

        named = [
            c.args[1]
            for c in sp.model.addConstr.call_args_list
            if len(c.args) > 1 and str(c.args[1]).startswith("fsf_")
        ]
        self.assertEqual(sorted(named), ["fsf_P1", "fsf_P2"])

    def test_fsf_max_floored_at_zero(self):
        sp = _stub_balance_planner({"P1": 0.8})
        _run_balance_pipeline(sp)
        self.assertEqual(sp.model.addVar.call_args.kwargs["lb"], 0.0)

    def test_shortfall_cap_uses_balance_slack(self):
        sp = _stub_balance_planner({"P1": 0.8})
        sp.config.read_string(
            "[semester.balance]\nglobal_shortfall_slack = 1.4\n"
        )
        sp.model.ObjVal = 100.0
        with patch.object(sp, "_objective_weighted_theta", return_value=0):
            _run_balance_pipeline(sp)

        capped = [
            c
            for c in sp.model.addConstr.call_args_list
            if len(c.args) > 1 and c.args[1] == "balance_shortfall_cap"
        ]
        self.assertEqual(len(capped), 1)


class TestBalanceIntegration(unittest.TestCase):
    def test_fsf_max_is_the_largest_gap(self):
        sp = SemesterPlanner("examples/hello_world/config_hello_world.ini")
        sp._constraint_fillfactor(max_fill=1.0)
        sp.model.setObjective(sp._objective_weighted_theta(), GRB.MINIMIZE)
        sp.optimize_model("shortfall")
        f_shortfall = pd.Series(sp.model.getAttr("X", sp.F))

        maff = pd.Series(1.0, index=f_shortfall.index)
        bottleneck = maff.index[0]
        maff[bottleneck] = 0.01

        sp._constraint_fillfactor(min_fill=f_shortfall, max_fill=1.0)
        fsf_max = sp.model.addVar(lb=0.0, name="fsf_max")
        for p in sp.F:
            sp.model.addConstr(fsf_max >= float(maff[p]) - sp.F[p], f"fsf_{p}")
        sp.model.setObjective(fsf_max, GRB.MINIMIZE)
        sp.optimize_model("balance")

        f_balance = pd.Series(sp.model.getAttr("X", sp.F))
        gaps = (maff - f_balance).clip(lower=0.0)
        self.assertAlmostEqual(fsf_max.X, gaps.max(), places=6)
        self.assertLess(gaps[bottleneck], fsf_max.X)

    def test_balance_redistributes_on_a_full_telescope(self):
        sp = SemesterPlanner(
            "examples/priorities/symmetric_toy_model/config_benchmark.ini"
        )
        sp._constraint_fillfactor(max_fill=1.0)
        sp.model.setObjective(sp._objective_weighted_theta(), GRB.MINIMIZE)
        sp.optimize_model("shortfall")
        theta_min = sp.model.ObjVal
        f_shortfall = pd.Series(sp.model.getAttr("X", sp.F))
        maff = sp.programs["max_feasible_fill"].reindex(f_shortfall.index)

        sp.model.addConstr(
            sp._objective_weighted_theta() <= theta_min * 1.1, "cap"
        )
        sp._constraint_fillfactor(max_fill=1.0)
        fsf_max = sp.model.addVar(lb=0.0, name="fsf_max")
        for p in sp.F:
            sp.model.addConstr(fsf_max >= float(maff[p]) - sp.F[p], f"fsf_{p}")
        sp.model.setObjective(fsf_max, GRB.MINIMIZE)
        sp.optimize_model("balance")

        f_balance = pd.Series(sp.model.getAttr("X", sp.F))
        gap_before = (maff - f_shortfall).max()
        gap_after = (maff - f_balance).max()
        self.assertLess(gap_after, gap_before)
        self.assertLess(
            f_balance.max() - f_balance.min(),
            f_shortfall.max() - f_shortfall.min(),
        )


class TestPastSlotsByProgram(unittest.TestCase):
    def test_includes_inactive_past(self):
        rfa = _throttle_requests()
        past = pd.DataFrame(
            {
                "unique_id": ["ACT", "ACT", "INACT", "INACT", "INACT", "INACT"],
                "timestamp": ["2026-02-01T10:00"] * 6,
            }
        )
        sp = _make_throttle_planner(rfa, past)
        self.assertEqual(sp.programs.loc["2026A_X001", "past_slots"], 18)

    def test_inactive_only_still_counts(self):
        rfa = _throttle_requests()
        past = pd.DataFrame(
            {
                "unique_id": ["INACT", "INACT"],
                "timestamp": ["2026-02-01T10:00", "2026-02-02T10:00"],
            }
        )
        sp = _make_throttle_planner(rfa, past)
        self.assertEqual(sp.programs.loc["2026A_X001", "past_slots"], 6)

    def test_empty_past(self):
        rfa = _throttle_requests()
        past = pd.DataFrame(columns=["unique_id", "timestamp"])
        sp = _make_throttle_planner(rfa, past)
        self.assertEqual(sp.programs.loc["2026A_X001", "past_slots"], 0)

    def test_integer_unique_id_join(self):
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

        sp = _make_throttle_planner(rfa, past)
        self.assertEqual(sp.programs.loc["2026A_X001", "past_slots"], 6)


if __name__ == "__main__":
    unittest.main()
