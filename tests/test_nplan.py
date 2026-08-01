"""Tests for astroq.nplan: TTP config resolution and nplan_weight handoff."""

import tempfile
import unittest
from configparser import ConfigParser
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd

from astroq.nplan import (
    _AUTO_ACS_TARGET_THRESHOLD,
    _resolve_nplan_weights,
    _resolve_night_ttp_config,
    _resolve_slew_fn,
)
from astroq.queue.hirescps.queue import HIRESCPS


def _night_cfg(**kwargs) -> ConfigParser:
    cfg = ConfigParser()
    cfg.add_section("night")
    for key, value in kwargs.items():
        cfg.set("night", key, str(value))
    return cfg


class TestResolveSlewFn(unittest.TestCase):
    def setUp(self):
        self.queue = HIRESCPS()

    def test_single_state(self):
        fn = _resolve_slew_fn(self.queue, 1)
        self.assertEqual(fn.__name__, "slew_fn")

    def test_two_state(self):
        fn = _resolve_slew_fn(self.queue, 2)
        self.assertEqual(fn.__name__, "slew_fn_state")

    def test_invalid_n_states(self):
        with self.assertRaises(ValueError):
            _resolve_slew_fn(self.queue, 3)


class TestResolveNightTtpConfig(unittest.TestCase):
    def setUp(self):
        self.queue = HIRESCPS()

    def test_milp_defaults(self):
        cfg = _resolve_night_ttp_config(self.queue, _night_cfg(), n_targets=10)
        self.assertEqual(cfg["method"], "milp")
        self.assertEqual(cfg["n_states"], 1)
        self.assertEqual(cfg["milp_time"], 300.0)
        self.assertEqual(cfg["warmstart_time"], 0.0)

    def test_auto_small_queue(self):
        cfg = _resolve_night_ttp_config(
            self.queue,
            _night_cfg(method="auto", max_solve_time=300),
            n_targets=_AUTO_ACS_TARGET_THRESHOLD - 1,
        )
        self.assertEqual(cfg["method"], "milp")

    def test_auto_large_queue(self):
        cfg = _resolve_night_ttp_config(
            self.queue,
            _night_cfg(method="auto", warmstart_time=30, max_solve_time=300),
            n_targets=_AUTO_ACS_TARGET_THRESHOLD,
        )
        self.assertEqual(cfg["method"], "acs8+milp")
        self.assertEqual(cfg["warmstart_time"], 30.0)
        self.assertEqual(cfg["milp_time"], 270.0)

    def test_acs8_milp_time_split(self):
        cfg = _resolve_night_ttp_config(
            self.queue,
            _night_cfg(
                method="acs8+milp",
                n_states=2,
                warmstart_time=30,
                max_solve_time=120,
            ),
            n_targets=48,
        )
        self.assertEqual(cfg["method"], "acs8+milp")
        self.assertEqual(cfg["n_states"], 2)
        self.assertEqual(cfg["acs_starts"], 8)
        self.assertEqual(cfg["warmstart_time"], 30.0)
        self.assertEqual(cfg["milp_time"], 90.0)

    def test_norel_milp(self):
        cfg = _resolve_night_ttp_config(
            self.queue,
            _night_cfg(method="norel+milp", warmstart_time=60, max_solve_time=300),
            n_targets=30,
        )
        self.assertEqual(cfg["method"], "norel+milp")
        self.assertEqual(cfg["warmstart_time"], 60.0)
        self.assertEqual(cfg["milp_time"], 240.0)

    def test_invalid_method(self):
        with self.assertRaises(ValueError):
            _resolve_night_ttp_config(
                self.queue, _night_cfg(method="bogus"), n_targets=5
            )

    def test_warmstart_exceeds_total(self):
        with self.assertRaises(ValueError):
            _resolve_night_ttp_config(
                self.queue,
                _night_cfg(method="acs8+milp", warmstart_time=120, max_solve_time=60),
                n_targets=30,
            )


class TestResolveNplanWeights(unittest.TestCase):
    def test_missing_column_defaults_to_one(self):
        df = pd.DataFrame({"unique_id": ["A", "B"]})
        weights = _resolve_nplan_weights(df)
        np.testing.assert_array_equal(weights, [1.0, 1.0])

    def test_parses_explicit_values(self):
        df = pd.DataFrame({"unique_id": ["A", "B", "C"], "nplan_weight": [1.0, 3.0, 2.5]})
        weights = _resolve_nplan_weights(df)
        np.testing.assert_array_equal(weights, [1.0, 3.0, 2.5])

    def test_coerces_bad_values_to_one(self):
        df = pd.DataFrame({"unique_id": ["A", "B"], "nplan_weight": ["bad", 2.0]})
        weights = _resolve_nplan_weights(df)
        np.testing.assert_array_equal(weights, [1.0, 2.0])


class TestWriteRequestSelected(unittest.TestCase):
    def test_includes_nplan_weight_default(self):
        from astroq.splan import SemesterPlanner

        tmpdir = tempfile.mkdtemp()
        sp = SemesterPlanner.__new__(SemesterPlanner)
        sp.config = MagicMock()
        sp.config.get.side_effect = lambda section, key, **kw: {
            ("global", "workdir"): tmpdir,
            ("global", "current_day"): "2026-06-22",
        }[(section, key)]
        sp.access_obj = MagicMock(current_night_index=0)
        sp.Yrds = {("T1", 0): MagicMock(x=1), ("T2", 0): MagicMock(x=0)}
        sp.requests = pd.DataFrame(
            {
                "unique_id": ["T1", "T2"],
                "r": ["T1", "T2"],
                "target": ["star1", "star2"],
                "priority": ["p1", "p2"],
                "inactive": [False, False],
            }
        )

        sp.write_request_selected()
        out = pd.read_csv(Path(tmpdir) / "outputs" / "request_selected.csv")
        self.assertEqual(len(out), 1)
        self.assertEqual(out.iloc[0]["unique_id"], "T1")
        self.assertIn("nplan_weight", out.columns)
        self.assertEqual(float(out.iloc[0]["nplan_weight"]), 1.0)


if __name__ == "__main__":
    unittest.main()
