"""Tests for astroq.io input contracts and shared prep writers."""

import os
import shutil
import tempfile
import unittest

import pandas as pd
import pytest
from astropy.time import Time

import astroq.io
from astroq.io import (
    CUSTOM_SCHEMA,
    DEFAULT_MAX_FILL,
    DEFAULT_MIN_FILL,
    PAST_SCHEMA,
    read_csv,
    validate_past_in_semester,
)
from astroq.queue.prep_common import PROGRAMS_COLS, write_programs_csv


def _write_csv(df):
    path = os.path.join(tempfile.mkdtemp(prefix="astroq_schema_"), "input.csv")
    df.to_csv(path, index=False)
    return path


def _request_row(**overrides):
    row = {
        "unique_id": "T1",
        "target": "T1",
        "program_code": "2026A_X001",
        "ra": 10.0,
        "dec": 20.0,
        "exptime": 600,
        "n_exp": 1,
        "n_inter_max": 5,
        "tau_inter": 1,
        "n_intra_min": 1,
        "n_intra_max": 1,
        "tau_intra": 0,
        "inactive": False,
        "splan_weight": 2,
    }
    row.update(overrides)
    return row


class TestReadCsv(unittest.TestCase):
    def test_request_coerces_dtypes(self):
        path = _write_csv(pd.DataFrame([_request_row(unique_id=101)]))
        df = read_csv(path, "request")
        self.assertEqual(df["unique_id"].tolist(), ["101"])
        self.assertEqual(df["exptime"].dtype.kind, "f")
        self.assertEqual(df["n_exp"].dtype.kind, "i")
        self.assertEqual(df["inactive"].dtype.kind, "b")
        self.assertEqual(df["splan_weight"].dtype.kind, "f")

    def test_missing_column_raises(self):
        row = _request_row()
        del row["splan_weight"]
        path = _write_csv(pd.DataFrame([row]))
        with self.assertRaisesRegex(ValueError, "splan_weight"):
            read_csv(path, "request")

    def test_extra_columns_pass_through(self):
        path = _write_csv(pd.DataFrame([_request_row(comments="hi")]))
        df = read_csv(path, "request")
        self.assertIn("comments", df.columns)

    def test_duplicate_active_unique_id_raises(self):
        rows = [
            _request_row(unique_id="dup", inactive=False),
            _request_row(unique_id="dup", inactive=False, target="T2"),
        ]
        path = _write_csv(pd.DataFrame(rows))
        with self.assertRaisesRegex(ValueError, "Duplicate unique_id among active"):
            read_csv(path, "request")

    def test_duplicate_key_raises(self):
        progs = pd.DataFrame(
            {"program": ["P1", "P1"], "hours": [10.0, 20.0]}
        )
        path = _write_csv(progs)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            read_csv(path, "programs")

    def test_programs_fill_defaults(self):
        path = _write_csv(pd.DataFrame({"program": ["P1"], "hours": [12.0]}))
        df = read_csv(path, "programs")
        self.assertEqual(df.loc["P1", "min_fill"], DEFAULT_MIN_FILL)
        self.assertEqual(df.loc["P1", "max_fill"], DEFAULT_MAX_FILL)

    def test_programs_custom_max_fill(self):
        path = _write_csv(
            pd.DataFrame(
                {
                    "program": ["P1"],
                    "hours": [12.0],
                    "max_fill": [2.0],
                }
            )
        )
        df = read_csv(path, "programs")
        self.assertEqual(df.loc["P1", "max_fill"], 2.0)

    def test_programs_max_feasible_fill_absent_is_nan(self):
        path = _write_csv(pd.DataFrame({"program": ["P1"], "hours": [12.0]}))
        df = read_csv(path, "programs")
        self.assertTrue(pd.isna(df.loc["P1", "max_feasible_fill"]))

    def test_programs_max_feasible_fill_empty_cell_is_nan(self):
        path = _write_csv(
            pd.DataFrame(
                {
                    "program": ["P1", "P2"],
                    "hours": [12.0, 12.0],
                    "max_feasible_fill": [None, 0.0],
                }
            )
        )
        df = read_csv(path, "programs")
        self.assertTrue(pd.isna(df.loc["P1", "max_feasible_fill"]))
        self.assertEqual(df.loc["P2", "max_feasible_fill"], 0.0)

    def test_missing_file_raises_without_empty_ok(self):
        with self.assertRaises(FileNotFoundError):
            read_csv("/nonexistent/request.csv", "request")

    def test_empty_ok_missing_file(self):
        df = read_csv("/nonexistent/past.csv", "past")
        self.assertTrue(df.empty)
        self.assertEqual(list(df.columns), list(PAST_SCHEMA))

    def test_empty_ok_header_only(self):
        path = _write_csv(pd.DataFrame(columns=list(CUSTOM_SCHEMA)))
        df = read_csv(path, "custom")
        self.assertTrue(df.empty)

    def test_time_columns_parsed(self):
        alloc = pd.DataFrame(
            {"start": ["2026-02-01T05:00"], "stop": ["2026-02-01T15:00"]}
        )
        path = _write_csv(alloc)
        df = read_csv(path, "allocation")
        self.assertIsInstance(df["start"].iloc[0], Time)
        self.assertIsInstance(df["stop"].iloc[0], Time)

    def test_validate_past_in_semester_accepts_in_range(self):
        past = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "target": ["R1"],
                "timestamp": ["2026-02-15 12:00"],
                "exposure_time": [100.0],
            }
        )
        validate_past_in_semester(past, "2026-02-01", "2026-07-31")

    def test_validate_past_in_semester_empty_ok(self):
        validate_past_in_semester(
            pd.DataFrame(columns=list(PAST_SCHEMA)),
            "2026-02-01",
            "2026-07-31",
        )

    def test_validate_past_in_semester_rejects_oos(self):
        past = pd.DataFrame(
            {
                "unique_id": ["R1"],
                "target": ["R1"],
                "timestamp": ["2025-12-01 12:00"],
                "exposure_time": [100.0],
            }
        )
        with self.assertRaisesRegex(ValueError, "outside semester date range"):
            validate_past_in_semester(past, "2026-02-01", "2026-07-31")


@pytest.mark.parametrize(
    "override,error_pattern",
    [
        ({"exptime": None}, "null values.*exptime"),
        ({"n_intra_max": "None"}, None),
        ({"n_inter_max": 1.5}, "non-integer"),
        ({"inactive": "maybe"}, "boolean"),
    ],
)
def test_request_validation_errors(override, error_pattern):
    row = _request_row(**override)
    path = _write_csv(pd.DataFrame([row]))
    with pytest.raises(ValueError):
        read_csv(path, "request")


class TestWriteProgramsCsv(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="astroq_programs_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def write(self, df):
        write_programs_csv(df, self.tmp)
        path = os.path.join(self.tmp, "programs.csv")
        return path, pd.read_csv(path)

    def test_fill_defaults_are_written(self):
        _, out = self.write(pd.DataFrame({"program": ["P1"], "hours": [10.0]}))
        self.assertEqual(list(out.columns), list(PROGRAMS_COLS))
        self.assertEqual(out.loc[0, "min_fill"], astroq.io.DEFAULT_MIN_FILL)
        self.assertEqual(out.loc[0, "max_fill"], astroq.io.DEFAULT_MAX_FILL)

    def test_max_feasible_fill_written_empty(self):
        path, out = self.write(pd.DataFrame({"program": ["P1"], "hours": [10.0]}))
        self.assertTrue(pd.isna(out.loc[0, "max_feasible_fill"]))
        with open(path, encoding="utf-8") as fh:
            lines = fh.read().splitlines()
        self.assertTrue(lines[1].endswith(","), lines[1])

    def test_caller_values_preserved(self):
        _, out = self.write(
            pd.DataFrame(
                {
                    "program": ["P1"],
                    "hours": [10.0],
                    "min_fill": [0.4],
                    "max_fill": [0.9],
                }
            )
        )
        self.assertEqual(out.loc[0, "min_fill"], 0.4)
        self.assertEqual(out.loc[0, "max_fill"], 0.9)

    def test_extra_columns_pass_through_after_known_ones(self):
        _, out = self.write(
            pd.DataFrame({"program": ["P1"], "hours": [10.0], "priority": [1.0]})
        )
        self.assertEqual(list(out.columns), list(PROGRAMS_COLS) + ["priority"])

    def test_round_trips_through_read_csv(self):
        path, _ = self.write(
            pd.DataFrame({"program": ["P1", "P2"], "hours": [10.0, 20.0]})
        )
        df = read_csv(path, "programs")
        self.assertEqual(df.loc["P1", "min_fill"], astroq.io.DEFAULT_MIN_FILL)
        self.assertEqual(df.loc["P1", "max_fill"], astroq.io.DEFAULT_MAX_FILL)
        self.assertTrue(df["max_feasible_fill"].isna().all())


if __name__ == "__main__":
    unittest.main()
