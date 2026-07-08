"""Tests for the astroq.io input contracts (schema dicts + read_csv).

read_csv validates and coerces only -- repair belongs to the prep stage --
so every nonconforming input must raise, and conforming inputs must come back
with pinned dtypes.
"""

import os
import tempfile
import unittest

import pandas as pd
from astropy.time import Time

from astroq.io import (
    ALLOCATION_SCHEMA,
    CUSTOM_SCHEMA,
    PAST_SCHEMA,
    PROGRAMS_SCHEMA,
    REQUEST_SCHEMA,
    read_csv,
)


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

    def test_null_value_raises(self):
        path = _write_csv(pd.DataFrame([_request_row(exptime=None)]))
        with self.assertRaisesRegex(ValueError, "null values.*exptime"):
            read_csv(path, "request")

    def test_legacy_none_string_raises(self):
        path = _write_csv(pd.DataFrame([_request_row(n_intra_max="None")]))
        with self.assertRaises(ValueError):
            read_csv(path, "request")

    def test_non_integer_raises(self):
        path = _write_csv(pd.DataFrame([_request_row(n_inter_max=1.5)]))
        with self.assertRaisesRegex(ValueError, "non-integer"):
            read_csv(path, "request")

    def test_non_boolean_inactive_raises(self):
        path = _write_csv(pd.DataFrame([_request_row(inactive="maybe")]))
        with self.assertRaisesRegex(ValueError, "boolean"):
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
            {"program": ["P1", "P1"], "hours": [10.0, 20.0], "nights": [1.0, 2.0]}
        )
        path = _write_csv(progs)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            read_csv(path, "programs")

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


if __name__ == "__main__":
    unittest.main()
