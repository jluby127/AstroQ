"""Tests for astroq.queue.hirescps.prep."""

import io
import contextlib
import tempfile
import unittest
from pathlib import Path
from unittest import mock
from unittest.mock import patch

import pandas as pd
import pytest

from astroq.queue.hirescps import prep as hirescps_prep
from astroq.queue.hirescps.prep import (
    PRIORITY_TO_SPLAN_WEIGHT,
    REQUEST_COLS,
    REQUEST_COLS_READ,
    SHEET_REQUEST_COLS,
    _canonicalize_sheet_columns,
    _fetch_sheet_dataframe,
    _human_header_labels,
    attach_splan_weight,
    jump_query_to_past,
    priority_token_to_splan_weight,
)

_H382_SHEET_CSV = """\
Program code,Star name,Unique ID,RA,Dec,Nominal time,Maximum time,Exposures per visit,Nights per semester,Minimum inter-night cadence,Max visits per night,Min visits per night,Min. intra-night cadence,Min. elevation,Min. moon seperation,Weather band 1,Weather band 2,Weather band 3,Gaia ID,Effective temperature,J-band mag,V-band mag,Proper motion in RA,Proper motion in Dec,Epoch,Exp meter threshold,Inactive,Decker,Iodine cell in or out?,Priority,start,end,comments
string,string,string,string,string,int,int,int,int,int,int,int,float,int,int,bool,bool,bool,string,int,float,float,float,float,int,string,bool,string,string,string,string,string
Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,,,,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,Pass,
program_code,starname,unique_id,ra,dec,exptime,maxtime,n_exp,n_inter_max,tau_inter,n_intra_max,n_intra_min,tau_intra,minimum_elevation,minimum_moon_separation,weather_band_1,weather_band_2,weather_band_3,gaia_id,teff,jmag,Vmag,pmra,pmdec,epoch,exp_meter_threshold,inactive,decker,cell in/out?,priority,start,stop,
2026A_H382,HIP94931,HIP94931,19 19 00.5,41 38 05.0,230,500,1,1,1,1,1,1,30,30,TRUE,TRUE,TRUE,,,,8.8,104.775,-68.016,2000,250k,FALSE,C2,out,p1,,,Kepler-444. Position the decker so that the companion does not contaminate the spectrum.
"""


class TestJumpPastRename(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def _write_request(self, rows):
        pd.DataFrame(rows).to_csv(self.tmp / "request.csv", index=False)

    def test_case_mismatch_renamed_to_request_id(self):
        self._write_request([
            {"unique_id": "GaiaDR3_123", "target": "GaiaDR3_123", "n_exp": 1},
        ])
        frames = pd.DataFrame([{
            "starname": "GAIADR3_123",
            "timestamp": "2026-05-29 06:46",
            "exposure_time": 1200,
            "decker": "C2",
            "iodine_in": False,
        }])
        out = jump_query_to_past(frames, str(self.tmp / "request.csv"))
        self.assertEqual(len(out), 1)
        self.assertEqual(out.iloc[0]["unique_id"], "GaiaDR3_123")

    def test_inactive_request_still_renames(self):
        self._write_request([
            {
                "unique_id": "GaiaDR3_999",
                "target": "GaiaDR3_999",
                "n_exp": 1,
                "inactive": True,
            },
        ])
        frames = pd.DataFrame([{
            "starname": "gaiadr3_999",
            "timestamp": "2026-05-29 06:46",
            "exposure_time": 1200,
            "decker": "C2",
            "iodine_in": False,
        }])
        out = jump_query_to_past(frames, str(self.tmp / "request.csv"))
        self.assertEqual(out.iloc[0]["unique_id"], "GaiaDR3_999")

    def test_unmatched_jump_name_unchanged(self):
        self._write_request([
            {"unique_id": "HR6410", "target": "HR6410", "n_exp": 1},
        ])
        frames = pd.DataFrame([{
            "starname": "UNKNOWN_CAL",
            "timestamp": "2026-05-29 06:46",
            "exposure_time": 1200,
            "decker": "C2",
            "iodine_in": False,
        }])
        out = jump_query_to_past(frames, str(self.tmp / "request.csv"))
        self.assertEqual(out.iloc[0]["unique_id"], "UNKNOWN_CAL")

    def test_unchanged_match_not_in_rename_printout(self):
        self._write_request([
            {"unique_id": "HR6410", "target": "HR6410", "n_exp": 1},
            {"unique_id": "GaiaDR3_123", "target": "GaiaDR3_123", "n_exp": 1},
        ])
        frames = pd.DataFrame([
            {
                "starname": "HR6410",
                "timestamp": "2026-05-29 06:46",
                "exposure_time": 1200,
                "decker": "C2",
                "iodine_in": False,
            },
            {
                "starname": "GAIADR3_123",
                "timestamp": "2026-05-29 07:00",
                "exposure_time": 1200,
                "decker": "C2",
                "iodine_in": False,
            },
        ])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            out = jump_query_to_past(frames, str(self.tmp / "request.csv"))
        printed = buf.getvalue()
        self.assertEqual(out.iloc[0]["unique_id"], "HR6410")
        self.assertEqual(out.iloc[1]["unique_id"], "GaiaDR3_123")
        self.assertNotIn("HR6410 ->", printed)
        self.assertIn("GAIADR3_123 -> GaiaDR3_123", printed)

    def test_semester_start_excludes_pre_semester_frames(self):
        self._write_request([
            {"unique_id": "T1", "target": "T1", "n_exp": 1},
        ])
        frames = pd.DataFrame([
            {
                "starname": "T1",
                "timestamp": "2026-01-31 14:15",
                "exposure_time": 1200,
                "decker": "C2",
                "iodine_in": False,
            },
            {
                "starname": "T1",
                "timestamp": "2026-02-01 06:46",
                "exposure_time": 1200,
                "decker": "C2",
                "iodine_in": False,
            },
        ])
        out = jump_query_to_past(
            frames,
            str(self.tmp / "request.csv"),
            semester_start_day="2026-02-01",
        )
        self.assertEqual(len(out), 1)
        self.assertEqual(out.iloc[0]["timestamp"], "2026-02-01 06:46")

    def test_current_day_frames_included(self):
        self._write_request([
            {"unique_id": "T1", "target": "T1", "n_exp": 1},
        ])
        frames = pd.DataFrame([
            {
                "starname": "T1",
                "timestamp": "2026-06-25 14:15",
                "exposure_time": 1200,
                "decker": "C2",
                "iodine_in": False,
            },
            {
                "starname": "T1",
                "timestamp": "2026-06-26 06:46",
                "exposure_time": 1200,
                "decker": "C2",
                "iodine_in": False,
            },
        ])
        out = jump_query_to_past(
            frames,
            str(self.tmp / "request.csv"),
            semester_start_day="2026-02-01",
        )
        self.assertEqual(len(out), 2)
        self.assertEqual(
            sorted(out["timestamp"].tolist()),
            ["2026-06-25 14:15", "2026-06-26 06:46"],
        )


class TestSheetCommentColumns(unittest.TestCase):
    def test_human_header_labels(self):
        labels = _human_header_labels(_H382_SHEET_CSV)
        self.assertEqual(labels[-1], "comments")

    def test_blank_machine_header_resolved_from_human_label(self):
        human = _human_header_labels(_H382_SHEET_CSV)
        df = pd.read_csv(io.StringIO(_H382_SHEET_CSV), skiprows=3, dtype=str)
        df.columns = [str(c).strip() for c in df.columns]
        self.assertNotIn("comments", df.columns)

        df = _canonicalize_sheet_columns(df, human)
        self.assertIn("comments", df.columns)
        hip = df[df["unique_id"] == "HIP94931"].iloc[0]
        self.assertIn("Kepler-444", hip["comments"])

    @mock.patch("astroq.queue.hirescps.prep.requests.get")
    def test_fetch_sheet_dataframe_preserves_comments(self, mock_get):
        mock_get.return_value = mock.Mock(
            status_code=200,
            text=_H382_SHEET_CSV,
            headers={},
            raise_for_status=mock.Mock(),
        )
        url = (
            "https://docs.google.com/spreadsheets/d/test123/edit"
            "?gid=0#gid=0"
        )
        df = _fetch_sheet_dataframe(url)
        self.assertEqual(list(df.columns), list(REQUEST_COLS_READ))
        hip = df[df["unique_id"] == "HIP94931"].iloc[0]
        self.assertIn("Kepler-444", hip["comments"])


@pytest.mark.parametrize(
    "token,expected",
    [
        ("p1", 1),
        ("P2", 2),
        (" p3 ", 3),
        ("p4", 3),
        ("p5", 3),
        ("", None),
        (None, None),
        (float("nan"), None),
        ("p9", 3),
    ],
)
def test_priority_token_to_splan_weight(token, expected):
    assert priority_token_to_splan_weight(token) == expected


class TestAttachSplanWeight(unittest.TestCase):
    def _frame(self, priorities):
        return pd.DataFrame(
            {
                "unique_id": [f"T{i}" for i in range(len(priorities))],
                "priority": priorities,
            }
        )

    def test_maps_all_priority_tokens(self):
        df = self._frame(["p1", "p2", "p3", "p4", "p5"])
        out = attach_splan_weight(df)
        self.assertListEqual(out["splan_weight"].tolist(), [1, 2, 3, 3, 3])

    def test_blank_priority_uses_least_favored_tier_present(self):
        df = self._frame(["p1", "", "p2", None])
        out = attach_splan_weight(df)
        self.assertListEqual(out["splan_weight"].tolist(), [1, 2, 2, 2])

    def test_blank_priority_when_only_p1_present(self):
        df = self._frame(["p1", "", None])
        out = attach_splan_weight(df)
        self.assertListEqual(out["splan_weight"].tolist(), [1, 1, 1])

    def test_all_blank_priorities_default_to_two(self):
        df = self._frame(["", None, "nan"])
        out = attach_splan_weight(df)
        self.assertListEqual(out["splan_weight"].tolist(), [2, 2, 2])

    def test_empty_dataframe(self):
        out = attach_splan_weight(pd.DataFrame())
        self.assertIn("splan_weight", out.columns)

    def test_mapping_table_matches_spec(self):
        self.assertEqual(PRIORITY_TO_SPLAN_WEIGHT["p1"], 1)
        self.assertEqual(PRIORITY_TO_SPLAN_WEIGHT["p2"], 2)
        self.assertEqual(PRIORITY_TO_SPLAN_WEIGHT["p3"], 3)
        self.assertEqual(PRIORITY_TO_SPLAN_WEIGHT["p4"], 3)
        self.assertEqual(PRIORITY_TO_SPLAN_WEIGHT["p5"], 3)


class TestSheetRequestCols(unittest.TestCase):
    def test_sheet_cols_selectable_before_attach(self):
        row = {col: "" for col in REQUEST_COLS_READ}
        row.update(
            {
                "unique_id": "T1",
                "target": "T1",
                "program_code": "2026A_X001",
                "priority": "p2",
                "comments": "",
            }
        )
        df = pd.DataFrame([row])
        sheet_slice = df[SHEET_REQUEST_COLS]
        self.assertNotIn("splan_weight", sheet_slice.columns)

        out = attach_splan_weight(sheet_slice)
        final = out[REQUEST_COLS]
        self.assertEqual(final["splan_weight"].iloc[0], 2)


class TestPullAllScheduled(unittest.TestCase):
    def test_queries_hiresr_and_kpfcc(self):
        hires = pd.DataFrame(
            {
                "Date": ["2026-02-01"],
                "Time": ["07:00 - 13:00 (100%)"],
                "Instrument": ["HIRESr"],
                "ProjCode": ["N012"],
            }
        )
        kpf = pd.DataFrame(
            {
                "Date": ["2026-02-02"],
                "Time": ["05:00 - 12:00 (100%)"],
                "Instrument": ["KPF-CC"],
                "ProjCode": ["N153"],
            }
        )

        with patch.object(
            hirescps_prep,
            "_query_keck_schedule_form",
            side_effect=[hires, kpf],
        ) as query:
            df = hirescps_prep.pull_all_scheduled("2026-02-01", "2026-07-31")

        self.assertEqual(
            query.call_args_list[0].args[:3],
            ("HIRESr", "2026-02-01", "2026-07-31"),
        )
        self.assertEqual(
            query.call_args_list[1].args[:3],
            ("KPF-CC", "2026-02-01", "2026-07-31"),
        )
        self.assertEqual(len(df), 2)
        self.assertEqual(
            sorted(df["Instrument"].tolist()),
            ["HIRESr", "KPF-CC"],
        )

    def test_skips_instrument_on_query_failure(self):
        hires = pd.DataFrame(
            {
                "Date": ["2026-02-01"],
                "Time": ["07:00 - 13:00 (100%)"],
                "Instrument": ["HIRESr"],
                "ProjCode": ["N012"],
            }
        )

        with patch.object(
            hirescps_prep,
            "_query_keck_schedule_form",
            side_effect=[
                hires,
                RuntimeError("Unexpected response from schedule form (KPF-CC): "),
            ],
        ):
            df = hirescps_prep.pull_all_scheduled("2026-08-01", "2027-01-31")

        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["Instrument"], "HIRESr")

    def test_returns_empty_when_all_instruments_fail(self):
        with patch.object(
            hirescps_prep,
            "_query_keck_schedule_form",
            side_effect=RuntimeError("Unexpected response from schedule form (HIRESr): "),
        ):
            df = hirescps_prep.pull_all_scheduled("2026-08-01", "2027-01-31")

        self.assertTrue(df.empty)


if __name__ == "__main__":
    unittest.main()
