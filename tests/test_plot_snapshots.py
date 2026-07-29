"""Characterization snapshots for the plot package (figures + tables)."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import plotly.io as pio

import astroq.plot as pl
import astroq.nplan as nplan
from astroq.splan import SemesterPlanner

HELLO_H5 = (
    "examples/hello_world/2018B/2018-08-05/band1/outputs/semester_planner.h5"
)
NIGHT_H5 = (
    "examples/hello_world/2018B/2018-08-05/band1/outputs/night_planner.h5"
)
BASELINE_DIR = Path(__file__).resolve().parent / "baseline" / "plot"
REGEN = os.environ.get("ASTROQ_REGEN_PLOT_BASELINES", "") == "1"


def _digest(obj) -> str:
    if hasattr(obj, "to_json"):
        payload = obj.to_json()
    elif isinstance(obj, str):
        payload = obj
    else:
        payload = json.dumps(obj, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()


def _write_or_compare(name: str, digest: str):
    path = BASELINE_DIR / f"{name}.sha256"
    if REGEN or not path.exists():
        BASELINE_DIR.mkdir(parents=True, exist_ok=True)
        path.write_text(digest + "\n")
        return
    expected = path.read_text().strip()
    assert digest == expected, f"{name}: expected {expected}, got {digest}"


class TestPlotSnapshots(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="astroq_plot_snap_")
        cls.sp = SemesterPlanner.from_hdf5(HELLO_H5)
        cls.pd = pl.build_plot_data(cls.sp)
        cls.sel_all = cls.pd.select_all()
        cls.sel_prog = cls.pd.select_all(aggregate_by_program=True)

    def test_semester_figures(self):
        all_programs = sorted(self.pd.program_table.index)
        cases = [
            ("birdseye", pl.get_birdseye, self.sel_prog, {}),
            ("football", pl.get_football, self.sel_all, {"use_program_colors": True}),
            ("tau_inter", pl.get_tau_inter_line, self.sel_all, {"use_program_colors": True}),
            ("timebar", pl.get_timebar, self.sel_all, {"use_program_colors": True}),
            ("timebar_by_program", pl.get_timebar_by_program, self.sel_all, {}),
            ("rawobs", pl.get_rawobs, self.sel_all, {"use_program_colors": True}),
            (
                "completion_hist",
                pl.get_completion_histogram_by_weight,
                self.sel_all,
                {},
            ),
        ]
        with patch.object(pl, "_football_cache_dir", lambda sp: Path(self.tmp)):
            _write_or_compare(
                "cof_visits",
                _digest(pl.get_cof(self.pd, programs=all_programs)),
            )
            _write_or_compare(
                "cof_time",
                _digest(pl.get_cof(self.pd, programs=all_programs, units="time")),
            )
            for name, fn, sel, kwargs in cases:
                fig = fn(self.pd, sel, **kwargs)
                _write_or_compare(name, _digest(fig))

    def test_request_table_html(self):
        df = pl.get_request_frame(self.pd, self.sel_all)
        html = pl.request_frame_to_html(df)
        _write_or_compare("request_table", _digest(html))

    def test_night_figures(self):
        if not os.path.exists(NIGHT_H5):
            self.skipTest("night planner fixture missing")
        night_planner = nplan.NightPlanner.from_hdf5(NIGHT_H5)
        data_ttp = night_planner.solution
        night_start, _ = nplan.get_nightly_times_from_allocation(
            night_planner.allocation_file, night_planner.current_day
        )
        ladder = pl.get_ladder(data_ttp, night_start)
        _write_or_compare("ladder", _digest(ladder))
        script_df = pl.get_script_plan(night_planner)
        html = pl.nightplan_table_to_html(script_df)
        _write_or_compare("script_table", _digest(html))


if __name__ == "__main__":
    unittest.main()
