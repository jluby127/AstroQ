import re
import tempfile
import unittest

import pandas as pd

from astroq.queue.hirescps.script_columns import (
    EXPOSURE_WIDTH,
    SECTION_HEADER_WIDTH,
    TARGET_NAME_WIDTH,
    VMAG_WIDTH,
    format_cell_token,
    format_exposure_token,
    format_meter_token,
    format_section_header,
    format_vmag_token,
)
from astroq.queue.hirescps.bstars import format_bstar_row
from astroq.queue.hirescps.starlist import format_hires_row, write_starlist


class TestHirescpsScriptColumns(unittest.TestCase):
    def test_vmag_token_examples(self):
        self.assertEqual(format_vmag_token(4.7), " vmag=4.7")
        self.assertEqual(format_vmag_token(11.7), "vmag=11.7")
        self.assertEqual(len(format_vmag_token(4.7)), VMAG_WIDTH)

    def test_exposure_token_examples(self):
        self.assertEqual(format_exposure_token(5, 500), "    5/500")
        self.assertEqual(format_exposure_token(96, 193), "   96/193")
        self.assertEqual(format_exposure_token(1800, 1800), "1800/1800")
        self.assertEqual(len(format_exposure_token(5, 500)), EXPOSURE_WIDTH)

    def test_exposure_overflow_emits_unpadded(self):
        self.assertEqual(format_exposure_token(12341, 24682), "12341/24682")

    def test_meter_token(self):
        self.assertEqual(format_meter_token("10k"), " 10k")
        self.assertEqual(format_meter_token("250k"), "250k")

    def test_cell_token(self):
        self.assertEqual(format_cell_token("in"), " in")
        self.assertEqual(format_cell_token("out"), "out")

    def _sample_row(self, **overrides):
        base = {
            "target": "TOI2399",
            "ra": 175.473,
            "dec": -21.7967,
            "exptime": 96,
            "maxtime": 193,
            "n_exp": 1,
            "n_intra_max": 1,
            "program_code": "2026A_H399",
            "priority": "p2",
            "decker": "C2",
            "cell in/out?": "out",
            "exp_meter_threshold": "10k",
            "Vmag": 11.7,
            "pmra": 0.0,
            "pmdec": 0.0,
        }
        base.update(overrides)
        return pd.DataFrame([base])

    def test_column_right_edges_align(self):
        rows = [
            self._sample_row(Vmag=11.7, exptime=96, maxtime=193, exp_meter_threshold="10k"),
            self._sample_row(
                target="TIC142381532",
                Vmag=10.7,
                exptime=661,
                maxtime=1058,
                exp_meter_threshold="125k",
                decker="B5",
                **{"cell in/out?": "in"},
            ),
            self._sample_row(Vmag=4.7, exptime=5, maxtime=500, exp_meter_threshold="250k", decker="B3", n_exp=3, priority="p1", **{"cell in/out?": "in"}),
            self._sample_row(Vmag=13.3, exptime=1200, maxtime=1200, exp_meter_threshold="30k"),
        ]
        edges = []
        for row in rows:
            line = format_hires_row(
                row,
                "05:57",
                "05:57",
                "07:09",
                "2026-06-12",
            )
            vmag_t = re.search(r"vmag=\S+", line).group(0)
            exp_t = re.search(r"\d+/\d+", line).group(0)
            meter_t = re.search(r"\d+/\d+\s+(\S+k)", line).group(1)
            edges.append(
                (
                    line.index(vmag_t) + len(vmag_t),
                    line.index(exp_t) + len(exp_t),
                    line.index(meter_t) + len(meter_t),
                )
            )
        self.assertEqual(len(set(edges)), 1)

    def test_bstar_column_right_edges_match_requests(self):
        import astropy.units as u
        from astropy.coordinates import SkyCoord

        coord = SkyCoord(ra=0 * u.deg, dec=0 * u.deg)
        bstar_line = format_bstar_row("hr9098", coord, 4.5, 0.0, 0.0, "2026-07-18")
        request_line = format_hires_row(
            self._sample_row(Vmag=4.5, exptime=5, maxtime=500, exp_meter_threshold="250k", decker="B5"),
            "10:29",
            "10:29",
            "11:51",
            "2026-07-18",
        )

        def edges(line):
            vmag_t = re.search(r"vmag=\S+", line).group(0)
            exp_t = re.search(r"\d+/\d+", line).group(0)
            meter_t = re.search(r"\d+/\d+\s+(\S+k)", line).group(1)
            return (
                line.index(vmag_t) + len(vmag_t),
                line.index(exp_t) + len(exp_t),
                line.index(meter_t) + len(meter_t),
            )

        self.assertEqual(edges(bstar_line), edges(request_line))

    def test_comments_appended_after_epoch_and_obs(self):
        row = self._sample_row(
            pmra=10.0,
            pmdec=-5.0,
            comments="do N comp 7",
        )
        line = format_hires_row(
            row,
            None,
            None,
            None,
            "2026-07-18",
            omit_timing=True,
            obs_token="obs=1/3",
        )
        self.assertIn("epoch=2026.5 obs=1/3 do N comp 7", line)

    def test_comments_after_epoch_when_no_obs(self):
        row = self._sample_row(
            pmra=10.0,
            pmdec=-5.0,
            comments="do N comp 7",
        )
        line = format_hires_row(row, "10:29", "10:29", "11:51", "2026-07-18")
        self.assertIn("epoch=2026.5 do N comp 7", line)
        self.assertNotIn("obs=", line)

    def test_comments_fallback_observing_notes(self):
        row = self._sample_row(
            pmra=10.0,
            pmdec=-5.0,
            **{"Observing Notes": "check focus"},
        )
        line = format_hires_row(
            row,
            None,
            None,
            None,
            "2026-07-18",
            omit_timing=True,
            obs_token="obs=0/1",
        )
        self.assertIn("epoch=2026.5 obs=0/1 check focus", line)

    def test_section_header_format(self):
        for label in ("EXTRAS", "2026A-Requests-All"):
            header = format_section_header(label)
            self.assertEqual(len(header), SECTION_HEADER_WIDTH)
            self.assertIn(f"__{label}__", header)
            self.assertNotIn(" ", header)

    def test_bstar_name_padding(self):
        import astropy.units as u
        from astropy.coordinates import SkyCoord

        coord = SkyCoord(ra=0 * u.deg, dec=0 * u.deg)
        line = format_bstar_row("hr9098", coord, 4.5, 0.0, 0.0, "2026-07-18")
        self.assertTrue(
            line.startswith(" " * (TARGET_NAME_WIDTH - len("HR9098")) + "HR9098")
        )

    def test_bstar_row_hr8976_golden(self):
        import astropy.units as u
        from astropy.coordinates import SkyCoord

        coord = SkyCoord("23 40 24.7 +44 20 02", unit=(u.hourangle, u.deg))
        line = format_bstar_row("hr8976", coord, 4.1, 0.0, 0.0, "2026-07-01")
        expected = (
            "          HR8976 23 40 24.7 +44 20 02 2000  vmag=4.1     5/500 "
            "250k B5 1x  in p3 XX epoch=2026.5"
        )
        self.assertEqual(line, expected)


class TestScriptSectionLayout(unittest.TestCase):
    def test_section_headers_no_blank_lines(self):
        from unittest import mock
        from astropy.time import Time

        active = pd.DataFrame(
            [
                {
                    "unique_id": "u1",
                    "target": "Star1",
                    "ra": 180.0,
                    "dec": 20.0,
                    "exptime": 60,
                    "maxtime": 120,
                    "n_exp": 1,
                    "program_code": "2026A_H001",
                    "priority": "p2",
                    "decker": "C2",
                    "cell in/out?": "out",
                    "exp_meter_threshold": "10k",
                    "Vmag": 6.0,
                    "pmra": 0.0,
                    "pmdec": 0.0,
                    "past_nights_observed": 0,
                    "n_inter_max": 1,
                }
            ]
        )
        empty_schedule = pd.DataFrame(
            [
                {
                    "is_anchor": True,
                    "scheduled": False,
                    "order": 0,
                    "unique_id": "anchor",
                    "t_start": 0.0,
                    "t_earliest_start": 0.0,
                    "t_latest_finish": 0.0,
                }
            ]
        )

        with tempfile.TemporaryDirectory() as outdir, mock.patch(
            "astroq.queue.hirescps.starlist.build_bstars_section",
            return_value=[format_section_header("B-Stars-Coordinates-Advanced"), "HR9098 stub"],
        ):
            lines = write_starlist(
                active.copy(),
                empty_schedule,
                Time("2026-06-22T05:00:00"),
                [],
                "2026-06-22",
                outdir,
                all_active_requests=active.copy(),
            )

        extras_header = format_section_header("EXTRAS")
        self.assertIn(extras_header, lines)
        extras_idx = lines.index(extras_header)
        self.assertNotEqual(lines[extras_idx + 1], "")

        backup_header = format_section_header("2026A-Requests-All")
        backup_idx = lines.index(backup_header)
        self.assertNotEqual(lines[backup_idx - 1], "")
        self.assertNotEqual(lines[backup_idx + 1], "")

        for i in range(len(lines) - 1):
            if lines[i] == "" and lines[i + 1] == "":
                self.fail(f"consecutive blank lines at index {i}")


class TestAccessibleAt(unittest.TestCase):
    """Access.accessible_at reuses the compute_altaz pointing gate."""

    def _build_access(self, request_frame):
        from astroq.access import Access
        from astroq.queue.hirescps.queue import HIRESCPS

        return Access(
            queue=HIRESCPS(),
            request_frame=request_frame,
            semester_start_date="2026-06-22",
            semester_length=1,
            slot_size=10.0,
        )

    def test_southern_target_never_accessible(self):
        import numpy as np
        import astropy.units as u
        from astropy.time import Time

        # dec = -80 deg is always below the horizon at Keck (lat ~ +19.8).
        rf = pd.DataFrame(
            {"unique_id": ["up", "down"], "ra": [180.0, 180.0], "dec": [20.0, -80.0]}
        )
        access = self._build_access(rf)
        day = Time("2026-06-22")
        times = day + np.linspace(0, 1, 49) * u.day  # 30-min cadence over 24h

        mask = access.accessible_at(times)
        self.assertEqual(mask.shape, (2, len(times)))
        # The dec=+20 target is accessible at some point in the day...
        self.assertTrue(mask[0].any())
        # ...but the dec=-80 target never is.
        self.assertFalse(mask[1].any())


class TestTwilightSections(unittest.TestCase):
    """Evening/Morning Twilight backup blocks filter to V<8 and the uid sets."""

    def _request_row(self, uid, target, vmag):
        return {
            "unique_id": uid,
            "target": target,
            "ra": 180.0,
            "dec": 20.0,
            "exptime": 60,
            "maxtime": 120,
            "n_exp": 1,
            "program_code": "2026A_H001",
            "priority": "p2",
            "decker": "C2",
            "cell in/out?": "out",
            "exp_meter_threshold": "10k",
            "Vmag": vmag,
            "pmra": 0.0,
            "pmdec": 0.0,
            "past_nights_observed": 0,
            "n_inter_max": 5,
        }

    def test_twilight_blocks_emitted_and_filtered(self):
        from unittest import mock
        from astropy.time import Time

        # bright + visible, bright + not-visible, faint + visible
        active = pd.DataFrame(
            [
                self._request_row("bright_vis", "BrightVis", 6.0),
                self._request_row("bright_novis", "BrightNoVis", 5.0),
                self._request_row("faint_vis", "FaintVis", 11.0),
            ]
        )
        # A single anchor row (filtered out by ~is_anchor) leaves an empty
        # on-sky schedule with correct boolean dtypes.
        empty_schedule = pd.DataFrame(
            [
                {
                    "is_anchor": True,
                    "scheduled": False,
                    "order": 0,
                    "unique_id": "anchor",
                    "t_start": 0.0,
                    "t_earliest_start": 0.0,
                    "t_latest_finish": 0.0,
                }
            ]
        )

        with tempfile.TemporaryDirectory() as outdir, mock.patch(
            "astroq.queue.hirescps.starlist.build_bstars_section", return_value=[]
        ):
            lines = write_starlist(
                active.copy(),
                empty_schedule,
                Time("2026-06-22T05:00:00"),
                [],
                "2026-06-22",
                outdir,
                all_active_requests=active.copy(),
                evening_twilight_uids={"bright_vis", "faint_vis"},
                morning_twilight_uids=set(),
            )

        text = "\n".join(lines)
        self.assertIn("__2026A-Evening-Twilight__", text)
        self.assertIn("__2026A-Morning-Twilight__", text)

        evening_block = text.split("__2026A-Evening-Twilight__", 1)[1].split(
            "__2026A-Morning-Twilight__", 1
        )[0]
        # Bright + in evening set -> present; faint filtered out; bright not in set absent.
        self.assertIn("BrightVis", evening_block)
        self.assertNotIn("FaintVis", evening_block)
        self.assertNotIn("BrightNoVis", evening_block)


if __name__ == "__main__":
    unittest.main()
