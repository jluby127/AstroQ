"""Tests for :mod:`astroq.queue.base`."""

import unittest

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

from astroq.queue.hirescps.queue import HIRESCPS


class TestVisitTiming(unittest.TestCase):
    def setUp(self):
        self.q = HIRESCPS()

    def test_visit_seconds_one_slew_scalar(self):
        # 300s x 2 shots + 30s readout + 30s slew = 660s (HIRESCPS constants)
        self.assertEqual(self.q.visit_seconds(300, 2), 660.0)

    def test_visit_seconds_one_slew_vector(self):
        exptime = np.array([300.0, 600.0])
        n_exp = np.array([1, 3])
        got = self.q.visit_seconds(exptime, n_exp)
        want = exptime * n_exp + self.q.readout_time * (n_exp - 1) + self.q.slew_overhead_mean
        np.testing.assert_allclose(got, want)

    def test_visit_duration_vs_visit_seconds(self):
        exptime, n_exp = 300, 2
        duration_min = self.q.visit_duration(exptime, n_exp)
        duration_s = duration_min * 60.0
        semester_s = self.q.visit_seconds(exptime, n_exp)
        self.assertAlmostEqual(semester_s - duration_s, self.q.slew_overhead_mean)


class TestIsAccessible(unittest.TestCase):
    def setUp(self):
        self.q = HIRESCPS()

    def test_below_elevation_clamp(self):
        self.assertFalse(self.q.is_accessible(27.9, 180.0))
        self.assertTrue(self.q.is_accessible(28.1, 180.0))

    def test_above_elevation_clamp(self):
        self.assertTrue(self.q.is_accessible(84.9, 180.0))
        self.assertFalse(self.q.is_accessible(85.0, 180.0))

    def test_nasmyth_deck_interior(self):
        self.assertFalse(self.q.is_accessible(20.0, 90.0))

    def test_nasmyth_deck_exterior(self):
        self.assertTrue(self.q.is_accessible(40.0, 200.0))


class TestSlewFn(unittest.TestCase):
    def setUp(self):
        self.q = HIRESCPS()

    def test_slew_fn_output_shape(self):
        coord_a = SkyCoord(ra=180.0 * u.deg, dec=30.0 * u.deg, frame="icrs")
        coord_b = SkyCoord(ra=200.0 * u.deg, dec=25.0 * u.deg, frame="icrs")
        window_start = Time(["2026-05-09T06:00:00", "2026-05-09T07:00:00"])
        window_end = Time(["2026-05-09T06:30:00", "2026-05-09T07:30:00"])
        tau = self.q.slew_fn(coord_a, coord_b, window_start, window_end)
        self.assertEqual(tau.shape, (1, 2))
        self.assertTrue(np.all(np.isfinite(tau)))
        self.assertTrue(np.all(tau >= 0))


if __name__ == "__main__":
    unittest.main()
