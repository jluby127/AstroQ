"""Tests for astroq.access: local-noon observing-day grid and night windows."""

import unittest
from datetime import timezone

import numpy as np
import pandas as pd
from astropy.time import Time

import astroq.access as ac
from astroq.queue.hirescps.queue import HIRESCPS


class TestObservingDayGrid(unittest.TestCase):
    def setUp(self):
        self.queue = HIRESCPS()

    def test_allocation_spanning_utc_midnight_single_night(self):
        """UTC-spanning allocation maps to one civil noon-start night d."""
        # Keck: civil 2026-06-01 noon HST = 2026-06-01 22:00 UTC through
        # 2026-06-02 22:00 UTC. Block below sits on morning of 2026-06-02 UTC.
        alloc = pd.DataFrame(
            {
                "start": [Time("2026-06-02T02:00:00", scale="utc")],
                "stop": [Time("2026-06-02T10:00:00", scale="utc")],
            }
        )
        rf = pd.DataFrame({"unique_id": ["t1"], "ra": [120.0], "dec": [20.0]})
        access = ac.Access(
            self.queue,
            rf,
            "2026-06-01",
            3,
            10,
            allocation=alloc,
        )
        allocated = access.compute_allocated()[0]
        nights_with_alloc = np.where(allocated.any(axis=1))[0]
        self.assertEqual(len(nights_with_alloc), 1)
        d = int(nights_with_alloc[0])
        self.assertEqual(access.all_dates_array[d], "2026-06-01")

    def test_local_midnight_near_center_slot(self):
        """Local midnight falls at row n_slots//2 on the noon-anchored grid."""
        access = ac.Access(
            self.queue,
            pd.DataFrame({"unique_id": ["t1"], "ra": [0.0], "dec": [0.0]}),
            "2026-06-01",
            2,
            10,
        )
        center = access.nslots // 2
        dt = access.slotmidpoints[0, center].utc.to_datetime()
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        local = dt.astimezone(ac._observatory_tz(self.queue.observatory))
        self.assertEqual(local.hour, 0)
        self.assertLess(local.minute, access.slot_size)

    def test_past_utc_morning_maps_to_previous_civil_night(self):
        """UTC morning after local midnight still belongs to prior noon-start label."""
        # 2026-06-02 08:00 UTC = 2026-06-01 22:00 HST → civil night 2026-06-01
        label = ac.civil_night_label(
            Time("2026-06-02T08:00:00", scale="utc"),
            self.queue.observatory,
        )
        self.assertEqual(label, "2026-06-01")

    def test_past_utc_evening_maps_to_same_civil_night(self):
        # 2026-06-02 06:00 UTC = 2026-06-01 20:00 HST → still night 2026-06-01
        label = ac.civil_night_label(
            Time("2026-06-02T06:00:00", scale="utc"),
            self.queue.observatory,
        )
        self.assertEqual(label, "2026-06-01")

    def test_night_window_contains_evening(self):
        access = ac.Access(
            self.queue,
            pd.DataFrame({"unique_id": ["t1"], "ra": [0.0], "dec": [0.0]}),
            "2026-06-01",
            2,
            10,
        )
        start, end = access.night_window("2026-06-01")
        evening = Time("2026-06-02T06:00:00", scale="utc")
        self.assertLessEqual(start, evening)
        self.assertLess(evening, end)


if __name__ == "__main__":
    unittest.main()
