"""Tests for the state-aware (cable-wrap) TTP slew model.

These are self-contained: they build small synthetic request sets and solve
the TTP directly, so they do not depend on the ``examples/`` sandbox.
"""

import unittest

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.time import Time
import astropy.units as u

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel


NIGHT_START = Time("2026-05-09T06:00:00", format="isot")
NIGHT_END = Time("2026-05-09T09:30:00", format="isot")


def _requests(ras, decs, queue, *, t_visit_min=5.0):
    """Build a QTable assuming full-window availability for each target."""
    n = len(ras)
    coords = SkyCoord(np.asarray(ras) * u.deg, np.asarray(decs) * u.deg, frame="icrs")
    return QTable(
        {
            "unique_id": np.array([f"T{i}" for i in range(n)], dtype=object),
            "coord": coords,
            "first_available": Time([NIGHT_START.isot] * n),
            "last_available": Time([NIGHT_END.isot] * n),
            "t_visit": np.full(n, t_visit_min) * u.min,
            "n_intra_max": np.ones(n, dtype=int),
            "tau_intra": np.zeros(n) * u.hr,
            "priority": np.full(n, 10.0),
        },
        copy=False,
    )


def _solve(queue, requests, n_states, *, runtime=60):
    slew_fn = queue.slew_fn_state if n_states > 1 else queue.slew_fn
    tm = TTPModel(
        requests=requests,
        night_start=NIGHT_START,
        night_end=NIGHT_END,
        slew_fn=slew_fn,
        n_slots=queue.nSlots,
        n_states=n_states,
    )
    tm.build_nodes()
    tm.build_arcs()
    tm.build_model()
    tm.model.params.OutputFlag = 0
    tm.model.params.TimeLimit = runtime
    tm.model.params.MIPGap = 0.01
    tm.model.update()
    tm.run_model()
    tm.build_schedule()
    return tm


class TestEncoderGeometry(unittest.TestCase):
    def setUp(self):
        self.q = HIRESCPS()

    def test_wrap_states_defined(self):
        self.assertEqual(self.q.n_states, 2)
        names = [s[0] for s in self.q.wrap_states]
        self.assertEqual(names, ["N", "S"])

    def test_encoder_az_overlap_two_states(self):
        # Western overlap region (sky az 215-315) reachable in BOTH windings.
        for az in (215.0, 250.0, 280.0, 315.0):
            n = self.q._encoder_az(az, -145.0, 90.0)
            s = self.q._encoder_az(az, 90.0, 315.0)
            self.assertTrue(np.isfinite(n), f"N-wrap should reach az={az}")
            self.assertTrue(np.isfinite(s), f"S-wrap should reach az={az}")

    def test_encoder_az_neutral_single_state(self):
        # az ~150 (south) only reachable in S-wrap; az ~30 (NE) only in N-wrap.
        self.assertTrue(np.isnan(self.q._encoder_az(150.0, -145.0, 90.0)))
        self.assertTrue(np.isfinite(self.q._encoder_az(150.0, 90.0, 315.0)))
        self.assertTrue(np.isfinite(self.q._encoder_az(30.0, -145.0, 90.0)))
        self.assertTrue(np.isnan(self.q._encoder_az(30.0, 90.0, 315.0)))

    def test_below_horizon_is_inaccessible(self):
        # Regression: the elevation clamp must exclude below-horizon points
        # (alt < 0), not just the 0-18 deg band, or availability windows let
        # the TTP schedule targets after they have set.
        alt = np.array([-90.0, -30.0, -0.1, 5.0, 17.9, 18.1, 45.0, 84.9, 85.1])
        az = np.full_like(alt, 200.0)  # away from the nasmyth deck az range
        ok = self.q.is_accessible(alt, az)
        self.assertFalse(ok[alt < 18.0].any(), "alt < 18 must be inaccessible")
        self.assertTrue(ok[(alt >= 18.0) & (alt <= 85.0)].all())

    def test_encoder_az_values(self):
        # N-wrap maps az>=215 to az-360 (continuous through north).
        self.assertAlmostEqual(self.q._encoder_az(235.0, -145.0, 90.0), -125.0)
        self.assertAlmostEqual(self.q._encoder_az(45.0, -145.0, 90.0), 45.0)
        self.assertAlmostEqual(self.q._encoder_az(270.0, 90.0, 315.0), 270.0)


class TestSlewFnState(unittest.TestCase):
    def setUp(self):
        self.q = HIRESCPS()

    def test_shape_and_nan_propagation(self):
        # One western (overlap) target, one southern (S-only) target.
        coords = SkyCoord([250.0, 165.0] * u.deg, [20.0, -10.0] * u.deg, frame="icrs")
        a = coords[[0, 0, 1, 1]]
        b = coords[[0, 1, 0, 1]]
        ws = Time([NIGHT_START.isot])
        we = Time([NIGHT_END.isot])
        tau = self.q.slew_fn_state(a, b, ws, we)
        self.assertEqual(tau.shape, (4, 1, 2, 2))
        self.assertTrue(np.all(tau[np.isfinite(tau)] >= 0.0))

    def test_self_pair_zero_diagonal(self):
        # Short window so the target stays in-wrap; same target+state -> 0 slew.
        coords = SkyCoord([250.0] * u.deg, [20.0] * u.deg, frame="icrs")
        ws = Time([NIGHT_START.isot])
        we = ws + 1.0 * u.min
        tau = self.q.slew_fn_state(coords, coords, ws, we)
        diag = np.array([tau[0, 0, 0, 0], tau[0, 0, 1, 1]])
        # at least one winding must be reachable, and any reachable self-pair
        # has zero slew.
        self.assertTrue(np.isfinite(diag).any())
        for v in diag[np.isfinite(diag)]:
            self.assertAlmostEqual(v, 0.0, places=6)


class TestTwoStateSolve(unittest.TestCase):
    def setUp(self):
        self.q = HIRESCPS()
        # Mix: western overlap, northern, southern targets.
        self.ras = [184.0, 186.0, 170.0, 176.0, 192.0, 198.0]
        self.decs = [13.0, 15.0, 5.0, 27.0, 51.0, -3.0]
        self.req = _requests(self.ras, self.decs, self.q)

    def test_two_state_schedules_and_assigns_wrap(self):
        tm = _solve(self.q, self.req, 2)
        self.assertIsNotNone(tm.schedule)
        sched = tm.schedule[tm.schedule["scheduled"]]
        self.assertGreater(len(sched), 0)
        self.assertIn("wrap_state", tm.schedule.columns)
        # Every scheduled node must have a defined wrap state in {0, 1}.
        self.assertTrue(sched["wrap_state"].isin([0, 1]).all())

    def test_no_selected_arc_uses_unreachable_winding(self):
        """Each scheduled node must be physically reachable in its wrap state."""
        tm = _solve(self.q, self.req, 2)
        sched = tm.schedule[tm.schedule["scheduled"]].copy()
        from astropy.time import TimeDelta

        t_obs = NIGHT_START + TimeDelta(sched["t_start"].to_numpy() * 60.0, format="sec")
        coords = SkyCoord(
            sched["ra"].to_numpy() * u.deg, sched["dec"].to_numpy() * u.deg, frame="icrs"
        )
        az = np.atleast_1d(self.q.observatory.altaz(t_obs, coords).az.deg)
        for k, state in enumerate(sched["wrap_state"].to_numpy().astype(int)):
            _, lo, hi = self.q.wrap_states[state]
            enc = self.q._encoder_az(az[k], lo, hi)
            self.assertTrue(
                np.isfinite(enc),
                f"node observed at az={az[k]:.1f} not reachable in state {state}",
            )

    def test_single_state_regression_still_solves(self):
        tm = _solve(self.q, self.req, 1)
        self.assertIsNotNone(tm.schedule)
        # Single-state is the degenerate S == 1 case: wrap_state is present and
        # always 0 (one unified, state-indexed code path).
        self.assertIn("wrap_state", tm.schedule.columns)
        sched = tm.schedule[tm.schedule["scheduled"]]
        self.assertTrue((sched["wrap_state"] == 0).all())
        self.assertEqual(tm.stats["n_requested"], len(self.ras))
        self.assertGreaterEqual(tm.stats["t_slew_sum"], 0.0)


if __name__ == "__main__":
    unittest.main()
