"""Tests for astroq.ttp: ACS warm-start and two-state cable-wrap model.

Self-contained synthetic request sets; ACS returns a tour dict
(:func:`astroq.ttp.acs.acs_warm_start`) that seeds the MILP for Gurobi.
"""

import unittest

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.time import Time
import astropy.units as u

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel
from astroq.ttp.acs import acs_warm_start


NIGHT_START = Time("2026-05-09T06:00:00", format="isot")
NIGHT_END = Time("2026-05-09T09:30:00", format="isot")


def _requests(ras, decs, *, t_visit_min=5.0):
    n = len(ras)
    coords = SkyCoord(np.asarray(ras) * u.deg, np.asarray(decs) * u.deg, frame="icrs")
    return QTable(
        {
            "unique_id": np.array([f"T{i}" for i in range(n)], dtype=object),
            "coord": coords,
            "time_earliest_start": Time([NIGHT_START.isot] * n),
            "time_latest_finish": Time([NIGHT_END.isot] * n),
            "t_visit": np.full(n, t_visit_min) * u.min,
            "n_intra_max": np.ones(n, dtype=int),
            "tau_intra": np.zeros(n) * u.hr,
            "weight": np.full(n, 1.0),
        },
        copy=False,
    )


def _build(requests, n_states):
    q = HIRESCPS()
    slew_fn = q.slew_fn_state if n_states > 1 else q.slew_fn
    tm = TTPModel(
        requests=requests,
        night_start=NIGHT_START,
        night_end=NIGHT_END,
        slew_fn=slew_fn,
        n_slots=q.nSlots,
        n_states=n_states,
    )
    tm.build_nodes()
    tm.build_arcs()
    return tm


def _uid_completions(tm, res):
    """Map ``{unique_id: [completion_minutes, ...]}`` from an ACS result."""
    out = {}
    for nid, ti in zip(res["order"], res["ti"]):
        uid = tm.nodes.at[nid, "unique_id"]
        out.setdefault(uid, []).append(float(ti))
    return out


RAS = [184.0, 186.0, 170.0, 176.0, 192.0, 198.0]
DECS = [13.0, 15.0, 5.0, 27.0, 51.0, -3.0]


class TestACSHeuristic(unittest.TestCase):
    def setUp(self):
        self.req = _requests(RAS, DECS)

    def test_single_state_feasible_and_windows(self):
        tm = _build(self.req, 1)
        res = acs_warm_start(tm, params={"time_limit_s": 1.0})
        self.assertTrue(res["feasible"])
        self.assertGreater(len(res["order"]), 0)
        # completion time within each node's window and the night budget.
        for nid, ti in zip(res["order"], res["ti"]):
            self.assertLessEqual(ti, float(tm.nodes.at[nid, "t_latest_finish"]) + 1e-6)
            self.assertLessEqual(ti, tm.dur_min + 1e-6)

    def test_matches_milp_optimum_small(self):
        """On a tiny instance the ACS should reach the MILP optimum."""
        tm_milp = _build(self.req, 1)
        tm_milp.build_model()
        tm_milp.model.params.OutputFlag = 0
        tm_milp.model.params.TimeLimit = 60
        tm_milp.model.params.MIPGap = 1e-4
        tm_milp.model.update()
        tm_milp.run_model()
        opt = float(tm_milp.model.ObjVal)

        tm = _build(self.req, 1)
        res = acs_warm_start(tm, params={"time_limit_s": 2.0})
        self.assertLessEqual(res["objective"], opt + 1e-6)
        self.assertGreaterEqual(res["objective"], opt - 1e-6)

    def test_seed_accepted_by_gurobi(self):
        """The MIPStart from the ACS tour must be loadable (feasible)."""
        tm = _build(self.req, 1)
        tm.build_model()
        tm.seed_from_tour(acs_warm_start(tm, params={"time_limit_s": 1.0}))
        tm.model.params.OutputFlag = 0
        tm.model.params.TimeLimit = 60
        tm.model.params.MIPGap = 1e-4
        tm.model.update()
        tm.run_model()
        tm.build_schedule()
        self.assertGreater(tm.model.SolCount, 0)
        self.assertGreater(len(tm.schedule[tm.schedule["scheduled"]]), 0)

    def test_two_state_assigns_valid_wrap(self):
        tm = _build(self.req, 2)
        res = acs_warm_start(tm, params={"time_limit_s": 1.0})
        self.assertTrue(res["feasible"])
        self.assertGreater(len(res["order"]), 0)
        # every kept node carries a valid wrap state.
        self.assertTrue(all(s in (0, 1) for s in res["states"].values()))

    def test_two_state_seed_accepted_by_gurobi(self):
        """A two-state ACS tour must seed the state-indexed MILP cleanly."""
        tm = _build(self.req, 2)
        tm.build_model()
        tm.seed_from_tour(acs_warm_start(tm, params={"time_limit_s": 1.0}))
        tm.model.params.OutputFlag = 0
        tm.model.params.TimeLimit = 60
        tm.model.params.MIPGap = 1e-4
        tm.model.update()
        tm.run_model()
        tm.build_schedule()
        self.assertGreater(tm.model.SolCount, 0)
        sched = tm.schedule[tm.schedule["scheduled"]]
        self.assertGreater(len(sched), 0)
        self.assertTrue(sched["wrap_state"].isin([0, 1]).all())

    def test_tau_intra_separation_respected(self):
        """Two visits of one target must be spaced by tau_intra."""
        req = QTable(
            {
                "unique_id": np.array(["A", "B"], dtype=object),
                "coord": SkyCoord([184.0, 186.0] * u.deg, [13.0, 15.0] * u.deg,
                                  frame="icrs"),
                "time_earliest_start": Time([NIGHT_START.isot] * 2),
                "time_latest_finish": Time([NIGHT_END.isot] * 2),
                "t_visit": np.full(2, 5.0) * u.min,
                "n_intra_max": np.array([2, 1], dtype=int),
                "tau_intra": np.array([0.5, 0.0]) * u.hr,  # 30 min for A
                "weight": np.full(2, 1.0),
            },
            copy=False,
        )
        tm = _build(req, 1)
        res = acs_warm_start(tm, params={"time_limit_s": 1.0})
        self.assertTrue(res["feasible"])
        a_comps = sorted(_uid_completions(tm, res).get("A", []))
        if len(a_comps) == 2:
            self.assertGreaterEqual(a_comps[1] - a_comps[0], 30.0 - 1e-6)


if __name__ == "__main__":
    unittest.main()


# --- two-state TTP model (from test_ttp_twostate.py) ---


def _twostate_requests(ras, decs, queue, *, t_visit_min=5.0):
    n = len(ras)
    coords = SkyCoord(np.asarray(ras) * u.deg, np.asarray(decs) * u.deg, frame="icrs")
    return QTable(
        {
            "unique_id": np.array([f"T{i}" for i in range(n)], dtype=object),
            "coord": coords,
            "time_earliest_start": Time([NIGHT_START.isot] * n),
            "time_latest_finish": Time([NIGHT_END.isot] * n),
            "t_visit": np.full(n, t_visit_min) * u.min,
            "n_intra_max": np.ones(n, dtype=int),
            "tau_intra": np.zeros(n) * u.hr,
            "weight": np.full(n, 1.0),
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
        for az in (215.0, 250.0, 280.0, 315.0):
            n = self.q._encoder_az(az, -145.0, 90.0)
            s = self.q._encoder_az(az, 90.0, 315.0)
            self.assertTrue(np.isfinite(n), f"N-wrap should reach az={az}")
            self.assertTrue(np.isfinite(s), f"S-wrap should reach az={az}")

    def test_encoder_az_neutral_single_state(self):
        self.assertTrue(np.isnan(self.q._encoder_az(150.0, -145.0, 90.0)))
        self.assertTrue(np.isfinite(self.q._encoder_az(150.0, 90.0, 315.0)))
        self.assertTrue(np.isfinite(self.q._encoder_az(30.0, -145.0, 90.0)))
        self.assertTrue(np.isnan(self.q._encoder_az(30.0, 90.0, 315.0)))

    def test_below_horizon_is_inaccessible(self):
        alt = np.array([-90.0, -30.0, -0.1, 5.0, 27.9, 28.1, 45.0, 84.9, 85.1])
        az = np.full_like(alt, 200.0)
        ok = self.q.is_accessible(alt, az)
        self.assertFalse(ok[alt < 28.0].any(), "alt < 28 must be inaccessible")
        self.assertTrue(ok[(alt >= 28.0) & (alt <= 85.0)].all())

    def test_encoder_az_values(self):
        self.assertAlmostEqual(self.q._encoder_az(235.0, -145.0, 90.0), -125.0)
        self.assertAlmostEqual(self.q._encoder_az(45.0, -145.0, 90.0), 45.0)
        self.assertAlmostEqual(self.q._encoder_az(270.0, 90.0, 315.0), 270.0)


class TestSlewFnState(unittest.TestCase):
    def setUp(self):
        self.q = HIRESCPS()

    def test_shape_and_nan_propagation(self):
        coords = SkyCoord([250.0, 165.0] * u.deg, [20.0, -10.0] * u.deg, frame="icrs")
        a = coords[[0, 0, 1, 1]]
        b = coords[[0, 1, 0, 1]]
        ws = Time([NIGHT_START.isot])
        we = Time([NIGHT_END.isot])
        tau = self.q.slew_fn_state(a, b, ws, we)
        self.assertEqual(tau.shape, (4, 1, 2, 2))
        self.assertTrue(np.all(tau[np.isfinite(tau)] >= 0.0))

    def test_self_pair_zero_diagonal(self):
        coords = SkyCoord([250.0] * u.deg, [20.0] * u.deg, frame="icrs")
        ws = Time([NIGHT_START.isot])
        we = ws + 1.0 * u.min
        tau = self.q.slew_fn_state(coords, coords, ws, we)
        diag = np.array([tau[0, 0, 0, 0], tau[0, 0, 1, 1]])
        self.assertTrue(np.isfinite(diag).any())
        for v in diag[np.isfinite(diag)]:
            self.assertAlmostEqual(v, 0.0, places=6)


class TestTwoStateSolve(unittest.TestCase):
    def setUp(self):
        self.q = HIRESCPS()
        self.ras = [184.0, 186.0, 170.0, 176.0, 192.0, 198.0]
        self.decs = [13.0, 15.0, 5.0, 27.0, 51.0, -3.0]
        self.req = _twostate_requests(self.ras, self.decs, self.q)

    def test_two_state_schedules_and_assigns_wrap(self):
        tm = _solve(self.q, self.req, 2)
        self.assertIsNotNone(tm.schedule)
        sched = tm.schedule[tm.schedule["scheduled"]]
        self.assertGreater(len(sched), 0)
        self.assertIn("wrap_state", tm.schedule.columns)
        self.assertTrue(sched["wrap_state"].isin([0, 1]).all())

    def test_no_selected_arc_uses_unreachable_winding(self):
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
        self.assertIn("wrap_state", tm.schedule.columns)
        sched = tm.schedule[tm.schedule["scheduled"]]
        self.assertTrue((sched["wrap_state"] == 0).all())
        self.assertEqual(tm.stats["n_requested"], len(self.ras))
        self.assertGreaterEqual(tm.stats["t_slew_sum"], 0.0)
