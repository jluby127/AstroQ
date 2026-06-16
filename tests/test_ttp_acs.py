"""Tests for the ACS heuristic (astroq.ttp.acs) and its model integration.

Self-contained: small synthetic request sets solved directly, mirroring
``tests/test_ttp_twostate.py``.
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


def _requests(ras, decs, *, t_visit_min=5.0):
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


RAS = [184.0, 186.0, 170.0, 176.0, 192.0, 198.0]
DECS = [13.0, 15.0, 5.0, 27.0, 51.0, -3.0]


class TestACSHeuristic(unittest.TestCase):
    def setUp(self):
        self.req = _requests(RAS, DECS)

    def test_single_state_feasible_and_windows(self):
        tm = _build(self.req, 1)
        res = tm.run_heuristic(params={"time_limit_s": 1.0})
        self.assertTrue(res["feasible"])
        sched = tm.schedule[tm.schedule["scheduled"]]
        self.assertGreater(len(sched), 0)
        # completion time within each node's window and the night budget.
        self.assertTrue((sched["t_end"] <= sched["t_late"] + 1e-6).all())
        self.assertTrue((sched["t_end"] <= tm.dur_min + 1e-6).all())
        self.assertGreaterEqual(tm.stats["t_slew_sum"], 0.0)

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
        res = tm.run_heuristic(params={"time_limit_s": 2.0})
        self.assertLessEqual(res["objective"], opt + 1e-6)
        self.assertGreaterEqual(res["objective"], opt - 1e-6)

    def test_seed_accepted_by_gurobi(self):
        """The MIPStart from the ACS tour must be loadable (feasible)."""
        tm = _build(self.req, 1)
        tm.build_model()
        tm.seed_from_tour(tm.run_heuristic(params={"time_limit_s": 1.0}))
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
        res = tm.run_heuristic(params={"time_limit_s": 1.0})
        self.assertTrue(res["feasible"])
        sched = tm.schedule[tm.schedule["scheduled"]]
        self.assertGreater(len(sched), 0)
        self.assertIn("wrap_state", tm.schedule.columns)
        self.assertTrue(sched["wrap_state"].isin([0, 1]).all())

    def test_tau_intra_separation_respected(self):
        """Two visits of one target must be spaced by tau_intra."""
        req = QTable(
            {
                "unique_id": np.array(["A", "B"], dtype=object),
                "coord": SkyCoord([184.0, 186.0] * u.deg, [13.0, 15.0] * u.deg,
                                  frame="icrs"),
                "first_available": Time([NIGHT_START.isot] * 2),
                "last_available": Time([NIGHT_END.isot] * 2),
                "t_visit": np.full(2, 5.0) * u.min,
                "n_intra_max": np.array([2, 1], dtype=int),
                "tau_intra": np.array([0.5, 0.0]) * u.hr,  # 30 min for A
                "priority": np.full(2, 10.0),
            },
            copy=False,
        )
        tm = _build(req, 1)
        res = tm.run_heuristic(params={"time_limit_s": 1.0})
        self.assertTrue(res["feasible"])
        sched = tm.schedule[tm.schedule["scheduled"]].sort_values("t_start")
        a_visits = sched[sched["unique_id"] == "A"].sort_values("t_end")
        if len(a_visits) == 2:
            gap = a_visits["t_end"].to_numpy()
            self.assertGreaterEqual(gap[1] - gap[0], 30.0 - 1e-6)


if __name__ == "__main__":
    unittest.main()
