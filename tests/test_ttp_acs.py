"""Tests for the ACS warm-start heuristic (astroq.ttp.acs) and its integration.

Self-contained: small synthetic request sets solved directly, mirroring
``tests/test_ttp_twostate.py``. The ACS now returns only a tour dict
(:func:`astroq.ttp.acs.acs_warm_start`); a schedule is produced by seeding the
MILP and running Gurobi (see ``test_seed_accepted_by_gurobi``).
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
            "first_available": Time([NIGHT_START.isot] * n),
            "last_available": Time([NIGHT_END.isot] * n),
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
            self.assertLessEqual(ti, float(tm.nodes.at[nid, "t_late"]) + 1e-6)
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
                "first_available": Time([NIGHT_START.isot] * 2),
                "last_available": Time([NIGHT_END.isot] * 2),
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
