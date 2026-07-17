"""Demonstration: single-state vs two-state (cable-wrap) TTP on Keck-I.

Solves the same set of requests two ways -- the legacy single-cut slew model
(``wrap_limit = 235``) and the state-aware two-wrap model -- then re-evaluates
*both* resulting visit orders under the physically correct encoder-azimuth slew
model. The legacy schedule tends to incur long "unwind" slews on western,
setting targets that straddle the 235 deg cut; the two-state schedule avoids
them by observing those targets in the south wrap.

Usage:

    python -m astroq.scripts.demo_twostate \
        [REQUESTS_CSV] [NIGHT_START] [NIGHT_END]

Defaults reproduce the committed demo (``demo_twostate_requests.csv``, a half
night on 2026-05-09).
"""

# Standard library imports
import argparse
import os

# Third-party imports
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.time import Time, TimeDelta
import astropy.units as u

# Local imports
from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel

_HERE = os.path.dirname(__file__)
_DEFAULT_CSV = os.path.join(_HERE, "demo_twostate_requests.csv")


def build_requests(df, queue, night_start, night_end, *, n_samples=120):
    """QTable for ``TTPModel`` with full-window accessibility per target."""
    times = Time(np.linspace(night_start.jd, night_end.jd, n_samples), format="jd")
    coords = SkyCoord(df.ra.values * u.deg, df.dec.values * u.deg, frame="icrs")
    aa = queue.observatory.altaz(times, coords, grid_times_targets=True)
    ok = queue.is_accessible(aa.alt.deg, aa.az.deg)
    jd = times.jd
    first_jd = np.where(ok, jd[None, :], np.inf).min(axis=1)
    last_jd = np.where(ok, jd[None, :], -np.inf).max(axis=1)
    no_good = ~ok.any(axis=1)
    first_jd[no_good] = night_start.jd
    last_jd[no_good] = night_end.jd
    visit_min = queue.visit_duration(df.exptime.to_numpy(), df.n_exp.to_numpy())
    return QTable(
        {
            "unique_id": df.unique_id.to_numpy(dtype=object),
            "coord": coords,
            "time_earliest_start": Time(first_jd, format="jd"),
            "time_latest_finish": Time(last_jd, format="jd"),
            "t_visit": np.asarray(visit_min, dtype=float) * u.min,
            "n_intra_max": df.n_intra_max.to_numpy(dtype=int),
            "tau_intra": df.tau_intra.to_numpy(dtype=float) * u.hr,
            "weight": np.ones(len(df), dtype=float),
        },
        copy=False,
    )


def solve(requests, queue, night_start, night_end, n_states, *, runtime, optgap):
    slew_fn = queue.slew_fn_state if n_states > 1 else queue.slew_fn
    tm = TTPModel(
        requests=requests,
        night_start=night_start,
        night_end=night_end,
        slew_fn=slew_fn,
        n_slots=queue.nSlots,
        n_states=n_states,
    )
    tm.build_nodes()
    tm.build_arcs()
    tm.build_model()
    tm.model.params.OutputFlag = 0
    tm.model.params.TimeLimit = runtime
    tm.model.params.MIPGap = optgap
    tm.model.update()
    tm.run_model()
    tm.build_schedule()
    return tm


def physical_slew_minutes(schedule, queue, night_start):
    """Min total *physical* slew (min) for a schedule's visit order.

    Re-evaluates the chosen visit order under the encoder-azimuth wrap model:
    each target is resolved into the wrap state(s) reachable at its observed
    time, and a shortest-path DP over those states returns the minimum
    achievable total slew respecting wrap continuity. ``inf`` means the order
    is physically impossible in any wrap assignment.
    """
    sc = schedule[schedule["scheduled"]].sort_values("order")
    if len(sc) == 0:
        return 0.0
    t_obs = night_start + TimeDelta(sc["t_start"].to_numpy() * 60.0, format="sec")
    coords = SkyCoord(sc.ra.values * u.deg, sc.dec.values * u.deg, frame="icrs")
    aa = queue.observatory.altaz(t_obs, coords)
    az = np.atleast_1d(aa.az.deg)
    alt = np.atleast_1d(aa.alt.deg)
    states = queue.wrap_states
    S = len(states)
    # enc[node, state]
    enc = np.stack(
        [queue._encoder_az(az, lo, hi) for _, lo, hi in states], axis=1
    )
    rate = 60.0 * float(queue.slew_rate)

    n = len(sc)
    INF = float("inf")
    # dp[s] = min slew to reach node k in state s
    dp = np.where(np.isfinite(enc[0]), 0.0, INF)
    for k in range(1, n):
        nxt = np.full(S, INF)
        for sj in range(S):
            if not np.isfinite(enc[k, sj]):
                continue
            best = INF
            for si in range(S):
                if not np.isfinite(dp[si]) or not np.isfinite(enc[k - 1, si]):
                    continue
                d = max(abs(enc[k - 1, si] - enc[k, sj]),
                        abs(alt[k - 1] - alt[k])) / rate
                best = min(best, dp[si] + d)
            nxt[sj] = best
        dp = nxt
    return float(np.min(dp))


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("requests_csv", nargs="?", default=_DEFAULT_CSV)
    ap.add_argument("night_start", nargs="?", default="2026-05-09T06:00:00")
    ap.add_argument("night_end", nargs="?", default="2026-05-09T09:30:00")
    ap.add_argument("--runtime", type=int, default=120)
    ap.add_argument("--optgap", type=float, default=0.01)
    return ap.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.requests_csv)
    night_start = Time(args.night_start, format="isot")
    night_end = Time(args.night_end, format="isot")
    queue = HIRESCPS()
    requests = build_requests(df, queue, night_start, night_end)

    tm1 = solve(requests, queue, night_start, night_end, 1,
                runtime=args.runtime, optgap=args.optgap)
    tm2 = solve(requests, queue, night_start, night_end, 2,
                runtime=args.runtime, optgap=args.optgap)

    phys1 = physical_slew_minutes(tm1.schedule, queue, night_start)
    phys2 = physical_slew_minutes(tm2.schedule, queue, night_start)

    print(f"Requests: {tm1.stats['n_requested']}   "
          f"night: {args.night_start} -> {args.night_end}\n")
    hdr = f"{'model':<14}{'scheduled':>10}{'modeled slew':>15}{'physical slew':>16}"
    print(hdr)
    print("-" * len(hdr))
    print(f"{'single-state':<14}{tm1.stats['n_scheduled']:>10}"
          f"{tm1.stats['t_slew_sum']:>14.2f}m{phys1:>15.2f}m")
    print(f"{'two-state':<14}{tm2.stats['n_scheduled']:>10}"
          f"{tm2.stats['t_slew_sum']:>14.2f}m{phys2:>15.2f}m")
    print("\n'physical slew' re-evaluates each visit order under the encoder "
          "wrap model.\nThe single-state order pays for long western unwinds "
          "the legacy metric hides.")

    sc = tm2.schedule[tm2.schedule.scheduled].sort_values("order")
    names = dict(zip(("N", "S"), ("N", "S")))
    wraps = [list(queue.wrap_states[int(s)][0:1])[0] for s in sc.wrap_state]
    print("\ntwo-state order (target : wrap):")
    print("  " + "  ".join(f"{u}:{w}" for u, w in zip(sc.unique_id, wraps)))


if __name__ == "__main__":
    main()
