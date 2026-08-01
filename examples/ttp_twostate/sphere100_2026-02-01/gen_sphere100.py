"""Generate a 100-target sphere benchmark for the TTP (Keck, 2026-02-01).

Targets are drawn uniformly on the celestial sphere, filtered to those with at
least one hour of Keck accessibility (``queue.is_accessible``) during nautical
twilight on 2026-02-01. Exposure times are set so that, if every target could be
observed with exactly 0.5 min slew between visits, the night fills exactly:

    night_minutes = n_targets * (t_visit + 0.5 min)

With ``n_exp = 1`` this implies ``exptime = 60 * (night_min / n - 0.5)`` seconds
(visit duration equals exposure time; no readout between shots).

Writes ``request_selected.csv`` and ``meta.json`` next to this script.

Usage (astroq-testing env, repo root):

    python examples/ttp_twostate/sphere100_2026-02-01/gen_sphere100.py
"""

import json
import os

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.utils.iers import conf
import astropy.units as u

from astroq.queue.hirescps.queue import HIRESCPS

conf.auto_max_age = None  # allow future dates without live IERS download

HERE = os.path.dirname(os.path.abspath(__file__))
NIGHT_DATE = "2026-02-01"
N_TARGETS = 100
MIN_OBS_MIN = 60.0
SLEW_BUDGET_MIN = 0.5
GRID_MIN = 2.0
SEED = 42


def keck_nautical_night(observatory, date_iso):
    """Nautical-twilight (-12 deg) bounds for one Keck night."""
    day = Time(date_iso, format="isot", scale="utc")
    start = observatory.twilight_evening_nautical(day, which="next")
    stop = observatory.twilight_morning_nautical(day, which="next")
    return start, stop


def uniform_sphere_radec(rng, n):
    """Uniform random directions on the sphere (ra deg, dec deg)."""
    u_cos = rng.uniform(-1.0, 1.0, n)
    dec = np.degrees(np.arcsin(u_cos))
    ra = rng.uniform(0.0, 360.0, n)
    return ra, dec


def observability_minutes(queue, coords, night_start, night_end, *, grid_min=2.0):
    """Per-target accessible span (minutes) over the night grid."""
    n_times = max(int((night_end.jd - night_start.jd) * 24 * 60 / grid_min), 10)
    times = Time(np.linspace(night_start.jd, night_end.jd, n_times), format="jd")
    aa = queue.observatory.altaz(times, coords, grid_times_targets=True)
    ok = queue.is_accessible(aa.alt.deg, aa.az.deg)  # (ntargets, ntimes)
    jd = times.jd
    first_jd = np.where(ok, jd[None, :], np.inf).min(axis=1)
    last_jd = np.where(ok, jd[None, :], -np.inf).max(axis=1)
    never = ~ok.any(axis=1)
    span_min = (last_jd - first_jd) * 24 * 60
    span_min[never] = 0.0
    return span_min, never


def sample_targets(queue, night_start, night_end, *, n_targets, min_obs_min, rng):
    """Reject/accept random sphere directions until ``n_targets`` qualify."""
    selected_ra, selected_dec, selected_span = [], [], []
    batch = max(n_targets * 4, 400)
    attempts = 0
    while len(selected_ra) < n_targets:
        attempts += 1
        ra, dec = uniform_sphere_radec(rng, batch)
        coords = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
        span, never = observability_minutes(
            queue, coords, night_start, night_end, grid_min=GRID_MIN
        )
        good = (~never) & (span >= min_obs_min)
        for r, d, s in zip(ra[good], dec[good], span[good]):
            if len(selected_ra) >= n_targets:
                break
            selected_ra.append(float(r))
            selected_dec.append(float(d))
            selected_span.append(float(s))
        if attempts > 200:
            raise RuntimeError(
                f"only found {len(selected_ra)}/{n_targets} targets after "
                f"{attempts} batches; relax MIN_OBS_MIN or increase batch size"
            )
    return np.array(selected_ra), np.array(selected_dec), np.array(selected_span)


def main():
    queue = HIRESCPS()
    night_start, night_end = keck_nautical_night(queue.observatory, NIGHT_DATE)
    night_min = (night_end.jd - night_start.jd) * 24 * 60

    slot_min = night_min / N_TARGETS
    visit_min = slot_min - SLEW_BUDGET_MIN
    if visit_min <= 0:
        raise ValueError("night too short for the requested target count / slew budget")
    exptime_s = visit_min * 60.0  # n_exp=1 => visit_min = exptime/60

    rng = np.random.default_rng(SEED)
    ra, dec, span = sample_targets(
        queue,
        night_start,
        night_end,
        n_targets=N_TARGETS,
        min_obs_min=MIN_OBS_MIN,
        rng=rng,
    )

    rows = []
    for k, (r, d, s) in enumerate(zip(ra, dec, span)):
        uid = f"S{k:03d}"
        rows.append(
            {
                "unique_id": uid,
                "target": f"sphere_{k:03d}",
                "ra": round(r, 6),
                "dec": round(d, 6),
                "exptime": round(exptime_s, 1),
                "n_exp": 1,
                "n_intra_max": 1,
                "tau_intra": 0.0,
                "priority": 10.0,
            }
        )
    df = pd.DataFrame(rows)
    csv_path = os.path.join(HERE, "request_selected.csv")
    df.to_csv(csv_path, index=False)

    meta = {
        "night_date": NIGHT_DATE,
        "night_start": night_start.isot,
        "night_end": night_end.isot,
        "night_minutes": round(night_min, 3),
        "n_targets": N_TARGETS,
        "min_observability_minutes": MIN_OBS_MIN,
        "slew_budget_minutes_per_target": SLEW_BUDGET_MIN,
        "slot_minutes_per_target": round(slot_min, 4),
        "visit_minutes_per_target": round(visit_min, 4),
        "exptime_seconds": round(exptime_s, 1),
        "packing_identity": (
            f"{N_TARGETS} * ({round(visit_min, 4)} + {SLEW_BUDGET_MIN}) "
            f"= {round(N_TARGETS * (visit_min + SLEW_BUDGET_MIN), 3)} "
            f"(night {round(night_min, 3)})"
        ),
        "seed": SEED,
        "observability_minutes": {
            "min": round(float(span.min()), 1),
            "median": round(float(np.median(span)), 1),
            "max": round(float(span.max()), 1),
        },
    }
    meta_path = os.path.join(HERE, "meta.json")
    with open(meta_path, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")

    print(f"Wrote {csv_path} ({len(df)} targets)")
    print(f"Wrote {meta_path}")
    print(f"  night : {night_start.isot} -> {night_end.isot}  ({night_min:.1f} min)")
    print(f"  slot  : {slot_min:.3f} min/target  (visit {visit_min:.3f} + slew {SLEW_BUDGET_MIN})")
    print(f"  exptime: {exptime_s:.1f} s  (n_exp=1)")
    print(
        f"  obs window (min): min={span.min():.0f}  "
        f"median={np.median(span):.0f}  max={span.max():.0f}"
    )


if __name__ == "__main__":
    main()
