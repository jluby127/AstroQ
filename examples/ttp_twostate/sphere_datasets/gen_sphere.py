"""Generate an N-target sphere benchmark for the TTP (Keck, 2026-02-01).

Targets are drawn uniformly on the celestial sphere, filtered to those with at
least one hour of Keck accessibility during nautical twilight. Exposure times
are set so the night fills exactly::

    night_minutes = n_targets * (t_visit + 0.5 min)

Usage (astroq-testing env, repo root)::

    python examples/ttp_twostate/sphere_datasets/gen_sphere.py --n 100
    python examples/ttp_twostate/sphere_datasets/gen_sphere.py --n 125 --out-dir custom/path
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.utils.iers import conf
import astropy.units as u

from astroq.queue.hirescps.queue import HIRESCPS

conf.auto_max_age = None

HERE = os.path.dirname(os.path.abspath(__file__))
NIGHT_DATE = "2026-02-01"
SLEW_BUDGET_MIN = 0.5
GRID_MIN = 2.0


def keck_nautical_night(observatory, date_iso):
    day = Time(date_iso, format="isot", scale="utc")
    start = observatory.twilight_evening_nautical(day, which="next")
    stop = observatory.twilight_morning_nautical(day, which="next")
    return start, stop


def uniform_sphere_radec(rng, n):
    u_cos = rng.uniform(-1.0, 1.0, n)
    dec = np.degrees(np.arcsin(u_cos))
    ra = rng.uniform(0.0, 360.0, n)
    return ra, dec


def observability_minutes(queue, coords, night_start, night_end, *, grid_min=2.0):
    n_times = max(int((night_end.jd - night_start.jd) * 24 * 60 / grid_min), 10)
    times = Time(np.linspace(night_start.jd, night_end.jd, n_times), format="jd")
    aa = queue.observatory.altaz(times, coords, grid_times_targets=True)
    ok = queue.is_accessible(aa.alt.deg, aa.az.deg)
    jd = times.jd
    first_jd = np.where(ok, jd[None, :], np.inf).min(axis=1)
    last_jd = np.where(ok, jd[None, :], -np.inf).max(axis=1)
    never = ~ok.any(axis=1)
    span_min = (last_jd - first_jd) * 24 * 60
    span_min[never] = 0.0
    return span_min, never


def sample_targets(queue, night_start, night_end, *, n_targets, min_obs_min, rng):
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
                f"{attempts} batches; relax min_obs_min or increase batch size"
            )
    return np.array(selected_ra), np.array(selected_dec), np.array(selected_span)


def generate(n_targets, *, seed=42, min_obs_min=60.0, out_dir=None):
    queue = HIRESCPS()
    night_start, night_end = keck_nautical_night(queue.observatory, NIGHT_DATE)
    night_min = (night_end.jd - night_start.jd) * 24 * 60

    slot_min = night_min / n_targets
    visit_min = slot_min - SLEW_BUDGET_MIN
    if visit_min <= 0:
        raise ValueError("night too short for the requested target count / slew budget")
    exptime_s = visit_min * 60.0

    rng = np.random.default_rng(seed)
    ra, dec, span = sample_targets(
        queue,
        night_start,
        night_end,
        n_targets=n_targets,
        min_obs_min=min_obs_min,
        rng=rng,
    )

    if out_dir is None:
        out_dir = os.path.join(HERE, f"sphere{n_targets}_2026-02-01")
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    for k, (r, d) in enumerate(zip(ra, dec)):
        rows.append(
            {
                "unique_id": f"S{k:03d}",
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
    csv_path = os.path.join(out_dir, "request_selected.csv")
    df.to_csv(csv_path, index=False)

    meta = {
        "night_date": NIGHT_DATE,
        "night_start": night_start.isot,
        "night_end": night_end.isot,
        "night_minutes": round(night_min, 3),
        "n_targets": n_targets,
        "min_observability_minutes": min_obs_min,
        "slew_budget_minutes_per_target": SLEW_BUDGET_MIN,
        "slot_minutes_per_target": round(slot_min, 4),
        "visit_minutes_per_target": round(visit_min, 4),
        "exptime_seconds": round(exptime_s, 1),
        "packing_identity": (
            f"{n_targets} * ({round(visit_min, 4)} + {SLEW_BUDGET_MIN}) "
            f"= {round(n_targets * (visit_min + SLEW_BUDGET_MIN), 3)} "
            f"(night {round(night_min, 3)})"
        ),
        "seed": seed,
        "observability_minutes": {
            "min": round(float(span.min()), 1),
            "median": round(float(np.median(span)), 1),
            "max": round(float(span.max()), 1),
        },
    }
    meta_path = os.path.join(out_dir, "meta.json")
    with open(meta_path, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")

    return out_dir, df, meta, night_start, night_end


def main():
    parser = argparse.ArgumentParser(description="Generate sphere N-target benchmark")
    parser.add_argument("--n", type=int, required=True, help="number of targets")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-obs-min", type=float, default=60.0)
    parser.add_argument("--out-dir", default=None)
    args = parser.parse_args()

    out_dir, df, meta, night_start, night_end = generate(
        args.n,
        seed=args.seed,
        min_obs_min=args.min_obs_min,
        out_dir=args.out_dir,
    )
    print(f"Wrote {out_dir}/request_selected.csv ({len(df)} targets)")
    print(f"  night : {night_start.isot} -> {night_end.isot}")
    print(f"  packing: {meta['packing_identity']}")


if __name__ == "__main__":
    main()
