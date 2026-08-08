"""Generate a staggered adversarial 58-target night (v4: random dec sampling).

Same staggered rise/set design as v2/v3, but the per-bucket declinations are
drawn *randomly* (uniform) within their range instead of placed on a regular
grid, so the per-slew costs sample a realistic spread of dec separations:

- South bucket: 29 targets with ``dec`` ~ Uniform[-40, -20].
- North bucket: 29 targets with ``dec`` ~ Uniform[+20, +40].

Within each bucket targets are split by timing class (even split):

- early-setting (~7/bucket): up at dusk, drops below access before t0 + 0.25 L.
- late-rising  (~7/bucket): not accessible until after t0 + 0.75 L.
- mid          (~15/bucket): broadly available, transiting mid-night.

Net over 58: ~14 early-setting (1/4) + ~14 late-rising (1/4) + ~30 mid (1/2).

RA is chosen empirically: for each fixed ``dec`` we scan RA, compute the real
accessibility window over the night (``queue.is_accessible``), and pick the RA
whose window matches the requested timing class. Every target uses a fixed
~10 min single-shot visit (``exptime = 600 s``). The dec draws use a fixed seed
so the set is reproducible.

Usage:
    python examples/ttp_twostate/adversarial-4_2026-06-12_runtime-600/gen_adversarial.py
"""

import os

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

from astroq.queue.hirescps.queue import HIRESCPS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_CSV = os.path.join(HERE, "request_selected.csv")

NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
NIGHT_END = Time("2026-06-12T14:48:00", format="isot")

N_PER = 29                 # per dec bucket -> 58 total
N_SET = 7                  # early-setting per bucket
N_RISE = 7                 # late-rising per bucket
N_MID = N_PER - N_SET - N_RISE   # 15
EXPTIME_S = 600            # 10 min single-shot visit

# v4: draw declinations randomly (uniform) within each bucket's range.
SEED = 42
SOUTH_DEC_RANGE = (-40.0, -20.0)
NORTH_DEC_RANGE = (20.0, 40.0)

# Night fractions: set before Q1 ends, rise after 3/4.
F_Q1 = 0.25
F_Q4 = 0.75

RA_STEP = 1.0              # deg
TIME_STEP_MIN = 5.0        # min


def ra_window_grid(queue, dec, t_grid):
    """Accessibility window fractions (first, last) for each RA on a grid.

    Returns (ra_grid, frac_first, frac_last) where the fracs are fractions of
    the night [0,1]; NaN where the target is never accessible at that RA.
    """
    ra_grid = np.arange(0.0, 360.0, RA_STEP)
    coords = SkyCoord(ra_grid * u.deg, np.full_like(ra_grid, dec) * u.deg, frame="icrs")
    aa = queue.observatory.altaz(t_grid, coords, grid_times_targets=True)
    ok = queue.is_accessible(aa.alt.deg, aa.az.deg)  # (nRA, nT)
    jd = t_grid.jd
    L = jd[-1] - jd[0]
    ff = np.where(ok, (jd[None, :] - jd[0]) / L, np.nan)
    frac_first = np.nanmin(ff, axis=1)
    frac_last = np.nanmax(ff, axis=1)
    return ra_grid, frac_first, frac_last


def pick_ra(ra_grid, frac_first, frac_last, timing, target_frac):
    """Pick the RA index for the requested timing class nearest ``target_frac``.

    ``target_frac`` is the desired night-fraction of the key event for the class
    (set time for "set", rise time for "rise", window center for "mid"). Varying
    it across the targets in a class spreads them out in RA while preserving the
    timing semantics.
    """
    dur = frac_last - frac_first
    valid = np.isfinite(frac_first) & (dur > 0.03)  # >~1 visit of coverage
    if timing == "set":
        # up near dusk, gone before end of Q1; choose set-time nearest target
        m = valid & (frac_first <= 0.06) & (frac_last <= F_Q1 + 0.02)
        if not m.any():
            m = valid & (frac_first <= 0.12) & (frac_last <= F_Q1 + 0.06)
        if not m.any():
            return None
        return int(np.argmin(np.where(m, np.abs(frac_last - target_frac), np.inf)))
    if timing == "rise":
        # rises after the 3/4 mark, up through dawn; rise-time nearest target
        m = valid & (frac_last >= 0.94) & (frac_first >= F_Q4 - 0.02)
        if not m.any():
            m = valid & (frac_last >= 0.88) & (frac_first >= F_Q4 - 0.08)
        if not m.any():
            return None
        return int(np.argmin(np.where(m, np.abs(frac_first - target_frac), np.inf)))
    # mid: broadly available; pick window-center nearest target_frac
    center = (frac_first + frac_last) / 2
    m = valid & (dur > 0.45)
    if not m.any():
        m = valid
    return int(np.argmin(np.where(m, np.abs(center - target_frac), np.inf)))


def az_alt_at(queue, ra, dec, t):
    aa = queue.observatory.altaz(Time([t]), SkyCoord(ra * u.deg, dec * u.deg, frame="icrs"))
    return float(aa.az.deg[0]), float(aa.alt.deg[0])


def build():
    queue = HIRESCPS()
    rng = np.random.default_rng(SEED)
    south_decs = rng.uniform(*SOUTH_DEC_RANGE, size=N_PER)
    north_decs = rng.uniform(*NORTH_DEC_RANGE, size=N_PER)

    n_t = max(int((NIGHT_END.jd - NIGHT_START.jd) * 24 * 60 / TIME_STEP_MIN), 10)
    t_grid = Time(np.linspace(NIGHT_START.jd, NIGHT_END.jd, n_t), format="jd")
    L_jd = NIGHT_END.jd - NIGHT_START.jd

    # Spread each class across a range of event-times so the targets fan out in
    # RA instead of stacking at one RA. set: set-time across Q1; rise: rise-time
    # across Q4; mid: window-center across the bulk of the night.
    set_fracs = np.linspace(0.08, F_Q1, N_SET)
    rise_fracs = np.linspace(F_Q4, 0.95, N_RISE)
    mid_fracs = np.linspace(0.20, 0.80, N_MID)
    timing_seq = (
        [("set", f) for f in set_fracs]
        + [("rise", f) for f in rise_fracs]
        + [("mid", f) for f in mid_fracs]
    )

    rows = []
    k = 0
    for bucket, decs in (("S", south_decs), ("N", north_decs)):
        for j, dec in enumerate(decs):
            timing, target_frac = timing_seq[j]
            ra_grid, ff, fl = ra_window_grid(queue, dec, t_grid)
            idx = pick_ra(ra_grid, ff, fl, timing, target_frac)
            if idx is None:  # fall back to widest window
                idx = pick_ra(ra_grid, ff, fl, "mid", 0.5)
                timing = timing + "?"
            ra = float(ra_grid[idx])
            frac_first, frac_last = float(ff[idx]), float(fl[idx])
            # az/alt at window midpoint
            t_mid = Time(NIGHT_START.jd + (frac_first + frac_last) / 2 * L_jd, format="jd")
            az, alt = az_alt_at(queue, ra, dec, t_mid)
            rows.append(
                {
                    "program_code": f"ADV_{bucket}",
                    "target": f"adv_{bucket.lower()}_{j:02d}",
                    "unique_id": f"ADV{k:03d}",
                    "ra": round(ra, 5),
                    "dec": round(float(dec), 5),
                    "exptime": EXPTIME_S,
                    "maxtime": EXPTIME_S,
                    "n_exp": 1,
                    "n_inter_max": 1,
                    "tau_inter": 1,
                    "n_intra_max": 1,
                    "n_intra_min": 1,
                    "tau_intra": 0.0,
                    "minimum_elevation": 18,
                    "minimum_moon_separation": 0,
                    "priority": 10,
                    "bucket": bucket,
                    "timing_class": timing,
                    "az_mid": round(az, 1),
                    "alt_mid": round(alt, 1),
                    "frac_first": round(frac_first, 3),
                    "frac_last": round(frac_last, 3),
                }
            )
            k += 1

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    n_set = (df.timing_class.str.startswith("set")).sum()
    n_rise = (df.timing_class.str.startswith("rise")).sum()
    n_mid = (df.timing_class.str.startswith("mid")).sum()
    print(f"wrote {OUT_CSV}  ({len(df)} targets: {n_set} early-set, "
          f"{n_rise} late-rise, {n_mid} mid)")
    print(df[["unique_id", "bucket", "timing_class", "dec", "ra",
              "az_mid", "alt_mid", "frac_first", "frac_last"]].to_string(index=False))


if __name__ == "__main__":
    build()
