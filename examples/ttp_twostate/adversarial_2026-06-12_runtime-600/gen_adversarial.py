"""Generate an adversarial 58-target night that forces south<->NW slews.

The legacy single-state slew model (``wrap_limit = 235``) charges any slew that
crosses sky az 235 deg the long way around, while the two-state model keeps the
whole western sky (az 90-315) continuous in the South wrap. To expose that, we
synthesize two interleaved buckets of targets:

- South bucket (29): ``dec`` below the Keck latitude, placed near az ~180 (due
  south) at their window center -- never west of ~225, so they never cross 235.
- NW bucket (29): ``dec`` > 20 (transits north), placed setting in the NW at
  az ~300 (kept <= 315 so they stay inside the South wrap).

Each target gets a tight window centered on a time ``t_k``; the two buckets are
interleaved in time (S, NW, S, NW, ...) so the feasible visit order is forced to
alternate -- which is cheap for two-state but pays repeated 235 crossings for
single-state.

Every target uses a fixed ~10 min single-shot visit (``exptime = 600 s``).

Usage:
    python examples/ttp_twostate/adversarial_2026-06-12_runtime-600/gen_adversarial.py
"""

import os

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
import astropy.units as u

from astroq.queue.hirescps.queue import HIRESCPS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_CSV = os.path.join(HERE, "request_selected.csv")

NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
NIGHT_END = Time("2026-06-12T14:48:00", format="isot")

N_PER = 29                 # targets per bucket -> 58 total
EXPTIME_S = 600            # 10 min single-shot visit
DELTA_MIN = 15.0           # half-window (min); full window = 2*DELTA_MIN

# Bucket geometry: (target sky az, allowed az band, dec values cycled).
SOUTH_AZ = 180.0
SOUTH_BAND = (120.0, 232.0)        # stay east of the 235 cut
SOUTH_DECS = np.linspace(-20.0, 10.0, N_PER)

NW_AZ = 300.0
NW_BAND = (272.0, 314.0)           # NW, inside South wrap (<315)
NW_DECS = np.linspace(22.0, 38.0, N_PER)

ALT_MIN, ALT_MAX = 20.0, 84.0      # comfortably inside the 18/85 clamps


def find_ra(queue, dec, t, target_az, az_lo, az_hi):
    """Scan RA (deg) for the pointing closest to ``target_az`` at time ``t``.

    Returns (ra_deg, az_deg, alt_deg) for the best alt-feasible RA whose azimuth
    falls in [az_lo, az_hi], or ``None`` if no RA satisfies the band + altitude.
    """
    ra_grid = np.arange(0.0, 360.0, 0.25)
    coords = SkyCoord(ra_grid * u.deg, np.full_like(ra_grid, dec) * u.deg, frame="icrs")
    aa = queue.observatory.altaz(Time([t]), coords, grid_times_targets=True)
    az = aa.az.deg[:, 0]
    alt = aa.alt.deg[:, 0]
    ok = (alt >= ALT_MIN) & (alt <= ALT_MAX) & (az >= az_lo) & (az <= az_hi)
    if not ok.any():
        return None
    idx = np.where(ok)[0]
    best = idx[np.argmin(np.abs(az[idx] - target_az))]
    return float(ra_grid[best]), float(az[best]), float(alt[best])


def build():
    queue = HIRESCPS()
    night_len_min = (NIGHT_END - NIGHT_START).to_value(u.min)

    rows = []
    # Interleave the two buckets in time: even slot -> South, odd slot -> NW.
    n_total = 2 * N_PER
    s_i = nw_i = 0
    for k in range(n_total):
        frac = (k + 0.5) / n_total
        t_k = NIGHT_START + TimeDelta(frac * night_len_min * 60.0, format="sec")
        if k % 2 == 0:
            bucket, dec = "S", SOUTH_DECS[s_i]
            target_az, (lo, hi) = SOUTH_AZ, SOUTH_BAND
            s_i += 1
        else:
            bucket, dec = "NW", NW_DECS[nw_i]
            target_az, (lo, hi) = NW_AZ, NW_BAND
            nw_i += 1

        hit = find_ra(queue, dec, t_k, target_az, lo, hi)
        if hit is None:
            # Fall back: widen the band to anything alt-feasible nearest target az.
            hit = find_ra(queue, dec, t_k, target_az, 0.0, 360.0)
        ra, az, alt = hit

        w0 = t_k - TimeDelta(DELTA_MIN * 60.0, format="sec")
        w1 = t_k + TimeDelta(DELTA_MIN * 60.0, format="sec")
        rows.append(
            {
                "program_code": f"ADV_{bucket}",
                "target": f"adv_{bucket.lower()}_{(s_i if bucket=='S' else nw_i):02d}",
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
                "az_center": round(az, 2),
                "alt_center": round(alt, 2),
                "window_start_ut": w0.isot,
                "window_end_ut": w1.isot,
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    n_s = (df.bucket == "S").sum()
    n_nw = (df.bucket == "NW").sum()
    print(f"wrote {OUT_CSV}  ({len(df)} targets: {n_s} south, {n_nw} NW)")
    print(df[["unique_id", "bucket", "dec", "az_center", "alt_center"]].to_string(index=False))


if __name__ == "__main__":
    build()
