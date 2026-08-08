"""Shared observability-window request builder for sphere benchmarks."""

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.time import Time
import astropy.units as u


def build_requests_access(df, queue, night_start, night_end, *, grid_min=2.0):
    """QTable whose availability is the real observability window."""
    coords = SkyCoord(df.ra.values * u.deg, df.dec.values * u.deg, frame="icrs")
    visit_min = queue.visit_duration(df.exptime.to_numpy(), df.n_exp.to_numpy())

    n = max(int((night_end.jd - night_start.jd) * 24 * 60 / grid_min), 10)
    t = Time(np.linspace(night_start.jd, night_end.jd, n), format="jd")
    aa = queue.observatory.altaz(t, coords, grid_times_targets=True)
    ok = queue.is_accessible(aa.alt.deg, aa.az.deg)

    jd = t.jd
    first_jd = np.where(ok, jd[None, :], np.inf).min(axis=1)
    last_jd = np.where(ok, jd[None, :], -np.inf).max(axis=1)
    no_good = ~ok.any(axis=1)
    first_jd[no_good] = night_start.jd
    last_jd[no_good] = night_end.jd

    return QTable(
        {
            "unique_id": df.unique_id.to_numpy(dtype=object),
            "coord": coords,
            "first_available": Time(first_jd, format="jd"),
            "last_available": Time(last_jd, format="jd"),
            "t_visit": np.asarray(visit_min, dtype=float) * u.min,
            "n_intra_max": df.n_intra_max.to_numpy(dtype=int),
            "tau_intra": df.tau_intra.to_numpy(dtype=float) * u.hr,
            "priority": df.priority.to_numpy(dtype=float),
        },
        copy=False,
    )
