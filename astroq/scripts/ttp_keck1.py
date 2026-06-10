"""Standalone TTP runner for HIRES-CPS on Keck-I.

Solve the Traveling Telescope Problem for one night, given a target CSV and
ISO night-start/stop times. Hard-wired to ``HIRESCPS`` so the recipe is
concrete; designed as a debug / planning-meeting tool, not a pipeline entry
point.

Usage:

    python -m astroq.scripts.ttp_keck1 REQUESTS_CSV NIGHT_START NIGHT_END [--outdir DIR]
        [--runtime SEC] [--optgap FRAC]

``REQUESTS_CSV`` must contain (at minimum):

    unique_id, ra, dec, exptime, n_exp, n_intra_max, tau_intra

``priority`` is optional (defaults to 10) and ``target`` improves hover text.
``NIGHT_START`` / ``NIGHT_END`` are ISO-8601 UTC, e.g. ``2026-05-09T05:30:00``.

How to adapt:

- For another telescope, swap ``HIRESCPS()`` for the queue you want (e.g.
  ``KPFCC()``); the rest of the pipeline is queue-agnostic.
- For richer accessibility (moon, past visits, custom windows), use the full
  ``astroq plan-night`` pipeline instead.
"""

# Standard library imports
import argparse
import os
import sys

# Third-party imports
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.time import Time
import astropy.units as u

# Local imports
from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp import plot as tplot
from astroq.ttp.model import TTPModel


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("requests_csv", help="Path to the request CSV.")
    ap.add_argument("night_start", help="ISO-8601 UTC start of the night.")
    ap.add_argument("night_end", help="ISO-8601 UTC end of the night.")
    ap.add_argument("--outdir", default="ttp_out",
                    help="Output directory (default: ./ttp_out).")
    ap.add_argument("--runtime", type=int, default=300,
                    help="Gurobi time limit, seconds (default: 300).")
    ap.add_argument("--optgap", type=float, default=0.01,
                    help="Gurobi MIP gap (default: 0.01).")
    return ap.parse_args()


def first_last_available(df, queue, night_start, night_end, *, n_samples=120):
    """Inline alt/az sweep; sets ``first_available`` / ``last_available`` ISO.

    Samples the night uniformly in JD, gates each (target, time) sample
    through ``queue.is_accessible``, and reduces with min/max to the first
    and last accessible JD per target. Targets with no accessible sample
    fall back to ``night_end``.
    """
    times = Time(np.linspace(night_start.jd, night_end.jd, n_samples), format="jd")
    coords = SkyCoord(df.ra.values * u.deg, df.dec.values * u.deg, frame="icrs")
    aa = queue.observatory.altaz(times, coords, grid_times_targets=True)
    ok = queue.is_accessible(aa.alt.deg, aa.az.deg)             # (n_targets, n_samples)

    jd = times.jd
    first_jd = np.where(ok, jd[None, :], np.inf).min(axis=1)
    last_jd = np.where(ok, jd[None, :], -np.inf).max(axis=1)
    no_good = ~ok.any(axis=1)
    first_jd[no_good] = night_end.jd
    last_jd[no_good] = night_end.jd

    df["first_available"] = [t[:16] for t in Time(first_jd, format="jd").iso]
    df["last_available"] = [t[:16] for t in Time(last_jd, format="jd").iso]
    return df


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.requests_csv)
    # Force numeric priority; the request CSV's ``p1``/``p2``/``p3`` strings
    # are user-facing labels and not multiplication-compatible with Gurobi
    # variables. Matches the default in ``nplan.NightPlanner``.
    df["priority"] = 10

    night_start = Time(args.night_start, format="isot")
    night_end = Time(args.night_end, format="isot")

    queue = HIRESCPS()                              # swap for other instruments (e.g. KPFCC())
    df = first_last_available(df, queue, night_start, night_end)

    # TTPModel expects astropy-typed columns inside a QTable (same recipe as
    # nplan.NightPlanner.run_ttp).
    df["n_intra_max"] = pd.to_numeric(df["n_intra_max"], errors="coerce").fillna(1)
    df["tau_intra"] = pd.to_numeric(
        df["tau_intra"].replace("None", np.nan), errors="coerce"
    ).fillna(0.0)
    visit_minutes = queue.visit_duration(
        df["exptime"].to_numpy(), df["n_exp"].to_numpy()
    )
    requests = QTable(
        {
            "unique_id": df["unique_id"].to_numpy(dtype=object),
            "coord": SkyCoord(
                df["ra"].to_numpy() * u.deg,
                df["dec"].to_numpy() * u.deg,
                frame="icrs",
            ),
            "first_available": Time(df["first_available"].tolist()),
            "last_available": Time(df["last_available"].tolist()),
            "t_visit": np.asarray(visit_minutes, dtype=float) * u.min,
            "n_intra_max": df["n_intra_max"].to_numpy(dtype=int),
            "tau_intra": df["tau_intra"].to_numpy(dtype=float) * u.hr,
            "priority": df["priority"].to_numpy(dtype=float),
        },
        copy=False,
    )

    tm = TTPModel(
        requests=requests,
        night_start=night_start,
        night_end=night_end,
        slew_fn=queue.slew_fn,
        n_slots=queue.nSlots,
    )
    tm.build_nodes()
    tm.build_arcs()
    tm.build_model()
    tm.model.params.TimeLimit = args.runtime
    tm.model.params.MIPGap = args.optgap
    tm.model.params.OutputFlag = 1
    tm.model.update()
    tm.run_model()
    if tm.model.SolCount == 0:
        sys.exit("TTP produced no schedule; exiting.")
    tm.build_schedule()
    print(tm.to_string())

    # Attach human-readable target (falls back to unique_id).
    names = dict(zip(df.unique_id, df.get("target", df.unique_id)))
    tm.schedule["target"] = tm.schedule.unique_id.map(names).fillna(tm.schedule.unique_id)

    schedule_csv = os.path.join(args.outdir, "schedule.csv")
    tm.schedule.to_csv(schedule_csv, index=False)
    print(f"Wrote {schedule_csv}")

    # Plot adapters expect observatory metadata on the model (set by plan-night / from_hdf5).
    tm.observer = queue.observatory
    tm.wrap_limit = queue.wrap_limit
    tm.slew_rate = queue.slew_rate
    tm.readout_time = queue.readout_time

    tplot.plot_path_2D_interactive(tm, night_start_time=night_start) \
        .write_html(os.path.join(args.outdir, "slew_path.html"))
    tplot.get_slew_animation_plotly(
        tm, args.requests_csv, inaccessible_zones=queue.inaccessible_zones,
    ).write_html(os.path.join(args.outdir, "slew_animation.html"))
    print(f"Wrote {args.outdir}/slew_path.html and {args.outdir}/slew_animation.html")


if __name__ == "__main__":
    main()
