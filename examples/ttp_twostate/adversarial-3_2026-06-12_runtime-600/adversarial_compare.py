"""Single- vs two-state TTP on the staggered adversarial target set.

Availability comes from real observability (``queue.is_accessible`` over the
whole night); the target RAs are staggered (see ``gen_adversarial.py``) so a
quarter set in the first quarter of the night and a quarter rise only in the
last quarter, breaking the trivial "all-south-first, cross once" order.

For each mode this persists ``schedule.csv`` + slew_path/slew_animation/ladder
plots and a Gurobi log, and writes ``comparison.json`` / ``comparison.txt`` with
modeled slew, physical slew, MIP gap, the count of sky-az 235 crossings, and how
many early-setting targets each model managed to grab.

Usage (astroq-testing env, from repo root):
    python examples/ttp_twostate/adversarial-2_2026-06-12_runtime-600/adversarial_compare.py        # both
    python examples/ttp_twostate/adversarial-2_2026-06-12_runtime-600/adversarial_compare.py 2      # two-state only
"""

import json
import os
import sys

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.time import Time, TimeDelta
import astropy.units as u

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel
import astroq.ttp.plot as tplot
import astroq.plot as aqplot
from astroq.scripts.demo_twostate import physical_slew_minutes

OUT = os.path.dirname(os.path.abspath(__file__))
REQUEST_CSV = os.path.join(OUT, "request_selected.csv")
NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
NIGHT_END = Time("2026-06-12T14:48:00", format="isot")
RUNTIME_S = int(os.environ.get("TTP_RUNTIME_S", "600"))
MIP_GAP = 0.005
AZ_CUT = 235.0

MODES = {
    1: {"label": "single-state", "outdir": "out_single", "log": "gurobi_single_state.log"},
    2: {"label": "two-state", "outdir": "out_two", "log": "gurobi_two_state.log"},
}


def build_requests_access(df, queue, night_start, night_end, *, grid_min=2.0):
    """QTable whose availability is the real observability window.

    For each target we sample the whole night at ``grid_min`` cadence and run
    ``queue.is_accessible`` (alt in [18,85], outside the Nasmyth deck zone) --
    the same access computation the production pipeline uses. ``first/last
    available`` are the first/last accessible sample. Targets never accessible
    keep the full night so the solver can simply drop them.
    """
    coords = SkyCoord(df.ra.values * u.deg, df.dec.values * u.deg, frame="icrs")
    visit_min = queue.visit_duration(df.exptime.to_numpy(), df.n_exp.to_numpy())

    n = max(int((night_end.jd - night_start.jd) * 24 * 60 / grid_min), 10)
    t = Time(np.linspace(night_start.jd, night_end.jd, n), format="jd")
    aa = queue.observatory.altaz(t, coords, grid_times_targets=True)
    ok = queue.is_accessible(aa.alt.deg, aa.az.deg)  # (ntargets, ntimes)

    jd = t.jd
    first_jd = np.where(ok, jd[None, :], np.inf).min(axis=1)
    last_jd = np.where(ok, jd[None, :], -np.inf).max(axis=1)
    no_good = ~ok.any(axis=1)
    first_jd[no_good] = night_start.jd
    last_jd[no_good] = night_end.jd

    avail_min = (last_jd - first_jd) * 24 * 60
    print(f"  observability windows (min): "
          f"min={avail_min.min():.0f} median={np.median(avail_min):.0f} "
          f"max={avail_min.max():.0f}; never-up={int(no_good.sum())}")

    return QTable(
        {
            "unique_id": df.unique_id.to_numpy(dtype=object),
            "coord": coords,
            "first_available": Time(first_jd, format="jd"),
            "last_available": Time(last_jd, format="jd"),
            "t_visit": np.asarray(visit_min, dtype=float) * u.min,
            "n_intra_max": df.n_intra_max.to_numpy(dtype=int),
            "tau_intra": df.tau_intra.to_numpy(dtype=float) * u.hr,
            "priority": np.full(len(df), 10.0),
        },
        copy=False,
    )


def count_az_cut_crossings(schedule, queue, night_start, az_cut=AZ_CUT):
    """Count consecutive scheduled pairs straddling ``az_cut`` (legacy-expensive)."""
    sc = schedule[schedule["scheduled"]].sort_values("order")
    if len(sc) < 2:
        return 0
    t_obs = night_start + TimeDelta(sc["t_start"].to_numpy() * 60.0, format="sec")
    coords = SkyCoord(sc.ra.values * u.deg, sc.dec.values * u.deg, frame="icrs")
    az = np.atleast_1d(queue.observatory.altaz(t_obs, coords).az.deg)
    side = az > az_cut
    return int(np.sum(side[1:] != side[:-1]))


def run(queue, requests, df, n_states):
    cfg = MODES[n_states]
    outdir = os.path.join(OUT, cfg["outdir"])
    logfile = os.path.join(OUT, cfg["log"])
    os.makedirs(outdir, exist_ok=True)

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
    if os.path.exists(logfile):
        os.remove(logfile)
    tm.model.params.LogFile = logfile
    tm.model.params.LogToConsole = 0
    tm.model.params.TimeLimit = RUNTIME_S
    tm.model.params.MIPGap = MIP_GAP
    tm.model.update()
    tm.run_model()
    tm.build_schedule()

    tm.schedule.to_csv(os.path.join(outdir, "schedule.csv"), index=False)
    tm.observer = queue.observatory
    tm.wrap_limit = queue.wrap_limit
    tm.wrap_states = queue.wrap_states if n_states > 1 else None
    tplot.plot_path_2D_interactive(tm, night_start_time=NIGHT_START).write_html(
        os.path.join(outdir, "slew_path.html")
    )
    tplot.get_slew_animation_plotly(
        tm, REQUEST_CSV, inaccessible_zones=queue.inaccessible_zones
    ).write_html(os.path.join(outdir, "slew_animation.html"))
    hover = [c for c in ("target", "exptime", "n_exp") if c not in tm.schedule.columns]
    if hover:
        tm.schedule = tm.schedule.merge(df[["unique_id", *hover]], on="unique_id", how="left")
    aqplot.get_ladder(tm, NIGHT_START).write_html(os.path.join(outdir, "ladder.html"))

    sc = tm.schedule[tm.schedule["scheduled"]]
    set_ids = set(df.loc[df.timing_class.str.startswith("set"), "unique_id"])
    n_set_total = len(set_ids)
    n_set_grabbed = int(sc["unique_id"].isin(set_ids).sum())
    return {
        "label": cfg["label"],
        "n_states": n_states,
        "scheduled": int(len(sc)),
        "n_requested": int(tm.stats["n_requested"]),
        "early_set_grabbed": n_set_grabbed,
        "early_set_total": n_set_total,
        "modeled_slew_min": round(float(tm.stats["t_slew_sum"]), 2),
        "physical_slew_min": round(float(physical_slew_minutes(tm.schedule, queue, NIGHT_START)), 2),
        "az235_crossings": count_az_cut_crossings(tm.schedule, queue, NIGHT_START),
        "gurobi_status": int(tm.model.Status),
        "mip_gap": round(float(tm.model.MIPGap), 5),
        "solve_runtime_s": round(float(tm.model.Runtime), 1),
        "objective": round(float(tm.model.ObjVal), 3),
        "objective_bound": round(float(tm.model.ObjBound), 3),
        "outdir": cfg["outdir"],
        "log_file": cfg["log"],
    }


def write_summary(results):
    by_label = {r["label"]: r for r in results}
    meta = {
        "request_csv": "request_selected.csv",
        "night_start": NIGHT_START.isot,
        "night_end": NIGHT_END.isot,
        "night_minutes": round((NIGHT_END - NIGHT_START).to_value("min"), 1),
        "time_limit_s": RUNTIME_S,
        "mip_gap_target": MIP_GAP,
        "results": [by_label[l] for l in ("single-state", "two-state") if l in by_label],
    }
    with open(os.path.join(OUT, "comparison.json"), "w") as f:
        json.dump(meta, f, indent=2)

    r0 = meta["results"][0]
    lines = [
        "Staggered adversarial TTP comparison (real observability, rise/set staggered)",
        f"  request_csv : request_selected.csv ({r0['n_requested']} targets; "
        f"{r0['early_set_total']} early-setting)",
        f"  night       : {NIGHT_START.isot} -> {NIGHT_END.isot} ({meta['night_minutes']:.0f} min)",
        f"  time limit  : {RUNTIME_S}s   MIPGap target: {MIP_GAP}",
        "",
        f"{'model':<14}{'sched':>7}{'set grab':>9}{'phys slew':>11}{'mod slew':>10}"
        f"{'az235x':>8}{'gap':>9}{'solve s':>9}",
        "-" * 77,
    ]
    for r in meta["results"]:
        lines.append(
            f"{r['label']:<14}{r['scheduled']:>4}/{r['n_requested']:<2}"
            f"{r['early_set_grabbed']:>5}/{r['early_set_total']:<3}"
            f"{r['physical_slew_min']:>10.2f}m{r['modeled_slew_min']:>9.2f}m"
            f"{r['az235_crossings']:>8}{r['mip_gap']*100:>8.2f}%{r['solve_runtime_s']:>9.1f}"
        )
    lines += [
        "",
        "Gurobi logs: " + ", ".join(r["log_file"] for r in meta["results"]),
        "'set grab' = early-setting targets scheduled; 'phys slew' re-evaluates each "
        "order under the encoder wrap model; 'az235x' counts az-235 straddling transitions.",
    ]
    txt = "\n".join(lines)
    with open(os.path.join(OUT, "comparison.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)


def main():
    states = [int(a) for a in sys.argv[1:]] or [1, 2]
    existing = {}
    cj = os.path.join(OUT, "comparison.json")
    if os.path.exists(cj):
        with open(cj) as f:
            for r in json.load(f).get("results", []):
                existing[r["label"]] = r

    queue = HIRESCPS()
    df = pd.read_csv(REQUEST_CSV)
    requests = build_requests_access(df, queue, NIGHT_START, NIGHT_END)

    for s in states:
        r = run(queue, requests, df, s)
        existing[r["label"]] = r

    write_summary(list(existing.values()))


if __name__ == "__main__":
    main()
