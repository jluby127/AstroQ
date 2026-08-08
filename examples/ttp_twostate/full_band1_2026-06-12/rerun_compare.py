"""Rerun the single- vs two-state TTP comparison on the static request_selected.

For each solved mode this persists everything needed to inspect/re-render
without re-solving:

- ``gurobi_<mode>.log``         -- full Gurobi solver log (in this folder)
- ``<outdir>/schedule.csv``     -- the solved schedule
- ``<outdir>/slew_path.html``   -- az/alt path plot
- ``<outdir>/slew_animation.html`` -- animated polar slew plot

and updates ``comparison.json`` / ``comparison.txt`` (merging per-mode rows so
running a single mode does not clobber the other).

Usage (from repo root, astroq-testing env):

    python examples/ttp_twostate/full_band1_2026-06-12/rerun_compare.py        # both modes
    python examples/ttp_twostate/full_band1_2026-06-12/rerun_compare.py 2      # two-state only
    python examples/ttp_twostate/full_band1_2026-06-12/rerun_compare.py 1      # single-state only
"""

import json
import os
import sys

import pandas as pd
from astropy.time import Time

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel
import astroq.ttp.plot as tplot
import astroq.plot as aqplot
from astroq.scripts.demo_twostate import build_requests, physical_slew_minutes

BASE = os.path.dirname(os.path.abspath(__file__))
# Output directory and solve budget are overridable so the same script can drive
# multiple comparison runs (e.g. a longer 3600 s budget into a sibling folder).
OUT = os.environ.get("TTP_OUTDIR", BASE)
os.makedirs(OUT, exist_ok=True)
# request_selected.csv lives in OUT if present, else fall back to the base folder.
REQUEST_CSV = os.path.join(OUT, "request_selected.csv")
if not os.path.exists(REQUEST_CSV):
    REQUEST_CSV = os.path.join(BASE, "request_selected.csv")
NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
NIGHT_END = Time("2026-06-12T14:48:00", format="isot")
RUNTIME_S = int(os.environ.get("TTP_RUNTIME_S", "600"))
MIP_GAP = 0.005

MODES = {
    1: {"label": "single-state", "outdir": "out_single", "log": "gurobi_single_state.log"},
    2: {"label": "two-state", "outdir": "out_two", "log": "gurobi_two_state.log"},
}


def run(queue, requests, n_states):
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

    # Persist the schedule and re-render the plots straight off this solve.
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
    # Ladder needs exptime/n_exp/target for hover; merge them onto the schedule.
    req = pd.read_csv(REQUEST_CSV)
    hover_cols = [c for c in ("target", "exptime", "n_exp") if c not in tm.schedule.columns]
    if hover_cols:
        tm.schedule = tm.schedule.merge(
            req[["unique_id", *hover_cols]], on="unique_id", how="left"
        )
    aqplot.get_ladder(tm, NIGHT_START).write_html(os.path.join(outdir, "ladder.html"))

    sc = tm.schedule[tm.schedule["scheduled"]]
    return {
        "label": cfg["label"],
        "n_states": n_states,
        "scheduled": int(len(sc)),
        "n_requested": int(tm.stats["n_requested"]),
        "modeled_slew_min": round(float(tm.stats["t_slew_sum"]), 2),
        "physical_slew_min": round(
            float(physical_slew_minutes(tm.schedule, queue, NIGHT_START)), 2
        ),
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
        "results": [by_label[lbl] for lbl in ("single-state", "two-state") if lbl in by_label],
    }
    with open(os.path.join(OUT, "comparison.json"), "w") as f:
        json.dump(meta, f, indent=2)

    lines = [
        "TTP single- vs two-state comparison",
        f"  request_csv : request_selected.csv ({meta['results'][0]['n_requested']} targets)",
        f"  night       : {NIGHT_START.isot} -> {NIGHT_END.isot} "
        f"({meta['night_minutes']:.0f} min)",
        f"  time limit  : {RUNTIME_S}s   MIPGap target: {MIP_GAP}",
        "",
        f"{'model':<14}{'sched':>7}{'phys slew':>11}{'mod slew':>10}{'gap':>9}{'solve s':>9}",
        "-" * 60,
    ]
    for r in meta["results"]:
        lines.append(
            f"{r['label']:<14}{r['scheduled']:>4}/{r['n_requested']:<2}"
            f"{r['physical_slew_min']:>10.2f}m{r['modeled_slew_min']:>9.2f}m"
            f"{r['mip_gap']*100:>8.2f}%{r['solve_runtime_s']:>9.1f}"
        )
    lines += [
        "",
        "Gurobi logs: " + ", ".join(r["log_file"] for r in meta["results"]),
        "'phys slew' re-evaluates each order under the encoder wrap model "
        "(the fair, common metric).",
    ]
    txt = "\n".join(lines)
    with open(os.path.join(OUT, "comparison.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)


def main():
    states = [int(a) for a in sys.argv[1:]] or [1, 2]

    # Start from any existing results so a single-mode run keeps the other row.
    existing = {}
    cj = os.path.join(OUT, "comparison.json")
    if os.path.exists(cj):
        with open(cj) as f:
            for r in json.load(f).get("results", []):
                existing[r["label"]] = r

    queue = HIRESCPS()
    df = pd.read_csv(REQUEST_CSV)
    requests = build_requests(df, queue, NIGHT_START, NIGHT_END)

    for s in states:
        r = run(queue, requests, s)
        existing[r["label"]] = r

    write_summary(list(existing.values()))


if __name__ == "__main__":
    main()
