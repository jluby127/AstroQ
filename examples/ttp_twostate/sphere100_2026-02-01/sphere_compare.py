"""Single- vs two-state TTP on the uniform-sphere 100-target benchmark.

Targets are randomly distributed on the sky, each observable for >= 1 hour on
2026-02-01 (Keck nautical twilight). Exposure times are tuned so the night
packs exactly if every visit uses 0.5 min slew (see ``meta.json``).

Usage (astroq-testing env, repo root):

    TTP_RUNTIME_S=600 python examples/ttp_twostate/sphere100_2026-02-01/sphere_compare.py
    TTP_RUNTIME_S=600 python examples/ttp_twostate/sphere100_2026-02-01/sphere_compare.py 2  # two-state only
"""

import importlib.util
import json
import os
import sys

import pandas as pd
from astropy.time import Time

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel
import astroq.ttp.plot as tplot
import astroq.plot as aqplot
from astroq.scripts.demo_twostate import physical_slew_minutes

OUT = os.path.dirname(os.path.abspath(__file__))
REQUEST_CSV = os.path.join(OUT, "request_selected.csv")
META_JSON = os.path.join(OUT, "meta.json")
RUNTIME_S = int(os.environ.get("TTP_RUNTIME_S", "600"))
MIP_GAP = 0.005

MODES = {
    1: {"label": "single-state", "outdir": "out_single", "log": "gurobi_single_state.log"},
    2: {"label": "two-state", "outdir": "out_two", "log": "gurobi_two_state.log"},
}


def _load_build_requests_access():
    path = os.path.join(OUT, "requests_access.py")
    spec = importlib.util.spec_from_file_location("requests_access", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.build_requests_access


def load_meta():
    with open(META_JSON, encoding="utf-8") as fh:
        return json.load(fh)


def run(queue, requests, df, night_start, night_end, n_states):
    cfg = MODES[n_states]
    outdir = os.path.join(OUT, cfg["outdir"])
    logfile = os.path.join(OUT, cfg["log"])
    os.makedirs(outdir, exist_ok=True)

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
    tplot.plot_path_2D_interactive(tm, night_start_time=night_start).write_html(
        os.path.join(outdir, "slew_path.html")
    )
    tplot.get_slew_animation_plotly(
        tm, REQUEST_CSV, inaccessible_zones=queue.inaccessible_zones
    ).write_html(os.path.join(outdir, "slew_animation.html"))
    hover = [c for c in ("target", "exptime", "n_exp") if c not in tm.schedule.columns]
    if hover:
        tm.schedule = tm.schedule.merge(df[["unique_id", *hover]], on="unique_id", how="left")
    aqplot.get_ladder(tm, night_start).write_html(os.path.join(outdir, "ladder.html"))

    sc = tm.schedule[tm.schedule["scheduled"]]
    return {
        "label": cfg["label"],
        "n_states": n_states,
        "scheduled": int(len(sc)),
        "n_requested": int(tm.stats["n_requested"]),
        "modeled_slew_min": round(float(tm.stats["t_slew_sum"]), 2),
        "physical_slew_min": round(
            float(physical_slew_minutes(tm.schedule, queue, night_start)), 2
        ),
        "idle_min": round(float(tm.stats["t_idle_sum"]), 2),
        "gurobi_status": int(tm.model.Status),
        "mip_gap": round(float(tm.model.MIPGap), 5),
        "solve_runtime_s": round(float(tm.model.Runtime), 1),
        "objective": round(float(tm.model.ObjVal), 3),
        "objective_bound": round(float(tm.model.ObjBound), 3),
        "outdir": cfg["outdir"],
        "log_file": cfg["log"],
    }


def main():
    meta = load_meta()
    night_start = Time(meta["night_start"], format="isot")
    night_end = Time(meta["night_end"], format="isot")

    which = [1, 2] if len(sys.argv) <= 1 else [int(sys.argv[1])]
    queue = HIRESCPS()
    df = pd.read_csv(REQUEST_CSV)
    build_requests_access = _load_build_requests_access()
    requests = build_requests_access(df, queue, night_start, night_end)

    results = [
        run(queue, requests, df, night_start, night_end, s)
        for s in which
    ]

    summary = {
        "request_csv": "request_selected.csv",
        "meta_json": "meta.json",
        "night_start": night_start.isot,
        "night_end": night_end.isot,
        "night_minutes": meta["night_minutes"],
        "packing_identity": meta["packing_identity"],
        "time_limit_s": RUNTIME_S,
        "mip_gap_target": MIP_GAP,
        "results": results,
    }
    with open(os.path.join(OUT, "comparison.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    lines = [
        "Uniform-sphere 100-target TTP comparison (Keck 2026-02-01)",
        f"  {meta['packing_identity']}",
        f"  night : {night_start.isot} -> {night_end.isot} ({meta['night_minutes']:.0f} min)",
        f"  limit : {RUNTIME_S}s   MIPGap: {MIP_GAP}",
        "",
        f"{'model':<14}{'sched':>7}{'phys slew':>11}{'mod slew':>10}{'idle':>8}"
        f"{'gap':>9}{'solve s':>9}",
        "-" * 68,
    ]
    for r in results:
        lines.append(
            f"{r['label']:<14}{r['scheduled']:>4}/{r['n_requested']:<2}"
            f"{r['physical_slew_min']:>10.1f}{r['modeled_slew_min']:>10.1f}"
            f"{r['idle_min']:>8.1f}{r['mip_gap']:>8.3f}{r['solve_runtime_s']:>9.1f}"
        )
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "comparison.txt"), "w", encoding="utf-8") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
