"""ACS heuristic vs MILP on the v2 adversarial night (58 targets).

Compares, under a fixed MILP time budget:
  - two-state plain          : run_model()
  - two-state acs_seed       : seed_from_tour(ACS) then run_model()
  - two-state acs_only       : run_heuristic() (no Gurobi)
  - single-state acs_only    : run_heuristic() (no Gurobi)

Captures scheduled count, modeled + physical slew, gap/bound/nodes, and ACS
wall time, into result.txt.

Usage (astroq-testing env, repo root):
    TTP_RUNTIME_S=120 python examples/ttp_twostate/adversarial-2_2026-06-12_acs/acs_compare.py
"""

import importlib.util
import os
import sys
import time

import pandas as pd
from astropy.time import Time

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel
from astroq.scripts.demo_twostate import physical_slew_minutes
import astroq.ttp.plot as tplot
import astroq.plot as aqplot

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.join(os.path.dirname(HERE), "adversarial-2_2026-06-12_runtime-600")
REQUEST_CSV = os.path.join(V2, "request_selected.csv")
NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
NIGHT_END = Time("2026-06-12T14:48:00", format="isot")
RUNTIME_S = int(os.environ.get("TTP_RUNTIME_S", "120"))
ACS_S = float(os.environ.get("ACS_S", "10"))
MIP_GAP = 0.005

spec = importlib.util.spec_from_file_location(
    "advc", os.path.join(V2, "adversarial_compare.py")
)
advc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(advc)


def _render(tm, queue, df, n_states, outdir):
    """Persist schedule.csv and render slew-path / animation / ladder plots."""
    outdir = os.path.join(HERE, outdir)
    os.makedirs(outdir, exist_ok=True)
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
        tm.schedule = tm.schedule.merge(
            df[["unique_id", *hover]], on="unique_id", how="left"
        )
    aqplot.get_ladder(tm, NIGHT_START).write_html(os.path.join(outdir, "ladder.html"))


def _new_model(queue, requests, n_states):
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
    return tm


def _milp_row(tm, queue, label, n_states):
    sc = tm.schedule[tm.schedule["scheduled"]]
    return {
        "label": label,
        "scheduled": int(len(sc)),
        "n_requested": int(tm.stats["n_requested"]),
        "modeled_slew": round(float(tm.stats["t_slew_sum"]), 2),
        "physical_slew": round(
            float(physical_slew_minutes(tm.schedule, queue, NIGHT_START)), 2
        ),
        "gap": round(float(tm.model.MIPGap) * 100, 2),
        "obj": round(float(tm.model.ObjVal), 2),
        "bound": round(float(tm.model.ObjBound), 2),
        "nodes": int(tm.model.NodeCount),
        "solve_s": round(float(tm.model.Runtime), 1),
    }


def _acs_row(tm, queue, label):
    sc = tm.schedule[tm.schedule["scheduled"]]
    return {
        "label": label,
        "scheduled": int(len(sc)),
        "n_requested": int(tm.stats["n_requested"]),
        "modeled_slew": round(float(tm.stats["t_slew_sum"]), 2),
        "physical_slew": round(
            float(physical_slew_minutes(tm.schedule, queue, NIGHT_START)), 2
        ),
        "gap": float("nan"),
        "obj": round(float(tm.acs_result["objective"]), 2),
        "bound": float("nan"),
        "nodes": 0,
        "solve_s": 0.0,
    }


def main():
    queue = HIRESCPS()
    df = pd.read_csv(REQUEST_CSV)
    requests = advc.build_requests_access(df, queue, NIGHT_START, NIGHT_END)
    rows = []

    # two-state plain
    tm = _new_model(queue, requests, 2)
    tm.build_model()
    tm.model.params.OutputFlag = 0
    tm.model.params.TimeLimit = RUNTIME_S
    tm.model.params.MIPGap = MIP_GAP
    tm.model.params.PreSolve = 2
    tm.model.params.MIPFocus = 1
    tm.model.update()
    tm.run_model()
    tm.build_schedule()
    rows.append(_milp_row(tm, queue, "two plain", 2))
    _render(tm, queue, df, 2, "out_two_plain")

    # two-state acs_seed
    tm = _new_model(queue, requests, 2)
    tm.build_model()
    t0 = time.time()
    seed = tm.run_heuristic(params={"time_limit_s": ACS_S})
    acs_wall = time.time() - t0
    tm.seed_from_tour(seed)
    tm.model.params.OutputFlag = 0
    tm.model.params.TimeLimit = RUNTIME_S
    tm.model.params.MIPGap = MIP_GAP
    tm.model.params.PreSolve = 2
    tm.model.params.MIPFocus = 1
    tm.model.update()
    tm.run_model()
    tm.build_schedule()
    r = _milp_row(tm, queue, "two acs_seed", 2)
    r["label"] = f"two acs_seed (acs {acs_wall:.0f}s, seed obj {seed['objective']:.1f})"
    rows.append(r)
    _render(tm, queue, df, 2, "out_two_acs_seed")

    # two-state acs only
    tm = _new_model(queue, requests, 2)
    tm.run_heuristic(params={"time_limit_s": ACS_S})
    rows.append(_acs_row(tm, queue, "two acs_only"))
    _render(tm, queue, df, 2, "out_two_acs_only")

    # single-state acs only
    tm = _new_model(queue, requests, 1)
    tm.run_heuristic(params={"time_limit_s": ACS_S})
    rows.append(_acs_row(tm, queue, "single acs_only"))
    _render(tm, queue, df, 1, "out_single_acs_only")

    lines = [
        "ACS vs MILP (v2 adversarial, 58 targets)",
        f"  night      : {NIGHT_START.isot} -> {NIGHT_END.isot}",
        f"  MILP limit : {RUNTIME_S}s   ACS limit: {ACS_S}s   MIPGap: {MIP_GAP}",
        "",
        f"{'run':<42}{'sched':>7}{'phys':>9}{'mod':>9}{'gap':>8}"
        f"{'obj':>9}{'bound':>9}{'nodes':>9}{'solve s':>9}",
        "-" * 109,
    ]
    for r in rows:
        gap = "   nan" if r["gap"] != r["gap"] else f"{r['gap']:>6.2f}%"
        bound = "     nan" if r["bound"] != r["bound"] else f"{r['bound']:>9.2f}"
        lines.append(
            f"{r['label']:<42}{r['scheduled']:>4}/{r['n_requested']:<2}"
            f"{r['physical_slew']:>8.2f}{r['modeled_slew']:>9.2f}{gap:>8}"
            f"{r['obj']:>9.2f}{bound}{r['nodes']:>9}{r['solve_s']:>9.1f}"
        )
    txt = "\n".join(lines)
    with open(os.path.join(HERE, "result.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
