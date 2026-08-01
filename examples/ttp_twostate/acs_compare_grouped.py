"""ACS heuristic vs MILP across target groups, with rendered figures.

Parametrized version of ``adversarial-2_*/acs_compare.py``. Selects a target
group via the ``ACS_GROUP`` env var and writes results + slew animations to
``<group>_acs/``.

Groups:
  adversarial-2 | adversarial-4  -> adversarial_compare.build_requests_access
  full-band1                     -> demo_twostate.build_requests
  sphere100                      -> requests_access (night from meta.json)

Each group runs four solves under a fixed MILP budget:
  two plain | two acs_seed | two acs_only | single acs_only

Usage (astroq-testing env, repo root):
    ACS_GROUP=adversarial-4 TTP_RUNTIME_S=120 ACS_S=8 \
        python examples/ttp_twostate/acs_compare_grouped.py
"""

import importlib.util
import json
import os
import time

import pandas as pd
from astropy.time import Time

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel
from astroq.ttp.acs import acs_warm_start
from astroq.scripts.demo_twostate import physical_slew_minutes, build_requests
import astroq.ttp.plot as tplot
import astroq.plot as aqplot

ROOT = os.path.dirname(os.path.abspath(__file__))
NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
NIGHT_END = Time("2026-06-12T14:48:00", format="isot")
RUNTIME_S = int(os.environ.get("TTP_RUNTIME_S", "120"))
ACS_S = float(os.environ.get("ACS_S", "30"))
ACS_STARTS = int(os.environ.get("ACS_STARTS", "3"))
MIP_GAP = 0.005

GROUPS = {
    "adversarial-2": ("adversarial-2_2026-06-12_runtime-600", "access"),
    "adversarial-4": ("adversarial-4_2026-06-12_runtime-600", "access"),
    "full-band1": ("full_band1_2026-06-12", "demo"),
    "sphere100": ("sphere100_2026-02-01", "access"),
}


def _resolve_group():
    group = os.environ.get("ACS_GROUP", "adversarial-4")
    src_name, builder_kind = GROUPS[group]
    src = os.path.join(ROOT, src_name)
    request_csv = os.path.join(src, "request_selected.csv")
    if builder_kind == "access":
        for modname in ("adversarial_compare.py", "requests_access.py"):
            path = os.path.join(src, modname)
            if os.path.isfile(path):
                spec = importlib.util.spec_from_file_location("access_mod", path)
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                builder = mod.build_requests_access
                break
        else:
            raise FileNotFoundError(f"no access builder in {src}")
    else:
        builder = build_requests
    out_dir = os.path.join(ROOT, f"{group}_acs")
    os.makedirs(out_dir, exist_ok=True)
    return group, request_csv, builder, out_dir


def _night_bounds(request_csv):
    """Use ``meta.json`` night times when present (e.g. sphere100 benchmark)."""
    meta_path = os.path.join(os.path.dirname(request_csv), "meta.json")
    if os.path.isfile(meta_path):
        with open(meta_path, encoding="utf-8") as fh:
            meta = json.load(fh)
        return (
            Time(meta["night_start"], format="isot"),
            Time(meta["night_end"], format="isot"),
        )
    return NIGHT_START, NIGHT_END


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


def _render(tm, queue, df, n_states, request_csv, out_dir, sub):
    d = os.path.join(out_dir, sub)
    os.makedirs(d, exist_ok=True)
    tm.schedule.to_csv(os.path.join(d, "schedule.csv"), index=False)
    tm.observer = queue.observatory
    tm.wrap_limit = queue.wrap_limit
    tm.wrap_states = queue.wrap_states if n_states > 1 else None
    tplot.plot_path_2D_interactive(tm, night_start_time=NIGHT_START).write_html(
        os.path.join(d, "slew_path.html")
    )
    tplot.get_slew_animation_plotly(
        tm, request_csv, inaccessible_zones=queue.inaccessible_zones
    ).write_html(os.path.join(d, "slew_animation.html"))
    hover = [c for c in ("target", "exptime", "n_exp") if c not in tm.schedule.columns]
    if hover:
        tm.schedule = tm.schedule.merge(
            df[["unique_id", *hover]], on="unique_id", how="left"
        )
    aqplot.get_ladder(tm, NIGHT_START).write_html(os.path.join(d, "ladder.html"))


def _milp_row(tm, queue, label):
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
        "obj": round(float(tm._seed_obj), 2),
        "bound": float("nan"),
        "nodes": 0,
        "solve_s": 0.0,
    }


def _set_params(tm):
    tm.model.params.OutputFlag = 0
    tm.model.params.TimeLimit = RUNTIME_S
    tm.model.params.MIPGap = MIP_GAP
    tm.model.params.PreSolve = 2
    tm.model.params.MIPFocus = 1
    tm.model.update()


def _acs_only(tm, acs_kw):
    """ACS warm-start + a 2 s Gurobi polish (no real MILP budget).

    There is no heuristic-only code path anymore; an "ACS-only" schedule is
    just the warm-start seed lightly cleaned up by a very short solve.
    """
    tm.build_model()
    seed = acs_warm_start(tm, **acs_kw)
    tm.seed_from_tour(seed)
    tm.model.params.OutputFlag = 0
    tm.model.params.TimeLimit = 2
    tm.model.params.MIPGap = 0.5
    tm.model.update()
    tm.run_model()
    tm.build_schedule()
    tm._seed_obj = seed["objective"]
    return tm


def main():
    global NIGHT_START, NIGHT_END
    group, request_csv, builder, out_dir = _resolve_group()
    NIGHT_START, NIGHT_END = _night_bounds(request_csv)
    queue = HIRESCPS()
    df = pd.read_csv(request_csv)
    requests = builder(df, queue, NIGHT_START, NIGHT_END)
    rows = []

    # two-state plain
    tm = _new_model(queue, requests, 2)
    tm.build_model()
    _set_params(tm)
    tm.run_model()
    tm.build_schedule()
    rows.append(_milp_row(tm, queue, "two plain"))
    _render(tm, queue, df, 2, request_csv, out_dir, "out_two_plain")

    acs_kw = dict(params={"time_limit_s": ACS_S}, n_starts=ACS_STARTS)

    # two-state acs_seed
    tm = _new_model(queue, requests, 2)
    tm.build_model()
    t0 = time.time()
    seed = acs_warm_start(tm, **acs_kw)
    acs_wall = time.time() - t0
    tm.seed_from_tour(seed)
    _set_params(tm)
    tm.run_model()
    tm.build_schedule()
    r = _milp_row(tm, queue, "two acs_seed")
    r["label"] = f"two acs_seed (acs {acs_wall:.0f}s, seed obj {seed['objective']:.1f})"
    rows.append(r)
    _render(tm, queue, df, 2, request_csv, out_dir, "out_two_acs_seed")

    # two-state acs only (warm-start + 2s polish)
    tm = _new_model(queue, requests, 2)
    _acs_only(tm, acs_kw)
    rows.append(_acs_row(tm, queue, "two acs_only"))
    _render(tm, queue, df, 2, request_csv, out_dir, "out_two_acs_only")

    # single-state acs only (warm-start + 2s polish)
    tm = _new_model(queue, requests, 1)
    _acs_only(tm, acs_kw)
    rows.append(_acs_row(tm, queue, "single acs_only"))
    _render(tm, queue, df, 1, request_csv, out_dir, "out_single_acs_only")

    lines = [
        f"ACS vs MILP ({group}, {rows[0]['n_requested']} targets)",
        f"  night      : {NIGHT_START.isot} -> {NIGHT_END.isot}",
        f"  MILP limit : {RUNTIME_S}s   ACS: {ACS_STARTS}x{ACS_S}s (parallel)"
        f"   MIPGap: {MIP_GAP}",
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
    with open(os.path.join(out_dir, "result.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
