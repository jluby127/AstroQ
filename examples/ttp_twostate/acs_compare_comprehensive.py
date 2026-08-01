"""Comprehensive ACS/MILP benchmark matrix for sphere100 and full-band1.

Runs scenarios per target group (see TESTS below), writes figures +
``result.txt`` under ``<group>_comprehensive/``.

Usage (astroq-testing env, repo root):

    ACS_GROUP=sphere100 python examples/ttp_twostate/acs_compare_comprehensive.py
    ACS_GROUP=full-band1 python examples/ttp_twostate/acs_compare_comprehensive.py

Run a subset (merges into existing ``results.json`` / ``result.txt``)::

    ACS_GROUP=sphere100 ACS_SUBSET=out_single_acs8_120,out_single_acs8_570 \\
        python examples/ttp_twostate/acs_compare_comprehensive.py
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
DEFAULT_NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
DEFAULT_NIGHT_END = Time("2026-06-12T14:48:00", format="isot")
MIP_GAP = 0.005
ACS_S = 30.0

GROUPS = {
    "full-band1": ("full_band1_2026-06-12", "demo"),
    "sphere100": ("sphere100_2026-02-01", "access"),
}

# (label, outdir, n_states, acs_starts, acs_nb_max, milp_s, acs_only_polish_s)
# acs_only_polish_s: if set, use ACS + this short polish (no full MILP row stats)
# milp_s: Gurobi TimeLimit after optional ACS seed
TESTS = [
    ("two plain", "out_two_plain", 2, 0, 25, 120, None),
    ("two acs_seed", "out_two_acs_seed", 2, 3, 25, 120, None),
    ("two acs_only", "out_two_acs_only", 2, 3, 25, None, 2),
    ("single acs_only", "out_single_acs_only", 1, 3, 25, None, 2),
    ("two acs8_120", "out_two_acs8_120", 2, 8, 25, 120, None),
    ("two acs8_570", "out_two_acs8_570", 2, 8, 25, 570, None),
    ("two acs8_120_nb100", "out_two_acs8_120_nb100", 2, 8, 100, 120, None),
    ("single acs_seed", "out_single_acs_seed", 1, 3, 25, 120, None),
    ("single acs8_120", "out_single_acs8_120", 1, 8, 25, 120, None),
    ("single acs8_570", "out_single_acs8_570", 1, 8, 25, 570, None),
]


def _resolve_group():
    group = os.environ.get("ACS_GROUP", "sphere100")
    if group not in GROUPS:
        raise ValueError(f"ACS_GROUP must be one of {list(GROUPS)}")
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
    out_dir = os.path.join(ROOT, f"{group}_comprehensive")
    os.makedirs(out_dir, exist_ok=True)
    return group, request_csv, builder, out_dir


def _night_bounds(request_csv):
    meta_path = os.path.join(os.path.dirname(request_csv), "meta.json")
    if os.path.isfile(meta_path):
        with open(meta_path, encoding="utf-8") as fh:
            meta = json.load(fh)
        return (
            Time(meta["night_start"], format="isot"),
            Time(meta["night_end"], format="isot"),
        )
    return DEFAULT_NIGHT_START, DEFAULT_NIGHT_END


def _new_model(queue, requests, night_start, night_end, n_states):
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
    return tm


def _render(tm, queue, df, n_states, request_csv, out_dir, sub, night_start):
    d = os.path.join(out_dir, sub)
    os.makedirs(d, exist_ok=True)
    tm.schedule.to_csv(os.path.join(d, "schedule.csv"), index=False)
    tm.observer = queue.observatory
    tm.wrap_limit = queue.wrap_limit
    tm.wrap_states = queue.wrap_states if n_states > 1 else None
    tplot.plot_path_2D_interactive(tm, night_start_time=night_start).write_html(
        os.path.join(d, "slew_path.html")
    )
    anim = tplot.get_slew_animation_plotly(
        tm, request_csv, inaccessible_zones=queue.inaccessible_zones
    )
    tplot.save_slew_animation(anim, os.path.join(d, "slew_animation.html"))
    hover = [c for c in ("target", "exptime", "n_exp") if c not in tm.schedule.columns]
    if hover:
        tm.schedule = tm.schedule.merge(
            df[["unique_id", *hover]], on="unique_id", how="left"
        )
    aqplot.get_ladder(tm, night_start).write_html(os.path.join(d, "ladder.html"))


def _milp_row(tm, queue, label, *, acs_wall=None, seed_obj=None):
    sc = tm.schedule[tm.schedule["scheduled"]]
    row = {
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
        "acs_wall_s": round(acs_wall, 1) if acs_wall is not None else float("nan"),
        "seed_obj": round(seed_obj, 2) if seed_obj is not None else float("nan"),
    }
    return row


def _acs_only_row(tm, queue, label, *, seed_obj, acs_wall):
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
        "obj": round(float(tm.model.ObjVal), 2) if tm.model.SolCount else round(seed_obj, 2),
        "bound": float("nan"),
        "nodes": int(tm.model.NodeCount) if tm.model.SolCount else 0,
        "solve_s": round(float(tm.model.Runtime), 1),
        "acs_wall_s": round(acs_wall, 1) if acs_wall is not None else float("nan"),
        "seed_obj": round(seed_obj, 2),
    }


def _run_test(
    queue,
    requests,
    df,
    night_start,
    night_end,
    request_csv,
    out_dir,
    label,
    sub,
    n_states,
    acs_starts,
    acs_nb_max,
    milp_s,
    polish_s,
):
    global NIGHT_START
    NIGHT_START = night_start

    tm = _new_model(queue, requests, night_start, night_end, n_states)
    acs_wall = None
    seed_obj = None
    acs_kw = dict(
        params={"time_limit_s": ACS_S, "nb_max": acs_nb_max},
        n_starts=acs_starts,
        parallel=True,
    )

    if acs_starts > 0:
        tm.build_model()
        t0 = time.time()
        seed = acs_warm_start(tm, **acs_kw)
        acs_wall = time.time() - t0
        seed_obj = seed["objective"]
        tm.seed_from_tour(seed)
    else:
        tm.build_model()

    if polish_s is not None:
        tm.model.params.OutputFlag = 0
        tm.model.params.TimeLimit = polish_s
        tm.model.params.MIPGap = 0.5
        tm.model.update()
        tm.run_model()
        tm.build_schedule()
        row = _acs_only_row(
            tm, queue, label, seed_obj=seed_obj, acs_wall=acs_wall
        )
    else:
        tm.model.params.OutputFlag = 0
        tm.model.params.TimeLimit = milp_s
        tm.model.params.MIPGap = MIP_GAP
        tm.model.params.PreSolve = 2
        tm.model.params.MIPFocus = 1
        tm.model.update()
        tm.run_model()
        tm.build_schedule()
        disp = label
        if acs_starts > 0:
            disp = f"{label} (acs {acs_wall:.0f}s, seed {seed_obj:.1f})"
        row = _milp_row(tm, queue, disp, acs_wall=acs_wall, seed_obj=seed_obj)

    _render(tm, queue, df, n_states, request_csv, out_dir, sub, night_start)
    print(f"  done {label}: {row['scheduled']}/{row['n_requested']} sched", flush=True)
    row["outdir"] = sub
    return row


def _load_existing_rows(out_dir):
    path = os.path.join(out_dir, "results.json")
    if not os.path.isfile(path):
        return {}
    with open(path, encoding="utf-8") as fh:
        rows = json.load(fh)
    label_to_outdir = {spec[0]: spec[1] for spec in TESTS}
    merged = {}
    for r in rows:
        outdir = r.get("outdir")
        if not outdir:
            base = r["label"].split(" (acs")[0]
            outdir = label_to_outdir.get(base)
        if outdir:
            r["outdir"] = outdir
            merged[outdir] = r
    return merged


def _merge_rows(existing, new_rows):
    merged = dict(existing)
    for row in new_rows:
        merged[row["outdir"]] = row
    return merged


def _format_results(group, n_requested, merged_by_outdir):
    rows = []
    for spec in TESTS:
        sub = spec[1]
        if sub in merged_by_outdir:
            rows.append(merged_by_outdir[sub])
    lines = [
        f"Comprehensive ACS/MILP ({group}, {n_requested} targets)",
        f"  night      : {NIGHT_START.isot} -> {NIGHT_END.isot}",
        f"  ACS        : {ACS_S}s per start (parallel)   MIPGap: {MIP_GAP}",
        "",
        f"{'run':<46}{'sched':>7}{'phys':>9}{'mod':>9}{'gap':>8}"
        f"{'obj':>9}{'acs s':>7}{'milp s':>8}",
        "-" * 103,
    ]
    for r in rows:
        gap = "   nan" if r["gap"] != r["gap"] else f"{r['gap']:>6.2f}%"
        acs = "    nan" if r["acs_wall_s"] != r["acs_wall_s"] else f"{r['acs_wall_s']:>6.1f}"
        lines.append(
            f"{r['label']:<46}{r['scheduled']:>4}/{r['n_requested']:<2}"
            f"{r['physical_slew']:>8.2f}{r['modeled_slew']:>9.2f}{gap:>8}"
            f"{r['obj']:>9.2f}{acs}{r['solve_s']:>8.1f}"
        )
    return rows, "\n".join(lines)


def main():
    global NIGHT_START, NIGHT_END
    group, request_csv, builder, out_dir = _resolve_group()
    NIGHT_START, NIGHT_END = _night_bounds(request_csv)
    queue = HIRESCPS()
    df = pd.read_csv(request_csv)
    requests = builder(df, queue, NIGHT_START, NIGHT_END)

    print(
        f"Comprehensive ACS/MILP ({group}, {len(df)} targets)\n"
        f"  night : {NIGHT_START.isot} -> {NIGHT_END.isot}\n"
        f"  out   : {out_dir}\n",
        flush=True,
    )

    subset = os.environ.get("ACS_SUBSET", "").strip()
    subset_outdirs = {s.strip() for s in subset.split(",") if s.strip()} if subset else None

    existing = _load_existing_rows(out_dir)
    new_rows = []
    for spec in TESTS:
        label, sub, n_states, acs_starts, nb_max, milp_s, polish_s = spec
        if subset_outdirs is not None and sub not in subset_outdirs:
            continue
        print(f"Running {label}...", flush=True)
        new_rows.append(
            _run_test(
                queue,
                requests,
                df,
                NIGHT_START,
                NIGHT_END,
                request_csv,
                out_dir,
                label,
                sub,
                n_states,
                acs_starts,
                nb_max,
                milp_s,
                polish_s,
            )
        )

    merged = _merge_rows(existing, new_rows)
    n_requested = new_rows[0]["n_requested"] if new_rows else next(iter(merged.values()))["n_requested"]
    rows, txt = _format_results(group, n_requested, merged)
    with open(os.path.join(out_dir, "result.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    with open(os.path.join(out_dir, "results.json"), "w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=2)
        fh.write("\n")
    print(txt, flush=True)


if __name__ == "__main__":
    main()
