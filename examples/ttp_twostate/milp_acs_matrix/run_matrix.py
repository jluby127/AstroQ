"""MILP/ACS dimensionality sweep — 8 datasets × 10 methods (80 runs).

All outputs live under this directory (``runs/``, ``all_results.csv``, ``summary/``).

Usage (astroq-testing env, repo root)::

    python examples/ttp_twostate/milp_acs_matrix/run_matrix.py
    MATRIX_DATASETS=sphere50,sphere100 MATRIX_METHODS=two_acs8_120 python ...
    MATRIX_SKIP_PLOTS=1 MATRIX_FORCE=1 python ...
"""

import importlib.util
import json
import os
import sys
import time

import pandas as pd
from astropy.time import Time

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel
from astroq.ttp.acs import acs_warm_start
from astroq.scripts.demo_twostate import physical_slew_minutes, build_requests
import astroq.ttp.plot as tplot
import astroq.plot as aqplot

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SPHERE_DS = os.path.join(ROOT, "sphere_datasets")
DEFAULT_NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
DEFAULT_NIGHT_END = Time("2026-06-12T14:48:00", format="isot")

MIP_GAP = 0.005
ACS_S = 30.0
ACS_STARTS = 8
ACS_NB_MAX = 25
NOREL_WARMSTART_S = 30.0  # match ops/HIRES/common/config_template.ini

# (method_key, n_states, acs_starts, milp_s, solve_method)
METHODS = [
    ("single_milp_120", 1, 0, 120, "milp"),
    ("single_milp_600", 1, 0, 600, "milp"),
    ("single_acs8_120", 1, ACS_STARTS, 120, "acs8+milp"),
    ("single_acs8_600", 1, ACS_STARTS, 600, "acs8+milp"),
    ("single_norel_120", 1, 0, 120, "norel+milp"),
    ("two_milp_120", 2, 0, 120, "milp"),
    ("two_milp_600", 2, 0, 600, "milp"),
    ("two_acs8_120", 2, ACS_STARTS, 120, "acs8+milp"),
    ("two_acs8_600", 2, ACS_STARTS, 600, "acs8+milp"),
    ("two_norel_120", 2, 0, 120, "norel+milp"),
]

DATASETS = {
    "full-band1": {
        "n_targets": 58,
        "request_dir": os.path.join(ROOT, "full_band1_2026-06-12"),
        "builder": "demo",
    },
    "sphere10": {"n_targets": 10, "request_dir": os.path.join(SPHERE_DS, "sphere10_2026-02-01"), "builder": "access"},
    "sphere20": {"n_targets": 20, "request_dir": os.path.join(SPHERE_DS, "sphere20_2026-02-01"), "builder": "access"},
    "sphere40": {"n_targets": 40, "request_dir": os.path.join(SPHERE_DS, "sphere40_2026-02-01"), "builder": "access"},
    "sphere50": {"n_targets": 50, "request_dir": os.path.join(SPHERE_DS, "sphere50_2026-02-01"), "builder": "access"},
    "sphere75": {"n_targets": 75, "request_dir": os.path.join(SPHERE_DS, "sphere75_2026-02-01"), "builder": "access"},
    "sphere100": {"n_targets": 100, "request_dir": os.path.join(SPHERE_DS, "sphere100_2026-02-01"), "builder": "access"},
    "sphere125": {"n_targets": 125, "request_dir": os.path.join(SPHERE_DS, "sphere125_2026-02-01"), "builder": "access"},
}


def _load_access_builder():
    path = os.path.join(SPHERE_DS, "requests_access.py")
    spec = importlib.util.spec_from_file_location("sphere_access", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.build_requests_access


def _night_bounds(request_dir):
    meta_path = os.path.join(request_dir, "meta.json")
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


def _acs_checkpoint(tm, seed, queue, night_start):
    """Metrics at end of ACS, before MILP."""
    tm.schedule_from_tour(seed)
    sc = tm.schedule[tm.schedule["scheduled"]]
    return {
        "acs_wall_s": None,  # filled by caller
        "acs_objective": round(float(seed["objective"]), 2),
        "acs_scheduled": int(len(sc)),
        "acs_modeled_slew": round(float(tm.stats["t_slew_sum"]), 2),
        "acs_physical_slew": round(
            float(physical_slew_minutes(tm.schedule, queue, night_start)), 2
        ),
    }


def _render(tm, queue, df, n_states, request_csv, run_dir, night_start):
    tm.schedule.to_csv(os.path.join(run_dir, "schedule.csv"), index=False)
    tm.observer = queue.observatory
    tm.wrap_limit = queue.wrap_limit
    tm.wrap_states = queue.wrap_states if n_states > 1 else None
    tplot.plot_path_2D_interactive(tm, night_start_time=night_start).write_html(
        os.path.join(run_dir, "slew_path.html")
    )
    anim = tplot.get_slew_animation_plotly(
        tm, request_csv, inaccessible_zones=queue.inaccessible_zones
    )
    tplot.save_slew_animation(anim, os.path.join(run_dir, "slew_animation.html"))
    hover = [c for c in ("target", "exptime", "n_exp") if c not in tm.schedule.columns]
    if hover:
        tm.schedule = tm.schedule.merge(
            df[["unique_id", *hover]], on="unique_id", how="left"
        )
    aqplot.get_ladder(tm, night_start).write_html(os.path.join(run_dir, "ladder.html"))


def _apply_gurobi_params(tm, *, solve_method, milp_s):
    tm.model.params.OutputFlag = 0
    tm.model.params.MIPGap = MIP_GAP
    tm.model.params.PreSolve = 2
    tm.model.params.MIPFocus = 1
    if solve_method == "norel+milp":
        tm.model.params.NoRelHeurTime = NOREL_WARMSTART_S
        tm.model.params.TimeLimit = milp_s - NOREL_WARMSTART_S
        tm.model.params.Heuristics = 0.2
    else:
        tm.model.params.NoRelHeurTime = 0.0
        tm.model.params.TimeLimit = milp_s
    tm.model.update()
    return (
        NOREL_WARMSTART_S if solve_method == "norel+milp" else None,
        milp_s - NOREL_WARMSTART_S if solve_method == "norel+milp" else milp_s,
    )


def _run_one(
    dataset_key,
    ds_info,
    method_key,
    n_states,
    acs_starts,
    milp_s,
    solve_method,
    *,
    access_builder,
    skip_plots,
):
    request_dir = ds_info["request_dir"]
    request_csv = os.path.join(request_dir, "request_selected.csv")
    night_start, night_end = _night_bounds(request_dir)
    n_targets = ds_info["n_targets"]

    run_dir = os.path.join(HERE, "runs", dataset_key, method_key)
    run_json = os.path.join(run_dir, "run.json")
    force = os.environ.get("MATRIX_FORCE", "").strip() in ("1", "true", "yes")
    if os.path.isfile(run_json) and not force:
        with open(run_json, encoding="utf-8") as fh:
            return json.load(fh)

    os.makedirs(run_dir, exist_ok=True)
    queue = HIRESCPS()
    df = pd.read_csv(request_csv)
    if ds_info["builder"] == "access":
        requests = access_builder(df, queue, night_start, night_end)
    else:
        requests = build_requests(df, queue, night_start, night_end)

    row = {
        "dataset": dataset_key,
        "n_targets": n_targets,
        "method": method_key,
        "solve_method": solve_method,
        "n_states": n_states,
        "milp_s": milp_s,
        "acs_starts": acs_starts,
        "norel_heur_time": None,
        "milp_time_limit": None,
        "acs_wall_s": None,
        "acs_objective": None,
        "acs_scheduled": None,
        "acs_physical_slew": None,
        "acs_modeled_slew": None,
    }

    tm = _new_model(queue, requests, night_start, night_end, n_states)
    seed = None

    if solve_method == "acs8+milp":
        tm.build_model()
        acs_kw = dict(
            params={"time_limit_s": ACS_S, "nb_max": ACS_NB_MAX},
            n_starts=acs_starts,
            parallel=True,
        )
        t0 = time.time()
        seed = acs_warm_start(tm, **acs_kw)
        acs_wall = time.time() - t0
        ckpt = _acs_checkpoint(tm, seed, queue, night_start)
        ckpt["acs_wall_s"] = round(acs_wall, 1)
        row.update(ckpt)
        tm.seed_from_tour(seed)
    else:
        tm.build_model()

    norel_heur, milp_time_limit = _apply_gurobi_params(
        tm, solve_method=solve_method, milp_s=milp_s
    )
    if norel_heur is not None:
        row["norel_heur_time"] = norel_heur
        row["milp_time_limit"] = milp_time_limit

    tm.run_model()
    tm.build_schedule()

    sc = tm.schedule[tm.schedule["scheduled"]]
    row.update(
        {
            "scheduled": int(len(sc)),
            "n_requested": int(tm.stats["n_requested"]),
            "physical_slew": round(
                float(physical_slew_minutes(tm.schedule, queue, night_start)), 2
            ),
            "modeled_slew": round(float(tm.stats["t_slew_sum"]), 2),
            "objective": round(float(tm.model.ObjVal), 2),
            "bound": round(float(tm.model.ObjBound), 2),
            "gap_pct": round(float(tm.model.MIPGap) * 100, 2),
            "milp_solve_s": round(float(tm.model.Runtime), 1),
            "nodes": int(tm.model.NodeCount),
        }
    )

    if not skip_plots:
        _render(tm, queue, df, n_states, request_csv, run_dir, night_start)

    with open(run_json, "w", encoding="utf-8") as fh:
        json.dump(row, fh, indent=2)
        fh.write("\n")

    print(
        f"  {dataset_key}/{method_key}: {row['scheduled']}/{row['n_requested']} "
        f"obj={row['objective']}",
        flush=True,
    )
    return row


def _write_aggregate(rows):
    csv_path = os.path.join(HERE, "all_results.csv")
    json_path = os.path.join(HERE, "all_results.json")
    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False)
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=2)
        fh.write("\n")


def _load_all_rows():
    rows_by_key = {}
    for dataset_key in DATASETS:
        for method_key, *_ in METHODS:
            run_json = os.path.join(HERE, "runs", dataset_key, method_key, "run.json")
            if os.path.isfile(run_json):
                with open(run_json, encoding="utf-8") as fh:
                    row = json.load(fh)
                rows_by_key[(dataset_key, method_key)] = row
    ordered = []
    for dataset_key in DATASETS:
        for method_key, *_ in METHODS:
            key = (dataset_key, method_key)
            if key in rows_by_key:
                ordered.append(rows_by_key[key])
    return ordered


def main():
    os.makedirs(os.path.join(HERE, "runs"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "datasets"), exist_ok=True)

    ds_filter = os.environ.get("MATRIX_DATASETS", "").strip()
    method_filter = os.environ.get("MATRIX_METHODS", "").strip()
    skip_plots = os.environ.get("MATRIX_SKIP_PLOTS", "").strip() in ("1", "true", "yes")

    datasets = list(DATASETS.keys())
    methods = METHODS
    if ds_filter:
        datasets = [d.strip() for d in ds_filter.split(",") if d.strip()]
    if method_filter:
        allowed = {m.strip() for m in method_filter.split(",") if m.strip()}
        methods = [m for m in METHODS if m[0] in allowed]

    access_builder = _load_access_builder()
    print(
        f"MILP/ACS matrix: {len(datasets)} datasets × {len(methods)} methods\n"
        f"  out   : {HERE}\n"
        f"  plots : {'skip' if skip_plots else 'on'}\n",
        flush=True,
    )

    for dataset_key in datasets:
        if dataset_key not in DATASETS:
            print(f"Unknown dataset {dataset_key!r}, skip", file=sys.stderr)
            continue
        ds_info = DATASETS[dataset_key]
        print(f"Dataset {dataset_key} (N={ds_info['n_targets']})", flush=True)
        for method_key, n_states, acs_starts, milp_s, solve_method in methods:
            _run_one(
                dataset_key,
                ds_info,
                method_key,
                n_states,
                acs_starts,
                milp_s,
                solve_method,
                access_builder=access_builder,
                skip_plots=skip_plots,
            )

    rows = _load_all_rows()
    _write_aggregate(rows)
    print(f"\nWrote {len(rows)} rows to all_results.csv", flush=True)


if __name__ == "__main__":
    main()
