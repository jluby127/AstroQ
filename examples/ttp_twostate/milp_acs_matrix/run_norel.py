"""Explore Gurobi NoRel heuristic on matrix benchmarks.

NoRelHeurTime limits seconds spent in Gurobi's No-Relaxation heuristic before
the root LP relaxation. Default is 0 (disabled).

Usage::

    python examples/ttp_twostate/milp_acs_matrix/run_norel.py
    NOREL_DATASET=sphere100 python examples/ttp_twostate/milp_acs_matrix/run_norel.py
    NOREL_TIMES=0,30,60,120 NOREL_METHODS=single_milp_120,two_milp_120 python ...
    MILP_S=120 MATRIX_FORCE=1 python ...
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

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SPHERE_DS = os.path.join(ROOT, "sphere_datasets")

MIP_GAP = 0.005
ACS_S = 30.0
ACS_STARTS = 8
ACS_NB_MAX = 25

DATASETS = {
    "full-band1": {
        "n_targets": 58,
        "request_dir": os.path.join(ROOT, "full_band1_2026-06-12"),
        "builder": "demo",
    },
    "sphere100": {
        "n_targets": 100,
        "request_dir": os.path.join(SPHERE_DS, "sphere100_2026-02-01"),
        "builder": "access",
    },
}

METHODS = {
    "single_milp_120": dict(n_states=1, acs_starts=0),
    "two_milp_120": dict(n_states=2, acs_starts=0),
    "single_acs8_120": dict(n_states=1, acs_starts=ACS_STARTS),
    "two_acs8_120": dict(n_states=2, acs_starts=ACS_STARTS),
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
        return Time(meta["night_start"], format="isot"), Time(meta["night_end"], format="isot")
    return (
        Time("2026-06-12T05:54:00", format="isot"),
        Time("2026-06-12T14:48:00", format="isot"),
    )


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


def _run_one(
    dataset_key,
    ds_info,
    method_key,
    norel_heur_time,
    *,
    milp_s,
    force,
    log_gurobi,
    access_builder,
):
    cfg = METHODS[method_key]
    n_states = cfg["n_states"]
    acs_starts = cfg["acs_starts"]
    request_dir = ds_info["request_dir"]
    run_name = f"{method_key}_norel{norel_heur_time:g}s"
    run_dir = os.path.join(HERE, "runs", dataset_key, "norel", run_name)
    run_json = os.path.join(run_dir, "run.json")
    if os.path.isfile(run_json) and not force:
        with open(run_json, encoding="utf-8") as fh:
            return json.load(fh)

    os.makedirs(run_dir, exist_ok=True)
    queue = HIRESCPS()
    request_csv = os.path.join(request_dir, "request_selected.csv")
    df = pd.read_csv(request_csv)
    night_start, night_end = _night_bounds(request_dir)
    if ds_info["builder"] == "access":
        requests = access_builder(df, queue, night_start, night_end)
    else:
        requests = build_requests(df, queue, night_start, night_end)

    row = {
        "dataset": dataset_key,
        "n_targets": ds_info["n_targets"],
        "method": method_key,
        "norel_heur_time": float(norel_heur_time),
        "milp_s": milp_s,
        "n_states": n_states,
        "acs_starts": acs_starts,
    }

    tm = _new_model(queue, requests, night_start, night_end, n_states)
    tm.build_model()

    if acs_starts > 0:
        t0 = time.time()
        seed = acs_warm_start(
            tm,
            params={"time_limit_s": ACS_S, "nb_max": ACS_NB_MAX},
            n_starts=acs_starts,
            parallel=True,
        )
        row["acs_wall_s"] = round(time.time() - t0, 1)
        tm.seed_from_tour(seed)

    tm.model.params.OutputFlag = 1 if log_gurobi else 0
    if log_gurobi:
        tm.model.params.LogFile = os.path.join(run_dir, "gurobi.log")
    tm.model.params.TimeLimit = milp_s
    tm.model.params.MIPGap = MIP_GAP
    tm.model.params.PreSolve = 2
    tm.model.params.MIPFocus = 1
    tm.model.params.NoRelHeurTime = float(norel_heur_time)
    tm.model.update()

    t0 = time.time()
    tm.run_model()
    wall = time.time() - t0
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
            "objective": round(float(tm.model.ObjVal), 2) if tm.model.SolCount else None,
            "bound": round(float(tm.model.ObjBound), 2) if tm.model.SolCount else None,
            "gap_pct": round(float(tm.model.MIPGap) * 100, 2) if tm.model.SolCount else None,
            "milp_solve_s": round(float(tm.model.Runtime), 1),
            "wall_s": round(wall, 1),
            "nodes": int(tm.model.NodeCount),
            "status": int(tm.model.Status),
        }
    )
    tm.schedule.to_csv(os.path.join(run_dir, "schedule.csv"), index=False)

    with open(run_json, "w", encoding="utf-8") as fh:
        json.dump(row, fh, indent=2)
        fh.write("\n")

    print(
        f"  {dataset_key}/{run_name}: {row['scheduled']}/{row['n_requested']} "
        f"obj={row['objective']} gap={row['gap_pct']}% "
        f"phys_slew={row['physical_slew']} nodes={row['nodes']}",
        flush=True,
    )
    return row


def main():
    dataset_key = os.environ.get("NOREL_DATASET", "full-band1").strip()
    if dataset_key not in DATASETS:
        raise SystemExit(f"Unknown dataset {dataset_key!r}; choose from {sorted(DATASETS)}")

    ds_info = DATASETS[dataset_key]
    access_builder = _load_access_builder() if ds_info["builder"] == "access" else None

    milp_s = int(os.environ.get("MILP_S", "120"))
    force = os.environ.get("MATRIX_FORCE", "").strip() in ("1", "true", "yes")
    log_gurobi = os.environ.get("NOREL_LOG", "").strip() in ("1", "true", "yes")

    methods = os.environ.get(
        "NOREL_METHODS", "single_milp_120,two_milp_120,two_acs8_120"
    ).split(",")
    methods = [m.strip() for m in methods if m.strip()]
    for m in methods:
        if m not in METHODS:
            raise SystemExit(f"Unknown method {m!r}; choose from {sorted(METHODS)}")

    norel_times = [
        float(x.strip())
        for x in os.environ.get("NOREL_TIMES", "0,30,60,120").split(",")
        if x.strip()
    ]

    results_csv = os.path.join(HERE, f"norel_{dataset_key}_results.csv")
    print(
        f"{dataset_key} NoRel sweep: methods={methods} "
        f"norel_times={norel_times} milp_s={milp_s}",
        flush=True,
    )

    rows = []
    for method_key in methods:
        for norel_t in norel_times:
            rows.append(_run_one(
                dataset_key,
                ds_info,
                method_key,
                norel_t,
                milp_s=milp_s,
                force=force,
                log_gurobi=log_gurobi,
                access_builder=access_builder,
            ))

    df = pd.DataFrame(rows)
    df.to_csv(results_csv, index=False)
    print(f"\nWrote {results_csv}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
