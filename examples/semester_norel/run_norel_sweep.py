"""Sweep Gurobi NoRelHeurTime on the semester (splan) MILP.

Usage::

    python examples/semester_norel/run_norel_sweep.py CONFIG.ini
    NOREL_TIMES=0,30,60,120 python examples/semester_norel/run_norel_sweep.py \\
        2026A/2026-06-19/band1/config.ini
"""

import configparser
import os
import sys
import time

import pandas as pd

from astroq import splan
from gurobipy import GRB

HERE = os.path.dirname(os.path.abspath(__file__))


def _run_one(config_file, norel_heur_time):
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(config_file)
    time_limit = cfg.getfloat("semester.default.gurobi", "TimeLimit", fallback=300.0)
    if not cfg.has_section("semester.shortfall.gurobi"):
        cfg.add_section("semester.shortfall.gurobi")
    cfg.set("semester.shortfall.gurobi", "NoRelHeurTime", str(norel_heur_time))
    cfg.set("semester.shortfall.gurobi", "OutputFlag", "0")

    t0 = time.time()
    sp = splan.SemesterPlanner(config_file)
    sp.config = cfg
    sp.config.optionxform = str
    sp.run_model_shortfall()
    wall = time.time() - t0

    sched = sp.schedule
    t_visit_slots = sp.requests_active.set_index("unique_id")["t_visit_slots"]
    slots_per_visit = sched["unique_id"].map(t_visit_slots).fillna(1)
    scheduled_starting = len(sched)
    reserved = int((slots_per_visit - 1).clip(lower=0).sum())
    total_scheduled = scheduled_starting + reserved
    allocated = int(sp.access_record["is_allocated"][0].sum())
    rf_slots = sp.requests_active["t_visit_slots"]
    total_requested = int(
        (
            rf_slots
            * sp.requests_active["n_intra_max"]
            * sp.requests_active["n_inter_max"]
        ).sum()
    )
    util_requested_pct = (
        100 * total_scheduled / total_requested if total_requested else 0.0
    )

    m = sp.model
    row = {
        "config": config_file,
        "norel_heur_time": float(norel_heur_time),
        "time_limit": time_limit,
        "objective": round(float(m.ObjVal), 2) if m.SolCount else None,
        "bound": round(float(m.ObjBound), 2) if m.SolCount else None,
        "gap_pct": round(float(m.MIPGap) * 100, 4) if m.SolCount else None,
        "milp_solve_s": round(float(m.Runtime), 1),
        "wall_s": round(wall, 1),
        "nodes": int(m.NodeCount),
        "status": int(m.Status),
        "slots_scheduled": int(total_scheduled),
        "slots_starting": int(scheduled_starting),
        "util_requested_pct": round(float(util_requested_pct), 2),
    }
    print(
        f"  norel={norel_heur_time:g}s: obj={row['objective']} gap={row['gap_pct']}% "
        f"slots={row['slots_scheduled']} wall={row['wall_s']}s nodes={row['nodes']}",
        flush=True,
    )
    return sp, row


def main():
    config_file = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
        "NOREL_CONFIG",
        os.path.join(HERE, "../../2026A/2026-06-19/band1/config.ini"),
    )
    config_file = os.path.abspath(config_file)
    if not os.path.isfile(config_file):
        raise SystemExit(f"Config not found: {config_file}")

    norel_times = [
        float(x.strip())
        for x in os.environ.get("NOREL_TIMES", "0,30,60,120").split(",")
        if x.strip()
    ]
    out_dir = os.path.dirname(config_file)
    results_csv = os.path.join(out_dir, "norel_semester_results.csv")

    print(
        f"Semester NoRel sweep: config={config_file} norel_times={norel_times}",
        flush=True,
    )

    rows = []
    best_sp = None
    best_obj = float("inf")
    for norel_t in norel_times:
        sp, row = _run_one(config_file, norel_t)
        rows.append(row)
        if row["objective"] is not None and row["objective"] < best_obj:
            best_obj = row["objective"]
            best_sp = sp

    df = pd.DataFrame(rows)
    df.to_csv(results_csv, index=False)
    print(f"\nWrote {results_csv}")
    print(df.to_string(index=False))

    if best_sp is not None:
        best_sp.write_request_selected()
        best_sp.to_hdf5()
        print(f"\nSaved outputs from best run (obj={best_obj}) to {out_dir}/outputs/")


if __name__ == "__main__":
    main()
