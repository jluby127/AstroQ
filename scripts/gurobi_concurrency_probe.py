#!/usr/bin/env python3
"""Probe whether two Gurobi optimizes can run concurrently."""

from __future__ import annotations

import argparse
import json
import random
import time

import gurobipy as gp


def run_probe(
    label: str,
    n_vars: int,
    time_limit: float,
    hold_seconds: float,
    sleep_in_optimize: float,
) -> dict:
    started = time.time()
    try:
        model = gp.Model(f"probe_{label}")
        model.setParam("OutputFlag", 0)
        model.setParam("TimeLimit", time_limit)
        # Hard knapsack: irrational-ish weights keep the LP bound loose so the
        # MIP actually takes seconds instead of milliseconds.
        rng = random.Random(0)
        weights = [rng.uniform(1.0, 2.0) for _ in range(n_vars)]
        values = [w + rng.uniform(-0.01, 0.01) for w in weights]
        x = model.addVars(n_vars, vtype=gp.GRB.BINARY, name="x")
        model.setObjective(
            gp.quicksum(values[j] * x[j] for j in range(n_vars)), gp.GRB.MAXIMIZE
        )
        model.addConstr(
            gp.quicksum(weights[j] * x[j] for j in range(n_vars))
            <= sum(weights) / 2 + 0.5
        )
        model.setParam("MIPGap", 0.0)
        if hold_seconds > 0:
            time.sleep(hold_seconds)

        slept = {"done": False}

        def callback(model_cb, where):
            if sleep_in_optimize <= 0 or slept["done"]:
                return
            if where in (gp.GRB.Callback.POLLING, gp.GRB.Callback.MIP):
                slept["done"] = True
                time.sleep(sleep_in_optimize)

        optimize_started = time.time()
        if sleep_in_optimize > 0:
            model.optimize(callback)
        else:
            model.optimize()
        optimize_elapsed = time.time() - optimize_started
        status = int(model.Status)
        model.dispose()
        return {
            "label": label,
            "ok": True,
            "status": status,
            "optimize_started_epoch": round(optimize_started, 3),
            "optimize_elapsed_s": round(optimize_elapsed, 3),
            "total_elapsed_s": round(time.time() - started, 3),
        }
    except gp.GurobiError as exc:
        return {
            "label": label,
            "ok": False,
            "error_type": "GurobiError",
            "error": str(exc),
            "total_elapsed_s": round(time.time() - started, 3),
        }
    except Exception as exc:  # pragma: no cover - diagnostic script
        return {
            "label": label,
            "ok": False,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "total_elapsed_s": round(time.time() - started, 3),
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--n-vars", type=int, default=800)
    parser.add_argument("--time-limit", type=float, default=45.0)
    parser.add_argument(
        "--hold-seconds",
        type=float,
        default=0.0,
        help="Sleep after model build and before optimize.",
    )
    parser.add_argument(
        "--sleep-in-optimize",
        type=float,
        default=0.0,
        help="Sleep once inside the MIP callback to keep optimize active.",
    )
    args = parser.parse_args()
    result = run_probe(
        args.label,
        args.n_vars,
        args.time_limit,
        args.hold_seconds,
        args.sleep_in_optimize,
    )
    print(json.dumps(result), flush=True)
    return 0 if result.get("ok") else 1


if __name__ == "__main__":
    raise SystemExit(main())
