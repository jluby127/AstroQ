"""Sweep ACS runtime + hyperparameters on the full-band1 target set.

Builds the two-state TTP model once, then runs the ACS heuristic (no Gurobi)
across several parameter axes, measuring objective and wall time. Each config is
repeated with different seeds to expose the heuristic's variance.

Axes:
  time_limit_s : single-start runtime curve
  n_starts     : parallel multi-start (per-start budget fixed)
  nb_max       : neighbor-list cap per node
  beta         : greedy (heuristic) weight
  max_ants     : ants per iteration
  swap_window  : local-search swap span

Writes ``sweep_results.csv`` and PNG figures next to this file.

Usage (astroq-testing env, repo root):
    python examples/ttp_twostate/full-band1_acs_sweep/acs_sweep.py
"""

import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.ttp.model import TTPModel
from astroq.ttp.acs import _ACS
from astroq.scripts.demo_twostate import build_requests

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), "full_band1_2026-06-12")
REQUEST_CSV = os.path.join(SRC, "request_selected.csv")
NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
NIGHT_END = Time("2026-06-12T14:48:00", format="isot")

BASE = {"nb_max": 25, "beta": 2.0, "max_ants": 10, "swap_window": 4}


def build_model():
    queue = HIRESCPS()
    df = pd.read_csv(REQUEST_CSV)
    requests = build_requests(df, queue, NIGHT_START, NIGHT_END)
    tm = TTPModel(
        requests=requests,
        night_start=NIGHT_START,
        night_end=NIGHT_END,
        slew_fn=queue.slew_fn_state,
        n_slots=queue.nSlots,
        n_states=2,
    )
    tm.build_nodes()
    tm.build_arcs()
    return tm


def run_acs(tm, params, n_starts, seed):
    t0 = time.time()
    acs = _ACS(tm, params=params, rng=np.random.default_rng(seed))
    if n_starts > 1:
        res = acs.solve_multistart(n_starts, base_seed=seed, parallel=True)
    else:
        res = acs.solve()
    return res["objective"], len(res["order"]), time.time() - t0


def sweep(tm, axis, values, *, time_s, n_starts, reps, extra=None):
    rows = []
    for val in values:
        params = {**BASE, "time_limit_s": time_s}
        if extra:
            params.update(extra)
        if axis in params or axis == "time_limit_s":
            params[axis] = val
        ns = val if axis == "n_starts" else n_starts
        ts = val if axis == "time_limit_s" else time_s
        params["time_limit_s"] = ts
        objs, walls, sched = [], [], []
        for rep in range(reps):
            obj, n_sched, wall = run_acs(tm, params, ns, seed=1000 + 31 * rep)
            objs.append(obj)
            walls.append(wall)
            sched.append(n_sched)
        rows.append(
            {
                "axis": axis,
                "value": val,
                "time_s": ts,
                "n_starts": ns,
                "obj_best": max(objs),
                "obj_mean": float(np.mean(objs)),
                "obj_std": float(np.std(objs)),
                "sched_best": max(sched),
                "wall_mean": float(np.mean(walls)),
                "reps": reps,
            }
        )
        print(
            f"  {axis}={val!s:<6} obj best/mean/std = "
            f"{max(objs):7.2f}/{np.mean(objs):7.2f}/{np.std(objs):5.2f}"
            f"   sched={max(sched)}  wall={np.mean(walls):5.1f}s"
        )
    return rows


def main():
    tm = build_model()
    n_real = tm.N - 2
    print(f"full-band1 two-state model: {n_real} real nodes\n")
    all_rows = []

    print("[time_limit_s] single start")
    all_rows += sweep(
        tm, "time_limit_s", [2, 5, 10, 20, 40], time_s=10, n_starts=1, reps=4
    )
    print("[n_starts] per-start 10s (parallel)")
    all_rows += sweep(tm, "n_starts", [1, 3, 5, 10], time_s=10, n_starts=1, reps=3)
    print("[nb_max] 10s single start")
    all_rows += sweep(tm, "nb_max", [10, 15, 25, 40], time_s=10, n_starts=1, reps=3)
    print("[beta] 10s single start")
    all_rows += sweep(tm, "beta", [1.0, 2.0, 4.0], time_s=10, n_starts=1, reps=3)
    print("[max_ants] 10s single start")
    all_rows += sweep(tm, "max_ants", [5, 10, 20], time_s=10, n_starts=1, reps=3)
    print("[swap_window] 10s single start")
    all_rows += sweep(tm, "swap_window", [2, 4, 8], time_s=10, n_starts=1, reps=3)

    res = pd.DataFrame(all_rows)
    res.to_csv(os.path.join(HERE, "sweep_results.csv"), index=False)
    _plot(res)
    print("\nwrote sweep_results.csv + figures to", HERE)


def _plot(res):
    # Figure 1: runtime curve (single start) + n_starts curve.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    t = res[res.axis == "time_limit_s"].sort_values("value")
    axes[0].fill_between(
        t.value, t.obj_mean - t.obj_std, t.obj_mean + t.obj_std,
        alpha=0.2, color="C0",
    )
    axes[0].plot(t.value, t.obj_mean, "o-", color="C0", label="mean")
    axes[0].plot(t.value, t.obj_best, "^--", color="C1", label="best")
    axes[0].set(xlabel="per-start runtime (s)", ylabel="ACS objective",
                title="Runtime curve (single start)")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    ns = res[res.axis == "n_starts"].sort_values("value")
    axes[1].plot(ns.value, ns.obj_best, "s-", color="C2", label="best-of-N")
    axes[1].plot(ns.value, ns.obj_mean, "o--", color="C3", label="run mean")
    axes[1].set(xlabel="n_starts (10s each, parallel)", ylabel="ACS objective",
                title="Multi-start benefit")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "runtime_curves.png"), dpi=130)
    plt.close(fig)

    # Figure 2: categorical hyperparameters (best + mean bars).
    cats = ["nb_max", "beta", "max_ants", "swap_window"]
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.8))
    for ax, c in zip(axes, cats):
        d = res[res.axis == c].sort_values("value")
        x = np.arange(len(d))
        ax.bar(x - 0.2, d.obj_best, 0.4, label="best", color="C0")
        ax.bar(x + 0.2, d.obj_mean, 0.4, label="mean", color="C1")
        ax.set_xticks(x)
        ax.set_xticklabels([str(v) for v in d.value])
        ax.set(xlabel=c, title=c)
        lo = float(min(d.obj_mean.min(), d.obj_best.min()))
        hi = float(max(d.obj_best.max(), d.obj_mean.max()))
        ax.set_ylim(lo - 0.05 * (hi - lo + 1), hi + 0.05 * (hi - lo + 1))
        ax.grid(alpha=0.3, axis="y")
    axes[0].set_ylabel("ACS objective (10s)")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "hyperparams.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
