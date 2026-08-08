"""Generate ``slew_animation.gif`` for existing matrix runs (HTML already present).

Rebuilds the Plotly animation from each run's ``schedule.csv`` without
re-solving. Skips runs whose GIF is newer than the schedule unless
``GIF_FORCE=1``.

Usage::

    python examples/ttp_twostate/milp_acs_matrix/render_gifs.py
    MATRIX_DATASETS=sphere100 GIF_FORCE=1 python .../render_gifs.py
"""

import json
import os
import sys

import pandas as pd
from astropy.time import Time

from astroq.queue.hirescps.queue import HIRESCPS
from astroq.scripts.demo_twostate import build_requests
import astroq.ttp.plot as tplot

from run_matrix import (
    DATASETS,
    METHODS,
    _load_access_builder,
    _new_model,
    _night_bounds,
)

HERE = os.path.dirname(os.path.abspath(__file__))


def _render_gif(run_dir, dataset_key, ds_info, access_builder):
    schedule_path = os.path.join(run_dir, "schedule.csv")
    gif_path = os.path.join(run_dir, "slew_animation.gif")
    run_json = os.path.join(run_dir, "run.json")

    if not os.path.isfile(schedule_path) or not os.path.isfile(run_json):
        return None

    force = os.environ.get("GIF_FORCE", "").strip() in ("1", "true", "yes")
    if (
        not force
        and os.path.isfile(gif_path)
        and os.path.getmtime(gif_path) >= os.path.getmtime(schedule_path)
    ):
        return gif_path

    with open(run_json, encoding="utf-8") as fh:
        run = json.load(fh)
    n_states = int(run["n_states"])

    request_dir = ds_info["request_dir"]
    request_csv = os.path.join(request_dir, "request_selected.csv")
    night_start, night_end = _night_bounds(request_dir)

    queue = HIRESCPS()
    df = pd.read_csv(request_csv)
    if ds_info["builder"] == "access":
        requests = access_builder(df, queue, night_start, night_end)
    else:
        requests = build_requests(df, queue, night_start, night_end)

    tm = _new_model(queue, requests, night_start, night_end, n_states)
    tm.schedule = pd.read_csv(schedule_path)
    tm.observer = queue.observatory
    tm.wrap_limit = queue.wrap_limit
    tm.wrap_states = queue.wrap_states if n_states > 1 else None

    anim = tplot.get_slew_animation_plotly(
        tm, request_csv, inaccessible_zones=queue.inaccessible_zones
    )
    tplot.write_slew_animation_gif(anim, gif_path)
    return gif_path


def main():
    ds_filter = os.environ.get("MATRIX_DATASETS", "").strip()
    method_filter = os.environ.get("MATRIX_METHODS", "").strip()
    datasets = list(DATASETS.keys())
    methods = METHODS
    if ds_filter:
        datasets = [d.strip() for d in ds_filter.split(",") if d.strip()]
    if method_filter:
        allowed = {m.strip() for m in method_filter.split(",") if m.strip()}
        methods = [m for m in METHODS if m[0] in allowed]

    access_builder = _load_access_builder()
    n_done = 0
    for dataset_key in datasets:
        if dataset_key not in DATASETS:
            print(f"Unknown dataset {dataset_key!r}", file=sys.stderr)
            continue
        ds_info = DATASETS[dataset_key]
        for method_key, *_ in methods:
            run_dir = os.path.join(HERE, "runs", dataset_key, method_key)
            out = _render_gif(run_dir, dataset_key, ds_info, access_builder)
            if out:
                n_done += 1
                print(f"  {dataset_key}/{method_key} -> {out}", flush=True)

    print(f"Wrote/verified {n_done} GIFs", flush=True)


if __name__ == "__main__":
    main()
