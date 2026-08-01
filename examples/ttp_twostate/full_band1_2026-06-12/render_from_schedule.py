"""Re-render the TTP slew plots from a saved ``schedule.csv`` (no solve).

The Plotly plot helpers only read attributes off the model object, so we can
feed them a lightweight namespace built from a persisted schedule and skip the
(expensive) Gurobi solve entirely. State-awareness is auto-detected from the
presence of a populated ``wrap_state`` column.

Usage:

    python examples/ttp_twostate/full_band1_2026-06-12/render_from_schedule.py
"""

import os
import types

import pandas as pd
from astropy.time import Time

from astroq.queue.hirescps.queue import HIRESCPS
import astroq.ttp.plot as tplot
import astroq.plot as aqplot

OUT = os.path.dirname(os.path.abspath(__file__))
NIGHT_START = Time("2026-06-12T05:54:00", format="isot")
NIGHT_END = Time("2026-06-12T14:48:00", format="isot")
REQUEST_CSV = os.path.join(OUT, "request_selected.csv")


def render(outdir, queue):
    sched = pd.read_csv(os.path.join(outdir, "schedule.csv"))
    # The ladder hover text reads exptime/n_exp/target, which are not persisted
    # in schedule.csv; merge them back in from request_selected.csv by unique_id.
    req = pd.read_csv(REQUEST_CSV)
    hover_cols = [
        c for c in ("target", "exptime", "n_exp") if c not in sched.columns
    ]
    if hover_cols:
        sched = sched.merge(
            req[["unique_id", *hover_cols]], on="unique_id", how="left"
        )
    state_aware = "wrap_state" in sched.columns and sched["wrap_state"].notna().any()
    model = types.SimpleNamespace(
        schedule=sched,
        observer=queue.observatory,
        night_start=NIGHT_START,
        night_end=NIGHT_END,
        wrap_limit=queue.wrap_limit,
        wrap_states=queue.wrap_states if state_aware else None,
    )
    tplot.plot_path_2D_interactive(model, night_start_time=NIGHT_START).write_html(
        os.path.join(outdir, "slew_path.html")
    )
    tplot.get_slew_animation_plotly(
        model, REQUEST_CSV, inaccessible_zones=queue.inaccessible_zones
    ).write_html(os.path.join(outdir, "slew_animation.html"))
    aqplot.get_ladder(model, NIGHT_START).write_html(
        os.path.join(outdir, "ladder.html")
    )
    n = int(sched["scheduled"].sum())
    print(f"{os.path.basename(outdir):<10} state_aware={state_aware!s:<5} "
          f"scheduled={n}  -> slew_path.html, slew_animation.html, ladder.html")


def main():
    queue = HIRESCPS()
    for sub in ("out_single", "out_two"):
        render(os.path.join(OUT, sub), queue)


if __name__ == "__main__":
    main()
