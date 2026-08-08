"""Export paper figures from a solved semester planner.

Produces the three semester-level figures used in the paper, as vector PDFs:

    plots/cof_hours.pdf   cumulative time charged, per program, as a
                          percentage of each program's awarded hours
    plots/football.pdf    all-sky map of requests over the observability
                          background
    plots/birdseye.pdf    night-by-slot map of the semester plan, colored by
                          program

All three are rendered from the plan as it stood on the first night of the
semester, so they show the queue as it was intended before any weather or
execution losses.

Usage:
    python tools/make_figures.py <run_dir>

<run_dir> must contain outputs/semester_planner.h5 and programs.csv. Note that
the figures are read out of the HDF5, but build_plot_data also re-reads
programs.csv from the workdir recorded *inside* the HDF5 config, so that
directory must still exist.
"""

import argparse
import os
import sys

import numpy as np
import plotly.io as pio

import astroq.plot as pl
from astroq.splan import SemesterPlanner

# Kaleido otherwise stamps a "Loading [MathJax]..." banner into the corner of
# every exported figure. Nothing here needs MathJax.
pio.kaleido.scope.mathjax = None

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
PLOTS = os.path.join(PAPER, "plots")

# Kaleido rasterizes at these dimensions; PDF output stays vector regardless,
# but the aspect ratio and font scaling follow from them.
COF_SIZE = dict(width=1400, height=900)
FOOTBALL_SIZE = dict(width=1400, height=800)
BIRDSEYE_SIZE = dict(width=1400, height=900)


def crop_to_night(fig, plot_data, pad=3):
    """Trim the birdseye y-axis to the slots that are ever allocated.

    The slot axis covers a full 24 hr day, but only the observing night is ever
    allocated, so the uncropped figure is mostly empty margin.
    """
    # nulltime is (n_slots, n_nights) and is 1 where the slot is unallocated.
    allocated = np.nonzero((plot_data.nulltime == 0).any(axis=1))[0]
    if allocated.size == 0:
        return
    lo = max(int(allocated.min()) - pad, 0)
    hi = min(int(allocated.max()) + pad, plot_data.n_slots - 1)
    fig.update_yaxes(range=[lo, hi], autorange=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", help="run directory holding outputs/semester_planner.h5")
    args = ap.parse_args()

    h5 = os.path.join(args.run_dir, "outputs", "semester_planner.h5")
    if not os.path.exists(h5):
        sys.exit(f"no semester_planner.h5 under {args.run_dir}")

    planner = SemesterPlanner.from_hdf5(h5)
    plot_data = pl.build_plot_data(planner)
    os.makedirs(PLOTS, exist_ok=True)

    # units="time" is the hours variant; the default "requests" counts visits.
    fig = pl.get_cof(
        plot_data,
        programs=sorted(plot_data.program_table.index),
        units="time",
    )
    out = os.path.join(PLOTS, "cof_hours.pdf")
    fig.write_image(out, engine="kaleido", **COF_SIZE)
    print(f"wrote {out}")

    fig = pl.get_football(plot_data, plot_data.select_all(), use_program_colors=True)
    out = os.path.join(PLOTS, "football.pdf")
    fig.write_image(out, engine="kaleido", **FOOTBALL_SIZE)
    print(f"wrote {out}")

    # Aggregating by program collapses the per-request starmaps into one trace
    # per program, so the mosaic is colored the same way as the other two.
    fig = pl.get_birdseye(plot_data, plot_data.select_all(aggregate_by_program=True))
    crop_to_night(fig, plot_data)
    # The webapp's "Today" marker is meaningless in print: the plan shown is by
    # construction the one built on the first night.
    fig.update_layout(shapes=[], annotations=[])
    out = os.path.join(PLOTS, "birdseye.pdf")
    fig.write_image(out, engine="kaleido", **BIRDSEYE_SIZE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
