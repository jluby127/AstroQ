"""Per-night allocation-versus-usage table for a semester plan.

This is the analysis that diagnosed the 2026B U258 seasonal mismatch: it shows,
for every allocated night, how many hours each program actually consumed. A
program that owns a night but books zero hours on it is donating that time to
the rest of the queue.

Usage:
    python night_usage.py <run_dir> [--owners allocation_hires_cps_2026B.csv]

<run_dir> must contain allocation.csv, request.csv, and
outputs/{semester_plan.csv,semester_planner.h5}.

Two details that are easy to get wrong:

1. Day indexing. semester_plan.csv indexes `d` from semester_start_day in
   *observatory-local* nights, while allocation.csv timestamps are UT. The
   pre-midnight part of local night d falls on UT date d+1, so the join needs
   `ut_date = semester_start + d + 1`. Without the shift every night appears
   empty and the totals silently disagree with the run report.

2. Visit length. semester_plan.csv has one row per visit, not per slot, so
   hours must be weighted by `t_visit_slots` from the planner HDF5 (key
   `requests`) times the slot size. Counting rows undercounts long visits.
"""

import argparse
import csv
import datetime
import os

import pandas as pd

SLOT_MINUTES = 2


def block_hours(start, stop):
    a = datetime.datetime.fromisoformat(start)
    b = datetime.datetime.fromisoformat(stop)
    return (b - a).total_seconds() / 3600.0


def load_plan(run_dir, semester_start):
    requests = pd.read_hdf(os.path.join(run_dir, "outputs", "semester_planner.h5"),
                           "requests")
    t_visit = requests.set_index("unique_id")["t_visit_slots"].to_dict()
    program = requests.set_index("unique_id")["program_code"].to_dict()

    plan = pd.read_csv(os.path.join(run_dir, "outputs", "semester_plan.csv"),
                       dtype={"unique_id": str})
    plan["hours"] = plan["unique_id"].map(t_visit).astype(float) * SLOT_MINUTES / 60.0
    plan["program"] = plan["unique_id"].map(program)
    # Local night index -> UT date of the block that contains its evening half.
    plan["ut_date"] = plan["d"].map(
        lambda d: semester_start + datetime.timedelta(days=int(d) + 1)
    )
    return plan


def load_allocation(run_dir, owners_file):
    """Return [(ut_date, hours, owner)]; owner is None without an owners file."""
    owner_by_date = {}
    if owners_file:
        path = os.path.join(run_dir, owners_file)
        if os.path.exists(path):
            for row in csv.DictReader(open(path)):
                start = datetime.datetime.fromisoformat(row["start"])
                owner_by_date[start.date()] = row["ProjCode"].strip(",")

    blocks = []
    for row in csv.DictReader(open(os.path.join(run_dir, "allocation.csv"))):
        date = datetime.datetime.fromisoformat(row["start"]).date()
        blocks.append((date, block_hours(row["start"], row["stop"]),
                       owner_by_date.get(date)))
    return sorted(blocks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--start", default="2026-08-01", help="semester_start_day")
    ap.add_argument("--owners", default="allocation_hires_cps_2026B.csv",
                    help="crossmatched Keck file supplying per-block ProjCode")
    args = ap.parse_args()

    semester_start = datetime.date.fromisoformat(args.start)
    plan = load_plan(args.run_dir, semester_start)
    blocks = load_allocation(args.run_dir, args.owners)

    usage = plan.pivot_table(index="ut_date", columns="program", values="hours",
                             aggfunc="sum").fillna(0.0)

    print(f"{'UT night':11} {'alloc':>6} {'owner':7} {'used':>6} {'idle':>6}  by program")
    total_alloc = total_idle = 0.0
    for date, hours, owner in blocks:
        row = usage.loc[date] if date in usage.index else None
        used = float(row.sum()) if row is not None else 0.0
        detail = ""
        if row is not None:
            detail = " ".join(
                f"{c.split('_')[-1]}:{row[c]:.1f}" for c in usage.columns if row[c] > 0.05
            )
        total_alloc += hours
        total_idle += hours - used
        print(f"{date}  {hours:6.2f} {owner or '-':7} {used:6.1f} {hours - used:6.1f}  {detail}")

    print(f"\ntotal allocated {total_alloc:.2f} h   "
          f"scheduled {usage.values.sum():.2f} h   idle {total_idle:.2f} h")
    print("\nscheduled hours by program:")
    for program, hours in usage.sum().sort_values(ascending=False).items():
        print(f"  {program:12} {hours:6.2f}")


if __name__ == "__main__":
    main()
