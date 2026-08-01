"""Sequential night-by-night semester replan under simulated weather losses.

Walks the allocated nights of a semester in order. On each night the semester
plan is rebuilt from scratch exactly as it would be in operations, then a coin
is flipped: with probability ``--loss-prob`` the night is wholly lost and no
data are recorded, otherwise every visit the plan assigned to that night is
appended to ``past.csv`` and carried into the next night's replan.

The point is to see whether the balance stage keeps programs from drifting
apart as nights are lost. Run it twice on the same seed, once per arm:

    python tools/weather_sim.py --arm balance   --seed 1 --out sim/s1-balance
    python tools/weather_sim.py --arm nobalance --seed 1 --out sim/s1-nobalance

The weather realization depends only on ``--seed`` and the night list, never on
the arm, so the two runs see an identical sequence of lost nights and any
difference between them is attributable to the balance stage.

Two readings of fill shortfall are recorded, because they answer different
questions:

  fsf         maff recomputed nightly, as operations does. Measures the gap to
              what a program could still achieve from here. Weather lowers maff
              too, so this stays small even as programs fall short.
  fsf_frozen  maff held at its night-1 value. Measures how far a program has
              fallen below what looked achievable at the start of the semester.
"""

import argparse
import configparser
import csv
import datetime as dt
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

ARM_MODES = {
    "balance": "shortfall,balance,prioritize,fill-empty,fill-current-day",
    "nobalance": "shortfall,prioritize,fill-empty,fill-current-day",
}

# Inputs that are fixed for the whole semester and simply copied each night.
STATIC_INPUTS = ("request.csv", "allocation.csv", "custom.csv")


def night_label(ts):
    """Observing-night date for a UT timestamp at Maunakea (UTC-10).

    The whole Maunakea night falls on a single UT date, one day after the
    civil date the night is named for, so the label is just UT date minus one.
    """
    return (ts - dt.timedelta(days=1)).date()


def allocated_nights(allocation_csv):
    """Sorted observing-night dates that have at least one allocated block."""
    nights = {}
    with open(allocation_csv) as fh:
        for row in csv.DictReader(fh):
            start = dt.datetime.fromisoformat(row["start"])
            stop = dt.datetime.fromisoformat(row["stop"])
            label = night_label(start)
            lo, hi = nights.get(label, (start, stop))
            nights[label] = (min(lo, start), max(hi, stop))
    return sorted(nights.items())


def write_config(template, dest, workdir, current_day, mode, gurobi):
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(template)
    cfg.set("global", "workdir", str(workdir))
    cfg.set("global", "current_day", current_day.isoformat())
    cfg.set("semester", "mode", mode)
    for section, key, value in gurobi:
        if not cfg.has_section(section):
            cfg.add_section(section)
        cfg.set(section, key, str(value))
    with open(dest, "w") as fh:
        cfg.write(fh)


def parse_final_report(log_path):
    """Per-program rows from the last run report in an astroq log."""
    with open(log_path) as fh:
        text = fh.read()
    marker = "Run report (fill-current-day):"
    if marker not in text:
        raise RuntimeError(f"no final run report in {log_path}")
    tail = text[text.rindex(marker):]
    rows = {}
    for line in tail.splitlines():
        parts = line.split()
        # Report lines are prefixed by the log timestamp, so the program code is
        # not necessarily the first field.
        code = next((p for p in parts[:8] if p.startswith("2026B_")), None)
        if code is None:
            continue
        i = parts.index(code)
        vals = parts[i + 1:i + 11]
        if len(vals) < 10:
            continue
        rows[code] = {
            "awarded_hr": float(vals[0]),
            "requested_hr": float(vals[1]),
            "past_hr": float(vals[2]),
            "proj_hr": float(vals[3]),
            "fill": _pct(vals[7]),
            "maff": _pct(vals[8]),
        }
    if not rows:
        raise RuntimeError(f"could not parse program rows from {log_path}")
    return rows


def _pct(token):
    token = token.rstrip("%")
    return np.nan if token in ("-", "") else float(token) / 100.0


def parse_solve_stats(log_path):
    """Final gap and runtime for each pipeline stage.

    Non-convergence of the shortfall stage is the main threat to interpreting
    an arm comparison, since unfinished stage-1 optimization gets picked up by
    later stages, so the gaps are recorded alongside the results rather than
    assumed away.
    """
    import re

    with open(log_path) as fh:
        text = fh.read()
    stats = {}
    for stage in ("shortfall", "balance", "prioritize", "fill-empty",
                  "fill-current-day"):
        found = re.findall(
            rf"{stage} solve: ([\w ]+), objective=([\d.e+-]+), "
            rf"bound=[\d.e+-]+, gap=([\d.]+)%, (\d+)s",
            text,
        )
        if not found:
            continue
        # compute-max-fill emits one "shortfall solve" per program before the
        # pipeline's own; the pipeline stage is always the last.
        status, obj, gap, secs = found[-1]
        stats[f"{stage}_gap"] = float(gap) / 100.0
        stats[f"{stage}_s"] = int(secs)
        if stage == "shortfall":
            stats["shortfall_objective"] = float(obj)
            stats["shortfall_optimal"] = status.strip() == "optimal"
    return stats


def run(cmd, cwd, log_path):
    """Run an astroq subcommand, appending its output to the run's log.

    The run report is split across both streams: the "Run report (stage)"
    marker is logged to stderr while the tables beneath it are printed to
    stdout. They have to be interleaved into one stream, as the production
    Makefile does with `2>&1 | tee`, or the tables lose their marker.
    """
    result = subprocess.run(cmd, cwd=cwd, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    with open(log_path, "a") as fh:
        fh.write(result.stdout)
    if result.returncode != 0:
        sys.stderr.write(result.stdout[-4000:])
        raise RuntimeError(f"failed ({result.returncode}): {' '.join(cmd)}")


def past_rows_for(plan_tonight, requests, block):
    """Past-observation rows for the visits scheduled on one night.

    Timestamps are spread across the night's allocated window. The planner maps
    each UT timestamp to an observing night via ``civil_night_label``, so they
    must fall inside the block rather than merely on the right date.
    """
    start, stop = block
    span = (stop - start).total_seconds()
    n = max(len(plan_tonight), 1)
    exptime = requests.set_index("unique_id")["exptime"].to_dict()
    rows = []
    for k, (_, visit) in enumerate(plan_tonight.iterrows()):
        uid = str(visit["unique_id"])
        ts = start + dt.timedelta(seconds=span * (k + 0.5) / n)
        rows.append({
            "unique_id": uid,
            "target": visit["target"],
            "timestamp": ts.strftime("%Y-%m-%dT%H:%M:%S"),
            "exposure_time": exptime.get(uid, 0),
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=sorted(ARM_MODES), required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", required=True, help="output directory for this run")
    ap.add_argument("--source", required=True,
                    help="prepped run directory supplying the fixed inputs")
    ap.add_argument("--loss-prob", type=float, default=0.25)
    ap.add_argument("--nights", type=int, default=0,
                    help="stop after this many nights (0 = whole semester)")
    ap.add_argument("--threads", type=int, default=0,
                    help="Gurobi Threads per solve (0 = Gurobi default)")
    ap.add_argument("--time-limit", type=int, default=900,
                    help="shortfall-stage TimeLimit; must be generous enough "
                         "for the stage to converge or the arm comparison is "
                         "polluted by unfinished stage-1 optimization")
    ap.add_argument("--balance-time-limit", type=int, default=300)
    ap.add_argument("--norel", type=int, default=300)
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args()

    source = os.path.abspath(args.source)
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)

    template = os.path.join(source, "config.ini")
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(template)
    semester_start = dt.date.fromisoformat(cfg.get("global", "semester_start_day"))

    nights = allocated_nights(os.path.join(source, "allocation.csv"))
    if args.nights:
        nights = nights[:args.nights]

    # Weather depends only on the seed and the night list, so both arms see the
    # same realization.
    clear = np.random.default_rng(args.seed).random(len(nights)) >= args.loss_prob
    print(f"{args.arm} seed={args.seed}: {len(nights)} allocated nights, "
          f"{int((~clear).sum())} lost ({args.loss_prob:.0%} target)")

    requests = pd.read_csv(os.path.join(source, "request.csv"),
                           dtype={"unique_id": str})
    programs = pd.read_csv(os.path.join(source, "programs.csv"))
    programs = programs.drop(columns=["max_feasible_fill"], errors="ignore")

    gurobi = [
        ("semester.shortfall.gurobi", "TimeLimit", args.time_limit),
        ("semester.shortfall.gurobi", "NoRelHeurTime", args.norel),
        ("semester.balance.gurobi", "TimeLimit", args.balance_time_limit),
    ]
    if args.threads:
        gurobi.append(("semester.default.gurobi", "Threads", args.threads))

    past = []
    frozen_maff = None
    records = []

    for i, (night, block) in enumerate(nights):
        d_index = (night - semester_start).days
        rundir = os.path.join(out, night.isoformat())
        os.makedirs(rundir, exist_ok=True)
        for name in STATIC_INPUTS:
            shutil.copy(os.path.join(source, name), rundir)
        # The sky-availability grids depend only on the semester, so reuse them
        # rather than recomputing for every night.
        cache = os.path.join(rundir, "cache")
        if not os.path.exists(cache):
            shutil.copytree(os.path.join(source, "cache"), cache)
        programs.to_csv(os.path.join(rundir, "programs.csv"), index=False)
        pd.DataFrame(
            past, columns=["unique_id", "target", "timestamp", "exposure_time"]
        ).to_csv(os.path.join(rundir, "past.csv"), index=False)
        write_config(template, os.path.join(rundir, "config.ini"), rundir,
                     night, ARM_MODES[args.arm], gurobi)

        log_path = os.path.join(rundir, "astroq.log")
        open(log_path, "w").close()
        started = dt.datetime.now()
        for sub in ("compute-max-fill", "plan-semester"):
            run([args.python, "-m", "astroq.cli", sub, "-cf", "config.ini"],
                rundir, log_path)
        elapsed = (dt.datetime.now() - started).total_seconds()

        report = parse_final_report(log_path)
        solve_stats = parse_solve_stats(log_path)
        if frozen_maff is None:
            frozen_maff = {p: r["maff"] for p, r in report.items()}

        plan = pd.read_csv(os.path.join(rundir, "outputs", "semester_plan.csv"),
                           dtype={"unique_id": str})
        tonight = plan[plan["d"] == d_index]
        if clear[i]:
            past.extend(past_rows_for(tonight, requests, block))

        for program, row in report.items():
            records.append({
                "night_index": i,
                "date": night.isoformat(),
                "clear": bool(clear[i]),
                "visits_tonight": int(len(tonight)),
                "elapsed_s": round(elapsed),
                **solve_stats,
                "program": program,
                **row,
                "fsf": max(row["maff"] - row["fill"], 0.0),
                "fsf_frozen": max(frozen_maff[program] - row["fill"], 0.0),
            })
        pd.DataFrame(records).to_csv(os.path.join(out, "metrics.csv"), index=False)

        science = [r for r in records[-len(report):]
                   if r["awarded_hr"] > 0 and r["program"] != "2026B_E475"]
        fills = [r["fill"] for r in science]
        converged = "" if solve_stats.get("shortfall_optimal") else " SHORTFALL-GAP=%.1f%%" % (
            100 * solve_stats.get("shortfall_gap", float("nan"))
        )
        print(f"  [{i + 1}/{len(nights)}] {night} "
              f"{'clear ' if clear[i] else 'LOST  '} "
              f"visits={len(tonight):3d} fill min={min(fills):.2f} "
              f"max={max(fills):.2f} worst_fsf={max(r['fsf'] for r in science):.3f} "
              f"({elapsed:.0f}s){converged}", flush=True)

    print(f"wrote {os.path.join(out, 'metrics.csv')}")


if __name__ == "__main__":
    main()
