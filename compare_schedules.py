#!/usr/bin/env python3
"""
Compare two semester_plan.csv (or similar) schedule files and report meaningful differences.

Ignores small slot shifts within the same night (e.g. same request on same night, different slot).
Reports:
- Per-star: which stars have more/fewer observations in one schedule vs the other
- Per-night: which nights have different stars selected, with specifics

Usage:
    python compare_schedules.py -s1 schedule1.csv -s2 schedule2.csv -o report.txt
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_schedule(path: Path) -> pd.DataFrame:
    """Load a schedule CSV with columns r, d, s, name (or similar)."""
    df = pd.read_csv(path)
    # Normalize column names to lowercase
    df.columns = df.columns.str.strip().str.lower()
    required = {"r", "d", "s", "name"}
    if not required.issubset(set(df.columns)):
        raise ValueError(
            f"Schedule {path} must have columns r, d, s, name (got {list(df.columns)})"
        )
    return df


def observations_by_request_night(df: pd.DataFrame) -> set:
    """
    Set of (r, d) for each request r observed on night d.
    Ignores slot s so that small slot shifts on the same night count as the same observation.
    """
    return set(zip(df["r"].astype(str), df["d"].astype(int)))


def observations_by_star(df: pd.DataFrame) -> dict:
    """Count of (r, d) observations per star name."""
    pairs = set(zip(df["r"].astype(str), df["d"].astype(int)))
    # Map r -> name (take first occurrence)
    r_to_name = df.drop_duplicates("r").set_index("r")["name"].to_dict()
    by_star = {}
    for r, d in pairs:
        name = r_to_name.get(r, r)
        by_star[name] = by_star.get(name, 0) + 1
    return by_star


def stars_per_night(df: pd.DataFrame) -> dict:
    """For each night d, set of star names scheduled that night."""
    out = {}
    for _, row in df.iterrows():
        d = int(row["d"])
        name = row["name"]
        out.setdefault(d, set()).add(name)
    return {d: out[d] for d in sorted(out)}


def run_comparison(s1_path: Path, s2_path: Path, out_path: Path) -> None:
    """Load both schedules, compute diffs, write report."""
    df1 = load_schedule(s1_path)
    df2 = load_schedule(s2_path)

    obs1 = observations_by_request_night(df1)
    obs2 = observations_by_request_night(df2)

    by_star1 = observations_by_star(df1)
    by_star2 = observations_by_star(df2)

    stars1_per_night = stars_per_night(df1)
    stars2_per_night = stars_per_night(df2)

    all_stars = sorted(set(by_star1) | set(by_star2))
    all_nights = sorted(set(stars1_per_night) | set(stars2_per_night))

    lines = []
    lines.append("=" * 70)
    lines.append("Schedule comparison report")
    lines.append("=" * 70)
    lines.append(f"Schedule 1: {s1_path}")
    lines.append(f"Schedule 2: {s2_path}")
    lines.append("")
    lines.append("(Differences ignore small slot shifts within the same night.)")
    lines.append("")

    # ---- Per-star stats ----
    lines.append("-" * 70)
    lines.append("1. Per-star observation counts (meaningful differences)")
    lines.append("-" * 70)

    only_in_1 = []
    only_in_2 = []
    diff_count = []

    for star in all_stars:
        c1 = by_star1.get(star, 0)
        c2 = by_star2.get(star, 0)
        if c1 == 0 and c2 > 0:
            only_in_2.append((star, c2))
        elif c1 > 0 and c2 == 0:
            only_in_1.append((star, c1))
        elif c1 != c2:
            diff_count.append((star, c1, c2))

    if only_in_1:
        lines.append("")
        lines.append("Stars with observations only in Schedule 1:")
        for star, c in sorted(only_in_1, key=lambda x: -x[1]):
            lines.append(f"  {star}: {c} obs")
        lines.append("")

    if only_in_2:
        lines.append("Stars with observations only in Schedule 2:")
        for star, c in sorted(only_in_2, key=lambda x: -x[1]):
            lines.append(f"  {star}: {c} obs")
        lines.append("")

    if diff_count:
        lines.append("Stars with different observation counts (s1 vs s2):")
        for star, c1, c2 in sorted(diff_count, key=lambda x: -abs(x[1] - x[2])):
            lines.append(f"  {star}: {c1} vs {c2} (diff: {c2 - c1:+d})")
        lines.append("")

    if not (only_in_1 or only_in_2 or diff_count):
        lines.append("No per-star count differences.")
        lines.append("")

    # Summary counts
    lines.append("Summary:")
    lines.append(f"  Total (r,d) observations in Schedule 1: {len(obs1)}")
    lines.append(f"  Total (r,d) observations in Schedule 2: {len(obs2)}")
    lines.append(
        f"  Stars only in Schedule 1: {len(only_in_1)}  |  only in Schedule 2: {len(only_in_2)}  |  count diff: {len(diff_count)}"
    )
    lines.append("")

    # ---- Per-night: nights with different star sets ----
    lines.append("-" * 70)
    lines.append("2. Nights with different stars selected")
    lines.append("-" * 70)

    nights_diff = []
    for d in all_nights:
        set1 = stars1_per_night.get(d, set())
        set2 = stars2_per_night.get(d, set())
        if set1 != set2:
            only1 = set1 - set2
            only2 = set2 - set1
            nights_diff.append((d, set1, set2, only1, only2))

    if not nights_diff:
        lines.append("")
        lines.append("No nights with different star sets.")
    else:
        for d, set1, set2, only1, only2 in nights_diff:
            lines.append("")
            lines.append(f"Night d={d}:")
            if only1:
                lines.append(f"  Only in Schedule 1: {sorted(only1)}")
            if only2:
                lines.append(f"  Only in Schedule 2: {sorted(only2)}")
            lines.append(f"  Schedule 1 stars ({len(set1)}): {sorted(set1)}")
            lines.append(f"  Schedule 2 stars ({len(set2)}): {sorted(set2)}")
        lines.append("")
        lines.append(f"Total nights with different star sets: {len(nights_diff)}")

    lines.append("")
    lines.append("=" * 70)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Report written to {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare two schedule CSVs and write a difference report."
    )
    parser.add_argument(
        "-s1",
        required=True,
        type=Path,
        dest="schedule1",
        help="Path to first schedule CSV (e.g. semester_plan.csv)",
    )
    parser.add_argument(
        "-s2",
        required=True,
        type=Path,
        dest="schedule2",
        help="Path to second schedule CSV",
    )
    parser.add_argument(
        "-o",
        required=True,
        type=Path,
        dest="output",
        help="Path to output report file (e.g. report.txt)",
    )
    args = parser.parse_args()

    if not args.schedule1.exists():
        print(f"Error: Schedule 1 not found: {args.schedule1}", file=sys.stderr)
        sys.exit(1)
    if not args.schedule2.exists():
        print(f"Error: Schedule 2 not found: {args.schedule2}", file=sys.stderr)
        sys.exit(1)

    run_comparison(args.schedule1, args.schedule2, args.output)


if __name__ == "__main__":
    main()
