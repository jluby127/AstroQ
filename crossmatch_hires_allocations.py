#!/usr/bin/env python3
"""Crossmatch the all-scheduled Keck file with requested PI names.

By default this reads allocation_hires_all_scheduled.csv and the current
semester request URL CSV, then writes matched rows to allocation_hires_cps_<semester>.csv.
A match only requires that the schedule PI appears at least once in the PI column
of the request CSV.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


TIME_RE = re.compile(r"(\d{1,2}:\d{2})\s*-\s*(\d{1,2}:\d{2})")


def infer_current_semester(today: date | None = None) -> str:
    """Infer the Keck semester from the current date."""
    today = today or date.today()
    if 2 <= today.month <= 7:
        return f"{today.year}A"
    if today.month >= 8:
        return f"{today.year}B"
    return f"{today.year - 1}B"


def normalize_name(value: str) -> str:
    """Normalize a PI name for forgiving string matching."""
    value = (value or "").strip().casefold()
    value = re.sub(r"[.,]", "", value)
    value = re.sub(r"\s+", " ", value)
    return value


def parse_time_cell(value: str) -> Tuple[str, str]:
    """Extract start and stop times from a Keck schedule Time cell."""
    text = (value or "").strip()
    match = TIME_RE.search(text)
    if not match:
        raise ValueError(f"Could not parse start/stop times from: {text!r}")
    return match.group(1), match.group(2)


def load_request_pis(path: Path) -> Dict[str, str]:
    """Load the PI names from the request CSV PI column."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Request CSV is empty: {path}")
        if "PI" not in reader.fieldnames:
            raise ValueError(f"Request CSV must contain a PI column: {path}")

        result: Dict[str, str] = {}
        for row in reader:
            pi = (row.get("PI") or "").strip()
            if pi:
                result[normalize_name(pi)] = pi
        return result


def collect_matches(schedule_path: Path, request_pis: Dict[str, str]) -> List[dict]:
    """Filter all-scheduled rows to only those whose PI appears in the request list."""
    matched: List[dict] = []
    with schedule_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pi = (row.get("PI") or row.get("Principal") or "").strip()
            if not pi or normalize_name(pi) not in request_pis:
                continue

            date_value = (row.get("Date") or "").strip()
            start_time, stop_time = parse_time_cell(row.get("Time") or "")
            matched.append(
                {
                    "Date": date_value,
                    "StartTime": start_time,
                    "EndTime": stop_time,
                    "start": f"{date_value}T{start_time}",
                    "stop": f"{date_value}T{stop_time}",
                    "PI": pi,
                    "Instrument": (row.get("Instrument") or "").strip(),
                    "ProjCode": (row.get("ProjCode") or "").strip(),
                    "comment": "",
                }
            )

    matched.sort(key=lambda row: (row["Date"], row["StartTime"], row["PI"], row["Instrument"]))
    return matched


def write_matches(rows: Iterable[dict], output_path: Path) -> None:
    """Write the matched rows to CSV."""
    fieldnames = ["Date", "StartTime", "EndTime", "start", "stop", "PI", "Instrument", "ProjCode", "comment"]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Crossmatch the all-scheduled Keck allocation file against request CSV PI names."
    )
    parser.add_argument(
        "-s",
        "--semester",
        default=infer_current_semester(),
        help="Semester label such as 2026A or 2026B. Default: current semester.",
    )
    parser.add_argument(
        "--scheduled",
        default="allocation_hires_all_scheduled.csv",
        help="Path to the all-scheduled allocation file.",
    )
    parser.add_argument(
        "--requests",
        default=None,
        help="Path to the request URL CSV with a PI column. Default: request_urls_<semester>.csv",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output path for the PI-matched allocation CSV. Default: allocation_hires_cps_<semester>.csv",
    )
    args = parser.parse_args()

    request_path = Path(args.requests or f"request_urls_{args.semester}.csv")
    schedule_path = Path(args.scheduled)
    output_path = Path(args.output or f"allocation_hires_cps_{args.semester}.csv")

    request_pis = load_request_pis(request_path)
    matched_rows = collect_matches(schedule_path, request_pis)
    write_matches(matched_rows, output_path)

    print(f"Loaded {len(request_pis)} unique PI names from {request_path}")
    print(f"Matched {len(matched_rows)} schedule rows from {schedule_path}")
    print(f"Wrote PI-matched allocation file: {output_path}")

    counts = Counter(row["PI"] for row in matched_rows)
    for pi, count in sorted(counts.items(), key=lambda item: (-item[1], item[0].casefold())):
        print(f"  {pi}: {count}")


if __name__ == "__main__":
    main()
