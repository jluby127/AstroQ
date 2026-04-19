#!/usr/bin/env python3
"""Query the Keck observing schedule form and list PIs by instrument.

Defaults to the current semester, which for the current date resolves to 2026A.
By default it queries the instruments KPF, KPF-CC, and HIRES.
"""

from __future__ import annotations

import argparse
import csv
import io
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Iterable, List, Tuple

import requests

QUERY_URL = "https://www2.keck.hawaii.edu/observing/keckSchedule/queryForm.php"
DEFAULT_INSTRUMENTS = ["KPF", "KPF-CC", "HIRES"]


def canonicalize_instrument(name: str) -> str:
    """Map Keck schedule instrument variants to a canonical name."""
    value = (name or "").strip()
    lower = value.lower()
    if lower.startswith("hires"):
        return "HIRES"
    if lower.startswith("kpf-cc"):
        return "KPF-CC"
    if lower == "kpf":
        return "KPF"
    return value


def infer_current_semester(today: date | None = None) -> str:
    """Infer the Keck semester from the current date."""
    today = today or date.today()
    # Semester A: Feb-Jul of the same year; Semester B: Aug-Jan spanning years.
    if 2 <= today.month <= 7:
        return f"{today.year}A"
    if today.month >= 8:
        return f"{today.year}B"
    return f"{today.year - 1}B"


def semester_date_range(semester: str) -> Tuple[str, str]:
    """Convert a semester like 2026A into start/end dates."""
    semester = semester.strip().upper()
    if len(semester) != 5 or semester[-1] not in {"A", "B"} or not semester[:4].isdigit():
        raise ValueError(f"Invalid semester format: {semester!r}. Expected forms like 2026A or 2026B.")

    year = int(semester[:4])
    half = semester[-1]
    if half == "A":
        return f"{year}-02-01", f"{year}-07-31"
    return f"{year}-08-01", f"{year + 1}-01-31"


def fetch_schedule_csv(semester: str, instrument: str, timeout: int = 60) -> List[dict]:
    """Submit the schedule query form and return CSV rows as dictionaries."""
    start_date, end_date = semester_date_range(semester)
    payload = {
        "doQuery": "1",
        "table": "schedule",
        "Date": f"between {start_date} and {end_date}",
        "Instrument": instrument,
        "cb_Date": "on",
        "cb_TelNr": "on",
        "cb_Instrument": "on",
        "cb_Account": "on",
        "cb_Principal": "on",
        "cb_Institution": "on",
        "cb_ProjCode": "on",
        "excel": "on",
        "sched": "Query Tel Schedule",
    }

    response = requests.post(QUERY_URL, data=payload, timeout=timeout)
    response.raise_for_status()
    text = response.text.strip()
    if not text.startswith("Date,"):
        snippet = text[:200].replace("\n", " ")
        raise RuntimeError(f"Unexpected response from schedule form: {snippet}")

    reader = csv.DictReader(io.StringIO(text))
    rows = [dict(row) for row in reader]

    # The site performs substring matches for Instrument. Normalize variants like
    # HIRESr/HIRESb/HIRES to a single canonical instrument when filtering.
    canonical_requested = canonicalize_instrument(instrument).lower()
    exact_rows = []
    for row in rows:
        normalized = canonicalize_instrument(row.get("Instrument") or "")
        if normalized.lower() == canonical_requested:
            new_row = dict(row)
            new_row["Instrument"] = normalized
            exact_rows.append(new_row)
    return exact_rows


def write_draft_allocation(rows: Iterable[dict], output_path: str) -> Path:
    """Write queried schedule rows to a draft allocation CSV."""
    fieldnames = ["Date", "Time", "Dark", "TelNr", "Instrument", "Account", "PI", "Institution", "ProjCode", ""]
    sorted_rows = sorted(
        rows,
        key=lambda row: (
            (row.get("Date") or "").strip(),
            ((row.get("Time") or "").split("-")[0]).strip(),
            canonicalize_instrument(row.get("Instrument") or ""),
            (row.get("ProjCode") or "").strip(),
        ),
    )
    out_path = Path(output_path)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in sorted_rows:
            out_row = {key: (row.get(key) or "").strip() for key in fieldnames if key}
            out_row["Instrument"] = canonicalize_instrument(out_row.get("Instrument", ""))
            out_row[""] = ""
            writer.writerow(out_row)
    return out_path


def summarize_pis(rows: Iterable[dict]) -> Counter:
    """Count schedule entries per PI."""
    counts: Counter = Counter()
    for row in rows:
        pi = (row.get("PI") or row.get("Principal") or "").strip()
        if pi:
            counts[pi] += 1
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Query the Keck schedule form for PIs in a semester and instrument set."
    )
    parser.add_argument(
        "-s",
        "--semester",
        default=infer_current_semester(),
        help="Semester to query, such as 2026A or 2026B. Default: current semester.",
    )
    parser.add_argument(
        "-i",
        "--instruments",
        nargs="+",
        default=DEFAULT_INSTRUMENTS,
        help="Instrument names to query. Default: KPF KPF-CC HIRES",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Print all returned schedule rows in addition to the PI summary.",
    )
    parser.add_argument(
        "--draft-allocation",
        action="store_true",
        help="Write a draft allocation CSV from the queried rows.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Optional output filename for the draft allocation CSV.",
    )
    args = parser.parse_args()

    print(f"Semester: {args.semester}")
    start_date, end_date = semester_date_range(args.semester)
    print(f"Date range: {start_date} to {end_date}\n")

    all_rows: List[dict] = []
    for instrument in args.instruments:
        rows = fetch_schedule_csv(args.semester, instrument)
        all_rows.extend(rows)
        pi_counts = summarize_pis(rows)
        print(f"=== {canonicalize_instrument(instrument)} ===")
        print(f"Rows returned: {len(rows)}")
        print(f"Unique PIs: {len(pi_counts)}")
        for pi, count in sorted(pi_counts.items(), key=lambda item: (-item[1], item[0].lower())):
            print(f"  {pi}: {count}")
        if not pi_counts:
            print("  No PI rows found")

        if args.full and rows:
            print("\nFull rows:")
            for row in rows:
                date_value = (row.get("Date") or "").strip()
                proj = (row.get("ProjCode") or "").strip()
                pi = (row.get("PI") or row.get("Principal") or "").strip()
                time_value = (row.get("Time") or "").strip()
                print(f"  {date_value} | {time_value} | {pi} | {proj}")
        print()

    if args.draft_allocation:
        output_name = args.output or f"allocation_hires_cps_{args.semester}.csv"
        out_path = write_draft_allocation(all_rows, output_name)
        print(f"Draft allocation CSV written to: {out_path}")


if __name__ == "__main__":
    main()
