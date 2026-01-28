"""
Update backup OB JSON files with a dictionary of keyword-value overrides.

Reads a JSON file containing a list of OB objects, verifies that each update
keyword exists in the OB structure, applies updates to every OB, and writes
the result to an output file.

Usage
-----
CLI with updates from a JSON file::

    python -m astroq.queue.update_backup_OBs path/to/obs.json -u path/to/updates.json -o path/to/out.json

CLI with inline key-value pairs (VALUE is parsed as JSON; use quotes for strings)::

    python -m astroq.queue.update_backup_OBs path/to/obs.json -s Semid '"2026A_E476"' -s Semester '"2026A"' -s Progid '"E475"' -o out.json

Add a key under a section only if missing (dotted path, e.g. schedule.weather_band_3)::

    python -m astroq.queue.update_backup_OBs path/to/obs.json -a schedule.weather_band_3 true -o out.json

Programmatic::

    from astroq.queue.update_backup_OBs import update_backup_obs
    updated = update_backup_obs("obs.json", {"Semid": "2026A_E476", "Semester": "2026A"}, output_path="out.json")
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _collect_keys(obj: dict, keys: set[str] | None = None) -> set[str]:
    """Recursively collect all keys in nested dicts. Lists of dicts are traversed."""
    if keys is None:
        keys = set()
    if isinstance(obj, dict):
        for k, v in obj.items():
            keys.add(k)
            _collect_keys(v, keys)
    elif isinstance(obj, list):
        for item in obj:
            _collect_keys(item, keys)
    return keys


def _update_nested(obj: dict, updates: dict[str, object]) -> None:
    """Recursively update matching keys in nested dicts (mutates in place)."""
    if isinstance(obj, dict):
        for k, v in list(obj.items()):
            if k in updates:
                obj[k] = updates[k]
            else:
                _update_nested(v, updates)
    elif isinstance(obj, list):
        for item in obj:
            if isinstance(item, dict):
                _update_nested(item, updates)


def _add_if_missing(ob: dict, dotted_path: str, value: object) -> bool:
    """
    Add a key-value under a dotted path (e.g. schedule.weather_band_3) if missing.

    Traverses/create intermediate dicts as needed. If the final key already exists,
    does nothing and returns False. Otherwise adds it and returns True.
    """
    parts = dotted_path.strip().split(".")
    if len(parts) < 2:
        raise ValueError(
            f"Dotted path must have at least 'section.key', got {dotted_path!r}"
        )
    parent: dict = ob
    for p in parts[:-1]:
        if p not in parent:
            parent[p] = {}
        child = parent[p]
        if not isinstance(child, dict):
            raise ValueError(
                f"Cannot add under {dotted_path!r}: {p!r} is not a dict"
            )
        parent = child
    key = parts[-1]
    if key in parent:
        return False
    parent[key] = value
    return True


def update_backup_obs(
    json_path: str | Path,
    updates: dict[str, object] | None = None,
    adds: list[tuple[str, object]] | None = None,
    output_path: str | Path | None = None,
) -> list[dict]:
    """
    Load OBs from a JSON file, apply keyword updates and/or adds, optionally write output.

    Parameters
    ----------
    json_path : str or Path
        Path to the input JSON file (list of OB objects).
    updates : dict, optional
        Keyword-value pairs to apply. Each key must exist in at least one OB.
    adds : list of (dotted_path, value), optional
        Pairs to add under dotted paths (e.g. schedule.weather_band_3). Added only
        if the key is missing in that OB; existing keys are left unchanged.
    output_path : str or Path, optional
        Where to write the updated JSON. If None, nothing is written.

    Returns
    -------
    list of dict
        The updated list of OB objects.

    Raises
    ------
    KeyError
        If any update keyword is not found in any OB.
    ValueError
        If a dotted path is invalid or a path component is not a dict.
    """
    json_path = Path(json_path)
    with open(json_path, "r", encoding="utf-8") as f:
        obs = json.load(f)

    if not isinstance(obs, list):
        raise ValueError(f"Expected a JSON array of OBs, got {type(obs).__name__}")

    updates = updates or {}
    adds = adds or []

    if updates:
        all_keys: set[str] = set()
        for ob in obs:
            _collect_keys(ob, all_keys)
        missing = set(updates.keys()) - all_keys
        if missing:
            raise KeyError(
                f"Update keywords not found in any OB: {sorted(missing)}. "
                f"Available keys include: {sorted(all_keys)[:30]}..."
            )
        for ob in obs:
            _update_nested(ob, updates)

    for dotted_path, value in adds:
        for ob in obs:
            _add_if_missing(ob, dotted_path, value)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(obs, f, indent=2, ensure_ascii=False)

    return obs


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Update backup OB JSON with keyword-value overrides.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "json_file",
        type=Path,
        help="Path to input JSON file (list of OB objects).",
    )
    parser.add_argument(
        "-u",
        "--updates-json",
        type=Path,
        metavar="PATH",
        help="Path to JSON file with update dict, e.g. {\"Semid\": \"2026A_E476\", \"Semester\": \"2026A\"}.",
    )
    parser.add_argument(
        "-s",
        "--set",
        nargs=2,
        action="append",
        metavar=("KEY", "VALUE"),
        dest="sets",
        help="Repeat to set KEY=VALUE. VALUE is parsed as JSON (strings need quotes).",
    )
    parser.add_argument(
        "-a",
        "--add",
        nargs=2,
        action="append",
        metavar=("SECTION.KEY", "VALUE"),
        dest="adds",
        help="Repeat to add KEY under SECTION if missing (e.g. schedule.weather_band_3). VALUE parsed as JSON.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="PATH",
        default=None,
        help="Output path for updated JSON. Default: input path with '_updated' before extension.",
    )
    args = parser.parse_args()

    if args.updates_json is None and not args.sets and not args.adds:
        parser.error(
            "Provide at least one of: --updates-json, --set KEY VALUE, or --add SECTION.KEY VALUE."
        )

    if args.updates_json is not None and args.sets:
        parser.error("Use either --updates-json or --set, not both.")

    return args


def _main() -> int:
    args = _parse_args()

    if args.updates_json is not None:
        with open(args.updates_json, "r", encoding="utf-8") as f:
            updates = json.load(f)
        if not isinstance(updates, dict):
            print("Error: --updates-json must contain a JSON object.", file=sys.stderr)
            return 1
    elif args.sets:
        updates = {}
        for k, v_str in args.sets:
            try:
                updates[k] = json.loads(v_str)
            except json.JSONDecodeError:
                updates[k] = v_str
    else:
        updates = {}

    adds: list[tuple[str, object]] = []
    for path, v_str in args.adds or []:
        try:
            val = json.loads(v_str)
        except json.JSONDecodeError:
            val = v_str
        adds.append((path, val))

    out = args.output
    if out is None:
        p = args.json_file
        out = p.parent / f"{p.stem}_updated{p.suffix}"

    try:
        update_backup_obs(args.json_file, updates=updates, adds=adds, output_path=out)
    except (KeyError, ValueError, OSError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"Updated {args.json_file} -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
