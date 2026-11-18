#!/usr/bin/env python3
"""Mark ``request.csv`` rows fainter than a V-band cutoff as inactive."""

from __future__ import annotations

import argparse
import configparser
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd

logs = logging.getLogger(__name__)


def _resolve_data_path(workdir: str | Path, raw: str) -> Path:
    p = Path(raw)
    return p if p.is_absolute() else Path(workdir) / p


def apply_vmag_max(request_df: pd.DataFrame, vmag_max: float) -> tuple[pd.DataFrame, int]:
    """Return a copy with ``inactive=True`` on rows where ``Vmag > vmag_max``.

    Rows with missing or non-numeric ``Vmag`` are left unchanged.
    """
    if "Vmag" not in request_df.columns:
        raise ValueError("request.csv missing required column 'Vmag'")

    out = request_df.copy()
    if "inactive" not in out.columns:
        out["inactive"] = False

    vmag = pd.to_numeric(
        out["Vmag"].replace("none", np.nan).replace("None", np.nan),
        errors="coerce",
    )
    faint = vmag.notna() & (vmag > float(vmag_max))
    out.loc[faint, "inactive"] = True
    return out, int(faint.sum())


def filter_request_file(request_path: os.PathLike | str, vmag_max: float) -> int:
    """Apply ``vmag_max`` in place; return count of rows marked inactive by cut."""
    path = Path(request_path)
    df = pd.read_csv(path)
    filtered, n_faint = apply_vmag_max(df, vmag_max)
    filtered.to_csv(path, index=False)
    logs.info(
        "Vmag cut V<=%.2f: marked %d row(s) inactive in %s",
        vmag_max,
        n_faint,
        path,
    )
    return n_faint


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-cf",
        "--config-file",
        type=Path,
        help="config.ini; request.csv resolved from [data] request_file",
    )
    parser.add_argument(
        "--request-csv",
        type=Path,
        help="request.csv path (overrides -cf)",
    )
    parser.add_argument(
        "--vmag-max",
        type=float,
        required=True,
        help="Mark targets with Vmag strictly greater than this value inactive",
    )
    args = parser.parse_args(argv)

    if args.request_csv is not None:
        request_path = args.request_csv
    elif args.config_file is not None:
        cfg = configparser.ConfigParser()
        cfg.read(args.config_file)
        workdir = cfg.get("global", "workdir")
        request_path = _resolve_data_path(workdir, cfg.get("data", "request_file"))
    else:
        parser.error("one of --config-file or --request-csv is required")

    n_faint = filter_request_file(request_path, args.vmag_max)
    print(
        f"Applied Vmag<={args.vmag_max:.2f} cut: {n_faint} target(s) marked inactive "
        f"in {request_path}"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
