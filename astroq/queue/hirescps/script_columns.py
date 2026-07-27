"""Fixed-width column formatters for HIRES nightly script lines."""

import logging

import numpy as np

logs = logging.getLogger(__name__)

TARGET_NAME_WIDTH = 16
VMAG_WIDTH = 9
EXPOSURE_WIDTH = 9  # max 4 digits + "/" + 4 digits
METER_WIDTH = 4
DECKER_WIDTH = 2
NEXP_WIDTH = 2
CELL_WIDTH = 3
PRIORITY_WIDTH = 2


def format_vmag_token(vmag: float) -> str:
    if not np.isfinite(vmag):
        vmag = 15.0
    return f"vmag={round(float(vmag), 1)}".rjust(VMAG_WIDTH)


def format_exposure_token(exptime: int, maxtime: int) -> str:
    token = f"{int(exptime)}/{int(maxtime)}"
    if len(token) > EXPOSURE_WIDTH:
        logs.warning(
            "exposure %r exceeds %d-char field (max 9999/9999); emitting unpadded",
            token,
            EXPOSURE_WIDTH,
        )
        return token
    return token.rjust(EXPOSURE_WIDTH)


def format_meter_token(threshold: str) -> str:
    return str(threshold).rjust(METER_WIDTH)


def format_decker_token(decker: str) -> str:
    return str(decker).rjust(DECKER_WIDTH)


def format_nexp_token(n_exp: int) -> str:
    return f"{int(n_exp)}x".rjust(NEXP_WIDTH)


def format_cell_token(cell: str) -> str:
    return str(cell).strip().rjust(CELL_WIDTH)


def format_priority_token(priority: str) -> str:
    return str(priority).rjust(PRIORITY_WIDTH)


SECTION_HEADER_WIDTH = 96


def format_section_header(label: str, total_width: int = SECTION_HEADER_WIDTH) -> str:
    """MAGIQ section divider: X-padding around ``__label__`` (label must have no spaces)."""
    token = f"__{label}__"
    pad = total_width - len(token)
    left = pad // 2
    right = pad - left
    return "X" * left + token + "X" * right
