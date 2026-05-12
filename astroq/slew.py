"""
Simplified Keck-I-style slew time model for the in-house night ILP.

Numbers are calibrated to reproduce typical real-world Keck-I slew + acquisition
times: ~30-90 s for nearby fields, up to ~2 min for max separations. They are
NOT a high-fidelity replacement for the TTP Keck1 model -- the goal is to give
the night-plan ILP a per-pair upper bound on the time the telescope is
unavailable between two consecutive science exposures.

Constants
---------
SLEW_RATE_AZ : float
    Effective azimuth slew rate in deg/s (dome-limited, ~0.7 deg/s for Keck-I).
SLEW_RATE_EL : float
    Effective elevation slew rate in deg/s (~0.5 deg/s mechanical for Keck-I).
SETTLE_SECONDS : float
    Fixed per-slew overhead for dome settle + acquisition + autoguider lock.
"""

import numpy as np
from astropy.time import Time
import astropy.units as u

SLEW_RATE_AZ = 0.85
SLEW_RATE_EL = 0.65
SETTLE_SECONDS = 15.0


def slew_seconds_at(alts, azs):
    """Per-pair slew time in seconds from N-vectors of alt/az.

    All inputs are 1-D arrays of length N: alt[i] and az[i] are the
    instantaneous alt/az of target i (deg). The returned (N, N) matrix
    entry [i, j] is the slew time from target i to target j at this instant.

    The model is intentionally simple: a fixed acquisition/settle overhead
    plus the larger of (delta-az / az_rate, delta-alt / el_rate). It is
    calibrated to reproduce typical Keck-I + HIRES slew + acquisition times
    (~30-90 s for nearby fields, up to ~2 min for max separations).

    Args:
        alts (np.ndarray): (N,) array of altitudes in degrees at one instant.
        azs (np.ndarray): (N,) array of azimuths in degrees at one instant.

    Returns:
        np.ndarray: (N, N) slew time matrix, seconds. Diagonal is zero.
    """
    daz = np.abs(azs[:, None] - azs[None, :])
    daz = np.minimum(daz, 360.0 - daz)
    dalt = np.abs(alts[:, None] - alts[None, :])
    slew = SETTLE_SECONDS + np.maximum(daz / SLEW_RATE_AZ, dalt / SLEW_RATE_EL)
    np.fill_diagonal(slew, 0.0)
    return slew
