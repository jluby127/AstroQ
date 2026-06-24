"""Magellan :class:`Queue` subclass: telescope/instrument descriptor."""

import astroplan as apl

from astroq.queue.base import Queue


class Magellan(Queue):
    """Magellan twin 6.5m telescopes at Las Campanas Observatory.

    Pointing policy: no telescope tracks below 30 deg elevation in any
    azimuth. Additional azimuth-dependent limits can be added later via
    :attr:`inaccessible_zones`.
    """

    slew_rate = 1.0  # deg/s; placeholder — tune per instrument/site measurements
    wrap_limit = None  # no legacy single-cut azimuth wrap model
    nSlots = 1
    readout_time = 30.0  # seconds; placeholder
    slew_overhead_mean = 60.0  # seconds; splan-only mean slew + acquisition estimate

    # Inaccessible (alt, az) boxes, degrees. (az_min, az_max, alt_min, alt_max).
    inaccessible_zones = [
        (0.0, 360.0, -90.0, 29.999999),  # below 30 deg elevation (30 deg is OK)
    ]

    access_constraints = (
        "altaz", "future", "moon", "custom", "inter", "allocated", "clear",
    )

    def __init__(self):
        self.observatory = apl.Observer.at_site(
            "Las Campanas Observatory",
            name="Magellan",
            timezone="Chile/Continental",
        )

    def write_starlist(self, *args, **kwargs):
        from astroq.queue.magellan.starlist import write_starlist

        return write_starlist(*args, **kwargs)
