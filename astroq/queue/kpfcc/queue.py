"""KPF-CC :class:`Queue` subclass: telescope/instrument descriptor."""

import astroplan as apl

from astroq.queue.base import Queue


class KPFCC(Queue):
    """KPF-CC on Keck-I.

    KPF is permanently installed on Keck-I, so the telescope/instrument
    pairing is unique and this single class fully describes the queue.

    Pointing geometry references the Keck-I limits page:
    https://www2.keck.hawaii.edu/inst/common/TelLimits.html

    Geometry constants are duplicated from :class:`HIRESCPS` rather than shared
    via a mixin: the two queues may legitimately diverge on elevation policy.
    """

    slew_rate = 0.6  # deg/s
    wrap_limit = 270.0  # deg azimuth
    nSlots = 1  # TTP slew-slot granularity
    # Same readout/slew numbers as HIRES-CPS in production today; revisit if
    # KPF's measured detector readout or acquisition time diverges.
    readout_time = 45.0
    slew_overhead_mean = 60.0

    # Inaccessible (alt, az) boxes, degrees. (az_min, az_max, alt_min, alt_max).
    # See Queue.is_accessible. Duplicated from HIRESCPS; may diverge over time.
    inaccessible_zones = [
        (5.3, 146.2, 0.0, 33.3),  # Nasmyth deck obstruction
        (0.0, 360.0, 0.0, 18.0),  # below 18 deg elevation clamp
        (0.0, 360.0, 85.0, 90.0),  # above 85 deg elevation clamp
    ]

    # Constraints `Access` should compute for KPF-CC. Full set including
    # weather; `run_weather_loss` still gates whether ``compute_clear`` actually
    # samples losses or returns all-True.
    access_constraints = (
        "altaz", "future", "moon", "custom", "inter", "allocated", "clear",
    )

    def __init__(self):
        self.observatory = apl.Observer.at_site(
            "Keck Observatory", name="Keck", timezone="US/Hawaii"
        )

    def write_starlist(self, *args, **kwargs):
        from astroq.queue.kpfcc.starlist import write_starlist

        return write_starlist(*args, **kwargs)
