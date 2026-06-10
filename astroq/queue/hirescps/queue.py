"""HIRES-CPS :class:`Queue` subclass: telescope/instrument descriptor."""

from astroq.queue.base import Queue


class HIRESCPS(Queue):
    """HIRES-CPS on Keck-I.

    HIRES is permanently installed on Keck-I, so the telescope/instrument
    pairing is unique and this single class fully describes the queue.

    Pointing geometry references the Keck-I limits page:
    https://www2.keck.hawaii.edu/inst/common/TelLimits.html

    The upper elevation limit (85 deg) matches AstroQ's accessibility policy
    rather than TTP's historical Keck1 zenith limit (84 deg); the physical
    zenith limit is 88.9 deg.
    """

    slew_rate = 0.6  # deg/s; matches TTP Keck1 (6./10.)
    wrap_limit = 235.0  # deg azimuth
    nSlots = 4  # TTP slew-slot granularity
    readout_time = 45.0  # seconds; per-shot detector readout
    slew_overhead_mean = (
        60.0  # seconds; mean per-visit slew + acquisition (splan-only estimate)
    )

    # Inaccessible (alt, az) boxes, degrees. (az_min, az_max, alt_min, alt_max).
    # A sky point is excluded iff it lies inside ANY box. See Queue.is_accessible.
    # Duplicated independently in HIRESCPS and KPFCC; the two queues may
    # legitimately diverge on elevation policy.
    inaccessible_zones = [
        (5.3, 146.2, 0.0, 33.3),  # Nasmyth deck obstruction
        (0.0, 360.0, 0.0, 18.0),  # below 18 deg elevation clamp
        (0.0, 360.0, 85.0, 90.0),  # above 85 deg elevation clamp
    ]

    # Constraints `Access` should compute for HIRES-CPS. Omits ``"clear"``;
    # weather loss is handled separately and the cube defaults to all-True.
    access_constraints = (
        "altaz", "future", "moon", "custom", "inter", "allocated",
    )

    def __init__(self):
        import astroplan as apl

        self.observatory = apl.Observer.at_site(
            "Keck Observatory", name="Keck", timezone="US/Hawaii"
        )

    def write_starlist(self, *args, **kwargs):
        # Lazy import to keep `queue.py` cheap to import and avoid pulling in
        # numpy/astropy at module-load time when only the descriptor is needed.
        from astroq.queue.hirescps.starlist import write_starlist

        return write_starlist(*args, **kwargs)
