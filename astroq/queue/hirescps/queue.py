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
    wrap_limit = 235.0  # deg azimuth (legacy single-cut plot line)
    # Cable-wrap states for the state-aware TTP slew model. The two windings
    # overlap in the west (sky az ~215-315 deg); see Queue.wrap_states.
    #   North wrap: encoder az [-145, 90]  -> sky az [0,90] U [215,360)
    #   South wrap: encoder az [90, 315]   -> sky az [90,315]
    wrap_states = [
        ("N", -145.0, 90.0),
        ("S", 90.0, 315.0),
    ]
    nSlots = 4  # TTP slew-slot granularity
    readout_time = 30.0  # seconds; per-shot detector readout
    slew_overhead_mean = (
        30.0  # seconds; mean per-visit slew + acquisition (splan-only estimate)
    )

    # Inaccessible (alt, az) boxes, degrees. (az_min, az_max, alt_min, alt_max).
    # A sky point is excluded iff it lies inside ANY box. See Queue.is_accessible.
    # Duplicated independently in HIRESCPS and KPFCC; the two queues may
    # legitimately diverge on elevation policy.

    # Note 2026A Because of the issues with the ropes on the Keck shutters, the Keck 1 bottom
    # shutter is only opening to a position of 10 degrees (normally 2 degrees). In this
    # configuration, vignetting begins at elevation 28 degrees. 
    inaccessible_zones = [
        (5.3, 146.2, 0.0, 33.3),  # Nasmyth deck obstruction
        (0.0, 360.0, -90.0, 28.0),  # below 28 degerees.
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
