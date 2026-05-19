"""
Queue base class.

A :class:`Queue` represents a single, specific telescope + instrument
combination (e.g. HIRES-CPS on Keck-I; KPF-CC on Keck-I). It is the single
source of truth for:

- Site geometry (astroplan ``Observer``, slew rate, wrap).
- Inaccessible alt/az regions (``inaccessible_zones``), the single declarative
  spec replacing the old scattered ``nays_*`` / ``tel_min`` / ``tel_max`` /
  ``deckAzLim*`` / ``vigLim`` / ``zenLim`` constants and the per-subclass
  ``is_accessible`` / ``pointing_limits`` overrides.
- Per-instrument timing overheads (readout time, slew overhead).
- Visit-duration math used by both the semester planner (slot accounting)
  and the night planner (TTP MILP).
- Instrument-specific I/O (``write_starlist``).

The same Queue instance is shared by ``SemesterPlanner``, ``NightPlanner``,
and ``Access``. The TTP MILP in ``astroq.ttp.*`` does NOT depend on Queue;
it consumes the queue's primitive fields (``observer``, ``slew_rate``,
``wrap_limit``, ``readout_time``, ``nSlots``, ``inaccessible_zones``) as
explicit ``TTPModel`` kwargs. This keeps ``astroq.ttp`` a leaf module.

Concrete subclasses live in :mod:`astroq.queue.hirescps` and
:mod:`astroq.queue.kpfcc`. The factory :func:`astroq.queue.from_config` selects
the right subclass from the ``[global] queue`` config field.
"""

# Standard library imports
from __future__ import annotations

# Third-party imports
import numpy as np


class Queue:
    """Abstract base class for a (telescope + instrument) queue.

    Subclasses are expected to populate the following attributes in
    ``__init__`` and (optionally) override :meth:`is_accessible` if the
    geometry can't be expressed as alt/az rectangles.

    Attributes:
        observer (astroplan.Observer): site-aware observer object.
        slew_rate (float): mean telescope slew rate, degrees/second.
        wrap_limit (float | None): azimuth wrap limit, degrees. ``None``
            means no wrap.
        nSlots (int): TTP slew-slot granularity (kept ``int`` and Pascal-cased
            because TTP consumes it as ``observatory.nSlots``).
        readout_time (float): detector readout time between successive shots
            of a single visit, seconds. Canonical source for both
            ``visit_duration`` and splan's slot accounting.
        slew_overhead_mean (float): mean per-visit slew overhead used by
            the semester ILP as a constant estimate (since splan cannot know
            target ordering). Superseded at night-plan time by TTP's per-arc
            ``tau_slew`` tensor computed from the alt/az pointing geometry.
            Seconds.
        inaccessible_zones (list[tuple]): Boxes in (alt, az) space where the
            telescope cannot point. Each entry is
            ``(az_min, az_max, alt_min, alt_max)`` in degrees. A sky point is
            excluded iff it lies inside ANY box. Replaces the old per-subclass
            ``is_accessible`` overrides; the base-class implementation is
            generic and reads from this list.
    """

    observer = None
    slew_rate: float
    wrap_limit: float | None = None
    nSlots: int = 1
    readout_time: float
    slew_overhead_mean: float
    inaccessible_zones: list[tuple[float, float, float, float]] = []

    def is_accessible(self, alt, az):
        """Boolean mask of telescope-accessible (alt, az) pairs.

        Encapsulates all hard telescope geometry by checking the input
        coordinates against :attr:`inaccessible_zones`. Used by
        :class:`astroq.access.Access` as the single per-cell pointing gate.

        The returned mask matches the broadcast shape of ``alt`` and ``az``.
        Subclasses may override for non-rectangular geometries (e.g. a
        polygonal nasmyth shadow), but the default loop over rectangular
        zones suffices for all currently-supported telescopes.
        """
        alt = np.asarray(alt)
        az = np.asarray(az)
        excluded = np.zeros(np.broadcast(alt, az).shape, dtype=bool)
        for az_min, az_max, alt_min, alt_max in self.inaccessible_zones:
            excluded |= (
                (az >= az_min) & (az <= az_max)
                & (alt >= alt_min) & (alt <= alt_max)
            )
        return ~excluded

    def visit_duration(self, exptime_s, n_shots):
        """Total duration of one visit (n_shots shots), in *minutes*.

        Canonical formula: ``(exptime_s * n_shots + readout_time * (n_shots - 1)) / 60``.
        Consumed by ``TTPModel._visit_duration`` and the TTP plotly summary.
        """
        return (exptime_s * n_shots + self.readout_time * (n_shots - 1)) / 60.0

    def visit_seconds(self, exptime_s, n_exp, n_intra_max):
        """Splan-canonical per-visit seconds.

        Per-visit elapsed time charged by the semester ILP. Includes the
        raw exposure time, between-shot readouts, and a slew-overhead term
        that uses ``n_intra_max`` as the multiplier (preserved bug-for-bug
        from the historical formula; see note below).

        This is the single source of truth for the splan-style "visit
        seconds" calculation. Callers (:meth:`visit_slots`,
        :meth:`astroq.splan.SemesterPlanner._build_slots_required_dictionary`)
        do their own slot conversion (round vs ceil) on top of the seconds.

        Note: ``slew_overhead_mean * n_intra_max`` is the legacy formula
        used by splan and :meth:`visit_slots`. A per-visit slew should
        arguably be a single ``slew_overhead_mean`` rather than
        ``n_intra_max`` of them; this is flagged as a follow-up.
        """
        return (
            exptime_s * n_exp
            + self.readout_time * (n_exp - 1)
            + self.slew_overhead_mean * n_intra_max
        )

    def visit_slots(self, exptime_s, n_exp, slot_size_min, n_intra_max):
        """Slots needed for one entry in ``slots_needed_for_exposure_dict``.

        Computes seconds via :meth:`visit_seconds` then rounds to slots.
        """
        total_s = self.visit_seconds(exptime_s, n_exp, n_intra_max)
        slots = int(np.round(total_s / (slot_size_min * 60.0)))
        return max(1, slots)

    def write_starlist(self, *args, **kwargs):
        """Write tonight's starlist in the instrument-specific format.

        Concrete subclasses bind their module-level ``write_starlist`` function
        as this method.
        """
        raise NotImplementedError
