"""
Queue base class.

A :class:`Queue` represents a single, specific telescope + instrument
combination (e.g. HIRES-CPS on Keck-I; KPF-CC on Keck-I). It is the single
source of truth for:

- Site geometry (astroplan ``Observer``, slew rate, wrap, pointing/elevation
  limits, nasmyth deck obstruction).
- Per-instrument timing overheads (readout time, slew overhead).
- Visit-duration math used by both the semester planner (slot accounting)
  and the night planner (TTP MILP).
- Instrument-specific I/O (``write_starlist``).

The same Queue instance is shared by ``SemesterPlanner``, ``NightPlanner``,
``Access``, and the TTP MILP. The TTP code in ``astroq.ttp.*`` consumes a
Queue via duck typing as its ``observatory`` argument; the back-compat
property aliases (:attr:`slewrate`, :attr:`wrapLimitAngle`, :attr:`readOutTime`)
exist exclusively for that purpose.

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
    ``__init__`` and override the abstract methods below.

    Attributes:
        observer (astroplan.Observer): site-aware observer object.
        slew_rate (float): mean telescope slew rate, degrees/second.
        wrap_limit (float | None): azimuth wrap limit, degrees. ``None``
            means no wrap.
        nSlots (int): TTP slew-slot granularity (kept ``int`` and Pascal-cased
            because TTP consumes it as ``observatory.nSlots``).
        tel_min (float): hard lower elevation limit, degrees.
        tel_max (float): hard upper elevation limit, degrees.
        readout_time (float): detector readout time between successive shots
            of a single visit, seconds. Canonical source for both
            ``visit_duration`` and splan's slot accounting.
        slew_overhead (float): average per-visit slew overhead, seconds. Used
            only by splan's ``visit_slots``; the TTP MILP computes slew time
            from the pointing tensor directly.
    """

    observer = None
    slew_rate: float
    wrap_limit: float | None = None
    nSlots: int = 1
    tel_min: float
    tel_max: float
    readout_time: float
    slew_overhead: float

    def pointing_limits(self, az, unvignetted=True):
        """Return ``[lower, upper]`` elevation limits at azimuth ``az`` (deg).

        Used by the legacy TTP back-compat path (``ttp.star`` falls back to
        this if a star lacks ``First Available`` / ``Last Available``).
        """
        raise NotImplementedError

    def is_up(self, alt, az, unvignetted=True):
        """1/0 array marking which (alt, az) pairs are on the sky.

        Legacy TTP signature. Default impl is a thin wrapper around
        :meth:`is_accessible`; subclasses can override if they need the
        ``unvignetted`` knob.
        """
        return self.is_accessible(alt, az).astype(int)

    def is_accessible(self, alt, az):
        """Boolean mask of telescope-accessible (alt, az) pairs.

        Encapsulates *all* hard telescope geometry: deck/nasmyth obstruction,
        elevation clamps (``tel_min`` / ``tel_max``), wrap, vignette. Used by
        :class:`astroq.access.Access` as the single per-cell pointing gate.

        The mask matches the shape of the broadcast of ``alt`` and ``az``.
        """
        raise NotImplementedError

    def visit_duration(self, exptime_s, n_shots):
        """Total duration of one visit (n_shots shots), in *minutes*.

        Canonical formula: ``(exptime_s * n_shots + readout_time * (n_shots - 1)) / 60``.
        Consumed by :class:`astroq.ttp.star.star` (``expwithreadout``) and the
        TTP plotly summary in :mod:`astroq.ttp.model`.
        """
        return (exptime_s * n_shots + self.readout_time * (n_shots - 1)) / 60.0

    def visit_slots(self, exptime_s, n_exp, slot_size_min, n_intra_max):
        """Slots needed for one entry in ``slots_needed_for_exposure_dict``.

        Bug-for-bug compatible with the formula at
        ``splan._build_slots_required_dictionary`` prior to the queue refactor::

            overhead_s = readout_time * (n_exp - 1) + slew_overhead * n_intra_max
            slots      = max(1, round((exptime_s * n_exp + overhead_s) / (slot_size_min * 60)))

        Note: the ``slew_overhead * n_intra_max`` term is suspect for a
        per-visit quantity (per-visit slew should be a single ``slew_overhead``,
        not ``n_intra_max`` of them). Preserved here to keep schedules
        identical across the refactor; flagged as a follow-up.
        """
        overhead_s = self.readout_time * (n_exp - 1) + self.slew_overhead * n_intra_max
        total_s = exptime_s * n_exp + overhead_s
        slots = int(np.round(total_s / (slot_size_min * 60.0)))
        return max(1, slots)

    def write_starlist(self, *args, **kwargs):
        """Write tonight's starlist in the instrument-specific format.

        Concrete subclasses bind their module-level ``write_starlist`` function
        as this method.
        """
        raise NotImplementedError

    # --- TTP back-compat aliases ---------------------------------------------
    # Read-only properties so the canonical attribute name on the Queue
    # (``slew_rate``, ``wrap_limit``, ``readout_time``) can never drift from
    # the TTP-facing name (``slewrate``, ``wrapLimitAngle``, ``readOutTime``).

    @property
    def slewrate(self):
        return self.slew_rate

    @property
    def wrapLimitAngle(self):
        return self.wrap_limit

    @property
    def readOutTime(self):
        return self.readout_time

    @property
    def vigLim(self):
        """TTP plotter alias for the unvignetted lower elevation limit."""
        return self.tel_min

    @property
    def zenLim(self):
        """TTP plotter alias for the upper (zenith-side) elevation limit."""
        return self.tel_max
