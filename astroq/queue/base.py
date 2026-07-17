"""Base :class:`Queue` for telescope+instrument model.

Subclasses live in :mod:`astroq.queue.*`
"""

# Standard library imports
from __future__ import annotations

# Third-party imports
import numpy as np
import astropy.units as u

class Queue:
    """Abstract base class for a (telescope + instrument) queue.

    Subclasses are expected to populate the following attributes in
    ``__init__``

    Visit timing has two entry points: :meth:`visit_duration` (minutes,
    exposure + readout only; night/TTP) and :meth:`visit_seconds` (seconds,
    adds one ``slew_overhead_mean``; semester slot accounting).

    Attributes:
        observatory (astroplan.Observer): site-aware observer object.
        slew_rate (float): mean telescope slew rate, degrees/second.
        wrap_limit (float | None): azimuth wrap limit, degrees. ``None``
            means no wrap.
        wrap_states (list[tuple] | None): cable-wrap states for the
            state-aware TTP slew model; see :meth:`slew_fn_state`.
        nSlots (int): TTP slew-slot granularity (kept ``int`` and Pascal-cased
            because TTP consumes it as ``observatory.nSlots``).
        readout_time (float): detector readout time between successive shots
            of a single visit, seconds.
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

    observatory = None
    slew_rate: float
    wrap_limit: float | None = None
    wrap_states: list[tuple[str, float, float]] | None = None
    nSlots: int = 1
    readout_time: float
    slew_overhead_mean: float
    inaccessible_zones: tuple[tuple[float, float, float, float], ...] = ()

    def is_accessible(self, alt, az):
        """Boolean mask of telescope-accessible (alt, az) pairs.

        Encapsulates all hard telescope geometry by checking the input
        coordinates against :attr:`inaccessible_zones`. Used by
        :class:`astroq.access.Access` as the single per-cell pointing gate.

        The returned mask matches the broadcast shape of ``alt`` and ``az``.
        """
        alt = np.asarray(alt)
        az = np.asarray(az)
        excluded = np.zeros(np.broadcast(alt, az).shape, dtype=bool)
        for az_min, az_max, alt_min, alt_max in self.inaccessible_zones:
            excluded |= (
                (az >= az_min) & (az <= az_max) & (alt >= alt_min) & (alt <= alt_max)
            )
        return ~excluded

    def visit_duration(self, exptime_s, n_shots):
        """Total duration of one visit (n_shots shots), in *minutes*.

        Canonical formula: ``(exptime_s * n_shots + readout_time * (n_shots - 1)) / 60``.
        Consumed by the night planner (which wraps the result in a Quantity for
        ``TTPModel``) and ladder/script plot adapters.
        """
        return (exptime_s * n_shots + self.readout_time * (n_shots - 1)) / 60.0

    def _wrap_az(self, angle_deg):
        """Vectorized wrap-frame shift. ``wrap_limit=None`` means no shift."""
        if self.wrap_limit is None:
            return np.asarray(angle_deg)
        a = np.asarray(angle_deg) + (360 - self.wrap_limit)
        return np.where(a > 360, a - 360, a)

    def _short_az_sep(self, az_sep):
        """If telescope has no wrap, az slews never exceed 180 deg."""
        az_sep = np.asarray(az_sep)
        if self.wrap_limit is not None:
            return az_sep
        return np.where(az_sep > 180, 360 - az_sep, az_sep)

    #: Per-window sampling policy for :meth:`slew_fn`.
    _SLEW_SAMPLE_CADENCE_MIN = 30
    _SLEW_SAMPLES_PER_WINDOW_FLOOR = 3

    def _slew_window_samples(self, window_start, window_end):
        """Sample each slew window uniformly; return ``(times, M, n_samples)``."""
        M = len(window_start)
        win_dur_min = (window_end[0] - window_start[0]).to_value(u.min)
        n_samples = int(max(
            win_dur_min / self._SLEW_SAMPLE_CADENCE_MIN,
            self._SLEW_SAMPLES_PER_WINDOW_FLOOR,
        ))
        fracs = np.linspace(0.0, 1.0, n_samples)
        delta = window_end - window_start
        times_grid = window_start[:, None] + delta[:, None] * fracs[None, :]
        return times_grid.ravel(), M, n_samples

    def _altaz_pairs(self, times, coord_a, coord_b):
        """AltAz frames for ``coord_a`` and ``coord_b`` at ``times``."""
        altaz_a = self.observatory.altaz(times, coord_a, grid_times_targets=True)
        altaz_b = self.observatory.altaz(times, coord_b, grid_times_targets=True)
        return altaz_a, altaz_b

    def slew_fn(self, coord_a, coord_b, window_start, window_end):
        """Worst-case slew minutes per (pair, window).

        Implements the ``TTPModel`` slew_fn contract: ``coord_a`` and
        ``coord_b`` are pair-aligned 1-D ``SkyCoord`` arrays of length
        ``P`` (pair ``k`` is ``(coord_a[k], coord_b[k])``);
        ``window_start`` and ``window_end`` are 1-D ``Time`` arrays of
        length ``M`` giving the bounds of each slew slot. Returns an
        ``ndarray`` of shape ``(P, M)`` whose ``[k, m]`` entry is the
        worst-case slew time (minutes) from ``coord_a[k]`` to
        ``coord_b[k]`` over ``[window_start[m], window_end[m]]``.

        Each window is sampled internally at a cadence of
        :attr:`_SLEW_SAMPLE_CADENCE_MIN` minutes (floor:
        :attr:`_SLEW_SAMPLES_PER_WINDOW_FLOOR`) and reduced with ``max``.

        Assumes all windows have equal duration (TTPModel splits the
        night uniformly). For unequal windows this would need per-row
        sampling.
        """
        times, M, n_samples = self._slew_window_samples(window_start, window_end)
        altaz_a, altaz_b = self._altaz_pairs(times, coord_a, coord_b)
        az_sep = self._short_az_sep(
            np.abs(self._wrap_az(altaz_a.az.deg) - self._wrap_az(altaz_b.az.deg))
        )
        alt_sep = np.abs(altaz_a.alt.deg - altaz_b.alt.deg)
        tau = np.maximum(az_sep, alt_sep) / (60.0 * float(self.slew_rate))
        # tau shape: (P, M*n_samples). Reduce per window.
        return tau.reshape(-1, M, n_samples).max(axis=2)

    @property
    def n_states(self):
        """Number of cable-wrap states (1 for the legacy single-cut model)."""
        return 1 if self.wrap_states is None else len(self.wrap_states)

    @staticmethod
    def _encoder_az(az_deg, enc_min, enc_max):
        """Map sky azimuth(s) to encoder azimuth for one wrap state.

        Returns ``A + k*360`` (``k in {-1, 0, +1}``) that lands inside
        ``[enc_min, enc_max]``, else ``NaN``. Each wrap range here is narrower
        than 360 deg, so at most one winding is valid per state.
        """
        az = np.asarray(az_deg, dtype=float)
        out = np.full(az.shape, np.nan)
        for k in (-360.0, 0.0, 360.0):
            cand = az + k
            in_range = (cand >= enc_min) & (cand <= enc_max)
            out = np.where(np.isnan(out) & in_range, cand, out)
        return out

    def slew_fn_state(self, coord_a, coord_b, window_start, window_end):
        """State-aware worst-case slew minutes per (pair, window, si, sj).

        Like :meth:`slew_fn`, but resolves each pointing into every cable-wrap
        state in :attr:`wrap_states` and returns an array of shape
        ``(P, M, S, S)`` whose ``[k, m, si, sj]`` entry is the worst-case slew
        (minutes) from ``coord_a[k]`` (observed in state ``si``) to
        ``coord_b[k]`` (observed in state ``sj``) over window ``m``. The entry
        is ``NaN`` when either pointing is unreachable in its state anywhere in
        the window, so ``TTPModel.build_arcs`` drops that arc.

        The slew distance is the straight-line *encoder* azimuth difference
        (no wrap discontinuity), maxed against the elevation difference and
        divided by ``slew_rate``.
        """
        if self.wrap_states is None:
            raise RuntimeError(
                "slew_fn_state requires `wrap_states` to be defined on the queue"
            )
        states = self.wrap_states
        S = len(states)
        times, M, n_samples = self._slew_window_samples(window_start, window_end)
        altaz_a, altaz_b = self._altaz_pairs(times, coord_a, coord_b)
        az_a = altaz_a.az.deg
        az_b = altaz_b.az.deg
        alt_a = altaz_a.alt.deg
        alt_b = altaz_b.alt.deg

        # Encoder azimuth per state: shape (S, P, M*n_samples).
        enc_a = np.stack([self._encoder_az(az_a, lo, hi) for _, lo, hi in states])
        enc_b = np.stack([self._encoder_az(az_b, lo, hi) for _, lo, hi in states])

        # (Si, Sj, P, T); NaN propagates from unreachable windings.
        az_sep = np.abs(enc_a[:, None, :, :] - enc_b[None, :, :, :])
        alt_sep = np.abs(alt_a - alt_b)[None, None, :, :]
        tau = np.maximum(az_sep, alt_sep) / (60.0 * float(self.slew_rate))

        # Reduce per window with plain max so a single unreachable sample makes
        # the whole window infeasible (NaN). Shape -> (Si, Sj, P, M).
        P = az_a.shape[0]
        tau = tau.reshape(S, S, P, M, n_samples).max(axis=4)
        # Return (P, M, Si, Sj).
        return np.transpose(tau, (2, 3, 0, 1))

    def visit_seconds(self, exptime_s, n_exp):
        """Splan-canonical per-visit seconds.

        Per-visit elapsed time charged by the semester ILP:

        ``exptime_s * n_exp + readout_time * (n_exp - 1) + slew_overhead_mean``

        One ``slew_overhead_mean`` is charged per on-sky visit block; total
        semester slew time is accumulated elsewhere via ``t_visit_slots *
        n_intra_max * n_inter_max``.

        This is the single source of truth for the splan-style "visit
        seconds" calculation. Callers (:meth:`visit_slots`,
        :meth:`astroq.splan.SemesterPlanner._add_request_columns`) do
        their own slot conversion (round vs ceil) on top of the seconds.
        """
        return (
            exptime_s * n_exp
            + self.readout_time * (n_exp - 1)
            + self.slew_overhead_mean
        )

    def visit_slots(self, exptime_s, n_exp, slot_size_min):
        """Slots needed for one visit (scalar version of ``t_visit_slots``).

        Computes seconds via :meth:`visit_seconds` then rounds to slots.
        """
        total_s = self.visit_seconds(exptime_s, n_exp)
        slots = int(np.round(total_s / (slot_size_min * 60.0)))
        return max(1, slots)

    def write_starlist(self, *args, **kwargs):
        """Write tonight's starlist in the instrument-specific format.

        Concrete subclasses bind their module-level ``write_starlist`` function
        as this method.
        """
        raise NotImplementedError
