"""
Module for computing per-target, per-(night, slot) accessibility maps.

See ``Access.SUPPORTED_CONSTRAINTS`` and the ``compute_<name>`` methods for the
authoritative list of constraints actually applied. The ``Access`` instance is
stored on ``SemesterPlanner`` and reused for plotting.
"""

import logging
import os
from datetime import datetime, timedelta, timezone
from importlib.resources import files
from zoneinfo import ZoneInfo

import astropy as apy
import astropy.units as u
import astroplan as apl
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
from astropy.utils.iers import conf

conf.auto_max_age = None

logs = logging.getLogger(__name__)


def build_date_dictionary(semester_start_date, semester_length):
    """Single source of truth for the semester date grid.

    Args:
        semester_start_date (str): ``'YYYY-MM-DD'`` local civil date of night 0.
        semester_length (int): number of nights in the semester.

    """
    start = datetime.strptime(semester_start_date, "%Y-%m-%d")
    all_dates_array = [
        (start + timedelta(days=i)).strftime("%Y-%m-%d") for i in range(semester_length)
    ]
    all_dates_dict = {d: i for i, d in enumerate(all_dates_array)}
    return all_dates_array, all_dates_dict


def _observer_timezone(observer):
    """Return a tzinfo for ``observer`` (astroplan ``Observer``)."""
    tz = observer.timezone
    if isinstance(tz, str):
        return ZoneInfo(tz)
    return tz


def _localize_civil_noon(local_date, observer):
    """Local civil noon on ``YYYY-MM-DD`` at the observatory."""
    tz = _observer_timezone(observer)
    naive = datetime.strptime(local_date, "%Y-%m-%d").replace(
        hour=12, minute=0, second=0, microsecond=0
    )
    if hasattr(tz, "localize"):
        return tz.localize(naive)
    return naive.replace(tzinfo=tz)


def parse_utc_time(value):
    """Parse a CSV timestamp as UTC ``astropy.time.Time``."""
    t = Time(value, scale="utc")
    return t.utc


def observing_day_start(local_date, observer):
    """Start of observing night ``local_date`` (local civil noon) as UTC ``Time``."""
    local_noon = _localize_civil_noon(local_date, observer)
    return Time(local_noon)


def observing_day_end(local_date, observer):
    """End of observing night ``local_date`` (local civil noon on the next day)."""
    start = datetime.strptime(local_date, "%Y-%m-%d")
    next_date = (start + timedelta(days=1)).strftime("%Y-%m-%d")
    return observing_day_start(next_date, observer)


def observing_day_window(local_date, observer):
    """Half-open UTC interval ``[noon on local_date, noon on local_date + 1 day)``."""
    return observing_day_start(local_date, observer), observing_day_end(
        local_date, observer
    )


def observing_day_label_for_utc_time(utc_time, observer):
    """Local civil date label for the observing night containing ``utc_time``."""
    tz = _observer_timezone(observer)
    local_dt = parse_utc_time(utc_time).to_datetime(timezone=tz)
    shifted = local_dt - timedelta(hours=12)
    return shifted.strftime("%Y-%m-%d")


def observing_day_index_containing_time(utc_time, all_dates_array, observer):
    """Night index ``d`` whose observing-day window contains ``utc_time``."""
    label = observing_day_label_for_utc_time(utc_time, observer)
    mapping = {d: i for i, d in enumerate(all_dates_array)}
    if label not in mapping:
        raise ValueError(
            f"UTC time {utc_time} maps to observing day {label!r}, "
            f"outside semester grid {all_dates_array[0]}..{all_dates_array[-1]}"
        )
    return mapping[label]


def utc_time_to_observing_day_label(utc_time, all_dates_array, observer):
    """Observing-day label for ``utc_time``; raises if outside the semester grid."""
    d = observing_day_index_containing_time(utc_time, all_dates_array, observer)
    return all_dates_array[d]


def utc_interval_overlaps_observing_day(utc_start, utc_stop, local_date, observer):
    """True when UTC ``[utc_start, utc_stop]`` overlaps observing night ``local_date``."""
    win_start, win_end = observing_day_window(local_date, observer)
    start = parse_utc_time(utc_start)
    stop = parse_utc_time(utc_stop)
    return start < win_end and stop > win_start


class Access:
    """Accessibility maps for a collection of targets across a semester.

    Optional inputs (``allocation_file``, ``custom_file``,
    ``slots_needed_for_exposure``, weather) are opt-in via keyword. Passing
    ``None`` (or omitting them) makes the corresponding constraint a no-op.

    Attributes set by :meth:`build_access` / :meth:`build_windows`:
        first_available (astropy.time.Time): shape ``(ntargets, nnights)``,
            slot midpoint of the earliest observable slot per (target, night).
            JD ``0.0`` is a sentinel where ``has_observable`` is False; callers
            must gate on ``has_observable``. Consumed by ``nplan.run_ttp``.
        last_available (astropy.time.Time): same shape and sentinel convention;
            slot midpoint of the latest observable slot. Consumed by
            ``nplan.run_ttp``.
        has_observable (np.ndarray[bool]): shape ``(ntargets, nnights)``, True
            where at least one slot is observable. Consumed by
            ``nplan.run_ttp`` and ``test_sample``.

    Args:
        queue (astroq.queue.base.Queue): instrument/telescope queue.
            Provides ``observer``, ``is_accessible``, ``access_constraints``.
        request_frame (pandas.DataFrame): target list. Required columns:
            ``unique_id``, ``ra`` (deg), ``dec`` (deg). Optional column
            ``t_visit_slots`` (int >= 1) drives the multi-slot exposure
            dilation in :meth:`build_access`; if absent, defaults to 1
            per target (no dilation).
        semester_start_date (str): ``'YYYY-MM-DD'`` local civil date of night 0.
        semester_length (int): number of nights in the semester.
        slot_size (int): slot length in minutes; must divide 1440 evenly.

    Keyword Args:
        current_day (str, optional): today's local civil ``'YYYY-MM-DD'`` for the
            ``compute_future`` mask. Defaults to ``semester_start_date``.
        allocation_file (str, optional): path to ``allocation.csv``. ``None``
            treats every slot as allocated.
        custom_file (str, optional): path to ``custom.csv`` (PI windows).
            ``None`` skips custom-window restriction.
        run_weather_loss (bool, optional): if True, ``compute_clear`` samples
            historical weather losses; otherwise the cube is all-True.
        weather_loss_file (str, optional): override CSV for historical losses.
            Defaults to Maunakea data shipped with the package.

    Example (standalone):

        >>> import pandas as pd
        >>> from astroq.queue.hirescps.queue import HIRESCPS
        >>> from astroq.access import Access
        >>> df = pd.DataFrame({
        ...     "unique_id": ["a", "b"],
        ...     "ra": [10.0, 200.0],
        ...     "dec": [20.0, -10.0],
        ... })
        >>> acc = Access(HIRESCPS(), df, "2026-02-01", 184, 5)
        >>> rec = acc.build_access()
    """

    #: Canonical schema of constraint cubes packed into the recarray returned
    #: by :meth:`build_access`. Each ``Queue`` subclass declares which of these
    #: ``Access`` actually computes via ``Queue.access_constraints``;
    #: unlisted names default to all-True cubes.
    SUPPORTED_CONSTRAINTS = (
        "altaz", "future", "moon", "night", "custom", "inter", "allocated", "clear",
    )

    def __init__(
        self,
        queue,
        request_frame,
        semester_start_date,
        semester_length,
        slot_size,
        *,
        current_day=None,
        allocation_file=None,
        custom_file=None,
        run_weather_loss=False,
        weather_loss_file=None,
    ):
        self.queue = queue
        self.observatory = queue.observatory

        slot_size = float(slot_size)
        if 1440 % slot_size != 0:
            raise ValueError(
                f"slot_size={slot_size} must evenly divide 1440 minutes/day."
            )

        self.semester_start_date = semester_start_date
        self.semester_length = int(semester_length)
        self.slot_size = slot_size
        self.current_day = (
            current_day if current_day is not None else semester_start_date
        )

        self.all_dates_array, self.all_dates_dict = build_date_dictionary(
            self.semester_start_date, self.semester_length
        )

        # Fill optional per-row columns so downstream compute_* code can assume
        # they exist. Copy to avoid mutating caller's frame.
        rf = request_frame.copy()
        for col, default in (
            ("minimum_elevation", 0.0),
            ("minimum_moon_separation", 0.0),
            ("tau_inter", 0),
            # Multi-shot dilation: per-target full-visit duration in slots.
            # Default 1 (no dilation) so a stand-alone caller can skip the
            # full splan-style exposure accounting.
            ("t_visit_slots", 1),
        ):
            if col not in rf.columns:
                rf[col] = default
        self.request_frame = rf

        self.ntargets = len(self.request_frame)
        self.nnights = self.semester_length
        self.nslots = int(1440 / self.slot_size)
        self._access_shape = (self.ntargets, self.nnights, self.nslots)

        # Opt-in constraint inputs. None == constraint is a no-op.
        self.allocation_file = allocation_file
        self.custom_file = custom_file
        self.run_weather_loss = run_weather_loss

        self.slot_size_time = TimeDelta(self.slot_size * u.min)
        coords = apy.coordinates.SkyCoord(
            self.request_frame.ra * u.deg, self.request_frame.dec * u.deg, frame="icrs"
        )
        self.targets = apl.FixedTarget(name=self.request_frame.unique_id, coord=coords)

        # Observing-day grid: each row d spans local civil noon → next local noon.
        night0 = self.all_dates_array[0]
        self.daily_start = observing_day_start(night0, self.observatory)
        self.daily_end = observing_day_end(night0, self.observatory)
        self.timegrid = Time(
            np.arange(self.daily_start.jd, self.daily_end.jd, self.slot_size_time.jd),
            format="jd",
            location=self.observatory.location,
        )
        self.timegrid = self.timegrid[np.argsort(self.timegrid.sidereal_time("mean"))]

        jd_rows = []
        for local_date in self.all_dates_array:
            day_start = observing_day_start(local_date, self.observatory)
            jd_rows.append(
                day_start.jd
                + (np.arange(self.nslots) + 0.5) * self.slot_size_time.jd
            )
        self.slotmidpoints = Time(
            np.array(jd_rows),
            format="jd",
            location=self.observatory.location,
        )
        self.slotmidpoints_oneday = self.slotmidpoints[0]

        # compute_clear reads weather_loss_file only when run_weather_loss
        # is True; otherwise the cube is unconditionally all-True.
        self.weather_loss_file = weather_loss_file

    def observing_night_bounds(self, local_date):
        """UTC ``Time`` bounds for observing night ``local_date`` (local civil label)."""
        return observing_day_window(local_date, self.observatory)

    # ------------------------------------------------------------------
    # Adapter for the planner pipeline. Wires SemesterPlanner attributes
    # into the standalone constructor.
    # ------------------------------------------------------------------

    @classmethod
    def from_planner(cls, planner):
        """Construct an ``Access`` from a :class:`SemesterPlanner` instance.

        The planner is consumed for its current state and is not retained,
        avoiding any circular references between planner and access. Trivial
        scalar fields are read straight from ``planner.config``; derived
        ones (``semester_length``) and path-resolved ones
        (``allocation_file``, ``custom_file``) come from planner properties.
        """
        cfg = planner.config
        weather_loss_file = cfg.get(
            "semester", "weather_loss_file", fallback=None
        ) or None
        return cls(
            queue=planner.queue,
            request_frame=planner.requests_frame,
            semester_start_date=cfg.get("global", "semester_start_day"),
            semester_length=planner.semester_length,
            slot_size=cfg.getfloat("semester", "slot_size"),
            current_day=cfg.get("global", "current_day"),
            allocation_file=planner.allocation_file,
            custom_file=planner.custom_file,
            run_weather_loss=cfg.getboolean("semester", "run_weather_loss"),
            weather_loss_file=weather_loss_file,
        )

    # ------------------------------------------------------------------
    # Each compute_<name> returns a freshly-allocated boolean array of shape
    # (ntargets, nnights, nslots)
    # ------------------------------------------------------------------

    def compute_altaz(self):
        """Per-slot telescope pointing accessibility.

        Hard geometry comes from ``self.queue.is_accessible``; the PI-supplied
        ``minimum_elevation`` overlay is applied here.

        Altitudes are computed on a single 24h LST-sorted time grid for night 0
        and back-mapped to every (night, slot) via sidereal-time lookup. This
        is correct for sidereal targets only.
        """
        altazes = self.observatory.altaz(
            self.timegrid, self.targets, grid_times_targets=True
        )
        alts = altazes.alt.deg
        is_altaz0 = self.queue.is_accessible(alts, altazes.az.deg)
        is_altaz0 &= alts >= self.request_frame["minimum_elevation"].values[:, np.newaxis]

        x = self.timegrid.sidereal_time("mean").value
        x_new = self.slotmidpoints.sidereal_time("mean").value
        idx = np.clip(np.searchsorted(x, x_new, side="left"), 0, len(x) - 1)
        return is_altaz0[:, idx]

    def accessible_at(self, times):
        """Per-target telescope accessibility at arbitrary ``times``.

        Same accessibility gate as :meth:`compute_altaz` (hard ``is_accessible``
        pointing geometry plus the PI ``minimum_elevation`` overlay), but
        evaluated directly at the supplied ``times`` rather than on the semester
        slot grid. Used for ad-hoc windows such as the twilight backup sections.

        Args:
            times (astropy.time.Time): scalar or array of evaluation times.

        Returns:
            np.ndarray: boolean mask shaped ``(ntargets, ntimes)`` (rows aligned
            with ``self.request_frame``).
        """
        altazes = self.observatory.altaz(times, self.targets, grid_times_targets=True)
        alts = altazes.alt.deg
        mask = self.queue.is_accessible(alts, altazes.az.deg)
        mask &= alts >= self.request_frame["minimum_elevation"].values[:, np.newaxis]
        return mask

    def compute_future(self):
        """Mask out nights before ``self.current_day`` (local civil date)."""
        cube = np.ones(self._access_shape, dtype=bool)
        cube[:, : self.all_dates_dict[self.current_day], :] = False
        return cube

    def compute_moon(self):
        """Per-target moon-separation gating, evaluated once per night at slot 0."""
        moon = apy.coordinates.get_moon(
            self.slotmidpoints[:, 0], self.observatory.location
        )
        ang_dist = apy.coordinates.angular_separation(
            self.targets.ra.reshape(-1, 1),
            self.targets.dec.reshape(-1, 1),
            moon.ra.reshape(1, -1),
            moon.dec.reshape(1, -1),
        )
        min_sep = self.request_frame["minimum_moon_separation"].values * u.deg
        ok_per_night = ang_dist.to(u.deg) > min_sep[:, np.newaxis]
        return np.broadcast_to(
            ok_per_night[:, :, np.newaxis], self._access_shape
        ).copy()

    def compute_night(self):
        """Per-slot dark mask (sun below -12 deg, nautical twilight)."""
        sun_below = self.observatory.is_night(
            self.slotmidpoints, horizon=-12 * u.deg
        )  # (nnights, nslots)
        return np.broadcast_to(
            sun_below[np.newaxis, :, :], self._access_shape
        ).copy()

    def compute_inter(self):
        """Block ``tau_inter`` nights after each target's last observation.

        Reads ``past_date_last_observed`` off ``self.request_frame`` (a
        local civil observing-day label, ``""`` if the target has no past
        observations). Falls back to all-True if the column is absent (e.g.
        standalone-Access use case).
        """
        cube = np.ones(self._access_shape, dtype=bool)
        rf = self.request_frame
        if "past_date_last_observed" not in rf.columns:
            return cube
        for itarget in range(self.ntargets):
            row = rf.iloc[itarget]
            last = row["past_date_last_observed"]
            if last and row["tau_inter"] > 1 and last in self.all_dates_dict:
                start = self.all_dates_dict[last]
                stop = min(start + int(row["tau_inter"]), self.nnights)
                cube[itarget, start:stop, :] = False
        return cube

    def compute_custom(self):
        """PI-supplied per-star observability windows.

        Targets not listed in ``custom.csv`` are unrestricted (all-True). For
        listed targets the first window replaces the all-True default and
        subsequent windows are OR-ed in.
        """
        cube = np.ones(self._access_shape, dtype=bool)
        if self.custom_file is None:
            return cube
        if not os.path.exists(self.custom_file):
            logs.warning(
                "Custom times file not found: %s. Using no custom constraints.",
                self.custom_file,
            )
            return cube

        custom = pd.read_csv(self.custom_file)
        if len(custom) == 0:
            return cube

        starid_to_index = {
            uid: idx for idx, uid in enumerate(self.request_frame["unique_id"])
        }
        custom["start"] = custom["start"].apply(parse_utc_time)
        custom["stop"] = custom["stop"].apply(parse_utc_time)
        for _, row in custom.iterrows():
            if row["unique_id"] not in starid_to_index:
                continue
            mask = (self.slotmidpoints >= row["start"]) & (
                self.slotmidpoints <= row["stop"]
            )
            i = starid_to_index[row["unique_id"]]
            # First window for this star: replace the all-True default. Sentinel
            # is "still all-True"; subsequent windows OR in.
            cube[i] = mask if np.all(cube[i]) else cube[i] | mask
        return cube

    def compute_allocated(self):
        """Per-night-per-slot allocation mask, broadcast to all targets.

        With ``allocation_file is None`` every slot is treated as allocated
        (standalone-Access use case).
        """
        per_night = np.ones(self._access_shape[1:], dtype=bool)
        if self.allocation_file is not None:
            alloc = pd.read_csv(self.allocation_file)
            alloc["start"] = alloc["start"].apply(parse_utc_time)
            alloc["stop"] = alloc["stop"].apply(parse_utc_time)
            per_night = np.zeros_like(per_night)
            for _, row in alloc.iterrows():
                per_night |= (self.slotmidpoints >= row["start"]) & (
                    self.slotmidpoints <= row["stop"]
                )
        return np.broadcast_to(
            per_night[np.newaxis, :, :], self._access_shape
        ).copy()

    def compute_clear(self, weather_loss_file=None):
        """Weather-loss gating.

        When ``run_weather_loss=False`` returns an all-True cube. Otherwise
        simulates per-night losses from historical data and tiles the
        per-night mask to every target.
        """
        if not self.run_weather_loss:
            logs.info("Pretending weather is always clear!")
            return np.ones(self._access_shape, dtype=bool)
        if self.weather_loss_file is None:
            raise ValueError(
                "run_weather_loss=True requires weather_loss_file to be set explicitly."
            )

        logs.info("Running weather loss model.")
        self.get_loss_stats(weather_loss_file or self.weather_loss_file)
        per_night = self.simulate_weather_losses(covariance=0.14)
        return np.broadcast_to(
            per_night[np.newaxis, :, :], self._access_shape
        ).copy()

    # ------------------------------------------------------------------
    # Self-mutating orchestrators. The build_ prefix marks side effects.
    # ------------------------------------------------------------------

    def build_access(self):
        """Build the access recarray and populate the TTP windowing attributes.

        Dispatches ``compute_<name>`` for every ``name`` in
        ``self.queue.access_constraints``; unlisted names default to all-True.
        Side effect: calls :meth:`build_windows`, setting
        ``self.first_available``, ``self.last_available``,
        ``self.has_observable``.

        Returns:
            np.recarray of shape ``(ntargets, nnights, nslots)`` per field, with
            fields ``is_<name>`` for ``name in SUPPORTED_CONSTRAINTS`` plus
            ``is_observable_now`` (slot-level clearance, AND-reduce of all
            constraint cubes) and ``is_observable`` (start-of-exposure mask
            narrowed so a multislot exposure of
            ``request_frame['t_visit_slots'][uid]`` slots fits before
            night-end).
        """
        cubes = {
            name: np.ones(self._access_shape, dtype=bool)
            for name in self.SUPPORTED_CONSTRAINTS
        }
        for name in self.queue.access_constraints:
            if name not in self.SUPPORTED_CONSTRAINTS:
                raise ValueError(f"Unsupported access constraint: {name!r}")
            cubes[name] = getattr(self, f"compute_{name}")()

        is_observable_now = np.logical_and.reduce(
            [cubes[n] for n in self.SUPPORTED_CONSTRAINTS]
        )

        # is_observable[t, d, s] = "an e_val-slot exposure can START at slot s
        # and fit before night-end". AND in shifted copies of is_observable_now,
        # then zero the last e_val - 1 slots (the shift loop never writes them).
        t_visit_slots = self.request_frame["t_visit_slots"].astype(int).to_numpy()
        is_observable = is_observable_now.copy()
        for itarget in range(self.ntargets):
            e_val = int(t_visit_slots[itarget])
            if e_val == 1:
                continue
            for shift in range(1, e_val):
                is_observable[itarget, :, :-shift] &= is_observable_now[
                    itarget, :, shift:
                ]
            is_observable[itarget, :, -(e_val - 1):] = False

        self.build_windows(is_observable)

        fields = {f"is_{n}": cubes[n] for n in self.SUPPORTED_CONSTRAINTS}
        fields["is_observable_now"] = is_observable_now
        fields["is_observable"] = is_observable
        return np.rec.fromarrays(list(fields.values()), names=list(fields))

    def build_windows(self, is_observable):
        """Populate per-(target, night) first/last observable slot midpoints.

        Sets, each shape ``(ntargets, nnights)``:

        - ``self.has_observable``: bool, True where at least one slot is
          observable.
        - ``self.first_available``: astropy ``Time`` at the earliest observable
          slot midpoint. JD is a sentinel ``0.0`` where ``has_observable`` is
          False; callers must gate on ``has_observable``.
        - ``self.last_available``: astropy ``Time`` at the latest observable
          slot midpoint, same sentinel convention.

        Args:
            is_observable: ``(ntargets, nnights, nslots)`` bool cube; the
                ``is_observable`` field of :meth:`build_access`'s return.
        """
        ntargets, nnights, nslots = self._access_shape
        self.has_observable = is_observable.any(axis=2)
        first_idx = np.argmax(is_observable, axis=2)
        last_idx = nslots - 1 - np.argmax(is_observable[..., ::-1], axis=2)

        # Sentinel JD 0.0; has_observable is the truth source for masking.
        # Time's location must match slotmidpoints' so item-assignment works.
        prefill = np.zeros((ntargets, nnights))
        self.first_available = Time(
            prefill, format="jd", scale="utc",
            location=self.observatory.location,
        )
        self.last_available = Time(
            prefill.copy(), format="jd", scale="utc",
            location=self.observatory.location,
        )

        mask = self.has_observable
        night_idx = np.broadcast_to(np.arange(nnights), (ntargets, nnights))
        self.first_available[mask] = self.slotmidpoints[
            night_idx[mask], first_idx[mask]
        ]
        self.last_available[mask] = self.slotmidpoints[night_idx[mask], last_idx[mask]]

    def observability(self, is_observable):
        """Long-form (unique_id, d, s) triples for every observable cell.

        Args:
            is_observable: bool cube of shape ``(ntargets, nnights, nslots)``
                aligned with ``self.request_frame`` row order. Pass
                ``access.is_observable`` from :meth:`build_access`, or any
                equivalently-shaped mask (e.g. the slot-clearance variant).

        Returns:
            pandas.DataFrame with columns ``unique_id``, ``d``, ``s``. One row
            per True cell, ordered ascending by ``(itarget, d, s)``.
        """
        itarget, d, s = np.nonzero(is_observable)
        uid = self.request_frame["unique_id"].to_numpy()[itarget]
        return pd.DataFrame({"unique_id": uid, "d": d, "s": s})

    def get_loss_stats(self, weather_loss_file):
        """
        Gather the loss probabilities for each night in the semester from the saved historical weather data.
        """
        # ``weather_loss_file`` is normally a bare filename shipped with the
        # package (resolved via ``astroq.data``). Absolute paths are honored so
        # callers can override with site-specific historical data.
        if os.path.isabs(weather_loss_file):
            weather_csv = weather_loss_file
        else:
            weather_csv = files("astroq.data").joinpath(weather_loss_file)
        historical_weather_data = pd.read_csv(weather_csv)
        loss_stats_this_semester = []
        for i, item in enumerate(self.all_dates_array):
            ind = historical_weather_data.index[
                historical_weather_data["Date"] == self.all_dates_array[i][5:]
            ].tolist()[0]
            loss_stats_this_semester.append(
                historical_weather_data["% Total Loss"][ind]
            )
        self.loss_stats_this_semester = loss_stats_this_semester

    def simulate_weather_losses(self, covariance=0.14):
        """
        Simulate nights totally lost to weather using historical data

        Args:
            covariance (float): the added percent chance that tomorrow will be
            lost if today is lost

        Returns:
            is_clear (array): Trues represent clear nights, Falses represent
            weathered nights
        """
        previous_day_was_lost = False
        is_clear = np.ones(self._access_shape[1:], dtype=bool)
        for i in range(len(self.loss_stats_this_semester)):
            value_to_beat = self.loss_stats_this_semester[i]
            if previous_day_was_lost:
                value_to_beat += covariance
            roll_the_dice = np.random.uniform(0.0, 1.0)

            if roll_the_dice < value_to_beat:
                # the night is simulated a total loss
                is_clear[i] = np.zeros(is_clear.shape[1])  # Set all slots to False
                previous_day_was_lost = True
            else:
                previous_day_was_lost = False
        logs.info(
            f"Total nights simulated as weathered out: {np.sum(~np.any(is_clear, axis=1))} of {len(is_clear)} nights remaining."
        )
        return is_clear
