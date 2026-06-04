"""
Module for computing the intersection of the various accessibility maps for all targets for the following constraints:
    - telescope pointing
    - telescope allocation
    - sky brightness
    - moon separation
    - internight cadence from past history
    - PI custom windows
    - simulated weather loss
    - enough time to complete the exposure tonight

The Access class is saved as an attribute of the splan object and used again in plotting.
"""

import logging
import os
from datetime import datetime, timedelta

import astropy as apy
import astropy.units as u
import astroplan as apl
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
from astropy.utils.iers import conf

conf.auto_max_age = None

DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

logs = logging.getLogger(__name__)


def build_date_dictionary(semester_start_date, semester_length):
    """Single source of truth for the semester date grid.

    Args:
        semester_start_date (str): ``'YYYY-MM-DD'`` ISO date of night 0.
        semester_length (int): number of nights in the semester.

    """
    start = datetime.strptime(semester_start_date, "%Y-%m-%d")
    all_dates_array = [
        (start + timedelta(days=i)).strftime("%Y-%m-%d") for i in range(semester_length)
    ]
    all_dates_dict = {d: i for i, d in enumerate(all_dates_array)}
    return all_dates_array, all_dates_dict


class Access:
    """Accessibility maps for a collection of targets across a semester.

    Optional inputs (``allocation_file``, ``custom_file``, ``past_history``,
    ``slots_needed_for_exposure``, weather) are opt-in via keyword. Passing
    ``None`` (or omitting them) makes the corresponding constraint a no-op,

    Args:
        queue (astroq.queue.base.Queue): instrument/telescope queue.
            Provides ``observer``, ``is_accessible``, ``access_constraints``.
        request_frame (pandas.DataFrame): target list. Required columns:
            ``unique_id``, ``ra`` (deg), ``dec`` (deg). 
        semester_start_date (str): ``'YYYY-MM-DD'`` ISO date of night 0 (UTC).
        semester_length (int): number of nights in the semester.
        slot_size (int): slot length in minutes; must divide 1440 evenly.

    Keyword Args:
        current_day (str, optional): today's ``'YYYY-MM-DD'`` for the
            ``compute_future`` mask. Defaults to ``semester_start_date``.
        slots_needed_for_exposure (dict[str, int], optional): per-``unique_id``
            slot count for the multislot exposure dilation. Defaults to 1
            per target (no dilation).
        allocation_file (str, optional): path to ``allocation.csv``. ``None``
            treats every slot as allocated.
        custom_file (str, optional): path to ``custom.csv`` (PI windows).
            ``None`` skips custom-window restriction.
        past_history (dict, optional): per-``unique_id`` history records used
            by ``compute_inter``. ``None`` means no internight cadence blocking.
        run_weather_loss (bool, optional): if True, ``compute_clear`` samples
            historical weather losses; otherwise the cube is all-True.
        weather_loss_file (str, optional): override CSV for historical losses.
            Defaults to Maunakea data shipped with the package.

    Example (standalone):

        >>> import pandas as pd
        >>> from astroq.queue.hirescps import HIRESCPS
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
        "altaz", "future", "moon", "custom", "inter", "allocated", "clear",
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
        slots_needed_for_exposure=None,
        allocation_file=None,
        custom_file=None,
        past_history=None,
        run_weather_loss=False,
        weather_loss_file=None,
    ):
        self.queue = queue
        self.observatory = queue.observer

        if 1440 % slot_size != 0:
            raise ValueError(
                f"slot_size={slot_size} must evenly divide 1440 minutes/day."
            )

        self.semester_start_date = semester_start_date
        self.semester_length = int(semester_length)
        self.slot_size = int(slot_size)
        self.current_day = (
            current_day if current_day is not None else semester_start_date
        )

        self.start_date = Time(self.semester_start_date, format="iso", scale="utc")
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
        self.past_history = past_history if past_history is not None else {}
        self.run_weather_loss = run_weather_loss

        # Multishot dilation: per-target slot counts. Default 1 per target
        # (no dilation) so a standalone caller can skip exposure accounting.
        if slots_needed_for_exposure is None:
            slots_needed_for_exposure = {
                uid: 1 for uid in self.request_frame["unique_id"]
            }
        self.slots_needed_for_exposure_dict = slots_needed_for_exposure

        self.slot_size_time = TimeDelta(self.slot_size * u.min)
        coords = apy.coordinates.SkyCoord(
            self.request_frame.ra * u.deg, self.request_frame.dec * u.deg, frame="icrs"
        )
        self.targets = apl.FixedTarget(name=self.request_frame.unique_id, coord=coords)

        # Time grid for one night, first night of the semester
        self.daily_start = Time(self.start_date, location=self.observatory.location)
        self.daily_end = self.daily_start + TimeDelta(1.0, format="jd")
        self.timegrid = Time(
            np.arange(self.daily_start.jd, self.daily_end.jd, self.slot_size_time.jd),
            format="jd",
            location=self.observatory.location,
        )
        self.timegrid = self.timegrid[np.argsort(self.timegrid.sidereal_time("mean"))]

        # Slot midpoint for all nights in semester 2D array (slots, nights)
        self.slotmidpoints_oneday = (
            self.daily_start + (np.arange(self.nslots) + 0.5) * self.slot_size * u.min
        )
        days = np.arange(self.nnights) * u.day
        self.slotmidpoints = (
            self.slotmidpoints_oneday[np.newaxis, :] + days[:, np.newaxis]
        )

        self.DATADIR = DATADIR
        # Default historical-loss CSV for Keck/Mauna Kea. compute_clear only
        # uses this when run_weather_loss is True; harmless otherwise.
        self.weather_loss_file = (
            weather_loss_file
            if weather_loss_file is not None
            else os.path.join(self.DATADIR, "maunakea_weather_loss_data.csv")
        )

    # ------------------------------------------------------------------
    # Adapter for the planner pipeline. Wires SemesterPlanner attributes
    # into the standalone constructor.
    # ------------------------------------------------------------------

    @classmethod
    def from_planner(cls, planner, *, queue=None):
        """Construct an ``Access`` from a :class:`SemesterPlanner` instance.

        The planner is consumed for its current state and is not retained;
        this avoids any circular references between planner and access.

        ``queue`` defaults to ``planner.queue``; pass it explicitly only when
        rehydrating from h5 before the planner's queue is populated.
        """
        return cls(
            queue=queue if queue is not None else planner.queue,
            request_frame=planner.requests_frame,
            semester_start_date=planner.semester_start_date,
            semester_length=planner.semester_length,
            slot_size=planner.slot_size,
            current_day=planner.current_day,
            slots_needed_for_exposure=planner.slots_needed_for_exposure_dict,
            allocation_file=planner.allocation_file,
            custom_file=planner.custom_file,
            past_history=planner.past_history,
            run_weather_loss=planner.run_weather_loss,
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

    def compute_future(self):
        """Mask out nights before ``self.current_day`` for every target."""
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

    def compute_inter(self):
        """Block ``tau_inter`` nights after each target's last observation."""
        cube = np.ones(self._access_shape, dtype=bool)
        for itarget in range(self.ntargets):
            row = self.request_frame.iloc[itarget]
            uid = row["unique_id"]
            if uid in self.past_history and row["tau_inter"] > 1:
                start = self.all_dates_dict[self.past_history[uid].date_last_observed]
                stop = min(start + row["tau_inter"], self.nnights)
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
        custom["start"] = custom["start"].apply(Time)
        custom["stop"] = custom["stop"].apply(Time)
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
            alloc["start"] = alloc["start"].apply(Time)
            alloc["stop"] = alloc["stop"].apply(Time)
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
        ``self.first_available``, ``self.last_available``, ``self.has_observable``.

        Returns:
            np.recarray of shape ``(ntargets, nnights, nslots)`` per field, with
            fields ``is_<name>`` for ``name in SUPPORTED_CONSTRAINTS`` plus
            ``is_observable_now`` (slot-level clearance, AND-reduce of all
            constraint cubes) and ``is_observable`` (start-of-exposure mask
            narrowed so a multislot exposure of ``slots_needed_for_exposure_dict[uid]``
            slots fits before night-end).
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
        is_observable = is_observable_now.copy()
        for itarget in range(self.ntargets):
            e_val = self.slots_needed_for_exposure_dict[
                self.request_frame.iloc[itarget]["unique_id"]
            ]
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

        - ``self.has_observable``: bool, True where at least one slot is observable.
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
        make_time = lambda jd: Time(
            jd, format="jd", scale="utc", location=self.observatory.location
        )
        self.first_available = make_time(prefill)
        self.last_available = make_time(prefill.copy())

        mask = self.has_observable
        night_idx = np.broadcast_to(np.arange(nnights), (ntargets, nnights))
        self.first_available[mask] = self.slotmidpoints[
            night_idx[mask], first_idx[mask]
        ]
        self.last_available[mask] = self.slotmidpoints[night_idx[mask], last_idx[mask]]

    def observability(self, access=None):
        """Long-form table of observable ``(unique_id, d, s)`` triples.

        Args:
            access: Optional recarray from :meth:`build_access` (built on the
                fly if ``None``).

        Returns:
            pandas.DataFrame with columns ``unique_id``, ``d`` (night index),
            ``s`` (slot index). One row per observable cell.
        """
        if access is None:
            access = self.build_access()
        ntargets, nnights, nslots = access.shape

        # specify indeces of 3D observability array
        itarget, inight, islot = np.mgrid[:ntargets, :nnights, :nslots]

        # define flat table to access maps
        df = pd.DataFrame(
            {
                "itarget": itarget.flatten(),
                "inight": inight.flatten(),
                "islot": islot.flatten(),
            }
        )
        df["is_observable"] = access.is_observable.flatten()
        df = pd.merge(
            self.request_frame[["unique_id"]].reset_index(drop=True),
            df,
            left_index=True,
            right_on="itarget",
        )
        namemap = {"starid": "unique_id", "inight": "d", "islot": "s"}
        df = df.query("is_observable").rename(columns=namemap)[namemap.values()]
        return df

    def get_loss_stats(self, weather_loss_file):
        """
        Gather the loss probabilities for each night in the semester from the saved historical weather data.
        """
        historical_weather_data = pd.read_csv(os.path.join(DATADIR, weather_loss_file))
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
            covariance (float): the added percent chance that tomorrow will be lost if today is lost

        Returns:
            is_clear (array): Trues represent clear nights, Falses represent weathered nights
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


def build_twilight_allocation_file(semester_planner):
    """
    Build an allocation.csv file where every night of the semester is allocated
    from evening to morning 12-degree twilight times.
    This is used exclusively by the football plot in the webapp.

    Args:
        semester_planner (SemesterPlanner): a semester planner object from splan.py

    Returns:
        twilight_file (str): Path to the created allocation.csv file
    """

    # Create the filename based on semester
    semester = (
        semester_planner.semester_start_date[:4] + semester_planner.semester_letter
    )
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    twilight_file = os.path.join(data_dir, f"{semester}_twilights.csv")

    # Check if file already exists
    if os.path.exists(twilight_file):
        return twilight_file

    # Create data directory if it doesn't exist
    os.makedirs(data_dir, exist_ok=True)

    # Use the planner's queue's observer (already an apl.Observer)
    observatory = semester_planner.queue.observer

    # Create allocation data
    allocation_data = []

    for date_str in semester_planner.all_dates_dict.keys():
        # Parse the date
        date = Time(date_str, format="iso", scale="utc")

        # Get 12-degree twilight times for this night
        evening_12 = observatory.twilight_evening_nautical(date, which="next")
        morning_12 = observatory.twilight_morning_nautical(date, which="next")

        # Add to allocation data
        allocation_data.append({"start": evening_12.isot, "stop": morning_12.isot})

    # Create DataFrame and save to CSV
    twilight_df = pd.DataFrame(allocation_data)
    twilight_df.to_csv(twilight_file, index=False)

    return twilight_file
