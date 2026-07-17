"""PlotData context: shared semester-plan state for all plot functions."""

from __future__ import annotations

import os
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Iterable

import numpy as np
import pandas as pd
import seaborn as sns

from astroq.plot._common import (
    cumulative_by_night,
    daily_visits_by_night,
    hours_per_night,
    _visit_counts_by_date,
    _visit_denominator,
)

MAP_NAMES = [
    "is_allocated",
    "is_custom",
    "is_altaz",
    "is_moon",
    "is_inter",
    "is_future",
    "is_observable_now",
]


@dataclass(frozen=True)
class PlotSelection:
    """Which requests/programs to include in a plot."""

    unique_ids: frozenset[str] | None = None
    program_codes: frozenset[str] | None = None
    aggregate_by_program: bool = False


@dataclass
class RequestView:
    """Per-request view for plot functions (replaces StarPlotter)."""

    unique_id: str
    target: str
    program: str
    inactive: bool
    ra: float
    dec: float
    exptime: int
    tau_inter: int
    total_observations_requested: int
    requested_visits: int
    total_requested_hours: float
    program_color_rgb: str
    star_color_rgb: str
    cume_observe: np.ndarray
    cume_observe_pct: np.ndarray
    observations_past: dict
    observations_future: dict
    starmap: np.ndarray
    maps: dict[str, np.ndarray]
    maps_names: list[str]
    allow_mapview: bool
    draw_lines: bool = False


@dataclass
class PlotData:
    """Shared plotting context built once per loaded semester planner."""

    semester_planner: Any
    requests: pd.DataFrame
    timeline: pd.DataFrame
    schedule: pd.DataFrame
    access: np.recarray
    uid_to_row: dict[str, int]
    nulltime: np.ndarray
    all_dates_array: list
    today_idx: int
    n_nights: int
    n_slots: int
    slot_size: float
    program_colors: dict[str, str]
    star_colors: dict[str, str]
    programs: pd.DataFrame
    _views_by_uid: dict[str, RequestView] = field(repr=False, default_factory=dict)
    _views_by_program: dict[str, RequestView] = field(repr=False, default_factory=dict)

    @property
    def program_dict(self) -> dict[str, list[RequestView]]:
        """Requests grouped by program code."""
        out: dict[str, list[RequestView]] = defaultdict(list)
        for view in self._views_by_uid.values():
            out[view.program].append(view)
        return dict(out)

    def select_all(self, *, aggregate_by_program: bool = False) -> PlotSelection:
        return PlotSelection(aggregate_by_program=aggregate_by_program)

    def select_program(self, program_code: str) -> PlotSelection:
        return PlotSelection(program_codes=frozenset({program_code}))

    def select_target(self, unique_id: str) -> PlotSelection:
        return PlotSelection(unique_ids=frozenset({str(unique_id)}))

    def starmap_for(self, uid: str) -> np.ndarray:
        return self._views_by_uid[str(uid)].starmap

    def maps_for(self, uid: str) -> dict[str, np.ndarray]:
        return self._views_by_uid[str(uid)].maps

    def views(self, selection: PlotSelection) -> list[RequestView]:
        if selection.aggregate_by_program:
            return self._program_views(selection)
        uids = self._selected_uids(selection)
        return [self._views_by_uid[uid] for uid in sorted(uids)]

    def _selected_uids(self, selection: PlotSelection) -> set[str]:
        all_uids = set(self._views_by_uid)
        if selection.unique_ids is not None:
            return all_uids & set(selection.unique_ids)
        if selection.program_codes is not None:
            return {
                uid
                for uid, view in self._views_by_uid.items()
                if view.program in selection.program_codes
            }
        return all_uids

    def _program_views(self, selection: PlotSelection) -> list[RequestView]:
        codes = sorted(self._selected_program_codes(selection))
        return [self._views_by_program[code] for code in codes if code in self._views_by_program]

    def _selected_program_codes(self, selection: PlotSelection) -> set[str]:
        if selection.program_codes is not None:
            return set(selection.program_codes)
        if selection.unique_ids is not None:
            return {
                self._views_by_uid[uid].program
                for uid in selection.unique_ids
                if uid in self._views_by_uid
            }
        return set(self._views_by_program)


def _build_starmap(
    semester_planner, forecast_df, uid: str, n_nights: int, n_slots: int
) -> np.ndarray:
    starmap = np.zeros((n_nights, n_slots), dtype=int)
    star_forecast = forecast_df[forecast_df["unique_id"] == str(uid)]
    if len(star_forecast) == 0:
        return starmap.T

    d_values = star_forecast["d"].values.astype(int)
    s_values = star_forecast["s"].values.astype(int)
    starmap[d_values, s_values] = 1

    rf = semester_planner.requests_active
    row = rf.loc[rf["unique_id"] == str(uid)]
    reserve_slots = int(row["t_visit_slots"].iloc[0]) if len(row) else 1
    for r in range(1, reserve_slots):
        starmap[d_values, s_values + r] = 1
    return starmap.T


def _request_stats(row, slot_size, queue):
    n_exp = int(row["n_exp"])
    n_intra_max = int(row["n_intra_max"])
    n_inter_max = int(row["n_inter_max"])
    exptime = int(row["exptime"])
    total_observations_requested = n_exp * n_intra_max * n_inter_max
    requested_visits = n_intra_max * n_inter_max
    total_requested_seconds = (
        total_observations_requested * exptime
        + queue.readout_time * (n_exp - 1) * n_inter_max
        + queue.slew_overhead_mean * n_intra_max * n_inter_max
    )
    total_requested_hours = total_requested_seconds / 3600
    return {
        "total_observations_requested": total_observations_requested,
        "requested_visits": requested_visits,
        "total_requested_hours": total_requested_hours,
        "total_requested_nights": total_requested_hours / hours_per_night,
    }


def _cume_pct(cume_observe, denominator, n_nights):
    if denominator > 0:
        return np.round(cume_observe / denominator * 100.0, 3)
    total_past_visits = int(cume_observe[-1]) if len(cume_observe) else 0
    if total_past_visits > 0:
        return np.round(cume_observe / total_past_visits * 100.0, 3)
    return np.zeros(n_nights)


def build_plot_data(semester_planner) -> PlotData:
    """Build shared plot context from a SemesterPlanner (replaces process_stars)."""
    access = semester_planner.access_record
    nulltime = np.array(1 - access["is_allocated"][0]).T

    forecast_df = semester_planner.schedule
    if forecast_df is None:
        forecast_df = pd.DataFrame(columns=["unique_id", "d", "s", "target"])
    else:
        forecast_df = forecast_df.copy()
        forecast_df["unique_id"] = forecast_df["unique_id"].astype(str)

    ps = semester_planner.timeline
    all_dates_array = semester_planner.access_obj.all_dates_array
    today_idx = semester_planner.access_obj.current_night_index
    n_nights = len(all_dates_array)
    slot_size = semester_planner.slot_size
    n_slots = int((24 * 60) / slot_size)

    programs = semester_planner.requests["program_code"].unique()
    colors = sns.color_palette("deep", len(programs))
    rgb_strings = [
        f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})" for r, g, b in colors
    ]
    program_colors = dict(zip(programs, rgb_strings))

    queue = semester_planner.queue
    views_by_uid: dict[str, RequestView] = {}
    star_colors: dict[str, str] = {}

    for _, row in semester_planner.requests.iterrows():
        uid = str(row["unique_id"])
        stats = _request_stats(row, slot_size, queue)
        daily_visits = daily_visits_by_night(
            ps, n_nights, group_col="unique_id", group_val=uid
        )
        cume_observe = np.cumsum(daily_visits)

        if row["inactive"]:
            stats["total_observations_requested"] = int(np.max(cume_observe))
            stats["requested_visits"] = int(np.max(cume_observe))
            stats["total_requested_hours"] = (
                stats["total_observations_requested"] * int(row["exptime"])
                + queue.slew_overhead_mean * stats["total_observations_requested"]
            ) / 3600

        denom = stats["requested_visits"]
        cume_observe_pct = _cume_pct(cume_observe, denom, n_nights)

        if len(rgb_strings) > 1:
            star_color = rgb_strings[np.random.randint(0, len(rgb_strings) - 1)]
        else:
            star_color = rgb_strings[0]
        star_colors[uid] = star_color

        try:
            target_idx = np.where(
                semester_planner.requests_active["unique_id"] == uid
            )[0][0]
            maps = {name: access[name][target_idx] for name in MAP_NAMES}
            allow_mapview = True
        except (IndexError, KeyError):
            maps = {
                name: np.zeros((n_nights, n_slots), dtype=bool) for name in MAP_NAMES
            }
            allow_mapview = False

        views_by_uid[uid] = RequestView(
            unique_id=uid,
            target=row["target"],
            program=str(row["program_code"]),
            inactive=bool(row["inactive"]),
            ra=float(row["ra"]),
            dec=float(row["dec"]),
            exptime=int(row["exptime"]),
            tau_inter=int(row["tau_inter"]),
            total_observations_requested=stats["total_observations_requested"],
            requested_visits=stats["requested_visits"],
            total_requested_hours=stats["total_requested_hours"],
            program_color_rgb=program_colors[str(row["program_code"])],
            star_color_rgb=star_color,
            cume_observe=cume_observe,
            cume_observe_pct=cume_observe_pct,
            observations_past=_visit_counts_by_date(
                ps,
                group_col="unique_id",
                group_val=uid,
                past=True,
                today_idx=today_idx,
                all_dates_array=all_dates_array,
            ),
            observations_future=_visit_counts_by_date(
                ps,
                group_col="unique_id",
                group_val=uid,
                past=False,
                today_idx=today_idx,
                all_dates_array=all_dates_array,
            ),
            starmap=_build_starmap(
                semester_planner, forecast_df, uid, n_nights, n_slots
            ),
            maps=maps,
            maps_names=list(MAP_NAMES),
            allow_mapview=allow_mapview,
        )

    views_by_program: dict[str, RequestView] = {}
    unique_programs = sorted({v.program for v in views_by_uid.values()})
    for prog_code in unique_programs:
        prog_views = [v for v in views_by_uid.values() if v.program == prog_code]
        prog_indices = [v.unique_id for v in prog_views]
        cume_observe = cumulative_by_night(
            ps, n_nights, group_col="program_code", group_val=prog_code, metric="visits"
        )
        max_value = sum(_visit_denominator(v) for v in prog_views)
        if max_value > 0:
            cume_observe_pct = np.round(cume_observe / max_value * 100, 2)
        else:
            total_past = int(cume_observe[-1])
            cume_observe_pct = (
                cume_observe / total_past * 100 if total_past > 0 else np.zeros(n_nights)
            )

        super_map = np.zeros(np.shape(prog_views[0].starmap))
        for v in prog_views:
            super_map += v.starmap

        combined_past: dict = {}
        for v in prog_views:
            for date, count in v.observations_past.items():
                combined_past[date] = combined_past.get(date, 0) + count

        ref = prog_views[0]
        views_by_program[prog_code] = RequestView(
            unique_id=prog_code,
            target=prog_code,
            program=prog_code,
            inactive=False,
            ra=0.0,
            dec=0.0,
            exptime=0,
            tau_inter=0,
            total_observations_requested=sum(_visit_denominator(v) for v in prog_views),
            requested_visits=sum(_visit_denominator(v) for v in prog_views),
            total_requested_hours=float(
                sum(v.total_requested_hours for v in prog_views)
            ),
            program_color_rgb=ref.program_color_rgb,
            star_color_rgb=ref.program_color_rgb,
            cume_observe=cume_observe,
            cume_observe_pct=cume_observe_pct,
            observations_past=combined_past,
            observations_future={},
            starmap=super_map,
            maps={},
            maps_names=list(MAP_NAMES),
            allow_mapview=False,
        )

    programs_df = pd.read_csv(
        os.path.join(semester_planner.config.get("global", "workdir"), "programs.csv")
    )

    uid_to_row = {
        str(uid): int(idx)
        for idx, uid in enumerate(semester_planner.requests_active["unique_id"])
    }

    return PlotData(
        semester_planner=semester_planner,
        requests=semester_planner.requests,
        timeline=ps,
        schedule=forecast_df,
        access=access,
        uid_to_row=uid_to_row,
        nulltime=nulltime,
        all_dates_array=all_dates_array,
        today_idx=today_idx,
        n_nights=n_nights,
        n_slots=n_slots,
        slot_size=slot_size,
        program_colors=program_colors,
        star_colors=star_colors,
        programs=programs_df,
        _views_by_uid=views_by_uid,
        _views_by_program=views_by_program,
    )
