"""PlotData context: shared semester-plan state for all plot functions.

``build_plot_data`` turns a solved ``SemesterPlanner`` into tidy per-request and
per-program DataFrames (plus the 3D access / forecast cubes that are genuinely
array-shaped). Plot functions consume these frames with groupby / merge / column
arithmetic instead of looping over per-object views.
"""

from __future__ import annotations

import os
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
import seaborn as sns

from astroq.plot._common import hours_per_night

MAP_NAMES = [
    "is_allocated",
    "is_custom",
    "is_altaz",
    "is_moon",
    "is_inter",
    "is_future",
    "is_observable_now",
]

# Columns shared by request_table and program_table so plots treat a per-request
# selection and a per-program (aggregate) selection uniformly.
_TABLE_COLUMNS = [
    "target",
    "program_code",
    "ra",
    "dec",
    "exptime",
    "tau_inter",
    "inactive",
    "total_observations_requested",
    "requested_visits",
    "total_requested_hours",
    "past_visits",
    "future_visits",
    "completion_pct",
    "splan_weight",
    "star_color",
    "program_color",
    "allow_mapview",
]


@dataclass(frozen=True)
class PlotSelection:
    """Which requests/programs to include in a plot."""

    unique_ids: frozenset[str] | None = None
    program_codes: frozenset[str] | None = None
    aggregate_by_program: bool = False


@dataclass
class RequestView:
    """Thin per-request handle for the webapp target index (id + display name)."""

    unique_id: str
    target: str
    program: str


@dataclass
class Selected:
    """Resolved selection: id list plus the frames sliced to those ids."""

    ids: list[str]
    is_program: bool
    table: pd.DataFrame
    cume_visits: pd.DataFrame
    cume_visits_pct: pd.DataFrame
    cume_slots: pd.DataFrame


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
    request_table: pd.DataFrame
    program_table: pd.DataFrame
    cume_visits: pd.DataFrame
    cume_visits_pct: pd.DataFrame
    cume_slots: pd.DataFrame
    cume_visits_program: pd.DataFrame
    cume_visits_pct_program: pd.DataFrame
    cume_slots_program: pd.DataFrame
    onsky_cadence: pd.DataFrame
    starmaps: dict[str, np.ndarray] = field(repr=False, default_factory=dict)
    program_starmaps: dict[str, np.ndarray] = field(repr=False, default_factory=dict)
    maps: dict[str, dict[str, np.ndarray]] = field(repr=False, default_factory=dict)

    @property
    def program_dict(self) -> dict[str, list[RequestView]]:
        """Requests grouped by program code (id + display name only)."""
        out: dict[str, list[RequestView]] = defaultdict(list)
        for prog, grp in self.request_table.groupby("program_code"):
            out[str(prog)] = [
                RequestView(unique_id=uid, target=target, program=str(prog))
                for uid, target in zip(grp.index, grp["target"])
            ]
        return dict(out)

    def select_all(self, *, aggregate_by_program: bool = False) -> PlotSelection:
        return PlotSelection(aggregate_by_program=aggregate_by_program)

    def select_program(self, program_code: str) -> PlotSelection:
        return PlotSelection(program_codes=frozenset({program_code}))

    def select_target(self, unique_id: str) -> PlotSelection:
        return PlotSelection(unique_ids=frozenset({str(unique_id)}))

    def starmap_for(self, id_: str, *, is_program: bool = False) -> np.ndarray:
        if is_program:
            return self.program_starmaps[str(id_)]
        return self.starmaps[str(id_)]

    def maps_for(self, uid: str) -> dict[str, np.ndarray]:
        return self.maps.get(
            str(uid),
            {n: np.zeros((self.n_nights, self.n_slots), dtype=bool) for n in MAP_NAMES},
        )

    def select(self, selection: PlotSelection) -> Selected:
        """Resolve a ``PlotSelection`` into id list + sliced tidy frames."""
        if selection.aggregate_by_program:
            ids = sorted(self._selected_program_codes(selection))
            ids = [c for c in ids if c in self.program_table.index]
            return Selected(
                ids=ids,
                is_program=True,
                table=self.program_table.loc[ids],
                cume_visits=self.cume_visits_program[ids],
                cume_visits_pct=self.cume_visits_pct_program[ids],
                cume_slots=self.cume_slots_program[ids],
            )
        ids = sorted(self._selected_uids(selection))
        return Selected(
            ids=ids,
            is_program=False,
            table=self.request_table.loc[ids],
            cume_visits=self.cume_visits[ids],
            cume_visits_pct=self.cume_visits_pct[ids],
            cume_slots=self.cume_slots[ids],
        )

    def _selected_uids(self, selection: PlotSelection) -> set[str]:
        all_uids = set(self.request_table.index)
        if selection.unique_ids is not None:
            return all_uids & {str(u) for u in selection.unique_ids}
        if selection.program_codes is not None:
            mask = self.request_table["program_code"].isin(selection.program_codes)
            return set(self.request_table.index[mask])
        return all_uids

    def _selected_program_codes(self, selection: PlotSelection) -> set[str]:
        if selection.program_codes is not None:
            return set(selection.program_codes)
        if selection.unique_ids is not None:
            uids = {str(u) for u in selection.unique_ids}
            present = [u for u in uids if u in self.request_table.index]
            return set(self.request_table.loc[present, "program_code"])
        return set(self.program_table.index)


def _program_colors(programs: np.ndarray) -> dict[str, str]:
    colors = sns.color_palette("deep", len(programs))
    rgb_strings = [
        f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})" for r, g, b in colors
    ]
    return dict(zip(programs, rgb_strings))


def _cume_matrix(ps, n_nights, group_col, columns, metric):
    """Cumulative visit count / slot sum per night, one column per group value."""
    tmp = ps.reset_index()  # index name 'd' -> column
    if metric == "visits":
        daily = tmp.groupby(["d", group_col]).size()
    else:
        daily = tmp.groupby(["d", group_col])["t_visit_slots"].sum()
    wide = (
        daily.unstack(group_col)
        .reindex(range(n_nights), fill_value=0)
        .reindex(columns=columns, fill_value=0)
        .fillna(0)
    )
    return wide.cumsum().astype(float)


def _cof_pct(cume, denom):
    """Cumulative % curve (COF): normalize by ``denom``, fall back to final value."""
    def curve(col):
        d = denom.get(col.name, 0)
        if d > 0:
            return (col / d * 100).round(2)
        last = col.iloc[-1]
        if last > 0:
            return (col / last * 100).round(2)
        return pd.Series(0.0, index=col.index)

    return cume.apply(curve)


def _completion_pct(last_cume, requested_visits):
    """Final completion % per request (round 3, matching legacy _cume_pct[-1])."""
    denom = requested_visits.where(requested_visits > 0)
    pct = (last_cume / denom * 100).round(3)
    pct = pct.where(requested_visits > 0, np.where(last_cume > 0, 100.0, 0.0))
    return pct


def _splan_weight(rf):
    """Numeric splan_weight per request; prefer 'splan_weight', fall back to 'weight'."""
    weight = pd.Series(np.nan, index=rf.index)
    if "splan_weight" in rf.columns:
        weight = pd.to_numeric(rf["splan_weight"], errors="coerce")
    if "weight" in rf.columns:
        fallback = pd.to_numeric(rf["weight"], errors="coerce")
        weight = weight.where(weight.notna(), fallback)
    return weight


def _build_request_table(semester_planner, cume_visits, program_colors, today_idx):
    """Vectorized per-request table (index = unique_id)."""
    rf = semester_planner.requests.copy()
    rf["unique_id"] = rf["unique_id"].astype(str)
    rf = rf.set_index("unique_id")
    queue = semester_planner.queue

    n_exp = rf["n_exp"].astype(int)
    n_intra_max = rf["n_intra_max"].astype(int)
    n_inter_max = rf["n_inter_max"].astype(int)
    exptime = rf["exptime"].astype(int)

    total_obs = (n_exp * n_intra_max * n_inter_max).astype(int)
    requested_visits = (n_intra_max * n_inter_max).astype(int)
    total_seconds = (
        total_obs * exptime
        + queue.readout_time * (n_exp - 1) * n_inter_max
        + queue.slew_overhead_mean * n_intra_max * n_inter_max
    )
    total_hours = total_seconds / 3600

    last_cume = cume_visits.iloc[-1].reindex(rf.index).fillna(0)
    inactive = rf["inactive"].astype(bool)
    total_obs = total_obs.mask(inactive, last_cume.astype(int))
    requested_visits = requested_visits.mask(inactive, last_cume.astype(int))
    inactive_hours = (
        total_obs * exptime + queue.slew_overhead_mean * total_obs
    ) / 3600
    total_hours = total_hours.mask(inactive, inactive_hours)

    if today_idx > 0:
        past_visits = cume_visits.iloc[today_idx - 1].reindex(rf.index).fillna(0)
    else:
        past_visits = pd.Series(0.0, index=rf.index)
    future_visits = last_cume - past_visits

    active_uids = set(semester_planner.requests_active["unique_id"].astype(str))
    program_code = rf["program_code"].astype(str)

    table = pd.DataFrame(index=rf.index)
    table["target"] = rf["target"]
    table["program_code"] = program_code
    table["ra"] = rf["ra"].astype(float)
    table["dec"] = rf["dec"].astype(float)
    table["exptime"] = exptime
    table["tau_inter"] = rf["tau_inter"].astype(int)
    table["inactive"] = inactive
    table["total_observations_requested"] = total_obs.astype(int)
    table["requested_visits"] = requested_visits.astype(int)
    table["total_requested_hours"] = total_hours.astype(float)
    table["past_visits"] = past_visits.astype(int)
    table["future_visits"] = future_visits.astype(int)
    table["completion_pct"] = _completion_pct(last_cume, requested_visits)
    table["splan_weight"] = _splan_weight(rf)
    table["program_color"] = program_code.map(program_colors)
    table["allow_mapview"] = table.index.to_series().isin(active_uids)
    return table


def _assign_star_colors(index, rgb_strings):
    """One random palette color per request (reproducible via module seed)."""
    n = len(rgb_strings)
    if n > 1:
        picks = np.random.randint(0, n, size=len(index))
    else:
        picks = np.zeros(len(index), dtype=int)
    return pd.Series([rgb_strings[i] for i in picks], index=index)


def _build_program_table(request_table, program_colors):
    """Aggregate request_table up to one row per program."""
    grouped = request_table.groupby("program_code")
    table = pd.DataFrame(index=grouped.size().index)
    table.index.name = None
    codes = table.index.astype(str)
    table["target"] = codes
    table["program_code"] = codes
    table["ra"] = 0.0
    table["dec"] = 0.0
    table["exptime"] = 0
    table["tau_inter"] = 0
    table["inactive"] = False
    table["total_observations_requested"] = grouped["requested_visits"].sum().astype(int)
    table["requested_visits"] = grouped["requested_visits"].sum().astype(int)
    table["total_requested_hours"] = grouped["total_requested_hours"].sum().astype(float)
    table["past_visits"] = grouped["past_visits"].sum().astype(int)
    table["future_visits"] = grouped["future_visits"].sum().astype(int)
    table["completion_pct"] = 0.0
    table["splan_weight"] = np.nan
    table["program_color"] = table.index.to_series().astype(str).map(program_colors)
    table["star_color"] = table["program_color"]
    table["allow_mapview"] = False
    return table


def _build_starmaps(forecast_df, request_table, uids, n_nights, n_slots):
    """Per-request forecast starmaps (n_slots, n_nights) via scattered assignment.

    Each scheduled (d, s) marks t_visit_slots consecutive slots. Built as one
    (n_uid, n_nights, n_slots) cube with vectorized assignment, then transposed
    per request.
    """
    row_of_uid = {uid: i for i, uid in enumerate(uids)}
    tv_of_uid = request_table["t_visit_slots"] if "t_visit_slots" in request_table else None
    cube = np.zeros((len(uids), n_nights, n_slots), dtype=int)

    if len(forecast_df):
        rows = forecast_df["unique_id"].map(row_of_uid)
        valid = rows.notna()
        r_idx = rows[valid].astype(int).to_numpy()
        d_idx = forecast_df.loc[valid, "d"].astype(int).to_numpy()
        s_idx = forecast_df.loc[valid, "s"].astype(int).to_numpy()
        if tv_of_uid is not None:
            tv = (
                forecast_df.loc[valid, "unique_id"]
                .map(tv_of_uid)
                .fillna(1)
                .astype(int)
                .to_numpy()
            )
        else:
            tv = np.ones(len(r_idx), dtype=int)
        max_tv = int(tv.max()) if len(tv) else 1
        for off in range(max_tv):
            m = off < tv
            cols = s_idx[m] + off
            inb = cols < n_slots
            cube[r_idx[m][inb], d_idx[m][inb], cols[inb]] = 1

    starmaps = {uid: cube[i].T for i, uid in enumerate(uids)}
    return starmaps, cube


def _build_maps(access, semester_planner, uids, n_nights, n_slots):
    """Per-request access cubes {map_name: (n_nights, n_slots)} for active requests."""
    active = semester_planner.requests_active["unique_id"].astype(str).tolist()
    row_of_active = {uid: i for i, uid in enumerate(active)}
    maps: dict[str, dict[str, np.ndarray]] = {}
    for uid in uids:
        row = row_of_active.get(uid)
        if row is None:
            continue
        maps[uid] = {name: access[name][row] for name in MAP_NAMES}
    return maps


def _build_onsky_cadence(cume_visits):
    """Tidy (unique_id, onsky_tau_inter): gaps between successive observed nights.

    Matches the legacy ``np.diff(np.where(np.diff(cume_observe) > 0)[0])``: gaps
    are computed from the nights where the cumulative visit count increases.
    """
    records: list[tuple[str, int]] = []

    def gaps(col):
        increases = np.where(np.diff(col.to_numpy()) > 0)[0]
        for gap in np.diff(increases):
            records.append((col.name, int(gap)))
        return None

    cume_visits.apply(gaps)
    return pd.DataFrame(records, columns=["unique_id", "onsky_tau_inter"])


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

    ps = semester_planner.timeline.copy()
    ps["unique_id"] = ps["unique_id"].astype(str)
    ps["program_code"] = ps["program_code"].astype(str)

    all_dates_array = semester_planner.access_obj.all_dates_array
    today_idx = semester_planner.access_obj.current_night_index
    n_nights = len(all_dates_array)
    slot_size = semester_planner.slot_size
    n_slots = int((24 * 60) / slot_size)
    slots_per_hour = 60 / slot_size

    uids = semester_planner.requests["unique_id"].astype(str).tolist()
    program_codes = semester_planner.requests["program_code"].astype(str).unique()
    program_colors = _program_colors(program_codes)

    colors = sns.color_palette("deep", len(program_codes))
    rgb_strings = [
        f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})" for r, g, b in colors
    ]

    # Per-request cumulative visits / slots (night x unique_id).
    cume_visits = _cume_matrix(ps, n_nights, "unique_id", uids, "visits")
    cume_slots = _cume_matrix(ps, n_nights, "unique_id", uids, "slots")

    request_table = _build_request_table(
        semester_planner, cume_visits, program_colors, today_idx
    )
    request_table["star_color"] = _assign_star_colors(request_table.index, rgb_strings)
    request_table["t_visit_slots"] = (
        semester_planner.requests.set_index(
            semester_planner.requests["unique_id"].astype(str)
        )["t_visit_slots"]
    )
    star_colors = request_table["star_color"].to_dict()

    requested_visits = request_table["requested_visits"]
    cume_visits_pct = _cof_pct(cume_visits, requested_visits.to_dict())

    program_table = _build_program_table(request_table, program_colors)
    program_list = program_table.index.tolist()
    cume_visits_program = _cume_matrix(
        ps, n_nights, "program_code", program_list, "visits"
    )
    cume_slots_program = _cume_matrix(
        ps, n_nights, "program_code", program_list, "slots"
    )
    cume_visits_pct_program = _cof_pct(
        cume_visits_program, program_table["requested_visits"].to_dict()
    )

    starmaps, cube = _build_starmaps(
        forecast_df, request_table, uids, n_nights, n_slots
    )
    program_starmaps = {}
    row_of_uid = {uid: i for i, uid in enumerate(uids)}
    for code, grp in request_table.groupby("program_code"):
        rows = [row_of_uid[u] for u in grp.index]
        # Float to match the legacy super-map (np.zeros defaults to float64),
        # which plotly serializes differently from an int heatmap.
        program_starmaps[str(code)] = cube[rows].sum(axis=0).T.astype(float)

    maps = _build_maps(access, semester_planner, uids, n_nights, n_slots)
    onsky_cadence = _build_onsky_cadence(cume_visits)

    programs_df = pd.read_csv(
        os.path.join(semester_planner.config.get("global", "workdir"), "programs.csv")
    )

    uid_to_row = {
        str(uid): int(idx)
        for idx, uid in enumerate(semester_planner.requests_active["unique_id"])
    }

    request_table = request_table.drop(columns=["t_visit_slots"])

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
        request_table=request_table,
        program_table=program_table,
        cume_visits=cume_visits,
        cume_visits_pct=cume_visits_pct,
        cume_slots=cume_slots,
        cume_visits_program=cume_visits_program,
        cume_visits_pct_program=cume_visits_pct_program,
        cume_slots_program=cume_slots_program,
        onsky_cadence=onsky_cadence,
        starmaps=starmaps,
        program_starmaps=program_starmaps,
        maps=maps,
    )
