"""
Module for night-level observation planning and optimization.
Uses the vendored TTP MILP solver (``astroq.ttp.model.TTPModel``) to optimize
nightly observation sequences.
"""

from __future__ import annotations

import json
import logging
import os
from configparser import ConfigParser
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time, TimeDelta

from astroq.splan import SemesterPlanner
from astroq.ttp import model

logs = logging.getLogger(__name__)

NIGHT_PLANNER_H5_SCHEMA = 5

_OBJECT_TYPED_COLUMNS = [
    "coord",
    "first_available",
    "last_available",
    "t_visit",
    "tau_intra",
]

_SOLUTION_ATTRS = (
    ("night_start_jd", "night_start", "time"),
    ("night_end_jd", "night_end", "time"),
    ("solution_stats_json", "stats", "dict_json"),
)


def _load_config_from_text(text: str) -> ConfigParser:
    cfg = ConfigParser()
    cfg.read_string(text)
    return cfg


def _resolve_data_path(workdir: str | Path, raw: str) -> str:
    p = Path(raw)
    return str(p if p.is_absolute() else Path(workdir) / p)


def _json_to_native(x):
    if isinstance(x, np.ndarray):
        return _json_to_native(x.tolist())
    if isinstance(x, (list, tuple)):
        return [_json_to_native(v) for v in x]
    if isinstance(x, dict):
        return {k: _json_to_native(v) for k, v in x.items()}
    if isinstance(x, (np.integer, np.int64, np.int32)):
        return int(x)
    if isinstance(x, (np.floating, np.float64, np.float32)):
        return float(x)
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    return x


def _quantity_series_to_value(series, target_unit):
    if len(series) == 0:
        return np.array([], dtype=float)
    sample = series.iloc[0]
    if isinstance(sample, u.Quantity):
        return u.Quantity(series.values).to_value(target_unit)
    return np.asarray(series, dtype=float)


def _scalarize_requests_frame(rf):
    out = rf.drop(columns=_OBJECT_TYPED_COLUMNS, errors="ignore").copy()
    if "coord" not in rf.columns:
        return out
    out["ra"] = rf["coord"].apply(
        lambda c: float(c.ra.deg) if c is not None else np.nan,
    )
    out["dec"] = rf["coord"].apply(
        lambda c: float(c.dec.deg) if c is not None else np.nan,
    )
    if "first_available" in rf.columns:
        out["first_available_jd"] = rf["first_available"].apply(
            lambda t: float(t.jd) if t is not None else np.nan,
        )
    if "last_available" in rf.columns:
        out["last_available_jd"] = rf["last_available"].apply(
            lambda t: float(t.jd) if t is not None else np.nan,
        )
    if "t_visit" in rf.columns:
        out["t_visit_min"] = _quantity_series_to_value(rf["t_visit"], u.min)
    if "tau_intra" in rf.columns:
        out["tau_intra_hr"] = _quantity_series_to_value(rf["tau_intra"], u.hr)
    return out


def _restore_requests_frame_from_scalar_df(df):
    if "ra" not in df.columns or "dec" not in df.columns:
        raise ValueError(
            "night_planner.h5 requests_frame lacks ra/dec scalar columns; "
            "re-run plan-night with a current AstroQ build."
        )
    priority = (
        df["priority"].to_numpy()
        if "priority" in df.columns
        else np.full(len(df), 10.0)
    )
    return pd.DataFrame(
        {
            "unique_id": df["unique_id"].to_numpy(),
            "n_intra_max": df["n_intra_max"].to_numpy(),
            "priority": priority,
            "coord": SkyCoord(
                df["ra"].to_numpy() * u.deg,
                df["dec"].to_numpy() * u.deg,
                frame="icrs",
            ),
            "first_available": Time(df["first_available_jd"].to_numpy(), format="jd"),
            "last_available": Time(df["last_available_jd"].to_numpy(), format="jd"),
            "t_visit": df["t_visit_min"].to_numpy() * u.min,
            "tau_intra": df["tau_intra_hr"].to_numpy() * u.hr,
        }
    )


def _visit_duration_minutes(series):
    if len(series) == 0:
        return np.array([], dtype=float)
    sample = series.iloc[0]
    if isinstance(sample, u.Quantity):
        return u.Quantity(series.values).to_value(u.min)
    return np.asarray(series, dtype=float)


def _enrich_schedule_with_positions(sched, requests_frame):
    sched = sched.copy()
    lookup = pd.DataFrame(
        {
            "unique_id": requests_frame["unique_id"],
            "ra": requests_frame["coord"].apply(lambda c: float(c.ra.deg)),
            "dec": requests_frame["coord"].apply(lambda c: float(c.dec.deg)),
        }
    ).drop_duplicates("unique_id")
    for col in ("ra", "dec"):
        if col in sched.columns:
            sched = sched.drop(columns=[col])
    return sched.merge(lookup, on="unique_id", how="left").reset_index(drop=True)


def _prepare_schedule_for_hdf5(sched):
    sched = sched.drop(columns=_OBJECT_TYPED_COLUMNS, errors="ignore").copy()
    if "t_visit" in sched.columns:
        sched["t_visit_min"] = _visit_duration_minutes(sched["t_visit"])
    return sched


def _schedule_request_metadata(output_directory, requests_frame):
    path = os.path.join(output_directory, "request_selected.csv")
    if os.path.exists(path):
        meta = pd.read_csv(path)
    else:
        meta = requests_frame.copy()
    if meta.empty or "unique_id" not in meta.columns:
        return pd.DataFrame(columns=["unique_id", "exptime", "n_exp"])
    meta = meta[["unique_id", "exptime", "n_exp"]].drop_duplicates("unique_id")
    meta["exptime"] = pd.to_numeric(meta["exptime"], errors="coerce") / 60.0
    meta["n_exp"] = pd.to_numeric(meta["n_exp"], errors="coerce")
    return meta


def _restore_schedule_for_plots(sched, requests_frame, *, output_directory=None):
    sched = sched.copy()
    if "t_visit" not in sched.columns:
        if "t_visit_min" in sched.columns:
            sched["t_visit"] = sched["t_visit_min"].to_numpy(dtype=float)
        else:
            lookup = pd.DataFrame(
                {
                    "unique_id": requests_frame["unique_id"],
                    "t_visit": _visit_duration_minutes(requests_frame["t_visit"]),
                }
            ).drop_duplicates("unique_id")
            sched = sched.merge(lookup, on="unique_id", how="left")
    if "exptime" not in sched.columns or "n_exp" not in sched.columns:
        meta = _schedule_request_metadata(output_directory, requests_frame)
        if not meta.empty:
            for col in ("exptime", "n_exp"):
                if col in sched.columns:
                    sched = sched.drop(columns=[col])
            sched = sched.merge(meta, on="unique_id", how="left")
    return sched.reset_index(drop=True)


def _hygiene_selected_df(selected_df: pd.DataFrame) -> pd.DataFrame:
    out = selected_df.copy()
    out["n_intra_max"] = out["n_intra_max"].replace("None", np.nan).fillna(1)
    if "n_intra_min" in out.columns:
        out["n_intra_min"] = out["n_intra_min"].replace("None", np.nan).fillna(1)
    out["tau_intra"] = out["tau_intra"].replace("None", np.nan).fillna(0.0)
    return out


class NightPlanner:
    """TTP night planner: requires a saved ``semester_planner.h5`` from plan-semester."""

    def __init__(self, config_file):
        self._config_ini_text = Path(config_file).read_text()
        self.config = _load_config_from_text(self._config_ini_text)
        out = Path(self.config.get("global", "workdir")) / "outputs"
        self.semester_planner = SemesterPlanner.from_hdf5(out / "semester_planner.h5")
        self.queue = self.semester_planner.queue

    @property
    def semester_directory(self) -> str:
        return str(self.config.get("global", "workdir"))

    @property
    def output_directory(self) -> str:
        return os.path.join(self.semester_directory, "outputs")

    @property
    def current_day(self) -> str:
        return str(self.config.get("global", "current_day"))

    @property
    def allocation_file(self) -> str:
        return _resolve_data_path(
            self.semester_directory,
            self.config.get("data", "allocation_file"),
        )

    def _night_index(self) -> int:
        return self.semester_planner.all_dates_dict[self.current_day]

    def _night_bounds(self):
        try:
            start, stop = get_nightly_times_from_allocation(
                self.allocation_file,
                self.current_day,
            )
        except ValueError:
            logs.info(
                "No allocation for %s; skipping TTP.",
                self.current_day,
            )
            return None
        hours = np.round((stop.jd - start.jd) * 24, 3)
        logs.info("Time in Night for Observations: %s hours.", hours)
        return start, stop

    def _load_selected(self) -> pd.DataFrame | None:
        path = os.path.join(self.output_directory, "request_selected.csv")
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path} not found. Please run the scheduler first."
            )
        selected_df = pd.read_csv(path)
        if selected_df.empty:
            logs.info("No targets in %s; skipping TTP.", path)
            return None
        return _hygiene_selected_df(selected_df)

    def _write_outputs(
        self,
        tm: model.TTPModel,
        selected_df: pd.DataFrame,
        observation_start_time: Time,
    ) -> None:
        id_to_name = dict(zip(selected_df["unique_id"], selected_df["starname"]))
        tm.schedule["starname"] = tm.schedule["unique_id"].map(
            lambda uid: id_to_name.get(uid, "NO MATCHING NAME")
        )
        self.solution = tm

        observers_path = self.output_directory
        on_sky = tm.schedule[~tm.schedule["is_anchor"]]
        scheduled_df = on_sky[on_sky["scheduled"]].sort_values("order")
        extras_df = on_sky[~on_sky["scheduled"]]

        observe_order_file = os.path.join(
            observers_path,
            f"ObserveOrder_{self.current_day}.txt",
        )
        rows = []
        for _, row in scheduled_df.iterrows():
            ts = TimeDelta(row["t_start"] * 60, format="sec") + observation_start_time
            rows.append(
                {
                    "unique_id": str(row["unique_id"]),
                    "Target": row["starname"],
                    "StartExposure": str(ts)[11:16],
                }
            )
        for _, row in extras_df.iterrows():
            rows.append(
                {
                    "unique_id": str(row["unique_id"]),
                    "Target": row["starname"],
                    "StartExposure": "24:00",
                }
            )
        pd.DataFrame(rows).to_csv(observe_order_file, index=False)

        self.queue.write_starlist(
            selected_df,
            tm.schedule,
            observation_start_time,
            [],
            str(self.current_day),
            observers_path,
            all_active_requests=self.semester_planner.requests_frame,
            past_history=self.semester_planner.past_history,
        )

    def run_ttp(self):
        """Run TTP for tonight's ``request_selected.csv`` targets."""
        os.makedirs(self.output_directory, exist_ok=True)
        bounds = self._night_bounds()
        observation_start_time, observation_stop_time = bounds

        # grab first and last available slots for each target
        sp = self.semester_planner
        obs = sp.access_obj.observability()
        obs = obs[obs.d == self._night_index()]
        times = pd.DataFrame(
            {
                "s": np.arange(sp.access_obj.slotmidpoints.shape[1]),
                "time": sp.access_obj.slotmidpoints[self._night_index()],
            }
        )
        obs = pd.merge(obs, times, on="s")
        obs = (
            obs.sort_values(["unique_id", "s"])
            .groupby("unique_id", as_index=False)
            .agg(
                first_available=("time", "first"),
                last_available=("time", "last"),
            )
        )

        selected_df = self._load_selected()
        merged = pd.merge(selected_df, obs, on="unique_id")
        visit_min = self.queue.visit_duration(
            merged["exptime"].to_numpy(),
            merged["n_exp"].to_numpy(),
        )
        reqs = pd.DataFrame(
            {
                "unique_id": merged["unique_id"].to_numpy(),
                "n_intra_max": merged["n_intra_max"].to_numpy(dtype=int),
                "first_available": merged["first_available"].tolist(),
                "last_available": merged["last_available"].tolist(),
                "priority": np.full(len(merged), 10.0),
            }
        )
        reqs["coord"] = SkyCoord(
            merged["ra"].to_numpy() * u.deg,
            merged["dec"].to_numpy() * u.deg,
            frame="icrs",
        )
        reqs["t_visit"] = visit_min * u.min
        reqs["tau_intra"] = merged["tau_intra"].to_numpy(dtype=float) * u.hr

        tm = model.TTPModel(
            requests_frame=reqs,
            night_start=observation_start_time,
            night_end=observation_stop_time,
            slew_fn=self.queue.slew_fn,
            n_slots=self.queue.nSlots,
        )
        tm.build_nodes()
        tm.build_arcs()
        tm.build_model()
        tm.model.params.TimeLimit = self.config.getint("night", "max_solve_time")
        tm.model.params.MIPGap = self.config.getfloat("night", "max_solve_gap")
        tm.model.params.OutputFlag = int(
            self.config.getboolean("night", "show_gurobi_output")
        )
        tm.model.params.PreSolve = 2
        tm.model.params.MIPFocus = 1
        tm.model.params.Heuristics = 0.2
        tm.model.update()
        tm.run_model()
        if tm.model.SolCount == 0:
            logs.warning("TTP produced no schedule; skipping night-plan outputs.")
            return None
        tm.build_schedule()
        logs.info("\n" + tm.to_string())
        del tm.model



        self._write_outputs(tm, selected_df, observation_start_time)
        return True

    def to_hdf5(self, hdf5_path=None):
        import h5py

        if hdf5_path is None:
            hdf5_path = os.path.join(self.output_directory, "night_planner.h5")
        if os.path.exists(hdf5_path):
            os.remove(hdf5_path)

        solution = self.solution
        rf = _scalarize_requests_frame(solution.requests_frame)
        sched = _enrich_schedule_with_positions(
            _prepare_schedule_for_hdf5(solution.schedule),
            solution.requests_frame,
        )
        fmt = "fixed" if sched.empty else "table"
        sched.to_hdf(hdf5_path, key="solution_schedule", mode="a", format=fmt)
        rf_fmt = "fixed" if rf.empty else "table"
        rf.to_hdf(hdf5_path, key="solution_requests_frame", mode="a", format=rf_fmt)

        with h5py.File(hdf5_path, "a") as f:
            f.attrs["schema_version"] = NIGHT_PLANNER_H5_SCHEMA
            f.attrs["config_ini_text"] = self._config_ini_text
            f.attrs["semester_planner_h5_path"] = os.path.join(
                self.output_directory,
                "semester_planner.h5",
            )
            for hdf5_key, attr_path, dtype in _SOLUTION_ATTRS:
                obj = solution
                for part in attr_path.split("."):
                    obj = getattr(obj, part)
                if dtype == "dict_json":
                    f.attrs[hdf5_key] = json.dumps(_json_to_native(obj))
                elif dtype == "time":
                    f.attrs[hdf5_key] = obj.jd

        return hdf5_path

    @classmethod
    def from_hdf5(cls, hdf5_path):
        import h5py

        with h5py.File(hdf5_path, "r") as f:
            schema = int(f.attrs.get("schema_version", 0))
            if schema != NIGHT_PLANNER_H5_SCHEMA:
                raise ValueError(
                    f"night_planner.h5 schema_version={schema} is unsupported "
                    f"(expected {NIGHT_PLANNER_H5_SCHEMA}). Re-run plan-night."
                )
            if "solution_schedule" not in f:
                raise AttributeError("solution_schedule not found in HDF5 file")
            config_ini_text = f.attrs["config_ini_text"]
            semester_planner_h5_path = f.attrs["semester_planner_h5_path"]
            solution_attrs = {
                key: (f.attrs[key], dtype)
                for key, _, dtype in _SOLUTION_ATTRS
            }

        instance = cls.__new__(cls)
        instance._config_ini_text = config_ini_text
        instance.config = _load_config_from_text(config_ini_text)

        if not os.path.exists(semester_planner_h5_path):
            raise FileNotFoundError(
                f"semester_planner.h5 not found at {semester_planner_h5_path}"
            )
        instance.semester_planner = SemesterPlanner.from_hdf5(semester_planner_h5_path)
        instance.queue = instance.semester_planner.queue

        solution_schedule = pd.read_hdf(hdf5_path, key="solution_schedule")
        solution_requests_disk = pd.read_hdf(hdf5_path, key="solution_requests_frame")

        solution = model.TTPModel.__new__(model.TTPModel)
        for hdf5_key, attr_name, dtype in _SOLUTION_ATTRS:
            raw, _ = solution_attrs[hdf5_key]
            if dtype == "dict_json":
                data = json.loads(raw)
                restored = {}
                for key, value in data.items():
                    restored[key] = np.array(value) if isinstance(value, list) else value
                setattr(solution, attr_name, restored)
            elif dtype == "time":
                setattr(solution, attr_name, Time(raw, format="jd"))

        solution.requests_frame = _restore_requests_frame_from_scalar_df(
            solution_requests_disk,
        )
        solution.schedule = _restore_schedule_for_plots(
            solution_schedule,
            solution.requests_frame,
            output_directory=instance.output_directory,
        )
        if "ra" not in solution.schedule.columns:
            solution.schedule = _enrich_schedule_with_positions(
                solution.schedule,
                solution.requests_frame,
            )

        queue = instance.queue
        solution.observer = queue.observer
        solution.slew_rate = queue.slew_rate
        solution.wrap_limit = queue.wrap_limit
        solution.readout_time = queue.readout_time
        solution.n_slots = queue.nSlots

        instance.solution = solution
        return instance


def get_nightly_times_from_allocation(allocation_file, current_day):
    """
    Extract start and stop times for a specific date from allocation.csv.

    Args:
        allocation_file (str): path to the allocation file
        current_day (str): the date to look for in YYYY-MM-DD format

    Returns:
       start_time (Time object): the start time of the allocation for the current day
       stop_time (Time object): the stop time of the allocation for the current day
    """
    allocated_times_frame = pd.read_csv(allocation_file)
    allocated_times_frame["start"] = allocated_times_frame["start"].apply(Time)
    allocated_times_frame["stop"] = allocated_times_frame["stop"].apply(Time)

    current_day_str = str(current_day)
    day_allocations = []
    for _, row in allocated_times_frame.iterrows():
        start_datetime = str(row["start"])[:10]
        if start_datetime == current_day_str:
            day_allocations.append(row)

    if not day_allocations:
        raise ValueError(f"No allocation found for date {current_day_str}")

    earliest_start = min(row["start"] for row in day_allocations)
    latest_stop = max(row["stop"] for row in day_allocations)
    return earliest_start, latest_stop
