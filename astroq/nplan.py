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
from astropy.table import QTable

from astroq.splan import SemesterPlanner
from astroq.ttp import model

logs = logging.getLogger(__name__)

NIGHT_PLANNER_H5_SCHEMA = 6

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


def _solution_to_requests_dataframe(solution: "model.TTPModel") -> pd.DataFrame:
    """Flatten ``solution.requests`` (QTable) into a scalar h5 frame."""
    r = solution.requests
    return pd.DataFrame(
        {
            "unique_id": np.asarray(r["unique_id"], dtype=object),
            "ra_deg": np.asarray(r["coord"].ra.deg, dtype=float),
            "dec_deg": np.asarray(r["coord"].dec.deg, dtype=float),
            "first_available_jd": np.asarray(r["first_available"].jd, dtype=float),
            "last_available_jd": np.asarray(r["last_available"].jd, dtype=float),
            "t_visit_min": np.asarray(r["t_visit"].to_value(u.min), dtype=float),
            "tau_intra_min": np.asarray(r["tau_intra"].to_value(u.min), dtype=float),
            "n_intra_max": np.asarray(r["n_intra_max"], dtype=int),
            "priority": np.asarray(r["priority"], dtype=float),
        }
    )


def _requests_dataframe_to_qtable(df: pd.DataFrame) -> QTable:
    """Rehydrate an h5 requests frame back to the QTable boundary shape."""
    return QTable(
        {
            "unique_id": np.asarray(df["unique_id"].to_numpy(), dtype=object),
            "coord": SkyCoord(
                df["ra_deg"].to_numpy() * u.deg,
                df["dec_deg"].to_numpy() * u.deg,
                frame="icrs",
            ),
            "first_available": Time(
                df["first_available_jd"].to_numpy(), format="jd", scale="utc"
            ),
            "last_available": Time(
                df["last_available_jd"].to_numpy(), format="jd", scale="utc"
            ),
            "t_visit": df["t_visit_min"].to_numpy() * u.min,
            "n_intra_max": df["n_intra_max"].to_numpy(dtype=int),
            "tau_intra": df["tau_intra_min"].to_numpy() * u.min,
            "priority": df["priority"].to_numpy(dtype=float),
        },
        copy=False,
    )


def _attach_request_metadata_to_schedule(
    sched: pd.DataFrame, output_directory: str | None
) -> pd.DataFrame:
    """Merge ``exptime`` (min) and ``n_exp`` into the schedule for plot adapters."""
    if output_directory is None:
        return sched.reset_index(drop=True)
    path = os.path.join(output_directory, "request_selected.csv")
    if not os.path.exists(path):
        return sched.reset_index(drop=True)
    meta = pd.read_csv(path)
    if meta.empty or "unique_id" not in meta.columns:
        return sched.reset_index(drop=True)
    meta = meta[["unique_id", "exptime", "n_exp"]].drop_duplicates("unique_id")
    meta["exptime"] = pd.to_numeric(meta["exptime"], errors="coerce") / 60.0
    meta["n_exp"] = pd.to_numeric(meta["n_exp"], errors="coerce")
    sched = sched.drop(
        columns=[c for c in ("exptime", "n_exp") if c in sched.columns]
    )
    return sched.merge(meta, on="unique_id", how="left").reset_index(drop=True)


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

        sp = self.semester_planner
        d = self._night_index()
        selected_df = self._load_selected()

        # SemesterPlanner.from_hdf5 only instantiates Access; populate the
        # accessibility cubes (and the derived first/last_available arrays)
        # before we slice them.
        if not hasattr(sp.access_obj, "first_available"):
            sp.access_obj.build_access()

        # Look up the row index of each selected target in the access
        # request_frame, then slice the precomputed (ntargets, nnights) Time
        # arrays for tonight.
        req_index = sp.access_obj.request_frame.set_index("unique_id").index
        row_idx = req_index.get_indexer(selected_df["unique_id"])
        first_available = sp.access_obj.first_available[row_idx, d]
        last_available = sp.access_obj.last_available[row_idx, d]

        visit_min = self.queue.visit_duration(
            selected_df["exptime"].to_numpy(),
            selected_df["n_exp"].to_numpy(),
        )
        requests = QTable(
            {
                "unique_id": selected_df["unique_id"].to_numpy(dtype=object),
                "coord": SkyCoord(
                    selected_df["ra"].to_numpy() * u.deg,
                    selected_df["dec"].to_numpy() * u.deg,
                    frame="icrs",
                ),
                "first_available": first_available,
                "last_available": last_available,
                "t_visit": np.asarray(visit_min, dtype=float) * u.min,
                "n_intra_max": selected_df["n_intra_max"].to_numpy(dtype=int),
                "tau_intra": selected_df["tau_intra"].to_numpy(dtype=float) * u.hr,
                "priority": np.full(len(selected_df), 10.0),
            },
            copy=False,
        )

        tm = model.TTPModel(
            requests=requests,
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
        # solution.schedule already has only scalar dtypes (no object cols).
        sched = solution.schedule.reset_index(drop=True)
        fmt = "fixed" if sched.empty else "table"
        sched.to_hdf(hdf5_path, key="solution_schedule", mode="a", format=fmt)

        rf = _solution_to_requests_dataframe(solution)
        rf_fmt = "fixed" if rf.empty else "table"
        rf.to_hdf(hdf5_path, key="solution_requests", mode="a", format=rf_fmt)

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
        solution_requests_disk = pd.read_hdf(hdf5_path, key="solution_requests")

        solution = model.TTPModel.__new__(model.TTPModel)
        for hdf5_key, attr_name, dtype in _SOLUTION_ATTRS:
            raw, _ = solution_attrs[hdf5_key]
            if dtype == "dict_json":
                data = json.loads(raw)
                restored = {}
                for key, value in data.items():
                    restored[key] = (
                        np.array(value) if isinstance(value, list) else value
                    )
                setattr(solution, attr_name, restored)
            elif dtype == "time":
                setattr(solution, attr_name, Time(raw, format="jd"))

        solution.requests = _requests_dataframe_to_qtable(solution_requests_disk)

        solution.schedule = _attach_request_metadata_to_schedule(
            solution_schedule, instance.output_directory
        )

        queue = instance.queue
        solution.observer = queue.observatory
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
