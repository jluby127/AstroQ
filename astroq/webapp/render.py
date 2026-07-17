"""
Shared HTML page builders for the AstroQ webapp and static archive export.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from datetime import datetime
from importlib.resources import files as _resource_files
from typing import Any, Optional

import jinja2
import numpy as np
import plotly.io as pio

import astroq.nplan as nplan
import astroq.plot as pl
import astroq.ttp.plot as tplot
from astroq.nplan import NightPlanner
from astroq.splan import SemesterPlanner

logs = logging.getLogger(__name__)

_TEMPLATE_DIR = str(_resource_files("astroq.webapp").joinpath("templates"))
_TEMPLATE_ENV = jinja2.Environment(
    loader=jinja2.FileSystemLoader(_TEMPLATE_DIR),
    autoescape=False,
    trim_blocks=True,
    lstrip_blocks=True,
)


@dataclass
class LoadedRun:
    """Planner data loaded from a run's outputs directory."""

    semester_planner: SemesterPlanner
    data_astroq: tuple
    semester_planner_timestamp: Optional[str]
    night_planner: Optional[NightPlanner] = None
    data_ttp: Any = None
    night_start_time: Any = None
    request_frame_path: Optional[str] = None


def _render(template_name: str, **context) -> str:
    template = _TEMPLATE_ENV.get_template(template_name)
    return template.render(**context)


def _fig_to_html(fig) -> str:
    return pio.to_html(fig, full_html=True, include_plotlyjs="cdn")


def load_planners_from_outputs(outputs_dir: str) -> LoadedRun:
    """
    Load semester and optional night planners from an outputs directory.

    Args:
        outputs_dir: Path to ``{workdir}/outputs/``.

    Returns:
        LoadedRun with semester planner (required) and night planner (optional).

    Raises:
        FileNotFoundError: if outputs_dir or semester_planner.h5 is missing.
        Exception: if semester planner fails to load.
    """
    if not os.path.isdir(outputs_dir):
        raise FileNotFoundError(f"Directory not found: {outputs_dir}")

    semester_planner_h5 = os.path.join(outputs_dir, "semester_planner.h5")
    night_planner_h5 = os.path.join(outputs_dir, "night_planner.h5")
    request_frame_path = os.path.join(outputs_dir, "request_selected.csv")

    if not os.path.exists(semester_planner_h5):
        raise FileNotFoundError(f"semester_planner.h5 not found in {outputs_dir}")

    semester_planner = SemesterPlanner.from_hdf5(semester_planner_h5)
    data_astroq = pl.process_stars(semester_planner)
    mtime = os.path.getmtime(semester_planner_h5)
    semester_planner_timestamp = datetime.fromtimestamp(mtime).strftime(
        "%Y-%m-%d %H:%M:%S"
    )

    night_planner = None
    data_ttp = None
    night_start_time = None
    try:
        night_planner = NightPlanner.from_hdf5(night_planner_h5)
        data_ttp = night_planner.solution
        night_start_time, _ = nplan.get_nightly_times_from_allocation(
            night_planner.allocation_file, night_planner.current_day
        )
    except Exception as e:
        logs.warning(
            "Failed to load night planner from %s: %s",
            night_planner_h5,
            e,
        )

    return LoadedRun(
        semester_planner=semester_planner,
        data_astroq=data_astroq,
        semester_planner_timestamp=semester_planner_timestamp,
        night_planner=night_planner,
        data_ttp=data_ttp,
        night_start_time=night_start_time,
        request_frame_path=request_frame_path,
    )


def load_planners_from_uptree(
    semester_code: str, date: str, band: str, uptree_path: str
) -> LoadedRun:
    """Load planners from ``{uptree}/{semester}/{date}/{band}/outputs/``."""
    outputs_dir = os.path.join(uptree_path, semester_code, date, band, "outputs")
    return load_planners_from_outputs(outputs_dir)


def resolve_outputs_dir(run_path: str, run_name: str | None = None) -> str:
    """Resolve a run directory to its ``outputs/`` path.

    Accepts:
    - ``run_path/outputs`` or ``run_path`` when it already contains
      ``semester_planner.h5`` (single-run mode; use flat ``/admin`` URLs).
    - ``run_path`` as a parent of multiple child runs, with ``run_name`` naming
      the child folder that contains ``outputs/`` (use ``/{run_name}/admin``).
    """
    run_path = os.path.abspath(run_path)

    def _outputs_with_h5(dir_path: str) -> str | None:
        h5 = os.path.join(dir_path, "semester_planner.h5")
        if os.path.isfile(h5):
            return dir_path
        return None

    if run_name is not None:
        for candidate in (
            os.path.join(run_path, run_name, "outputs"),
            os.path.join(run_path, run_name),
        ):
            resolved = _outputs_with_h5(candidate)
            if resolved is not None:
                return resolved
        raise FileNotFoundError(
            f"No semester_planner.h5 under {run_path!r}/{run_name!r} "
            f"(tried child/outputs and child folder itself)."
        )

    resolved = _outputs_with_h5(os.path.join(run_path, "outputs"))
    if resolved is not None:
        return resolved
    resolved = _outputs_with_h5(run_path)
    if resolved is not None:
        return resolved

    children = list_child_runs(run_path)
    if children:
        example = children[0]
        raise FileNotFoundError(
            f"{run_path!r} contains multiple runs ({', '.join(children)}). "
            f"Open /{{run_name}}/admin, e.g. /{example}/admin."
        )
    raise FileNotFoundError(
        f"No outputs directory found at {os.path.join(run_path, 'outputs')!r} "
        f"and no semester_planner.h5 in {run_path!r}. "
        f"Pass a run folder, its outputs/ directory, or a parent of run folders."
    )


def list_child_runs(parent: str) -> list[str]:
    """Child folder names under ``parent`` that contain ``outputs/semester_planner.h5``."""
    parent = os.path.abspath(parent)
    if not os.path.isdir(parent):
        return []
    runs = []
    for name in sorted(os.listdir(parent)):
        sub = os.path.join(parent, name)
        if not os.path.isdir(sub):
            continue
        if os.path.isfile(os.path.join(sub, "outputs", "semester_planner.h5")):
            runs.append(name)
    return runs


def is_parent_run_path(run_path: str) -> bool:
    """True when ``run_path`` is a parent of child runs, not a single run root."""
    run_path = os.path.abspath(run_path)
    if os.path.isfile(os.path.join(run_path, "outputs", "semester_planner.h5")):
        return False
    if os.path.isfile(os.path.join(run_path, "semester_planner.h5")):
        return False
    return bool(list_child_runs(run_path))


def route_context_from_planner(semester_planner: SemesterPlanner) -> tuple[str, str, str]:
    """Derive semester/date/band URL parts from a loaded planner config."""
    semester_code = semester_planner.config.get("global", "semester")
    date = semester_planner.config.get("global", "current_day")
    workdir = semester_planner.config.get("global", "workdir")
    band = os.path.basename(os.path.normpath(workdir))
    return semester_code, date, band


def build_admin_html(
    loaded: LoadedRun,
    semester_code: str,
    date: str,
    band: str,
    *,
    link_targets: bool = True,
) -> str:
    """Render the admin dashboard page."""
    all_stars = np.concatenate(list(loaded.data_astroq[0].values()))
    request_df = pl.get_request_frame(loaded.semester_planner, all_stars)
    if link_targets:
        request_table_html = pl.request_frame_to_html(
            request_df, semester_code, date, band
        )
    else:
        request_table_html = pl.request_frame_to_html(request_df)

    fig_cof1 = pl.get_cof(loaded.semester_planner, list(loaded.data_astroq[1].values()))
    fig_cof2 = pl.get_cof(
        loaded.semester_planner, list(loaded.data_astroq[1].values()), use_time=True
    )
    fig_completion_hist = pl.get_completion_histogram_by_weight(
        loaded.semester_planner, all_stars
    )
    fig_completion_scatter = pl.get_completion_vs_target_name(
        loaded.semester_planner, all_stars
    )
    fig_birdseye = pl.get_birdseye(
        loaded.semester_planner, loaded.data_astroq[2], list(loaded.data_astroq[1].values())
    )
    fig_football = pl.get_football(
        loaded.semester_planner, all_stars, use_program_colors=True
    )
    fig_tau_inter_line = pl.get_tau_inter_line(
        loaded.semester_planner, all_stars, use_program_colors=True
    )
    fig_timebar = pl.get_timebar(
        loaded.semester_planner, all_stars, use_program_colors=True
    )
    fig_timebar_by_program = pl.get_timebar_by_program(
        loaded.semester_planner, loaded.data_astroq[0]
    )
    fig_rawobs = pl.get_rawobs(
        loaded.semester_planner, all_stars, use_program_colors=True
    )

    figures_html = [
        _fig_to_html(fig_timebar),
        _fig_to_html(fig_timebar_by_program),
        _fig_to_html(fig_cof1),
        _fig_to_html(fig_cof2),
        _fig_to_html(fig_completion_hist),
        _fig_to_html(fig_completion_scatter),
        _fig_to_html(fig_birdseye),
        _fig_to_html(fig_rawobs),
        _fig_to_html(fig_tau_inter_line),
        _fig_to_html(fig_football),
    ]

    return _render(
        "admin.html",
        tables_html=[request_table_html],
        figures_html=figures_html,
        timestamp=loaded.semester_planner_timestamp,
    )


def build_nightplan_html(loaded: LoadedRun, band: str) -> str:
    """Render the night plan page. Raises ValueError if night planner is missing."""
    if loaded.data_ttp is None or loaded.night_planner is None:
        raise ValueError("No night planner data available")

    script_table_df = pl.get_script_plan(loaded.night_planner)
    ladder_fig = pl.get_ladder(loaded.data_ttp, loaded.night_start_time)
    slew_animation_fig = tplot.get_slew_animation_plotly(
        loaded.data_ttp,
        loaded.request_frame_path,
        animationStep=120,
        inaccessible_zones=loaded.night_planner.queue.inaccessible_zones,
    )
    slew_path_fig = tplot.plot_path_2D_interactive(
        loaded.data_ttp, night_start_time=loaded.night_start_time
    )

    figure_html_list = [
        pl.nightplan_table_to_html(
            script_table_df, table_id="script-table", page_size=100
        ),
        _fig_to_html(ladder_fig),
        _fig_to_html(slew_animation_fig),
        _fig_to_html(slew_path_fig),
    ]

    return _render(
        "nightplan.html",
        figure_html_list=figure_html_list,
        band=band,
    )


def build_program_html(
    loaded: LoadedRun,
    semester_code: str,
    date: str,
    band: str,
    program_code: str,
    *,
    link_targets: bool = True,
) -> str:
    """Render a program overview page."""
    if program_code not in loaded.data_astroq[0]:
        raise KeyError(f"Program {program_code} not found")

    program_stars = loaded.data_astroq[0][program_code]
    request_df = pl.get_request_frame(loaded.semester_planner, program_stars)
    if link_targets:
        request_table_html = pl.request_frame_to_html(
            request_df, semester_code, date, band
        )
    else:
        request_table_html = pl.request_frame_to_html(request_df)

    fig_cof = pl.get_cof(loaded.semester_planner, program_stars)
    fig_completion_hist = pl.get_completion_histogram_by_weight(
        loaded.semester_planner, program_stars
    )
    fig_birdseye = pl.get_birdseye(
        loaded.semester_planner, loaded.data_astroq[2], program_stars
    )
    fig_tau_inter_line = pl.get_tau_inter_line(loaded.semester_planner, program_stars)
    fig_football = pl.get_football(loaded.semester_planner, program_stars)
    fig_timebar = pl.get_timebar(
        loaded.semester_planner, program_stars, use_program_colors=True
    )
    fig_rawobs = pl.get_rawobs(loaded.semester_planner, program_stars)

    figures_html = [
        _fig_to_html(fig_timebar),
        _fig_to_html(fig_cof),
        _fig_to_html(fig_completion_hist),
        _fig_to_html(fig_birdseye),
        _fig_to_html(fig_rawobs),
        _fig_to_html(fig_tau_inter_line),
        _fig_to_html(fig_football),
    ]

    return _render(
        "semesterplan.html",
        programname=program_code,
        tables_html=[request_table_html],
        figures_html=figures_html,
        programs=[program_code],
        timestamp=loaded.semester_planner_timestamp,
    )


def build_star_html(
    loaded: LoadedRun, target: str, program_code: Optional[str] = None
) -> str:
    """Render a single-target page."""
    compare_target = target.lower().replace(" ", "")
    programs_to_search = (
        [program_code]
        if program_code and program_code in loaded.data_astroq[0]
        else loaded.data_astroq[0].keys()
    )

    for program in programs_to_search:
        for star_ind in range(len(loaded.data_astroq[0][program])):
            star_obj = loaded.data_astroq[0][program][star_ind]
            true_target = star_obj.target
            if true_target.lower().replace(" ", "") != compare_target:
                continue

            request_df = pl.get_request_frame(loaded.semester_planner, [star_obj])
            request_table_html = pl.request_frame_to_html(request_df)

            fig_cof = pl.get_cof(
                loaded.semester_planner, [loaded.data_astroq[0][program][star_ind]]
            )
            fig_birdseye = pl.get_birdseye(
                loaded.semester_planner, loaded.data_astroq[2], [star_obj]
            )
            fig_tau_inter_line = pl.get_tau_inter_line(
                loaded.semester_planner, [star_obj]
            )
            fig_football = pl.get_football(loaded.semester_planner, [star_obj])
            fig_rawobs = pl.get_rawobs(loaded.semester_planner, [star_obj])

            figures_html = [
                _fig_to_html(fig_cof),
                _fig_to_html(fig_birdseye),
                _fig_to_html(fig_rawobs),
                _fig_to_html(fig_tau_inter_line),
                _fig_to_html(fig_football),
            ]

            return _render(
                "star.html",
                target=true_target,
                tables_html=[request_table_html],
                figures_html=figures_html,
                timestamp=loaded.semester_planner_timestamp,
            )

    raise KeyError(
        f"Target {target} not found in programs {list(programs_to_search)}"
    )
