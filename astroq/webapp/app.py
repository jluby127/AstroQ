"""
Web application module for AstroQ.
"""

# Standard library imports
import logging
import os

logs = logging.getLogger(__name__)

# Third-party imports
from flask import Flask, render_template, abort
from socket import gethostname

# Local imports
from astroq.webapp import render as webrender

running_on_keck_machines = False

app = Flask(__name__)

# Global variables to store loaded data
loaded_run = None
_uptree_path = None
_run_path = None


def load_data_for_path(semester_code, date, band, uptree_path_param):
    """
    Load data for a specific semester_code/date/band combination

    Args:
        semester_code (str): the semester code
        date (str): the date in YYYY-MM-DD format
        band (str): the band
        uptree_path_param (str): the path to the uptree directory

    Returns:
        tuple: (success, message)
    """
    global loaded_run
    try:
        loaded_run = webrender.load_planners_from_uptree(
            semester_code, date, band, uptree_path_param
        )
    except FileNotFoundError as e:
        loaded_run = None
        return False, str(e)
    except Exception as e:
        loaded_run = None
        return False, f"Error loading semester planner: {str(e)}"
    return True, "Data loaded successfully"


def load_data_for_run_path(run_path_param, run_name=None):
    """Load data from a run directory (single run or child of a parent -rp path)."""
    global loaded_run
    try:
        outputs_dir = webrender.resolve_outputs_dir(run_path_param, run_name=run_name)
        loaded_run = webrender.load_planners_from_outputs(outputs_dir)
    except FileNotFoundError as e:
        loaded_run = None
        return False, str(e)
    except Exception as e:
        loaded_run = None
        return False, f"Error loading semester planner: {str(e)}"
    return True, "Data loaded successfully"


def _require_run_path():
    if not _run_path:
        abort(
            404,
            description="Run-path routes require launching with -rp/--run_path.",
        )


def _single_run_mode():
    """Flat /admin URLs apply only when -rp points at one run, not a parent folder."""
    return _run_path and not webrender.is_parent_run_path(_run_path)


def _parent_run_mode():
    return _run_path and webrender.is_parent_run_path(_run_path)


VALID_BANDS = [
    "band1",
    "band2",
    "band3",
    "full-band1",
    "full-band2",
    "full-band3",
]


def _validate_band(band):
    if band not in VALID_BANDS:
        abort(
            400,
            description=(
                "Band must be 'band1', 'band2', 'band3', "
                "'full-band1', 'full-band2', or 'full-band3'"
            ),
        )


def _homepage_navigation_text():
    if _run_path:
        if _parent_run_mode():
            runs = webrender.list_child_runs(_run_path)
            run_lines = "\n".join(f"    - /{name}/admin" for name in runs)
            return f"""
    This server was launched with -rp on a parent directory. Available runs:

{run_lines}

    Other pages per run:
    - /{{run_name}}/nightplan
    - /{{run_name}}/{{program_code}}
    - /{{run_name}}/{{program_code}}/{{target}}
    """
        return """
    This server was launched with -rp (single run path). Use flat URLs:

    - /admin
    - /nightplan
    - /{program_code}
    - /{program_code}/{target}

    Example:
    - http://localhost:50001/admin
    """
    if _uptree_path:
        return """
    This server was launched with -up (uptree path). Append to the URL:

    url/{semester_code}/{date}/{band}/{page}

    where:
    - semester_code is the four digit year and one letter semester
    - date is in format YYYY-MM-DD
    - band is either band1, band2, or band3 (or full-band1, full-band2, or full-band3)
    - page is one of: {program_code}, {program_code}/{target}, nightplan, or admin

    Examples:
    - /2025B/2025-01-15/band1/admin
    - /2025B/2025-01-15/band3/nightplan
    - /2025B/2025-01-15/band1/2025B_N001
    - /2025B/2025-01-15/band1/2025B_N001/HD4614
    """
    return """
    Launch the webapp with a data path:

    - astroq webapp -up /path/to/uptree
      then open /{semester}/{date}/{band}/admin

    - astroq webapp -rp /path/to/run
      then open /admin (single run) or /{run_name}/admin (parent of runs)

    The run path (-rp) can be a run folder with outputs/, outputs/ itself,
    or a parent directory containing multiple run folders.
    """


# New homepage with navigation instructions
@app.route("/", methods=["GET"])
def index():
    return render_template(
        "homepage.html", navigation_text=_homepage_navigation_text()
    )


# ------------------------------------------------------------------
# Parent -rp routes: /{run_name}/...
# ------------------------------------------------------------------


@app.route("/<run_name>/admin", methods=["GET"])
def run_parent_admin(run_name):
    _require_run_path()
    if not _parent_run_mode():
        abort(404)
    success, message = load_data_for_run_path(_run_path, run_name=run_name)
    if not success:
        return f"Error: {message}", 404
    semester_code, date, band = webrender.route_context_from_planner(
        loaded_run.semester_planner
    )
    return webrender.build_admin_html(
        loaded_run, semester_code, date, band, link_targets=True
    )


@app.route("/<run_name>/nightplan", methods=["GET"])
def run_parent_nightplan(run_name):
    _require_run_path()
    if not _parent_run_mode():
        abort(404)
    success, message = load_data_for_run_path(_run_path, run_name=run_name)
    if not success:
        return f"Error: {message}", 404
    _, _, band = webrender.route_context_from_planner(loaded_run.semester_planner)
    try:
        return webrender.build_nightplan_html(loaded_run, band)
    except ValueError as e:
        return f"Error: {e}", 404


@app.route("/<run_name>/<program_code>/<target>", methods=["GET"])
def run_parent_star_page(run_name, program_code, target):
    _require_run_path()
    if not _parent_run_mode():
        abort(404)
    if program_code in ("admin", "nightplan"):
        abort(404)
    success, message = load_data_for_run_path(_run_path, run_name=run_name)
    if not success:
        return f"Error: {message}", 404
    try:
        return webrender.build_star_html(loaded_run, target, program_code)
    except KeyError as e:
        return f"Error: {e}", 404


@app.route("/<run_name>/<page>", methods=["GET"])
def run_parent_dynamic_page(run_name, page):
    _require_run_path()
    if not _parent_run_mode():
        abort(404)
    if page in ("admin", "nightplan"):
        abort(404)
    success, message = load_data_for_run_path(_run_path, run_name=run_name)
    if not success:
        return f"Error: {message}", 404

    semester_code, date, band = webrender.route_context_from_planner(
        loaded_run.semester_planner
    )
    if page in loaded_run.plot_data.program_dict:
        try:
            return webrender.build_program_html(
                loaded_run, semester_code, date, band, page
            )
        except KeyError as e:
            return f"Error: {e}", 404
    abort(404, description=f"Page '{page}' not found for run '{run_name}'")


# ------------------------------------------------------------------
# Flat routes for -rp (single run path)
# ------------------------------------------------------------------


@app.route("/admin", methods=["GET"])
def run_path_admin():
    _require_run_path()
    if _parent_run_mode():
        runs = webrender.list_child_runs(_run_path)
        example = runs[0] if runs else "run_name"
        abort(
            404,
            description=(
                f"Use /{{run_name}}/admin for a child run, e.g. /{example}/admin."
            ),
        )
    success, message = load_data_for_run_path(_run_path)
    if not success:
        return f"Error: {message}", 404
    semester_code, date, band = webrender.route_context_from_planner(
        loaded_run.semester_planner
    )
    return webrender.build_admin_html(
        loaded_run, semester_code, date, band, link_targets=True
    )


@app.route("/nightplan", methods=["GET"])
def run_path_nightplan():
    _require_run_path()
    if _parent_run_mode():
        abort(404, description="Use /{run_name}/nightplan for a child run.")
    success, message = load_data_for_run_path(_run_path)
    if not success:
        return f"Error: {message}", 404
    _, _, band = webrender.route_context_from_planner(loaded_run.semester_planner)
    try:
        return webrender.build_nightplan_html(loaded_run, band)
    except ValueError as e:
        return f"Error: {e}", 404


@app.route("/<program_code>/<target>", methods=["GET"])
def run_path_star_page(program_code, target):
    if not _single_run_mode():
        abort(404)
    if program_code in ("admin", "nightplan"):
        abort(404)
    success, message = load_data_for_run_path(_run_path)
    if not success:
        return f"Error: {message}", 404
    try:
        return webrender.build_star_html(loaded_run, target, program_code)
    except KeyError as e:
        return f"Error: {e}", 404


@app.route("/<page>", methods=["GET"])
def run_path_dynamic_page(page):
    if not _single_run_mode():
        abort(404)
    success, message = load_data_for_run_path(_run_path)
    if not success:
        return f"Error: {message}", 404

    semester_code, date, band = webrender.route_context_from_planner(
        loaded_run.semester_planner
    )
    if page in loaded_run.plot_data.program_dict:
        try:
            return webrender.build_program_html(
                loaded_run, semester_code, date, band, page
            )
        except KeyError as e:
            return f"Error: {e}", 404
    abort(404, description=f"Page '{page}' not found")


# ------------------------------------------------------------------
# Uptree routes for -up
# ------------------------------------------------------------------


# Star page: /semester/date/band/program_code/target (star under program)
@app.route("/<semester_code>/<date>/<band>/<program_code>/<target>")
def star_page(semester_code, date, band, program_code, target):
    """Handle star page route: star is under program in URL."""
    global _uptree_path
    if not _uptree_path:
        abort(
            404,
            description="Uptree routes require launching with -up/--uptree_path.",
        )
    _validate_band(band)
    success, message = load_data_for_path(semester_code, date, band, _uptree_path)
    if not success:
        return f"Error: {message}", 404
    try:
        return webrender.build_star_html(loaded_run, target, program_code)
    except KeyError as e:
        return f"Error: {e}", 404


# Dynamic route for program, admin, nightplan
@app.route("/<semester_code>/<date>/<band>/<page>")
def dynamic_page(semester_code, date, band, page):
    """Handle program, admin, and nightplan routes."""
    global _uptree_path
    if not _uptree_path:
        abort(
            404,
            description="Uptree routes require launching with -up/--uptree_path.",
        )
    _validate_band(band)
    success, message = load_data_for_path(semester_code, date, band, _uptree_path)
    if not success:
        return f"Error: {message}", 404
    if page == "admin":
        return webrender.build_admin_html(
            loaded_run, semester_code, date, band, link_targets=True
        )
    if page == "nightplan":
        try:
            return webrender.build_nightplan_html(loaded_run, band)
        except ValueError as e:
            return f"Error: {e}", 404
    if page in loaded_run.plot_data.program_dict:
        try:
            return webrender.build_program_html(
                loaded_run, semester_code, date, band, page
            )
        except KeyError as e:
            return f"Error: {e}", 404
    abort(404, description=f"Page '{page}' not found")


def launch_app(uptree_path=None, run_path=None, port=50001):
    """Launch the Flask app."""
    global _uptree_path, _run_path
    _uptree_path = uptree_path
    _run_path = run_path

    if running_on_keck_machines:
        app.run(host=gethostname(), debug=False, use_reloader=False, port=port)
    else:
        app.run(debug=True, use_reloader=True, port=port)


if __name__ == "__main__":
    launch_app(uptree_path=".")
