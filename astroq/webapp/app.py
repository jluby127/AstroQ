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
uptree_path = None


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


# New homepage with navigation instructions
@app.route("/", methods=["GET"])
def index():
    navigation_text = """
    To navigate, append to the URL in the following way:
    url/{semester_code}/{date}/{band}/{page}

    where:
    - semester_code is the four digit year and one letter semester
    - date is in format YYYY-MM-DD
    - band is either band1, band2, or band3 (or full-band1, full-band2, or full-band3)
    - page is one of: {program_code}, {program_code}/{target}, nightplan, or admin

    Examples:
    - /2025B/2025-01-15/band1/admin
    - /2025B/2025-01-15/band3/nightplan
    - /2025B/2025-01-15/band1/2025B_N001                          (program overview)
    - /2025B/2025-01-15/band1/2025B_N001/HD4614                   (star under program)

    Note: program_code contains the semester information. Correct: 2025B_N001, Incorrect: N001
    Note: You only have access to the programs and stars for which you are a PI or named Co-I on the proposal coversheet.
    Note: Access to nightplan pages is for observers.
    Note: Access to admin pages is for the queue manager and observatory staff.
    """
    return render_template("homepage.html", navigation_text=navigation_text)


# Star page: /semester/date/band/program_code/target (star under program)
@app.route("/<semester_code>/<date>/<band>/<program_code>/<target>")
def star_page(semester_code, date, band, program_code, target):
    """Handle star page route: star is under program in URL."""
    global uptree_path
    _validate_band(band)
    success, message = load_data_for_path(semester_code, date, band, uptree_path)
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
    global uptree_path
    _validate_band(band)
    success, message = load_data_for_path(semester_code, date, band, uptree_path)
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
    if page in loaded_run.data_astroq[0]:
        try:
            return webrender.build_program_html(
                loaded_run, semester_code, date, band, page
            )
        except KeyError as e:
            return f"Error: {e}", 404
    abort(404, description=f"Page '{page}' not found")


def launch_app(uptree_path_param, port=50001):
    """Launch the Flask app"""
    global uptree_path
    uptree_path = uptree_path_param

    if running_on_keck_machines:
        app.run(host=gethostname(), debug=False, use_reloader=False, port=port)
    else:
        app.run(debug=True, use_reloader=True, port=port)


if __name__ == "__main__":
    launch_app(".")
