"""
Web application module for AstroQ.
"""

# Standard library imports
import base64
import os
import pickle
import threading
from configparser import ConfigParser
from io import BytesIO
        
# Third-party imports
import imageio.v3 as iio
import numpy as np
import pandas as pd
import plotly.io as pio
from flask import Flask, render_template, request, abort
from socket import gethostname 

# Local imports
import astroq.nplan as nplan
import astroq.plot as pl
import astroq.splan as splan
from astroq.splan import SemesterPlanner
from astroq.nplan import NightPlanner
from astroq.nplan import get_nightly_times_from_allocation

running_on_keck_machines = False

app = Flask(__name__, template_folder="../templates")

# Global variables to store loaded data
data_astroq = None
data_ttp = None
semester_planner = None
night_planner = None
uptree_path = None
semester_planner_timestamp = None

def load_data_for_path(semester_code, date, band, uptree_path):
    """
    Load data for a specific semester_code/date/band combination
    
    Args:
        semester_code (str): the semester code
        date (str): the date in YYYY-MM-DD format
        band (str): the band
        uptree_path (str): the path to the uptree directory

    Returns:
        success (bool): True if data loaded successfully, False otherwise
    """
    global data_astroq, data_ttp, semester_planner, night_planner, request_frame_path, night_start_time
    
    # Construct the workdir path based on URL parameters
    workdir = os.path.join(uptree_path, semester_code, date, band, "outputs")
    request_frame_path = os.path.join(workdir, 'request_selected.csv')
    
    # Check if the directory exists
    if not os.path.exists(workdir):
        return False, f"Directory not found: {workdir}"
    
    semester_planner_h5 = os.path.join(workdir, 'semester_planner.h5')
    night_planner_h5 = os.path.join(workdir, 'night_planner.h5')

    # Load semester planner
    global semester_planner_timestamp
    try:
        semester_planner = SemesterPlanner.from_hdf5(semester_planner_h5)
        data_astroq = pl.process_stars(semester_planner)
        # Get file modification time
        if os.path.exists(semester_planner_h5):
            from datetime import datetime
            mtime = os.path.getmtime(semester_planner_h5)
            semester_planner_timestamp = datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')
        else:
            semester_planner_timestamp = None
    except Exception as e:
        semester_planner = None
        data_astroq = None
        semester_planner_timestamp = None
        return False, f"Error loading semester planner: {str(e)}"
    
    # Load night planner (optional)
    try:
        print(night_planner_h5)
        night_planner = NightPlanner.from_hdf5(night_planner_h5)
        data_ttp = night_planner.solution

        # Get the night start time from allocation file (this is "Minute 0")
        night_start_time, _ = nplan.get_nightly_times_from_allocation(
            night_planner.allocation_file, 
            night_planner.current_day
        )
    except Exception as e:
        print(f"No night planner found")
        # import traceback
        # traceback.print_exc()
        night_planner = None
        data_ttp = None
    
    return True, "Data loaded successfully"

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
    - page is one of: {program_code}, {program_code}/{starname}, nightplan, or admin

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

# Star page: /semester/date/band/program_code/starname (star under program)
@app.route("/<semester_code>/<date>/<band>/<program_code>/<starname>")
def star_page(semester_code, date, band, program_code, starname):
    """Handle star page route: star is under program in URL."""
    global uptree_path
    if band not in ['band1', 'band2', 'band3', 'full-band1', 'full-band2', 'full-band3']:
        abort(400, description="Band must be 'band1', 'band2', 'band3', 'full-band1', 'full-band2', or 'full-band3'")
    success, message = load_data_for_path(semester_code, date, band, uptree_path)
    if not success:
        return f"Error: {message}", 404
    return render_star_page(starname, program_code)

# Dynamic route for program, admin, nightplan
@app.route("/<semester_code>/<date>/<band>/<page>")
def dynamic_page(semester_code, date, band, page):
    """Handle program, admin, and nightplan routes."""
    global uptree_path
    if band not in ['band1', 'band2', 'band3', 'full-band1', 'full-band2', 'full-band3']:
        abort(400, description="Band must be 'band1', 'band2', 'band3', 'full-band1', 'full-band2', or 'full-band3'")
    success, message = load_data_for_path(semester_code, date, band, uptree_path)
    if not success:
        return f"Error: {message}", 404
    if page == "admin":
        return render_admin_page(semester_code, date, band)
    elif page == "nightplan":
        return render_nightplan_page(band)
    elif page in data_astroq[0]:
        return render_program_page(semester_code, date, band, page)
    else:
        abort(404, description=f"Page '{page}' not found")

def render_admin_page(semester_code, date, band):
    """Render the admin page"""
    if data_astroq is None:
        return "Error: No data available", 404

    all_stars_from_all_programs = np.concatenate(list(data_astroq[0].values()))

    # Get request frame table for all stars, with starname as links under program
    request_df = pl.get_request_frame(semester_planner, all_stars_from_all_programs)
    request_table_html = pl.request_frame_to_html(request_df, semester_code, date, band)
    
    fig_cof1 = pl.get_cof(semester_planner, list(data_astroq[1].values()))
    fig_cof2 = pl.get_cof(semester_planner, list(data_astroq[1].values()), use_time=True)

    fig_birdseye = pl.get_birdseye(semester_planner, data_astroq[2], list(data_astroq[1].values()))
    fig_football = pl.get_football(semester_planner, all_stars_from_all_programs, use_program_colors=True)
    fig_tau_inter_line = pl.get_tau_inter_line(semester_planner, all_stars_from_all_programs, use_program_colors=True)
    fig_timebar = pl.get_timebar(semester_planner, all_stars_from_all_programs, use_program_colors=True)
    fig_timebar_by_program = pl.get_timebar_by_program(semester_planner, data_astroq[0])
    fig_rawobs = pl.get_rawobs(semester_planner, all_stars_from_all_programs, use_program_colors=True)

    fig_cof_html1 = pio.to_html(fig_cof1, full_html=True, include_plotlyjs='cdn')
    fig_cof_html2 = pio.to_html(fig_cof2, full_html=True, include_plotlyjs='cdn')
    fig_birdseye_html = pio.to_html(fig_birdseye, full_html=True, include_plotlyjs='cdn')
    fig_football_html = pio.to_html(fig_football, full_html=True, include_plotlyjs='cdn')
    fig_tau_inter_line_html = pio.to_html(fig_tau_inter_line, full_html=True, include_plotlyjs='cdn')
    fig_timebar_html = pio.to_html(fig_timebar, full_html=True, include_plotlyjs='cdn')
    fig_timebar_by_program_html = pio.to_html(fig_timebar_by_program, full_html=True, include_plotlyjs='cdn')
    fig_rawobs_html = pio.to_html(fig_rawobs, full_html=True, include_plotlyjs='cdn')

    figures_html = [fig_timebar_html, fig_timebar_by_program_html, fig_cof_html1, fig_cof_html2, fig_birdseye_html, fig_rawobs_html, fig_tau_inter_line_html, fig_football_html]

    return render_template("admin.html", tables_html=[request_table_html], figures_html=figures_html, timestamp=semester_planner_timestamp)

def render_program_page(semester_code, date, band, program_code):
    """Render the program overview page for a specific program"""
    if data_astroq is None:
        return "Error: No data available", 404
    
    # Get all stars in the specified program
    if program_code not in data_astroq[0]:
        return f"Error: Program {program_code} not found", 404
    
    program_stars = data_astroq[0][program_code]
    
    # Get request frame table for this program's stars, with starname as links
    request_df = pl.get_request_frame(semester_planner, program_stars)
    request_table_html = pl.request_frame_to_html(request_df, semester_code, date, band)
    
    # Create overview figures for this program
    fig_cof = pl.get_cof(semester_planner, program_stars)
    fig_birdseye = pl.get_birdseye(semester_planner, data_astroq[2], program_stars)
    fig_tau_inter_line = pl.get_tau_inter_line(semester_planner, program_stars)
    fig_football = pl.get_football(semester_planner, program_stars)
    fig_timebar = pl.get_timebar(semester_planner, program_stars, use_program_colors=True)
    fig_rawobs = pl.get_rawobs(semester_planner, program_stars)

    fig_cof_html = pio.to_html(fig_cof, full_html=True, include_plotlyjs='cdn')
    fig_birdseye_html = pio.to_html(fig_birdseye, full_html=True, include_plotlyjs='cdn')
    fig_tau_inter_line_html = pio.to_html(fig_tau_inter_line, full_html=True, include_plotlyjs='cdn')
    fig_football_html = pio.to_html(fig_football, full_html=True, include_plotlyjs='cdn')
    fig_timebar_html = pio.to_html(fig_timebar, full_html=True, include_plotlyjs='cdn')
    fig_rawobs_html = pio.to_html(fig_rawobs, full_html=True, include_plotlyjs='cdn')

    figures_html = [fig_timebar_html, fig_cof_html, fig_birdseye_html, fig_rawobs_html, fig_tau_inter_line_html, fig_football_html]
    
    return render_template("semesterplan.html", 
                         programname=program_code, 
                         tables_html=[request_table_html], 
                         figures_html=figures_html, 
                         programs=[program_code],
                         timestamp=semester_planner_timestamp)

def render_star_page(starname, program_code=None):
    """Render a specific star page. If program_code is given, only look in that program."""
    if data_astroq is None:
        return "Error: No data available", 404

    compare_starname = starname.lower().replace(' ', '')  # Lower case and remove all spaces
    programs_to_search = [program_code] if program_code and program_code in data_astroq[0] else data_astroq[0].keys()

    for program in programs_to_search:
        for star_ind in range(len(data_astroq[0][program])):
            star_obj = data_astroq[0][program][star_ind]

            true_starname = star_obj.starname
            object_compare_starname = true_starname.lower().replace(' ', '')

            if object_compare_starname == compare_starname:
                # Get request frame table for this specific star (no star links needed)
                request_df = pl.get_request_frame(semester_planner, [star_obj])
                request_table_html = pl.request_frame_to_html(request_df)

                fig_cof = pl.get_cof(semester_planner, [data_astroq[0][program][star_ind]])
                fig_birdseye = pl.get_birdseye(semester_planner, data_astroq[2], [star_obj])
                fig_tau_inter_line = pl.get_tau_inter_line(semester_planner, [star_obj])
                fig_football = pl.get_football(semester_planner, [star_obj])
                fig_rawobs = pl.get_rawobs(semester_planner, [star_obj])

                fig_cof_html = pio.to_html(fig_cof, full_html=True, include_plotlyjs='cdn')
                fig_birdseye_html = pio.to_html(fig_birdseye, full_html=True, include_plotlyjs='cdn')
                fig_tau_inter_line_html = pio.to_html(fig_tau_inter_line, full_html=True, include_plotlyjs='cdn')
                fig_football_html = pio.to_html(fig_football, full_html=True, include_plotlyjs='cdn')
                fig_rawobs_html = pio.to_html(fig_rawobs, full_html=True, include_plotlyjs='cdn')

                tables_html = [request_table_html]
                figures_html = [fig_cof_html, fig_birdseye_html, fig_rawobs_html, fig_tau_inter_line_html, fig_football_html]

                return render_template("star.html", starname=true_starname, tables_html=tables_html, figures_html=figures_html, timestamp=semester_planner_timestamp)
    
    return f"Error, star {starname} not found in programs {list(program_names)}"

def render_nightplan_page(band):
    """Render the night plan page"""
    if data_ttp is None:
        return "Error: No night planner data available", 404
    
    plots = ['script_table', 'slewgif', 'ladder', 'slewpath']
    
    script_table_df = pl.get_script_plan(night_planner)
    ladder_fig = pl.get_ladder(data_ttp, night_start_time)
    slew_animation_fig = pl.get_slew_animation_plotly(data_ttp, request_frame_path, animationStep=120)
    slew_path_fig = pl.plot_path_2D_interactive(data_ttp, night_start_time=night_start_time)
    
    script_table_html = pl.nightplan_table_to_html(script_table_df, table_id='script-table', page_size=100)
    # Convert figures to HTML
    ladder_html = pio.to_html(ladder_fig, full_html=True, include_plotlyjs='cdn')
    slew_animation_html = pio.to_html(slew_animation_fig, full_html=True, include_plotlyjs='cdn')
    slew_path_html = pio.to_html(slew_path_fig, full_html=True, include_plotlyjs='cdn')
    
    figure_html_list = [script_table_html, ladder_html, slew_animation_html, slew_path_html]

    return render_template("nightplan.html", starname=None, figure_html_list=figure_html_list, 
                         semester_planner=semester_planner, night_planner=night_planner, band=band)

@app.route("/<semester_code>/<date>/<band>/download_nightplan")
def download_nightplan(semester_code, date, band):
    """Download the Magiq formatted night plan file"""
    global uptree_path, semester_planner, night_planner
    
    # Validate parameters
    if band not in ['band1', 'band3']:
        abort(400, description="Band must be 'band1' or 'band3'")
    
    # Load data for this path
    success, message = load_data_for_path(semester_code, date, band, uptree_path)
    if not success:
        return f"Error: {message}", 404
    
    if semester_planner is None or night_planner is None:
        return "Error: No planner data available", 404
    
    try:
        # Construct the path to the script file
        script_file_path = os.path.join(semester_planner.output_directory, 
                                      f'script_{night_planner.current_day}_nominal.txt')
        
        if not os.path.exists(script_file_path):
            return "Error: Night plan file not found", 404
        
        # Return the file for download
        from flask import send_file
        return send_file(script_file_path, 
                        as_attachment=True,
                        download_name=f'script_{night_planner.current_day}_nominal.txt',
                        mimetype='text/plain')
        
    except Exception as e:
        return f"Error downloading file: {str(e)}", 500

def launch_app(uptree_path_param):
    """Launch the Flask app"""
    global uptree_path
    uptree_path = uptree_path_param
    
    if running_on_keck_machines:
        app.run(host=gethostname(), debug=False, use_reloader=False, port=50001)
    else:
        app.run(debug=True, use_reloader=True, port=50001)

if __name__ == "__main__":
    launch_app(".")
