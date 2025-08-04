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
import plotly.io as pio
from flask import Flask, render_template, request, abort
from socket import gethostname 

# Local imports
import astroq.dynamic as dn
import astroq.nplan as nplan
import astroq.plot as pl
import astroq.splan as splan

running_on_keck_machines = False

app = Flask(__name__, template_folder="../templates")

# Global variables to store loaded data
data_astroq = None
data_ttp = None
semester_planner = None
night_planner = None
uptree_path = None

def load_data_for_path(semester_code, date, band, uptree_path):
    """Load data for a specific semester_code/date/band combination"""
    global data_astroq, data_ttp, semester_planner, night_planner
    
    # Construct the workdir path based on URL parameters
    workdir = f"{uptree_path}/{semester_code}/{date}/{band}/outputs/"
    
    # Check if the directory exists
    if not os.path.exists(workdir):
        return False, f"Directory not found: {workdir}"
    
    semester_planner_pkl = os.path.join(workdir, 'semester_planner.pkl')
    night_planner_pkl = os.path.join(workdir, 'night_planner.pkl')
    
    # Load semester planner
    try:
        with open(semester_planner_pkl, 'rb') as f:
            semester_planner = pickle.load(f)
        data_astroq = pl.process_stars(semester_planner)
    except Exception as e:
        semester_planner = None
        data_astroq = None
        return False, f"Error loading semester planner: {str(e)}"
    
    # Load night planner (optional)
    try:
        with open(night_planner_pkl, 'rb') as f:
            night_planner = pickle.load(f)
            data_ttp = night_planner.solution
    except:
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
    - semester_code is 2025B (for example)
    - date is in format YYYY-MM-DD
    - band is either band1 or band3
    - page is one of the following: admin, semester_id, star/{starname}, or nightplan
    
    Examples:
    - /2025B/2025-01-15/band1/admin
    - /2025B/2025-01-15/band3/semester_id
    - /2025B/2025-01-15/band1/star/HD4614
    - /2025B/2025-01-15/band3/nightplan
    """
    return render_template("homepage.html", navigation_text=navigation_text)

# Dynamic route for all pages
@app.route("/<semester_code>/<date>/<band>/<page>")
@app.route("/<semester_code>/<date>/<band>/star/<starname>")
def dynamic_page(semester_code, date, band, page, starname=None):
    """Handle all dynamic routes based on URL parameters"""
    global uptree_path
    
    # Validate parameters
    if band not in ['band1', 'band3']:
        abort(400, description="Band must be 'band1' or 'band3'")
    
    # Load data for this path
    success, message = load_data_for_path(semester_code, date, band, uptree_path)
    if not success:
        return f"Error: {message}", 404
    
    # Route to appropriate page based on 'page' parameter
    if page == "admin":
        return render_admin_page()
    elif page == "semester_id":
        return render_semester_overview_page()
    elif page == "star" and starname:
        return render_star_page(starname)
    elif page == "nightplan":
        return render_nightplan_page()
    else:
        abort(404, description=f"Page '{page}' not found")

def render_admin_page():
    """Render the admin page"""
    if data_astroq is None:
        return "Error: No data available", 404
    
    all_stars_from_all_programs = np.concatenate(list(data_astroq[0].values()))
    
    tables_html = []
    fig_cof = dn.get_cof(semester_planner, list(data_astroq[1].values()))
    fig_birdseye = dn.get_birdseye(semester_planner, data_astroq[2], list(data_astroq[1].values()))
    fig_cof_html = pio.to_html(fig_cof, full_html=True, include_plotlyjs='cdn')
    fig_birdseye_html = pio.to_html(fig_birdseye, full_html=True, include_plotlyjs='cdn')
    figures_html = [fig_cof_html, fig_birdseye_html]

    return render_template("admin.html", tables_html=tables_html, figures_html=figures_html)

def render_semester_overview_page():
    """Render the semester overview page"""
    if data_astroq is None:
        return "Error: No data available", 404
    
    # Get all programs
    programs = list(data_astroq[0].keys())
    
    # Create overview figures
    all_stars_list = list(data_astroq[1].values())
    fig_cof = dn.get_cof(semester_planner, all_stars_list)
    fig_birdseye = dn.get_birdseye(semester_planner, data_astroq[2], all_stars_list)
    fig_cof_html = pio.to_html(fig_cof, full_html=True, include_plotlyjs='cdn')
    fig_birdseye_html = pio.to_html(fig_birdseye, full_html=True, include_plotlyjs='cdn')
    figures_html = [fig_cof_html, fig_birdseye_html]
    
    return render_template("semesterplan.html", programname=None, tables_html=[], figures_html=figures_html, programs=programs)

def render_star_page(starname):
    """Render a specific star page"""
    if data_astroq is None:
        return "Error: No data available", 404
    
    compare_starname = starname.lower().replace(' ', '') # Lower case and remove all spaces
    program_names = data_astroq[0].keys()

    for program in program_names:
        for star_ind in range(len(data_astroq[0][program])):
            star_obj = data_astroq[0][program][star_ind]

            true_starname = star_obj.starname
            object_compare_starname = true_starname.lower().replace(' ', '')
            
            if object_compare_starname == compare_starname:
                table_reqframe_html = dn.get_requests_frame(semester_planner, filter_condition=f"starname=='{true_starname}'")

                fig_cof = dn.get_cof(semester_planner, [data_astroq[0][program][star_ind]])
                fig_birdseye = dn.get_birdseye(semester_planner, data_astroq[2], [star_obj])
                fig_cof_html = pio.to_html(fig_cof, full_html=True, include_plotlyjs='cdn')
                fig_birdseye_html = pio.to_html(fig_birdseye, full_html=True, include_plotlyjs='cdn')

                tables_html = [table_reqframe_html]
                figures_html = [fig_cof_html, fig_birdseye_html]

                return render_template("star.html", starname=true_starname, tables_html=tables_html, figures_html=figures_html)
    
    return f"Error, star {starname} not found in programs {list(program_names)}"

def render_nightplan_page():
    """Render the night plan page"""
    if data_ttp is None:
        return "Error: No night planner data available", 404
    
    plots = ['script_table', 'slewgif', 'ladder', 'slewpath']

    script_table_df = dn.get_script_plan(semester_planner, data_ttp)
    ladder_fig = dn.get_ladder(data_ttp)
    slew_animation_figures = dn.get_slew_animation(data_ttp, animationStep=120)
    slew_path_fig = dn.plot_path_2D_interactive(data_ttp)
    
    # Convert dataframe to HTML
    script_table_html = dn.dataframe_to_html(script_table_df, sort_column=11)
    # Convert figures to HTML
    ladder_html = pio.to_html(ladder_fig, full_html=True, include_plotlyjs='cdn')
    slew_path_html = pio.to_html(slew_path_fig, full_html=True, include_plotlyjs='cdn')
    
    # Convert matplotlib figures to GIF and then to HTML
    gif_frames = []
    for fig in slew_animation_figures:
        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=100)
        buf.seek(0)
        gif_frames.append(iio.imread(buf))
        buf.close()
    
    gif_buf = BytesIO()
    iio.imwrite(gif_buf, gif_frames, format='gif', loop=0, duration=0.3)
    gif_buf.seek(0)
    
    gif_base64 = base64.b64encode(gif_buf.getvalue()).decode('utf-8')
    slew_animation_html = f'<img src="data:image/gif;base64,{gif_base64}" alt="Observing Animation"/>'
    gif_buf.close()
    
    figure_html_list = [script_table_html, ladder_html, slew_animation_html, slew_path_html]

    return render_template("nightplan.html", starname=None, figure_html_list=figure_html_list)

def launch_app(uptree_path):
    """Launch the Flask app"""
    global uptree_path
    uptree_path = uptree_path
    
    if running_on_keck_machines:
        app.run(host=gethostname(), debug=False, use_reloader=False, port=50001)
    else:
        app.run(debug=True, use_reloader=True, port=50001)

