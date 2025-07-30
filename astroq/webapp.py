import os
from flask import Flask, render_template, request
import numpy as np
import plotly.io as pio
import astroq.plot as pl
import astroq.dynamic as dn
import astroq.splan as splan
import astroq.nplan as nplan
import threading

app = Flask(__name__, template_folder="../templates")

# Shared homepage view
@app.route("/", methods=["GET", "POST"])
def index():
    pages = ["Overview", "Single Star"]
    selected_page = request.form.get("page", pages[0])  # Default to "Overview"

    figures_html = []

    if selected_page == "Overview":
        labels = ["All Programs COF", "Bird's Eye"]
        for label in labels:
            fig_html = dn.get_cof(semester_planner, all_stars_list) if label == "All Programs COF" else dn.get_birdseye(semester_planner, data_astroq[2], all_stars_list)
            figures_html.append(fig_html)

    return render_template("index.html", figures_html=figures_html, pages=pages, selected_page=selected_page)

# Admin page route
@app.route("/admin")
def admin():
    all_stars_from_all_programs = np.concatenate(list(data_astroq[0].values()))
    
    table_reqframe_html = dn.get_requests_frame(semester_planner, filter_condition=None)
    tables_html = [table_reqframe_html]
    fig_cof_html = dn.get_cof(semester_planner, list(data_astroq[1].values()))
    fig_birdseye_html = dn.get_birdseye(semester_planner, data_astroq[2], list(data_astroq[1].values()))
    figures_html = [fig_cof_html, fig_birdseye_html]

    return render_template("admin.html", tables_html=tables_html, figures_html=figures_html)


# /semesterplan landing page
@app.route("/semesterplan")
def semesterplan_home():
    return render_template("semesterplan.html", starname=None, figure_html=None)

# Page route for the semesterplan of a specific program
@app.route("/semesterplan/<programname>")
def single_program(programname):

    if programname not in data_astroq[1]:
        return f"Error: program '{programname}' not found."

    star_obj_list = list(data_astroq[0][programname])

    table_reqframe_html = dn.get_requests_frame(semester_planner, filter_condition=f"program_code=='{programname}'")
    tables_html = [table_reqframe_html]

    fig_cof_html = dn.get_cof(semester_planner, star_obj_list)
    fig_birdseye_html = dn.get_birdseye(semester_planner, data_astroq[2], star_obj_list)
    figures_html = [fig_cof_html, fig_birdseye_html]

    return render_template("semesterplan.html", programname=programname, tables_html=tables_html, figures_html=figures_html)

# /star landing page
@app.route("/star")
def star_home():
    return render_template("star.html", starname=None, figure_html=None)

# Star detail page route
@app.route("/star/<starname>")
def single_star(starname):

    compare_starname = starname.lower().replace(' ', '') # Lower case and remove all spaces
    program_names = data_astroq[0].keys()

    for program in program_names:
        for star_ind in range(len(data_astroq[0][program])):
            star_obj = data_astroq[0][program][star_ind]

            true_starname = star_obj.starname
            object_compare_starname = true_starname.lower().replace(' ', '')
            
            if object_compare_starname == compare_starname:
                table_reqframe_html = dn.get_requests_frame(semester_planner, filter_condition=f"starname=='{true_starname}'")

                fig_cof_html = dn.get_cof(semester_planner, [data_astroq[0][program][star_ind]])
                fig_birdseye_html = dn.get_birdseye(semester_planner, data_astroq[2], [star_obj])

                tables_html = [table_reqframe_html]
                figures_html = [fig_cof_html, fig_birdseye_html]

                return render_template("star.html", starname=true_starname, tables_html=tables_html, figures_html=figures_html)
    return f"Error, star {starname} not found in programs {list(program_names)}"


# Page route for single night plan
@app.route("/nightplan")
def nightplan():
    plots = ['script_table', 'slewgif', 'ladder', 'slewpath']

    if data_ttp is not None:
        script_table_html = dn.get_script_plan(config_file, data_ttp)
        ladder_html = dn.get_ladder(data_ttp)
        slew_animation_html = dn.get_slew_animation(data_ttp, animationStep=120)
        slew_path_html = dn.plot_path_2D_interactive(data_ttp)
        figure_html_list = [script_table_html, ladder_html, slew_animation_html, slew_path_html]

    else:
        figure_html_list = []
    return render_template("nightplan.html", starname=None, figure_html_list=figure_html_list)

def run_server():
    app.run(debug=False, use_reloader=False)

def run_server_coverage():
    # Run the Flask server in a background thread
    server_thread = threading.Thread(target=lambda: app.run(debug=False, use_reloader=False), daemon=True)
    server_thread.start()
    admin()
    single_program('U001')
    single_star('191939')
    nightplan()

def launch_app(config_file, flag=False):
    global data_astroq, data_ttp, semester_planner, night_planner, all_stars_list

    # Create SemesterPlanner and NightPlanner instead of manager
    semester_planner = splan.SemesterPlanner(config_file)
    night_planner = nplan.NightPlanner(config_file)

    data_astroq = pl.read_star_objects(night_planner.reports_directory + "star_objects.pkl")
    all_stars_list = [star_obj for star_obj_list in data_astroq[0].values() for star_obj in star_obj_list]

    ttp_path = os.path.join(night_planner.reports_directory, "ttp_data.pkl")
    if os.path.exists(ttp_path):
        data_ttp = pl.read_star_objects(ttp_path)
    else:
        data_ttp = None
    if flag:
        run_server_coverage()
    else:
        run_server()

@app.route('/shutdown')
def shutdown():
    func = request.environ.get('werkzeug.server.shutdown')
    if func:
        func()
    return 'Shutting down...'

if __name__=="__main__":

    cf =  'examples/hello_world/config_hello_world.ini'
    launch_app(cf)
    admin()

