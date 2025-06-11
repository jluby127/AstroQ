# webapp.py

from flask import Flask, render_template, request
import plotly.io as pio
import kpfcc.plot as pl
import kpfcc.management as mn

app = Flask(__name__, template_folder="../templates")

# Plot-building helper functions
def make_cof(star_obj_list):
    return pl.cof_builder(star_obj_list, manager, flag=True)

def make_birdseye(star_obj_list):
    return pl.generate_birds_eye(manager, data[2], star_obj_list)

# Shared homepage view
@app.route("/", methods=["GET", "POST"])
def index():
    pages = ["Overview", "Single Star"]
    selected_page = request.form.get("page", pages[0])  # Default to "Overview"

    figures_html = []

    if selected_page == "Overview":
        labels = ["All Programs COF", "Bird's Eye"]
        for label in labels:
            fig = make_cof(all_stars_list) if label == "All Programs COF" else make_birdseye(all_stars_list)
            fig_html = pio.to_html(fig, include_plotlyjs=False, full_html=False)
            figures_html.append(fig_html)

    return render_template("index.html", figures_html=figures_html, pages=pages, selected_page=selected_page)

# Admin page route
@app.route("/admin")
def admin():
    fig_cof = make_cof(list(data[1].values()))
    fig_birdseye = make_birdseye(list(data[1].values()))
    figures_html = [
        pio.to_html(fig_cof, include_plotlyjs=False, full_html=False),
        pio.to_html(fig_birdseye, include_plotlyjs=False, full_html=False)
    ]
    return render_template("admin.html", figures_html=figures_html)
    
    

# /semesterplan landing page
@app.route("/semesterplan")
def semesterplan_home():
    return render_template("semesterplan.html", starname=None, figure_html=None)
    
        
# Page route for the semesterplan of a specific program
@app.route("/semesterplan/<programname>")
def single_program(programname):
    
    if programname not in data[1]:
        return f"Error: program '{programname}' not found."
    
    star_obj_list = list(data[0][programname])
    
    fig_cof = make_cof(star_obj_list)
    fig_birdseye = make_birdseye(star_obj_list)
    
    # fig_cof.update_layout(autosize=True)
    # fig_birdseye.update_layout(autosize=True)
    # fig_cof.update_layout(xaxis=dict(scaleanchor="y", scaleratio=2))
    # fig_birdseye.update_layout(xaxis=dict(scaleanchor="y", scaleratio=2))
    
    figure_html = [pio.to_html(fig_cof, include_plotlyjs=False, full_html=False),
                   pio.to_html(fig_birdseye, include_plotlyjs=False, full_html=False)
               ]
                   
    return render_template("semesterplan.html", figure_html=figure_html)
    
    # for program in program_names:
    #     for star_ind in range(len(data[0][program])):
    #         star_obj = data[0][program][star_ind]
    #         if star_obj.starname == starname:
    #             # herb = data[0][program][star_ind]
    #             # herbie = list(herb)
    #             fig = make_cof([data[0][program][star_ind]])
    #             figure_html = [pio.to_html(fig, include_plotlyjs=False, full_html=False)]
    #             return render_template("star.html", figure_html=figure_html)
    # return f"Error, star {starname} not found in programs {list(program_names)}"
    
    
# /star landing page
@app.route("/star")
def star_home():
    return render_template("star.html", starname=None, figure_html=None)
    
        
# Star detail page route
@app.route("/star/<starname>")
def single_star(starname):
    
    program_names = data[0].keys()
    
    for program in program_names:
        for star_ind in range(len(data[0][program])):
            star_obj = data[0][program][star_ind]
            if star_obj.starname == starname:
                # herb = data[0][program][star_ind]
                # herbie = list(herb)
                fig_cof = make_cof([data[0][program][star_ind]])
                fig_birdseye = make_birdseye([data[0][program][star_ind]])
                
                figure_html = [pio.to_html(fig_cof, include_plotlyjs=False, full_html=False),
                               pio.to_html(fig_birdseye, include_plotlyjs=False, full_html=False)
                           ]
                
                return render_template("star.html", figure_html=figure_html)
    return f"Error, star {starname} not found in programs {list(program_names)}"

# # Observer page route
# @app.route("/observer")
# def observer():
#     fig = make_birdseye()
#     figure_html = [pio.to_html(fig, include_plotlyjs=False, full_html=False)]
#     return render_template("observer.html", figure_html=figure_html)
#
# # Program page route
# @app.route("/program")
# def program():
#     fig = make_cof()
#     figure_html = [pio.to_html(fig, include_plotlyjs=False, full_html=False)]
#     return render_template("program.html", figure_html=figure_html)

# Entry point for launching the app
def launch_app(pkl_path, config_file):
    global data, manager, all_stars_list

    data = pl.read_star_objects(pkl_path)
    manager = mn.data_admin(config_file)
    manager.run_admin()
    all_stars_list = [star_obj for star_obj_list in data[0].values() for star_obj in star_obj_list]

    app.run(debug=True)
    
if __name__=="__main__":
    
    # pk = 'examples/hello_world/reports/admin/2024-08-01/'
    # cf =  'examples/hello_world/config_hello_world.ini'
    
    pk = 'examples/bench/reports/admin/2024-08-01/'
    cf =  'examples/bench/config_benchmark.ini'
    
    launch_app(pk, cf)
    single_program('U001')
    # import pdb; pdb.set_trace()
    
