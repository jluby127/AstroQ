from flask import Flask, render_template, request
import plotly.io as pio
import plotly.graph_objs as go
import kpfcc.plot as pl
import kpfcc.management as mn
import os

app = Flask(__name__, template_folder="../templates")

def cof_birdseye(program):
    # file_root = "/Users/judahvz/research/code/GitHub/AstroQ/examples/hello_world/"
    # config_file = os.path.join(file_root, "config_hello_world.ini")
    # pkl_path = os.path.join(file_root, "reports/admin/2024-08-01/")
    # data = pl.read_star_objects(pkl_path)
    #
    # manager = mn.data_admin(config_file)
    # manager.run_admin()
    #
    # all_stars_list = [obj for obj_list in data[0].values() for obj in obj_list]

    if program == 'All Stars COF':
        return pl.cof_builder(all_stars_list, manager, flag=True)
    elif program == 'All Programs COF':
        return pl.cof_builder(list(data[1].values()), manager, flag=True)
    elif program == 'Birds Eye':
        return pl.generate_birds_eye(manager, data[2], all_stars_list)

@app.route("/", methods=['GET', 'POST'])
def index():
    labels = ["All Stars COF", "All Programs COF", "Birds Eye"]
    selected_label = request.form.get("plot_type", labels[0])  # Default to first
    fig = cof_birdseye(selected_label)
    fig_html = pio.to_html(fig, include_plotlyjs=False, full_html=False, div_id="main_plot")
    
    return render_template("index2.html", figure_html=fig_html, labels=labels, selected_label=selected_label)
    
    
def launch_app(pkl_path, config_file):
    """
    Wrapper for index() functions that
    accepts arguments and defines global
    variables for use in index()
    """
    
    global data, manager, all_stars_list
    
    data = pl.read_star_objects(pkl_path)
    
    manager = mn.data_admin(config_file)
    manager.run_admin()
    
    all_stars_list = [obj for obj_list in data[0].values() for obj in obj_list]
    
    app.run(debug=True)
    
    return

if __name__ == "__main__":
    
    from kpfcc import _ROOT

    pkl_path = os.path.join(_ROOT, "../examples/hello_world/reports/admin/2024-08-01/")
    config_file = os.path.join(_ROOT, "../examples/hello_world/config_hello_world.ini")
    
    launch_app(pkl_path, config_file)

