# webapp.py

from flask import Flask, render_template, request
import plotly.io as pio
import plotly.graph_objs as go
import kpfcc.plot as pl
import kpfcc.management as mn
import os

app = Flask(__name__, template_folder="../templates")

def cof_birdseye(program):
    # if program == 'All Stars COF':
    #     return pl.cof_builder(all_stars_list, manager, flag=True)
    if program == 'All Programs COF':
        return pl.cof_builder(list(data[1].values()), manager, flag=True)
    elif program == 'Birds Eye':
        return pl.generate_birds_eye(manager, data[2], all_stars_list)

@app.route("/", methods=['GET', 'POST'])
def index():
    pages = ["Overview", "Single Star"]
    selected_page = request.form.get("page", pages[0])  # Default to "Overview"

    figures_html = []

    if selected_page == "Overview":
        labels = ["All Programs COF", "Birds Eye"]
        for label in labels:
            fig = cof_birdseye(label)
            fig_html = pio.to_html(fig, include_plotlyjs=False, full_html=False)
            figures_html.append(fig_html)

    return render_template("index1.html", figures_html=figures_html, pages=pages, selected_page=selected_page)

def launch_app(pkl_path, config_file):
    global data, manager, all_stars_list

    data = pl.read_star_objects(pkl_path)
    manager = mn.data_admin(config_file)
    manager.run_admin()
    all_stars_list = [obj for obj_list in data[0].values() for obj in obj_list]

    app.run(debug=True)
