# webapp.py

from flask import Flask, render_template, request
import plotly.io as pio
import kpfcc.plot as pl
import kpfcc.management as mn

app = Flask(__name__, template_folder="../templates")

# Plot-building helper functions
def make_cof():
    return pl.cof_builder(list(data[1].values()), manager, flag=True)

def make_birdseye():
    return pl.generate_birds_eye(manager, data[2], all_stars_list)

# Shared homepage view
@app.route("/", methods=["GET", "POST"])
def index():
    pages = ["Overview", "Single Star"]
    selected_page = request.form.get("page", pages[0])  # Default to "Overview"

    figures_html = []

    if selected_page == "Overview":
        labels = ["All Programs COF", "Bird's Eye"]
        for label in labels:
            fig = make_cof() if label == "All Programs COF" else make_birdseye()
            fig_html = pio.to_html(fig, include_plotlyjs=False, full_html=False)
            figures_html.append(fig_html)

    return render_template("index.html", figures_html=figures_html, pages=pages, selected_page=selected_page)

# Admin page route
@app.route("/admin")
def admin():
    fig1 = make_cof()
    fig2 = make_birdseye()
    figures_html = [
        pio.to_html(fig1, include_plotlyjs=False, full_html=False),
        pio.to_html(fig2, include_plotlyjs=False, full_html=False)
    ]
    return render_template("admin.html", figures_html=figures_html)

# Observer page route
@app.route("/observer")
def observer():
    fig = make_birdseye()
    figure_html = [pio.to_html(fig, include_plotlyjs=False, full_html=False)]
    return render_template("observer.html", figure_html=figure_html)

# Program page route
@app.route("/program")
def program():
    fig = make_cof()
    figure_html = [pio.to_html(fig, include_plotlyjs=False, full_html=False)]
    return render_template("program.html", figure_html=figure_html)

# Entry point for launching the app
def launch_app(pkl_path, config_file):
    global data, manager, all_stars_list

    data = pl.read_star_objects(pkl_path)
    manager = mn.data_admin(config_file)
    manager.run_admin()
    all_stars_list = [obj for obj_list in data[0].values() for obj in obj_list]

    app.run(debug=True)
