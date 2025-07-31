import astroq.driver as dr
import astroq.benchmarking as bn
import argparse
from configparser import ConfigParser
import os
import astroq.management as mn
import astroq.splan as splan
import astroq.plot as pl
import astroq.dynamic as dn
import astroq.nplan as nplan
import unittest
from pathlib import Path

import multiprocessing
import time
import requests
import os
from astroq.webapp import launch_app, app
import threading

class TestClass(unittest.TestCase):

    def test_cli(self):
        dr.kpfcc(argparse.Namespace(kpfcc_subcommand=None))

    def test_helloworld(self):
        dr.kpfcc_plan_semester(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    # def test_round2_weather(self):
    #     dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world_yesbonus.ini'))
    #     dr.kpfcc_plan_semester(argparse.Namespace(config_file='examples/hello_world/config_hello_world_yesbonus.ini'))

    # def test_bench(self):
    #     dr.bench(argparse.Namespace(config_file='examples/bench/config_benchmark.ini', number_slots=12, thin=10))

    # def test_prep(self):
    #     dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/test_config.ini', allo_source='db', past_source='db'))

    def test_plot(self):
        #dr.plot_pkl(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
        #dr.plot_static(argparse.Namespace(path_to_semester_planner='examples/hello_world/outputs/semester_planner.pkl'))
        cf = 'examples/hello_world/config_hello_world.ini'
        semester_planner = splan.SemesterPlanner(cf)
        semester_planner.run_model()
        
        # Add the future_forecast attribute that the plot functions expect
        semester_planner.future_forecast = os.path.join(semester_planner.output_directory, "semester_plan.csv")
        
        data_astroq = pl.process_stars(semester_planner)
        saveout = semester_planner.output_directory + "/static_plots/"
        os.makedirs(saveout, exist_ok=True)

        # build the interactive plots
        fig_cof = dn.get_cof(semester_planner, list(data_astroq[1].values()))
        fig_birdseye = dn.get_birdseye(semester_planner, data_astroq[2], list(data_astroq[1].values()))
        # write the static versions to the reports directory
        fig_cof.write_image(os.path.join(saveout, "all_programs_COF.png"), width=1200, height=800)
        fig_birdseye.write_image(os.path.join(saveout, "all_stars_birdseye.png"), width=1200, height=800)

        # Use the new NightPlanner class for object-oriented night planning
        night_planner = nplan.NightPlanner(cf)
        night_planner.run_ttp()
        
        # Read the TTP data that was created
        ttp_path = os.path.join(night_planner.reports_directory, "ttp_data.pkl")
        if os.path.exists(ttp_path):
            data_ttp = pl.read_star_objects(ttp_path)
            # build the interactive plots
            script_table_df = dn.get_script_plan(cf, data_ttp)
            ladder_fig = dn.get_ladder(data_ttp)
            slew_animation_html = dn.get_slew_animation(data_ttp, animationStep=120)
            slew_path_fig = dn.plot_path_2D_interactive(data_ttp)
            # write the static versions to the reports directory
            script_table_df.to_csv(os.path.join(saveout, "script_table.csv"), index=False)
            ladder_fig.write_image(os.path.join(saveout, "ladder_plot.png"), width=1200, height=800)
            slew_path_fig.write_image(os.path.join(saveout, "slew_path_plot.png"), width=1200, height=800)
        


    # def test_ob_database_pull(self):
    #     dr.kpfcc_data(argparse.Namespace(pull_file='examples/pull_file.json', database_file='examples/recreate_paper/'))

    # def test_ttp(self):
    #     dr.ttp(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    # def test_history(self):
    #     dr.get_history(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    # def test_dynamic_plotting(self):
    #     dr.get_dynamics(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    # def test_webapp(self):
    #     config_path = 'examples/hello_world/config_hello_world.ini'

    #     # Launch Flask in a thread
    #     def run():
    #         launch_app(config_path, flag=True)

    #     thread = threading.Thread(target=run, daemon=True)
    #     thread.start()

    #     time.sleep(3)  # Wait for server to be ready
    #     try:
    #         response = requests.get("http://127.0.0.1:5000")
    #         assert response.status_code == 200
    #     finally:
    #         # Gracefully shut down the Flask app
    #         try:
    #             requests.get("http://127.0.0.1:5000/shutdown")
    #         except requests.exceptions.RequestException:
    #             pass  # The server might already be down
    #         thread.join(timeout=5)

    # def test_requests_vs_schedule(self):
    #     sch = 'examples/hello_world/outputs/semester_plan.csv'
    #     dr.requests_vs_schedule(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini', schedule_file=sch))

    # this is not working right now.
    # def test_simulate_history(self):
    #     dr.make_simulated_history(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    # we don't care about OIA yet
    # def test_oia(self):
    #     dr.kpfcc_prep(argparse.Namespace(config_file='examples/recreate_paper/oia1/config_oia1.ini'))
    #     dr.kpfcc_build(argparse.Namespace(config_file='examples/recreate_paper/oia1/config_oia1.ini'))
    #     dr.schedule(argparse.Namespace(request_file="examples/recreate_paper/oia1/outputs/2024-08-02/request_set.json", config_file='examples/recreate_paper/oia1/config_oia1.ini'))

if __name__=="__main__":
    unittest.main()
