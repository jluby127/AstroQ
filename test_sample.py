import astroq.driver as dr
import astroq.benchmarking as bn
import argparse
from configparser import ConfigParser
import os
import astroq.management as mn
import astroq.splan as splan
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

    def test_round2_weather(self):
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world_yesbonus.ini'))
        dr.kpfcc_plan_semester(argparse.Namespace(config_file='examples/hello_world/config_hello_world_yesbonus.ini'))

    def test_bench(self):
        dr.bench(argparse.Namespace(config_file='examples/bench/config_benchmark.ini', number_slots=12, thin=10))

    def test_prep(self):
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/test_config.ini', allo_source='db', past_source='db'))

    def test_plot(self):
        dr.plot_pkl(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
        dr.plot_static(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    def test_ob_database_pull(self):
        dr.kpfcc_data(argparse.Namespace(pull_file='examples/pull_file.json', database_file='examples/recreate_paper/'))

    def test_ttp(self):
        dr.ttp(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    def test_history(self):
        dr.get_history(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    def test_dynamic_plotting(self):
        dr.get_dynamics(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    def test_webapp(self):
        config_path = 'examples/hello_world/config_hello_world.ini'

        # Launch Flask in a thread
        def run():
            launch_app(config_path, flag=True)

        thread = threading.Thread(target=run, daemon=True)
        thread.start()

        time.sleep(3)  # Wait for server to be ready
        try:
            response = requests.get("http://127.0.0.1:5000")
            assert response.status_code == 200
        finally:
            # Gracefully shut down the Flask app
            try:
                requests.get("http://127.0.0.1:5000/shutdown")
            except requests.exceptions.RequestException:
                pass  # The server might already be down
            thread.join(timeout=5)

    def test_requests_vs_schedule(self):
        sch = 'examples/hello_world/outputs/semester_plan.csv'
        dr.requests_vs_schedule(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini', schedule_file=sch))

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
