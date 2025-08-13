import astroq.driver as dr
import astroq.benchmarking as bn
import argparse
from configparser import ConfigParser
import os
import astroq.splan as splan
import astroq.plot as pl
import astroq.nplan as nplan
import unittest
from pathlib import Path
from io import BytesIO
import imageio.v3 as iio

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
        dr.kpfcc_plan_semester(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini', run_band3=False))

    def test_round2_weather(self):
        dr.kpfcc_plan_semester(argparse.Namespace(config_file='examples/hello_world/config_hello_world_bonus_weather.ini', run_band3=False))

    def test_ttp(self):
        dr.ttp(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    def test_bench(self):
        dr.bench(argparse.Namespace(config_file='examples/bench/config_benchmark.ini', number_slots=12, thin=10))
        dr.ttp(argparse.Namespace(config_file='examples/bench/config_benchmark.ini'))
        dr.plot(argparse.Namespace(config_file='examples/bench/config_benchmark.ini'))

    def test_prep(self):
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world_prep.ini', allo_source='db', past_source='db', band3_program_code='2025B_E473'))

    def test_plot(self):
        dr.plot(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

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
