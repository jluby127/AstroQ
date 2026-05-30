import astroq.driver as dr
import argparse
import astroq.splan as splan
import astroq.nplan as nplan
from astroq.ttp.model import TTPModel
import unittest
from astropy.coordinates import SkyCoord
import os
from astroq.webapp import launch_app, app

class TestClass(unittest.TestCase):

    def test01_helloworld(self):
        dr.plan_semester(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini', run_band3=False))

    def test02_round2_weather(self):
        dr.plan_semester(argparse.Namespace(config_file='examples/hello_world/config_hello_world_bonus_weather.ini', run_band3=False))

    def test03_plan_night(self):
        dr.plan_night(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    def test04_bench(self):
        dr.bench(argparse.Namespace(config_file='examples/bench/config_benchmark.ini', number_slots=12, thin=10))
        dr.plan_night(argparse.Namespace(config_file='examples/bench/config_benchmark.ini'))
        dr.plot(argparse.Namespace(config_file='examples/bench/config_benchmark.ini'))

    def test05_generic_prep(self):
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world_prep.ini', allo_source='examples/hello_world/prepped/observatory_schedule.csv', past_source='examples/hello_world/prepped/jump_past_history.csv', request_source='examples/hello_world/prepped/request.csv', filler_programs='2025B_E473', band_number=1, is_full_band=False))
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world_prep.ini', allo_source='examples/hello_world/prepped/observatory_schedule.csv', past_source='examples/hello_world/prepped/jump_past_history.csv', request_source='examples/hello_world/prepped/request.csv', filler_programs='2025B_E473', band_number=3, is_full_band=True))

    def test06_kpfcc_prep(self):
       dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world_prep.ini', allo_source='db', past_source='db', request_source='db', filler_programs='2025B_E473', band_number=1, is_full_band=True))

    def test07_plot(self):
        dr.plot(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    def test08_requests_vs_schedule(self):
        sch = 'examples/hello_world/2018B/2018-08-05/band1/outputs/semester_plan.csv'
        dr.requests_vs_schedule(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini', schedule_file=sch))

if __name__=="__main__":
    unittest.main()
