import kpfcc.driver as dr
import kpfcc.request as rq
import kpfcc.benchmarking as bn
import argparse
from configparser import ConfigParser
import os
import kpfcc.management as mn
import kpfcc.scheduler as sch

class TestClass:

    def test_round2_weather(self):
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
        dr.kpfcc_build(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
        dr.schedule(argparse.Namespace(request_file="examples/hello_world/outputs/2024-08-01/request_set.json", config_file='examples/hello_world/config_hello_world.ini'))

    def test_bench(self):
        # override the tests of "do these files exist" and just run it all
        # but run a shortened version, so it doesn't take long, see shortcut=5


        print("Running benchmark test.")

        nR = 5
        nS = 5
        cf = 'examples/bench/config_benchmark.ini'

        # Initialize manager and compute request set on the fly
        # This is a hacky workaround. run_admin needs this file to exist. This can
        # lead to race conditions if benchmarking is run in parallel.
        config = ConfigParser()
        config.read(cf)
        upstream_path = eval(config.get('required', 'folder'), {"os": os})
        semester_directory = upstream_path
        requests_frame = bn.build_toy_model_from_paper(nS, hours_per_program = 100)
        if nR is not None:
            requests_frame = requests_frame.iloc[:nR][::10] # short benchmark
        requests_frame.to_csv(os.path.join(semester_directory, "inputs/Requests.csv"))
        manager = mn.data_admin(cf)
        manager.run_admin()

        # Build observability maps and request set
        print("Building valid indices.")
        strategy, observable = rq.define_indices_for_requests(manager)
        meta = rq.build_meta(cf)
        request_set = rq.RequestSet(meta, strategy, observable)
        current_day = str(config.get('required', 'current_day'))

        # Run the schedule
        schedule = sch.Scheduler(request_set, cf)
        schedule.run_model()

    """
    def test_prep(self):
        # test the creation of all upstream files, including the allocation map
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
    """


    def test_plot(self):
        dr.plot(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    # need environment variable to run
    """
    def test_ob_database_pull(self):
        dr.kpfcc_data(argparse.Namespace(pull_file='examples/pull_file.json', database_file='examples/recreate_paper/'))
    """


    # following five tests has not been tested
    """
    def test_ttp(self):
        dr.ttp(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
    def test_history(self):
        dr.get_history(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

    def test_backups(self):
        dr.backups(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
    def test_oia(self):
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/recreate_paper/oia1/config_oia1.ini'))
        dr.kpfcc_build(argparse.Namespace(config_file='examples/recreate_paper/oia1/config_oia1.ini'))
        dr.schedule(argparse.Namespace(request_file="examples/recreate_paper/oia1/outputs/2024-08-01/request_set.json", config_file='examples/recreate_paper/oia1/config_oia1.ini'))
    """
    
    def test_requests_vs_schedule(self):
        req = 'examples/hello_world/outputs/2024-08-01/request_set.json'
        sch = 'examples/hello_world/outputs/2024-08-01/serialized_outputs_sparse.csv'
        
        dr.requests_vs_schedule(argparse.Namespace(request_file=req, schedule_file=sch))
            
    
        


if __name__=="__main__":
    # import pdb; pdb.set_trace()
    tc = TestClass()
    tc.test_requests_vs_schedule()

































    
    