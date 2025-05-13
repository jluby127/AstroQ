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
        requests_frame = bn.build_toy_model_from_paper(hours_per_program = 100)
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
        import kpfcc.request as req
        import pandas as pd
        import numpy as np
        req = req.read_json('examples/hello_world/outputs/2024-08-01/request_set.json').strategy
        sch = pd.read_csv('examples/hello_world/outputs/2024-08-01/serialized_outputs_dense.csv')
        
        
        # First, ensure no repeated day/slot pairs (weakness: could allow missing pairs)
        assert sch.groupby(['d','s']).size().max()<=1
        #import pdb; pdb.set_trace()
        
        for star in req.id:
            star_request = req.query(f"id=='{star}'")
            star_schedule = sch.query(f"r=='{star}'") # Only the slots with the star listed
            
        # 1) t-visit: No stars scheduled during another star's slot
            ## - This is already covered by first check above. Overlap isn't possible if all day/slot pairs are unique
            
            
            # import pdb; pdb.set_trace()
        # 2) n_inter_max: Total number of nights a target is scheduled in the semester is less than n_inter_max
            n_inter_max = star_request['n_inter_max'].values[0]
            n_inter_sch = len(set(star_schedule.d)) # All unique nights with scheduled obs

            # Now make sure the number of visits is less than the limit
            n_inter_max_err = ("n_inter_max violated: "
                              f"{star} is scheduled too many times in the semester "
                              f"({n_inter_sch} > {n_inter_max})")
            assert n_inter_sch <= n_inter_max, nim_err
            
            
        # 3) n_intra_min, n_intra_max: N obs per day is between n_intra_min and n_intra_max
            
            # Upper/lower limits on N obs per day
            n_intra_min, n_intra_max = star_request[['n_intra_min', 'n_intra_max']].values[0]
            
            # Scheduled min/max number of obs per day
            n_intra_groupby = star_schedule.groupby(['d']).size()
            n_intra_min_sch, n_intra_max_sch = n_intra_groupby.min(), n_intra_groupby.max()
            
            
            # Ensure the target is never scheduled too few/many times in one night
            n_intra_min_err = ("n_intra_min violated: "
                              f"{star} is scheduled too few times in one night "
                              f"({n_intra_min_sch} obs vs {n_intra_min} obs)")
            assert n_intra_min <= n_intra_min_sch, n_intra_min_err
            
            n_intra_max_err = ("n_intra_max violated: "
                              f"{star} is scheduled too many times in one night "
                              f"({n_intra_max_sch} obs vs {n_intra_max} obs)")
            assert n_intra_max_sch <= n_intra_max, n_intra_max_err
            
            
        # 4) tau_inter: There must be at least tau_inter nights between successive observations of a target over the semester
            tau_inter = star_request[['tau_inter']].values[0] # min num of nights before another obs

            unique_days = np.sort(np.array(list(set(star_schedule.d))))
            max_day_gaps = np.max(unique_days[1:] - unique_days[:-1])
            
            
            if n_inter_max <= 1: # If only 1 obs per semester, no risk of spacing obs too closely
                pass
            else:
                # All gaps are greater than the min gap 
                tau_inter_err = ("tau_inter violated: "
                                f"two obs of {star} are not spaced by enough days "
                                f"({max_day_gaps} days vs {tau_inter} days)")
                assert max_day_gaps >= tau_inter, tau_inter_err
        
        
        # 5) tau_intra: There must be at least tau_intra slots between successive observations of a target in a single night
            tau_intra = star_request[['tau_intra']].values[0] # min num of slots before another obs
            
            max_slot_diffs = star_schedule.groupby('d').s.diff() # Group by day, then find successive differences between slot numbers in the same day. Differences are not computed between the last slot of one night and the first slot of the next night (those values are NaN)
            
            if n_intra_max <= 1: # If only 1 obs per night, no risk of spacing obs too closely
                pass
            else:
                tau_intra_err = ("tau_intra_violated: "
                                f"two obs of {star} are not spaced by enough slots "
                                f"({max_slot_diffs} vs {tau_intra})")
                assert max_slot_diffs >= tau_intra, tau_intra_err
            
    
        


if __name__=="__main__":
    # import pdb; pdb.set_trace()
    tc = TestClass()
    tc.test_requests_vs_schedule()

































    
    