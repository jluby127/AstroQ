import kpfcc.driver as dr
import kpfcc.request as rq
import kpfcc.benchmarking as bn
import argparse
from configparser import ConfigParser
import os
import kpfcc.management as mn
import kpfcc.scheduler as sch

class TestClass:

    # def test_round2_weather(self):
    #     dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
    #     dr.kpfcc_build(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
    #     dr.schedule(argparse.Namespace(request_file="examples/hello_world/outputs/2024-08-01/request_set.json", config_file='examples/hello_world/config_hello_world.ini'))
    #
    # def test_bench(self):
    #     # override the tests of "do these files exist" and just run it all
    #     # but run a shortened version, so it doesn't take long, see shortcut=5
    #
    #
    #     print("Running benchmark test.")
    #
    #     nR = 5
    #     nS = 5
    #     cf = 'examples/bench/config_benchmark.ini'
    #
    #     # Initialize manager and compute request set on the fly
    #     # This is a hacky workaround. run_admin needs this file to exist. This can
    #     # lead to race conditions if benchmarking is run in parallel.
    #     config = ConfigParser()
    #     config.read(cf)
    #     upstream_path = eval(config.get('required', 'folder'), {"os": os})
    #     semester_directory = upstream_path
    #     requests_frame = bn.build_toy_model_from_paper(hours_per_program = 100)
    #     if nR is not None:
    #         requests_frame = requests_frame.iloc[:nR][::10] # short benchmark
    #     requests_frame.to_csv(os.path.join(semester_directory, "inputs/Requests.csv"))
    #     manager = mn.data_admin(cf)
    #     manager.run_admin()
    #
    #     # Build observability maps and request set
    #     print("Building valid indices.")
    #     strategy, observable = rq.define_indices_for_requests(manager)
    #     meta = rq.build_meta(cf)
    #     request_set = rq.RequestSet(meta, strategy, observable)
    #     current_day = str(config.get('required', 'current_day'))
    #
    #     # Run the schedule
    #     schedule = sch.Scheduler(request_set, cf)
    #     schedule.run_model()

    """
    def test_prep(self):
        # test the creation of all upstream files, including the allocation map
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
    """


    # def test_plot(self):
    #     dr.plot(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))

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
        req = req.read_json('examples/hello_world/outputs/2024-08-01/request_set.json').strategy
        sch = pd.read_csv('examples/hello_world/outputs/2024-08-01/serialized_outputs_dense.csv')
        
        
        # First, ensure no repeated day/slot pairs (weakness: could allow missing pairs)
        assert sch.groupby(['d','s']).size().max()<=1
        #import pdb; pdb.set_trace()
        
        for star in req.id:
            star_schedule = req.query(f"id=='{star}'") # Only the slots with the star listed
            
        # 1) No stars scheduled during another star's slot
            ## - This is already covered by first check above. Overlap isn't possible if all day/slot pairs are unique
            
            
        # 2) Total number of nights in the semester is less than n_inter_max
            n_inter_max = star_schedule['n_inter_max']

            n_inter_sch = len(star_schedule)/t_visit # Should be evenly divisible
            
            assert np.isclose(n_inter_sch, int(n_inter_sch)) # First, require that n_inter_sch is an integer. If not, then it is not the case that every visit takes up t_visit slots
            assert n_inter_sch <= n_inter_max # Now make sure the number of visits is less than the limit
            
            
            
        # 3) N obs per day is between n_intra_min and n_intra_max
            
            # Upper/lower limits on N obs per day
            n_intra_min, n_intra_max = star_schedule[['n_intra_min', 'n_intra_max']].values[0]
            
            # Scheduled min/max number of obs per day
            n_intra_groupby = sch.query(f"id=='{star}'").groupby(['d']).size()
            n_intra_min_sch, n_intra_max_sch = n_intra_groupby.min(), n_intra_groupby.max()
            
            import pdb; pdb.set_trace()
            assert (n_intra_min <= n_intra_min_sch) & (n_intra_max_sch <= n_intra_max)
            
        
            
            
        
        
        # def visit_identifier(starname, schedule):
        #     """
        #     Find and group the slot indices of every
        #     visit to a star in a schedule csv
        #     """
        #
        #     #import pdb; pdb.set_trace()
        #     all_slots = schedule.query(f"r=='{starname}'") # All slots assigned to this star
        #
        #     all_visits = []
        #     for day in all_slots.d.values: # Each separate day
        #         same_day_slots = all_slots.query(f"d=={day}").index # Indices of all slots in that day
        #
        #         same_day_visits = [[same_day_slots[0]]] # The first visit is the first filled slot
        #
        #         for slot in same_day_slots[1:]:
        #             if slot == visits[-1][-1] + 1:
        #                 same_day_visits[-1].append(slot) # If it's consecutive, add to the latest visit
        #             else:
        #                 same_day_visits.append([slot]) # If not, create new visit
        #
        #         all_visits.append([day, same_day_visits])
        #
        #     return all_visits
        #
        # # This block assumes each star is requested with the same requirements (t_visit, etc.). This need not be true. To modify, include in schedule file an index to link each obs to a specific request
        # for starname in req.id:
        #
        #     all_visits = visit_identifier(starname, sch)
        #     import pdb; pdb.set_trace()
    
        


if __name__=="__main__":
    # import pdb; pdb.set_trace()
    tc = TestClass()
    tc.test_requests_vs_schedule()

































    
    