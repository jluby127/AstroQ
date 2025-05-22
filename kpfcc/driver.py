import os
import json
import numpy as np
import pandas as pd
from configparser import ConfigParser
from argparse import Namespace

import kpfcc.scheduler as sch
import kpfcc.request as rq
import kpfcc.management as mn
import kpfcc.benchmarking as bn
import kpfcc.blocks as ob
import kpfcc.plot as pl
import kpfcc.onsky as sk
import kpfcc.history as hs

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.ERROR)


def bench(args):
    print("Running benchmark test.")

    nR = args.number_requests
    nS = args.number_slots
    cf = args.config_file

    # Initialize manager and compute request set on the fly
    # This is a hacky workaround. run_admin needs this file to exist. This can
    # lead to race conditions if benchmarking is run in parallel.
    config = ConfigParser()
    config.read(cf)
    upstream_path = eval(config.get('required', 'folder'), {"os": os})
    semester_directory = upstream_path
    requests_frame = bn.build_toy_model_from_paper(nS)
    requests_frame.to_csv(os.path.join(semester_directory, "inputs/Requests.csv"))
    manager = mn.data_admin(cf)
    manager.run_admin()

    # Build observability maps and request set
    print("Building valid indices.")
    strategy, observable = rq.define_indices_for_requests(manager)
    meta = rq.build_meta(cf)
    request_set = rq.RequestSet(meta, strategy, observable)
    # print out the last rows of strategy to ensure the size of the model looks right
    request_set.to_json(os.path.join(manager.output_directory, "request_set.json"))
    # Run the schedule
    schedule = sch.Scheduler(request_set, cf)
    schedule.run_model()
    print("Done solving the schedule.")
    return

def kpfcc(args):
    print('    Entering kpfcc function in driver.py')
    print("this function doesn't do anything yet.")
    return

def kpfcc_build(args):

    cf = args.config_file
    print(f'    kpfcc_schedule function: config_file is {cf}')

    manager = mn.data_admin(cf)
    manager.run_admin()
    print("Building valid indices.")
    strategy, observable = rq.define_indices_for_requests(manager)
    meta = rq.build_meta(cf)
    request_set = rq.RequestSet(meta, strategy, observable)
    request_set.to_json(manager.output_directory + "request_set.json")
    return

def kpfcc_prep(args):
    cf = args.config_file
    print(f'    kpfcc_schedule function: config_file is {cf}')

    mn.prepare_new_semester(cf)
    return

def kpfcc_data(args):

    pull_file = args.pull_file
    print(f'    kpfcc_data function: using pull info from {pull_file}')

    savepath = args.database_file
    print(f'    kpfcc_data function: saving to {savepath}')

    with open(pull_file, "r") as f:
        data = json.load(f)
    semester = data["semester"]
    awarded_programs = data["awarded_programs"]

    if not os.path.exists(savepath):
        os.makedirs(savepath)

    OBs = ob.refresh_local_data(semester)
    good_obs, bad_obs_values, bad_obs_hasFields = ob.get_request_sheet(OBs, awarded_programs, savepath + "/Requests.csv")

    send_emails_with = []
    for i in range(len(bad_obs_values)):
        if bad_obs_values['metadata.semid'][i] in awarded_programs:
            send_emails_with.append(ob.inspect_row(bad_obs_hasFields, bad_obs_values, i))

    '''
    this is where code to automatically send emails will go.
    '''

    return

def schedule(args):
    
    log.debug("Debug message")
    log.info("Info message")
    log.warning("Warning message here")
    log.error("Error message")
    log.critical("Critical message")

    rf = args.request_file
    print(f'    kpfcc_schedule function: request_file is {rf}')
    cf = args.config_file
    print(f'    kpfcc_schedule function: config_file is {cf}')

    request_set = rq.read_json(rf)
    schedule = sch.Scheduler(request_set, cf)
    schedule.run_model()
    print("Done solving the schedule.")
    return

def plot(args):

    # so = args.schedule_object
    # tp = type(so)
    # print(f'    kpfcc_plot function: schedule object is {so} and type is {tp}')
    # print("this function doesn't do anything yet.")

    cf = args.config_file
    print(f'    kpfcc_schedule function: config_file is {cf}')
    pl.run_plot_suite(cf)

    return

def ttp(args):

    cf = args.config_file
    print(f'    kpfcc_schedule function: config_file is {cf}')

    manager = mn.data_admin(cf)
    manager.run_admin()

    sk.run_ttp(manager)
    sk.produce_bright_backups(manager)

def get_history(args):

    cf = args.config_file
    print(f'    kpfcc_schedule function: config_file is {cf}')

    manager = mn.data_admin(cf)
    database_info_dict = hs.build_past_history(manager.past_database_file, manager.requests_frame, manager.twilight_frame)


def requests_vs_schedule(args):
    
    rf = args.request_file
    sf = args.schedule_file
    
    req = rq.read_json(rf).strategy
    sch = pd.read_csv(sf)
    sch = sch.sort_values(by=['d', 's']).reset_index(drop=True) # Re-order into the real schedule
    
    # First, ensure no repeated day/slot pairs (does allow missing pairs)
    no_duplicate_slot_err = ("'No duplicate slot' condition violated: "
                             "At least one pair of rows corresponds to "
                             "the same day and slot.")
    assert sch.groupby(['d','s']).size().max()<=1, no_duplicate_slot_err

    for star in req.id:
        star_request = req.query(f"id=='{star}'")
        star_schedule = sch.query(f"r=='{star}'") # Only the slots with the star listed
        
    # 1) t_visit: No stars scheduled during another star's slot

        t_visit = star_request.t_visit.values[0] # Number of slots needed to complete observation
        star_inds = star_schedule.index
        day_slot = sch[['d', 's']]
        
        # Check the number of slots between consecutive obs. If they're on the same day, demand a minimum separation
        if star_inds.max() == day_slot.index.max(): # Special case: if this star includes the last obs in the whole schedule
            star_inds = star_inds[:-1] # Exclude the very last observation to avoid index err. That obs can't be overlapped by a later target anyway
            
        day_slot_diffs = day_slot.iloc[star_inds+1].reset_index() - day_slot.iloc[star_inds].reset_index()

        if len(day_slot_diffs.query('d==0'))==0: # If the target is always the last obs of the night, pass
            pass
        else:
            assert day_slot_diffs.query('d==0').s.min() >= t_visit-1, f"{star}"
        

    # 2) n_inter_max: Total number of nights a target is scheduled in the semester is less than n_inter_max
        n_inter_max = star_request['n_inter_max'].values[0]
        n_inter_sch = len(set(star_schedule.d)) # All unique nights with scheduled obs

        # Now make sure the number of visits is less than the limit
        n_inter_max_err = ("n_inter_max violated: "
                          f"{star} is scheduled too many times in the semester "
                          f"({n_inter_sch} > {n_inter_max})")
        assert n_inter_sch <= n_inter_max, n_inter_max_err
        
        
    # 3) n_intra_min, n_intra_max: N obs per day is between n_intra_min and n_intra_max
    
        # t_visit, the number of slots required to complete a single observation (aka visit)
        t_visit = req[req.id==star].t_visit.values
        
        # Upper/lower limits on N obs per day
        n_intra_min, n_intra_max = star_request[['n_intra_min', 'n_intra_max']].values[0]
        
        # Scheduled min/max number of obs per day
        n_intra_groupby = star_schedule.groupby(['d']).size() # The numerator gives the sum of all starting slots in which the target is observed in a day.
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
        min_day_gaps = np.min(unique_days[1:] - unique_days[:-1])
        
        
        if n_inter_max <= 1: # If only 1 obs per semester, no risk of spacing obs too closely
            pass
        else:
            # Require that all gaps are greater than the min gap 
            tau_inter_err = ("tau_inter violated: "
                            f"two obs of {star} are not spaced by enough days "
                            f"({min_day_gaps} days vs {tau_inter} days)")
            assert min_day_gaps >= tau_inter, tau_inter_err
    
    
    # 5) tau_intra: There must be at least tau_intra slots between successive observations of a target in a single night
        tau_intra = star_request[['tau_intra']].values[0] # min num of slots before another obs
        
        min_slot_diffs = star_schedule.groupby('d').s.diff().min() # Group by day, then find successive differences between slot numbers in the same day. Differences are not computed between the last slot of one night and the first slot of the next night (those values are NaN). The differences must all be AT LEAST tau_intra.
        
        if n_intra_max <= 1: # If only 1 obs per night, no risk of spacing obs too closely
            pass
        else:
            tau_intra_err = ("tau_intra_violated: "
                            f"two obs of {star} are not spaced by enough slots "
                            f"({min_slot_diffs} vs {tau_intra})")

            assert min_slot_diffs >= tau_intra, tau_intra_err
            
            

