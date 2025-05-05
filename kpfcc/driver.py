import os
from configparser import ConfigParser
from argparse import Namespace

import kpfcc.scheduler as sch
import kpfcc.request as rq
import kpfcc.management as mn
import kpfcc.benchmarking as bn

def bench(args):
    print("Running benchmark test.")

    nR = args.number_requests
    nS = args.number_slots
    cf = args.config_file

    # Initialize manager and compute request set on the fly
    # This is a hacky workaround. run_admin needs this file to exist. This can
    # lead to race conditions if benchmarking is run in parallel.
    from configparser import ConfigParser
    config = ConfigParser()
    config.read(cf)
    upstream_path = eval(config.get('required', 'folder'), {"os": os})
    semester_directory = upstream_path
    requests_frame = bn.build_toy_model_from_paper(hours_per_program = 100)
    requests_frame = requests_frame.iloc[:nR]#[::10] # downsample to 10x faster
    requests_frame.to_csv(os.path.join(semester_directory, "inputs/Requests.csv"))
    manager = mn.data_admin(cf)
    manager.run_admin()
    
    # Build observability maps and request set
    print("Building valid indices.")
    strategy, observable = rq.define_indices_for_requests(manager)
    meta = rq.build_meta(cf)
    request_set = rq.RequestSet(meta, strategy, observable)
    
    # Limit number of slots if needed
    if nS is not None:
        request_set = bn.set_nSlots_singles(nS, request_set)

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

def schedule(args):

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

    so = args.schedule_object
    tp = type(so)
    print(f'    kpfcc_plot function: schedule object is {so} and type is {tp}')
    print("this function doesn't do anything yet.")
    return
