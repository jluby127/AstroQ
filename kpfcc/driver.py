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
    requests_frame = bn.build_toy_model_from_paper(hours_per_program = 100)
    if nR is not None:
        requests_frame = requests_frame.iloc[:nR]#[::10]
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

def backups(args):

    cf = args.config_file
    print(f'    kpfcc_schedule function: config_file is {cf}')

    manager = mn.data_admin(cf)
    manager.run_admin()

    sk.produce_bright_backups(manager)

def get_history(args):

    cf = args.config_file
    print(f'    kpfcc_schedule function: config_file is {cf}')

    manager = mn.data_admin(cf)
    database_info_dict = hs.build_past_history(manager.past_database_file, manager.requests_frame, manager.twilight_frame)
