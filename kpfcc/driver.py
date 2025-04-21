import kpfcc.scheduler as sch
import kpfcc.request as rq
import kpfcc.management as mn
import kpfcc.benchmarking as bn


def bench(args):
    print("Running benchmark test.")

    nR = args.number_requests
    nS = args.number_slots
    cf = args.config_file

    print("Checking for toy model files.")
    rf = bn.do_benchmark_files_exist(cf)
    request_set = rq.read_json(rf)
    print("Parsing down size of model for desired test.")
    request_set = bn.firstN_Requests(nR, request_set)
    request_set = bn.set_nSlots_singles(nS, request_set)

    print("Running solver.")
    schedule = sch.Scheduler(request_set, cf)

    return    

def kpfcc(args):

    print('    Entering kpfcc function in driver.py')


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
    request_set.to_json(manager.output_directory)
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

    return
