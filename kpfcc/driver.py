import kpfcc.scheduler as sch
import kpfcc.request as rq
import kpfcc.management as mn


def bench(args):

    print("  Benchmark function in driver.py")

    if args.benchmark1:
        print("    Conducting benchmark 1")

    if args.benchmark2:
        print("    Conducting benchmark 2")

    return

def kpfcc(args):

    print('    Entering kpfcc function in driver.py')


    return

def kpfcc_build(args):

    print(f'    kpfcc_build function: boolarg is {args.boolean_argument}')


    return

def kpfcc_build(args):

    print(f'    kpfcc_build function: boolarg is {args.boolean_argument}')


    return

def kpfcc_schedule(args):

    rf = args.request_file
    print(f'    kpfcc_schedule function: request_file is {rf}')

    cf = args.config_file
    print(f'    kpfcc_schedule function: config_file is {cf}')

    # request_set = rq.RequestSet(manager)
    # request_set.read_from_json(rf)
    request_set = rq.read_json(rf)

    schedule = sch.Scheduler(request_set, cf)
    schedule.run_model()
    print("Done solving the schedule.")
    return

def kpfcc_plot(args):

    so = args.schedule_object
    tp = type(so)
    print(f'    kpfcc_plot function: schedule object is {so} and type is {tp}')

    return
