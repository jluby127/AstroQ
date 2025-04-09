import kpfcc.scheduler as sch
import kpfcc.requests as rq
import kpfcc.management as mn


def bench(args):

    print("  Benchmark function in driver.py")
<<<<<<< HEAD

    if args.benchmark1:
        print("    Conducting benchmark 1")

    if args.benchmark2:
        print("    Conducting benchmark 2")

=======
    print(f"    Conducting benchmark: {args.benchmark}")

    
    return
    
def schedule(args):
    
    rf = args.request_file
    cf = args.config_file
    type_rf = type(rf)
    type_cf = type(cf)
    print(f'    kpfcc_schedule function: request_file is {rf} and type is {type_rf}')
    print(f'    kpfcc_schedule function: config_file is {cf} and type is {type_cf}')
    
>>>>>>> d467d827c6e6bad210706b5278830a9c98655720
    return

def kpfcc(args):

    print('    Entering kpfcc function in driver.py')

<<<<<<< HEAD
    return
=======
def kpfcc_build(args):
    
    print(f'    kpfcc_build function: boolarg is {args.boolean_argument}')
    
    
    return
    
def kpfcc_plot(args):
    
    so = args.schedule_object
    tp = type(so)
    print(f'    kpfcc_plot function: schedule object is {so} and type is {tp}')
    
    return
    


































>>>>>>> d467d827c6e6bad210706b5278830a9c98655720

def kpfcc_build(args):

    print(f'    kpfcc_build function: boolarg is {args.boolean_argument}')


    return

def kpfcc_schedule(args):

    rf = args.request_file
    print(f'    kpfcc_schedule function: request_file is {rf}')

    cf = args.config_file
    print(f'    kpfcc_schedule function: config_file is {cf}')

    manager = mn.data_admin(cf)
    manager.run_admin()

    request_set = rq.RequestSet(manager)
    request_set.read_from_json(rf)

    schedule = sch.Scheduler(request_set, manager)
    schedule.run_model()
    print("Done solving the schedule.")
    return

def kpfcc_plot(args):

    so = args.schedule_object
    tp = type(so)
    print(f'    kpfcc_plot function: schedule object is {so} and type is {tp}')

    return
