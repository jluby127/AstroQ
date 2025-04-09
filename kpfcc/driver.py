


def bench(args):
    
    print("  Benchmark function in driver.py")
    print(f"    Conducting benchmark: {args.benchmark}")

    
    return
    
def schedule(args):
    
    rf = args.request_file
    cf = args.config_file
    type_rf = type(rf)
    type_cf = type(cf)
    print(f'    kpfcc_schedule function: request_file is {rf} and type is {type_rf}')
    print(f'    kpfcc_schedule function: config_file is {cf} and type is {type_cf}')
    
    return

def kpfcc(args):
    
    print('    Entering kpfcc function in driver.py')
    
    return

def kpfcc_build(args):
    
    print(f'    kpfcc_build function: boolarg is {args.boolean_argument}')
    
    
    return
    
def kpfcc_plot(args):
    
    so = args.schedule_object
    tp = type(so)
    print(f'    kpfcc_plot function: schedule object is {so} and type is {tp}')
    
    return
    















































