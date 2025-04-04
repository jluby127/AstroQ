


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
    
def kpfcc_schedule(args):
    
    rf = args.request_file
    tp = type(rf)
    print(f'    kpfcc_schedule function: request_file is {rf} and type is {tp}')
    
    return
    
def kpfcc_plot(args):
    
    so = args.schedule_object
    tp = type(so)
    print(f'    kpfcc_plot function: schedule object is {so} and type is {tp}')
    
    return















































