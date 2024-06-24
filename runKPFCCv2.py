import argparse

parser = argparse.ArgumentParser(description='Generate schedules with KPF-CC v2')

parser.add_argument('-d','--schedule_dates',action='append',help='Date(s) to be scheduled as strings in format YYYY-MM-DD. Must be in allocated_nights')
parser.add_argument('-f','--folder', help='Folder to save generated scripts and plots', default='/Users/jack/Desktop/')
parser.add_argument('-l','--time_limit', help='Max time spent optimizing (s)',type=int, default=300)
parser.add_argument('-s','--slot_size', help='The slot size (min)',type=int, default=10)
parser.add_argument('-b','--enforcedNO', help='The quarters that must not be chosen in optimal allocation', type=str, default='nofilename.csv')
parser.add_argument('-w','--enforcedYES', help='The quarters that yes must be chosen in optimal allocation', type=str, default='nofilename.csv')
parser.add_argument('-nq','--enforcedNonQueue', help='The slots that must not be chosen because they are allocated to Non Queue observations', type=str, default='nofilename.csv')


parser.add_argument('-g','--gurobi_output',action='store_true', help='Activate Gurobi console outputs',default=True)
parser.add_argument('-p','--plot_results',action='store_true', help='Turn on plotting for semester/upcoming quarter nights', default=True)
parser.add_argument('-o','--run_optimal_allocation',action='store_true', help='Turn on the optimal allocation determination', default=False)
parser.add_argument('-r','--run_extra_rounds',action='store_true', help='Turn on the additional rounds of scheduling', default=False)

args = parser.parse_args()

past = True
if past == True:
    # -------- 2024A Info ---------
    deskpath = "/Users/jack/Desktop/"
    dirpath = "/Users/jack/Desktop/" #'/Users/jack/Documents/Github/optimalAllocation/'
    pastDatabase = deskpath + "testJump_allObs2024A.csv" #'nofilename.csv' #
    starmap_template_filename = dirpath + "kpfcc_v2_files/KPF_2024A_Template_all.csv"
    request_sheet = dirpath + "kpfcc_v2_files/2024A_KPFCC_Requests.csv"
    allocated_nights =  dirpath + "kpfcc_v2_files/2024A_Binary_Schedule.txt"
    twilight_times =  dirpath + "kpfcc_v2_files/twilight_times_2024A.csv"
    turnFile = dirpath + "kpfcc_v2_files/turnOnOffDates_2024A.csv"
    access_map = dirpath + "kpfcc_v2_files/2024A_AccessibilityMaps__" + str(args.slot_size) + "minSlot_14Hr.pkl"
    enforcedNoDatesfile = 'nofilename.csv'
    enforcedYesDatesfile = 'nofilename.csv'
    nonqueueMap = 'nofilename.csv'
    nonqueueMap_str = 'nofilename.csv'
    nonqueueObs_info = 'nofilename.csv'

if past == False:
    # -------- 2024B Info ---------
    deskpath = "/Users/jack/Desktop/"
    dirpath = '/Users/jack/Documents/Github/optimalAllocation/'
    pastDatabase = 'nofilename.csv'
    starmap_template_filename = dirpath + "kpfcc_v2_files/KPF_2024B_Template_all.csv"
    request_sheet = dirpath + "kpfcc_v2_files/2024B_KPFCC_Requests.csv"
    allocated_nights =  deskpath + 'V9_Binary_Map.txt' #dirpath + 'nofilename.csv'
    twilight_times =  dirpath + "kpfcc_v2_files/twilight_times_2024B.csv"
    turnFile = "kpfcc_v2_files/turnOnOffDates_2024B.csv"
    access_map = dirpath + "kpfcc_v2_files/2024B_AccessibilityMaps__" + str(args.slot_size) + "minSlot_14Hr.pkl"
    enforcedNoDatesfile = 'nofilename.csv' #deskpath + "enforcedNODates_V8.csv"
    enforcedYesDatesfile = 'nofilename.csv' #deskpath + "enforcedYESDates_V8.csv"
    nonqueueMap = 'nofilename.csv'
    nonqueueMap_str = 'nofilename.csv'
    nonqueueObs_info = 'nofilename.csv'

# import sys
# sys.path.append(dirpath)
import kpfcc_v2_parse
kpfcc_v2_parse.runKPFCCv2(args.schedule_dates, request_sheet, allocated_nights, access_map, twilight_times,
                           args.folder + 'KPFCC_' + str(args.schedule_dates[0]) + '_Outputs/', args.slot_size,
                           args.run_optimal_allocation, args.run_extra_rounds,
                           pastDatabase, starmap_template_filename, turnFile,
                           enforcedNoDatesfile, enforcedYesDatesfile, nonqueueMap, nonqueueMap_str, nonqueueObs_info,
                           args.gurobi_output, args.plot_results, args.time_limit)
