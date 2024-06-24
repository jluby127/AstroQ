import argparse

parser = argparse.ArgumentParser(description='Generate schedules with KPF-CC v2')

parser.add_argument('-d','--schedule_dates',action='append',help='Date(s) to be scheduled as strings in format YYYY-MM-DD. Must be in allocated_nights')
parser.add_argument('-f','--folder', help='Folder to save generated scripts and plots', default='/Users/jack/Desktop/')
parser.add_argument('-t','--time_limit', help='Max time spent optimizing (s)',type=int, default=900)
parser.add_argument('-s','--slot_size', help='The slot size (min)',type=int, default=10)
parser.add_argument('-g','--gurobi_output',action='store_true', help='Activate Gurobi console outputs',default=True)
parser.add_argument('-p','--plot_results',action='store_true', help='Turn on plotting for semester/upcoming quarter nights', default=True)

args = parser.parse_args()

# -------- 2024A Info ---------
deskpath = "/Users/jack/Desktop/"
dirpath = "/Users/jack/Documents/KPF_CC/semesters/2024A/v2_testing/" #'/Users/jack/Documents/Github/optimalAllocation/'
pastDatabase = 'nofilename.csv' #dirpath + "testJump_allObs2024A.csv"
starmap_template_filename = dirpath + "kpfcc_v2_files/KPF_2024A_Template_all.csv"
request_sheet = dirpath + "kpfcc_v2_files/2024A_KPFCC_Requests.csv"
allocated_nights =  dirpath + "kpfcc_v2_files/2024A_Binary_Schedule.txt"
twilight_times =  dirpath + "kpfcc_v2_files/twilight_times_2024A.csv"
turnFile = dirpath + "kpfcc_v2_files/turnOnOffDates_2024A.csv"
access_map = dirpath + "kpfcc_v2_files/2024A_AccessibilityMaps__" + str(args.slot_size) + "minSlot_14Hr.pkl"
nonqueueMap = 'nofilename.csv' #deskpath + 'NonQueueMap_str.txt' #

import sys
sys.path.append(dirpath + "autoschedulerV2/")
import optimalAllocation as oa
oa.runOptimalAllocation(args.schedule_dates,
               request_sheet,
               access_map,
               twilight_times,
               args.folder + 'KPFCC_OptimalAllocation_' + str(args.schedule_dates[0]) + '_Outputs/',
               args.slot_size,
               starmap_template_filename,
               turnFile,
               'nofilename.csv',
               'nofilename.csv',
               nonqueueMap,
               args.gurobi_output,
               args.plot_results,
               args.time_limit)
