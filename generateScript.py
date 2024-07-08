import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta
import argparse

parser = argparse.ArgumentParser(description='Generate schedules with KPF-CC v2')

parser.add_argument('-d','--schedule_dates',action='append',help='Date(s) to be scheduled as strings in format YYYY-MM-DD. Must be in allocated_nights')
parser.add_argument('-f','--folder', help='Folder to save generated scripts and plots', default='/Users/jack/Desktop/')
parser.add_argument('-r','--run_extra_rounds',action='store_true', help='Turn on the additional rounds of scheduling', default=False)
parser.add_argument('-t','--time_limit', help='Max time spent optimizing (s)',type=int, default=300)
parser.add_argument('-s','--slot_size', help='The slot size (min)',type=int, default=10)
parser.add_argument('-g','--gurobi_output',action='store_true', help='Activate Gurobi console outputs',default=True)
parser.add_argument('-p','--plot_results',action='store_true', help='Turn on plotting for semester/upcoming quarter nights', default=True)

args = parser.parse_args()

# -------- 2024A Info ---------
deskpath = "/Users/jack/Desktop/"
dirpath = "/Users/jack/Documents/Github/optimalAllocation/"
pastDatabase = args.folder + 'KPFCC_' + str(args.schedule_dates[0]) + '_Outputs/' + "queryJumpDatabase.csv" #deskpath + "Jump_allObs2024A.csv" # 'nofilename.csv'
starmap_template_filename = dirpath + "semesterFiles/2024A/KPF_2024A_Template_all.csv"
request_sheet = dirpath + "semesterFiles/2024A/2024A_KPFCC_Requests.csv"
allocated_nights =  dirpath + "semesterFiles/2024A/2024A_Binary_Schedule.txt"
twilight_times =  dirpath + "semesterFiles/2024A/twilight_times_2024A.csv"
turnFile = dirpath + "semesterFiles/2024A/turnOnOffDates_2024A.csv"
access_map = dirpath + "semesterFiles/2024A/2024A_AccessibilityMaps__" + str(args.slot_size) + "minSlot_14Hr.pkl"
nonqueueMap =  'nofilename.csv'

import sys
sys.path.append(dirpath + "autoschedulerV2/")
import solveSemester as ssm
ssm.runKPFCCv2(args.schedule_dates,
                          request_sheet,
                          allocated_nights,
                          access_map,
                          twilight_times,
                          args.folder + 'KPFCC_' + str(args.schedule_dates[0]) + '_Outputs/',
                          args.slot_size,
                          args.run_extra_rounds,
                          pastDatabase,
                          starmap_template_filename,
                          turnFile,
                          nonqueueMap,
                          args.gurobi_output,
                          args.plot_results,
                          args.time_limit)


import processingFunctions as pf
import helperFunctions as hf
sys.path.append('/Users/jack/Documents/Github/ttp/')
import formatting
import telescope
import plotting
import model
savepath = args.folder + 'KPFCC_' + str(args.schedule_dates[0]) + '_Outputs/'
print("Prepare schedule for the TTP.")
tel = telescope.Keck1()
startstoptimes = pd.read_csv(dirpath + 'semesterFiles/2024A/Nightly_StartStop_Times.csv')
semester_start_date, semester_end_date, semesterLength, semesterYear, semesterLetter = hf.getSemesterInfo(args.schedule_dates[0])
all_dates_dict = hf.buildDayDateDictionary(semester_start_date, semesterLength)
the_schedule = np.loadtxt(savepath + 'raw_combined_semester_schedule_Round2.txt', delimiter=',', dtype=str)
for n in range(len(args.schedule_dates)):

    dayInSemester = all_dates_dict[args.schedule_dates[n]]
    idx = startstoptimes.index[startstoptimes['Date']==str(args.schedule_dates[n])][0]

    filltargets = np.loadtxt(savepath + 'gapFillerTargets.txt', dtype=str)
    toTTP = pf.prepareTTP(request_sheet, the_schedule[dayInSemester], filltargets)
    filename = savepath + '/Selected_' + str(args.schedule_dates[n]) + ".txt"
    toTTP.to_csv(filename, index=False)
    targlist = formatting.theTTP(filename)

    starttime = startstoptimes['Start'][idx]
    stoptime = startstoptimes['Stop'][idx]
    startObs = Time(str(args.schedule_dates[n]) + "T" + str(starttime), format='isot')
    endObs = Time(str(args.schedule_dates[n]) + "T" + str(stoptime), format='isot')
    total_time = np.round((endObs.jd-startObs.jd)*24,3)
    print("Time in Night for Observations: " + str(total_time) + " hours.")

    solution = model.TTPModel(startObs, endObs, targlist, tel)
    solutionDict = solution.schedule

    orderedList = []
    for i in range(len(solutionDict['Starname'])):
        orderedList.append(solutionDict['Starname'][i])

    plotting.writeStarList(solution.plotly, startObs, args.schedule_dates[n], outputdir=savepath)
    plotting.plot_path_2D(solution,outputdir=args.folder + 'KPFCC_' + str(args.schedule_dates[0]) + '_Outputs')
    plotting.nightPlan(solution.plotly, args.schedule_dates[n], outputdir=savepath)

    obs_and_times = pd.read_csv(savepath + 'ObserveOrder_' + str(args.schedule_dates[0]) + ".txt")
    all_targets_frame = pd.read_csv(request_sheet)
    gapFillers = np.loadtxt(savepath + 'gapFillerTargets.txt', delimiter=',', dtype=str)
    pf.write_starlist(all_targets_frame, obs_and_times, solution.extras, gapFillers, 'nominal', str(args.schedule_dates[0]), savepath)
