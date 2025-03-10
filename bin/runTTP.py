import sys
import pandas as pd
import numpy as np
from astropy.time import Time
from astropy.time import TimeDelta
sys.path.append('/Users/jack/Documents/Github/ttp/')
sys.path.append('/Users/jack/Documents/Github/optimalAllocation/autoschedulerV2/')
import processingFunctions as pf
import helperFunctions as hf
import formatting
import telescope
import plotting
import model

# NOTE: in order to run this file, you must have already run the generateScript.py file and produced a folder with all the output files

import argparse
parser = argparse.ArgumentParser(description='Generate schedules with KPF-CC v2')
parser.add_argument('-d','--schedule_dates',action='append',help='Date(s) to be scheduled as strings in format YYYY-MM-DD. Must be in allocated_nights')
parser.add_argument('-f','--folder', help='Folder to save generated scripts and plots', default=os.environ["KPFCC_SAVE_PATH"])
args = parser.parse_args()

dirpath = "/Users/jack/Documents/Github/optimalAllocation/"
savepath = args.folder + 'outputs/' + str(current_day[0]) + '_ttpOnly/'
request_sheet = args.folder + "inputs/Requests.csv"
startstop_sheet = args.folder + "inputs/2024B_NightlyStartStopTimes.csv"

# note to self: fix all the paths and environment variables
# add environment variable to the autoscheduler modules 

print("Prepare schedule for the TTP.")
tel = telescope.Keck1()
startstoptimes = pd.read_csv(startstop_sheet)
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

    solution = model.TTPModel(startObs, endObs, targlist, tel, savepath)
    solutionDict = solution.schedule

    orderedList = []
    for i in range(len(solutionDict['Starname'])):
        orderedList.append(solutionDict['Starname'][i])

    plotting.writeStarList(solution.plotly, startObs, args.schedule_dates[n], outputdir=savepath)
    plotting.plot_path_2D(solution, outputdir=savepath)
    plotting.nightPlan(solution.plotly, args.schedule_dates[n], outputdir=savepath)

    obs_and_times = pd.read_csv(savepath + 'ObserveOrder_' + str(args.schedule_dates[0]) + ".txt")
    all_targets_frame = pd.read_csv(request_sheet)
    gapFillers = np.loadtxt(savepath + 'gapFillerTargets.txt', delimiter=',', dtype=str)
    pf.write_starlist(all_targets_frame, obs_and_times, solution.extras, gapFillers, 'nominal', str(args.schedule_dates[0]), savepath)
