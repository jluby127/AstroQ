import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta
import argparse
import os

parser = argparse.ArgumentParser(description='Generate schedules with KPF-CC v2')
parser.add_argument('-d','--schedule_dates',action='append',help='Date(s) to be scheduled as strings in format YYYY-MM-DD. Must be in allocated_nights')
parser.add_argument('-f','--folder', help='Folder to save generated scripts and plots', default=os.environ["KPFCC_SAVE_PATH"])
parser.add_argument('-r','--run_extra_rounds',action='store_true', help='Turn on the additional rounds of scheduling', default=False)
parser.add_argument('-t','--time_limit', help='Max time spent optimizing (s)',type=int, default=300)
parser.add_argument('-s','--slot_size', help='The slot size (min)',type=int, default=10)
parser.add_argument('-g','--gurobi_output',action='store_true', help='Activate Gurobi console outputs',default=True)
parser.add_argument('-p','--plot_results',action='store_true', help='Turn on plotting for semester/upcoming quarter nights', default=True)
# hard code these in just for now
parser.add_argument('-b1','--backupStarList',action='store_true', help="File containing the backup star names/RAs/Decs/exptimes", default='/Users/jack/Documents/KPF_CC/StandardStarInfo/testBackupScriptGenerationFiles/backupStarlist.csv')
parser.add_argument('-b2','--backupStarObservability',action='store_true', help="File containing the backup star accessibilties", default='/Users/jack/Documents/KPF_CC/StandardStarInfo/testBackupScriptGenerationFiles/backupObservability_justJan.csv')
# parser.add_argument('-b1','--backupStarList',action='store_true', help="File containing the backup star names/RAs/Decs/exptimes", default='')
# parser.add_argument('-b2','--backupStarObservability',action='store_true', help="File containing the backup star accessibilties", default='')

args = parser.parse_args()

request_sheet = args.folder + "inputs/Requests.csv"
allocated_nights = args.folder + "inputs/2024B_Binary_Schedule.txt"
pastDatabase = 'nofilename.txt'#args.folder + "inputs/queryJumpDatabase.csv"
twilight_times = args.folder + "inputs/2024B_twilight_times.csv"
access_map = args.folder + "inputs/2024B_AccessMaps_" + str(args.slot_size) + "minSlots.txt"
turnFile = args.folder + "inputs/2024B_turnOnOffDates.csv"
starmap_template_filename = args.folder + "inputs/2024B_cadenceTemplateFile.csv"
nonqueueMap =  args.folder + 'inputs/2024B_NonQueueMap'  + str(args.slot_size) + '.txt'
specialMaps = args.folder + 'inputs/2024B_specialMaps_' + str(args.slot_size) + 'minSlots.txt'
zeroOutFile = 'nofilename.txt'
startstoptimes = pd.read_csv(args.folder + "inputs/2024B_NightlyStartStopTimes.csv")
backupStarList = args.backupStarList
backupStarObservability = args.backupStarObservability

import sys
path2modules = os.path.dirname(os.path.abspath(__file__))[:-3] #the -3 cuts of the "bin/" of the path to this current file
sys.path.append(path2modules + "kpfcc/")
import solveSemester as ssm
import processingFunctions as pf

ssm.runKPFCCv2(args.schedule_dates,
                          request_sheet,
                          allocated_nights,
                          access_map,
                          twilight_times,
                          args.folder + 'outputs/' + str(args.schedule_dates[0]) + '/',
                          args.slot_size,
                          args.run_extra_rounds,
                          pastDatabase,
                          starmap_template_filename,
                          turnFile,
                          nonqueueMap,
                          specialMaps,
                          zeroOutFile,
                          args.gurobi_output,
                          args.plot_results,
                          args.time_limit)
#
# import helperFunctions as hf
# sys.path.append(os.environ["TTP_PATH"] + "ttp/")
# import formatting
# import telescope
# import plotting
# import model
# savepath = args.folder + 'outputs/' + str(args.schedule_dates[0]) + '/'
# print("Prepare schedule for the TTP.")
# tel = telescope.Keck1()
# semester_start_date, semester_end_date, semesterLength, semesterYear, semesterLetter = hf.getSemesterInfo(args.schedule_dates[0])
# all_dates_dict = hf.buildDayDateDictionary(semester_start_date, semesterLength)
# the_schedule = np.loadtxt(savepath + 'raw_combined_semester_schedule_Round2.txt', delimiter=',', dtype=str)
# for n in range(len(args.schedule_dates)):
#
#     dayInSemester = all_dates_dict[args.schedule_dates[n]]
#     # starttime = startstoptimes['Start'][idx]
#     # stoptime = startstoptimes['Stop'][idx]
#     starttime = '05:00'
#     stoptime = '10:00'
#     startObs = Time(str(args.schedule_dates[n]) + "T" + str(starttime), format='isot')
#     endObs = Time(str(args.schedule_dates[n]) + "T" + str(stoptime), format='isot')
#     total_time = np.round((endObs.jd-startObs.jd)*24,3)
#     print("Time in Night for Observations: " + str(total_time) + " hours.")
#
#     idx = startstoptimes.index[startstoptimes['Date']==str(args.schedule_dates[n])][0]
#
#     Round2_Requests = np.loadtxt(savepath + 'Round2_Requests.txt', dtype=str)
#     toTTP = pf.prepareTTP(request_sheet, the_schedule[dayInSemester], Round2_Requests)
#
#     filename = savepath + '/Selected_' + str(args.schedule_dates[n]) + ".txt"
#     if os.path.exists(args.backupStarList) and os.path.exists(args.backupStarObservability):
#         print("Generating filler bright stars and backup weather script")
#         check = os.path.isdir(savepath + 'Backups/')
#         if not check:
#             os.makedirs(savepath + 'Backups/')
#         import backupStarFunctions as bsf
#         backup_starlist = pd.read_csv(args.backupStarList)
#         backUpObservability = pd.read_csv(args.backupStarObservability)
#         backups = bsf.getStarsForTonight(backup_starlist, backUpObservability, args.schedule_dates[0][5:], minimumUpTime=4.0)
#         backups.drop(columns='index', inplace=True)
#         subsetBackups = bsf.getTimesWorth(backups, 1.0) # get targets for 1 hours' worth of exposures
#
#         toTTP_all = pd.concat([toTTP, subsetBackups], ignore_index=True)
#         toTTP_all.to_csv(filename, index=False)
#         filename_backups = savepath + '/Backups_' + str(args.schedule_dates[n]) + ".txt"
#         backups_full = bsf.getTimesWorth(backups, total_time+1) # get targets for 1 hours' worth of exposures
#         backups_full.to_csv(filename_backups, index=False)
#         # print("median: ", np.median(backups['Exposure Time']))
#         # print("mean: ", np.mean(backups['Exposure Time']))
#         # print("max: ", max(backups['Exposure Time']))
#         # print("min: ", min(backups['Exposure Time']))
#         # totalopenshutterwithreadout = []
#         # for w in range(len(backups)):
#         #     if backup_starlist['ExpTime'][w] <= 150.0 and backup_starlist['ExpTime'][w] > 72.0:
#         #         nshot = 2
#         #     elif backup_starlist['ExpTime'][w] <= 72.0 and backup_starlist['ExpTime'][w] > 45.0:
#         #         nshot = 3
#         #     elif backup_starlist['ExpTime'][w] <= 45.0:
#         #         nshot = 5
#         #     else:
#         #         nshot = 1
#         #     totalexptime = backup_starlist['ExpTime'][w]*nshot + 45*(nshot-1)
#         #     totalopenshutterwithreadout.append(totalexptime)
#         # print("--------------------------------------")
#         # print("Now with nshots and read time accounted")
#         # print("--------------------------------------")
#         # print("median: ", np.median(totalopenshutterwithreadout))
#         # print("mean: ", np.mean(totalopenshutterwithreadout))
#         # print("max: ", max(totalopenshutterwithreadout))
#         # print("min: ", min(totalopenshutterwithreadout))
#     else:
#         toTTP.to_csv(filename, index=False)
#
#     targlist = formatting.theTTP(filename)
#     solution = model.TTPModel(startObs, endObs, targlist, tel, savepath)
#
#     plotting.writeStarList(solution.plotly, startObs, args.schedule_dates[n], outputdir=savepath)
#     plotting.plot_path_2D(solution,outputdir=savepath)
#     plotting.nightPlan(solution.plotly, args.schedule_dates[n], outputdir=savepath)
#
#     obs_and_times = pd.read_csv(savepath + 'ObserveOrder_' + str(args.schedule_dates[0]) + ".txt")
#     all_targets_frame = pd.read_csv(request_sheet)
#     # pf.write_starlist(all_targets_frame, obs_and_times, solution.extras, Round2_Requests, 'nominal', str(args.schedule_dates[0]), savepath)
#
#     if os.path.exists(args.backupStarList) and os.path.exists(args.backupStarObservability):
#         print('Generating independent backup weather script.')
#         targlist_backups = formatting.theTTP(filename_backups)
#
#         solution_b = model.TTPModel(startObs, endObs, targlist_backups, tel, savepath + "Backups/", runtime=600, optgap=0.05)
#
#         plotting.writeStarList(solution_b.plotly, startObs, args.schedule_dates[n], outputdir=savepath + "Backups/")
#         plotting.plot_path_2D(solution_b,outputdir=savepath+ "Backups/")
#         plotting.nightPlan(solution_b.plotly, args.schedule_dates[n], outputdir=savepath + "Backups/")
#
#         obs_and_times_b = pd.read_csv(savepath + 'Backups/ObserveOrder_' + str(args.schedule_dates[0]) + ".txt")
#         # pf.write_starlist(backup_starlist, obs_and_times_b, solution_b.extras, [], 'Backups', str(args.schedule_dates[0]), savepath + "Backups/")
#
# print("Semester is scheduled and script is generated. Complete.")
