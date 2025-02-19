import sys
import argparse
import os

import numpy as np
import pandas as pd

# the -3 cuts of the "bin/" of the path to this current file
path2modules = os.path.dirname(os.path.abspath(__file__))[:-3]
sys.path.append(path2modules  + "/kpfcc/")
import solve_semester as ss
import helper_functions as hf
import plotting_functions as ptf
import processing_functions as pf
import ttp_functions as ttp

parser = argparse.ArgumentParser(description='Generate schedules with KPF-CC v2')
parser.add_argument('-d','--schedule_dates',action='append',help='Date(s) to be scheduled; \
                            string in format YYYY-MM-DD. Use a -d flag between each date.')
parser.add_argument('-f','--folder', help='Folder to save all outputs',
                            default=os.environ["KPFCC_SAVE_PATH"])

parser.add_argument('-a','--run_scheduler', help='Turn off the autoscheduler', action='store_false')
parser.add_argument('-p','--run_plots', help='Turn off the plotting', action='store_false')
parser.add_argument('-ttp','--run_ttp', help='Turn off the TTP.', action='store_false')
parser.add_argument('-e','--exclude_history', help='Do not parse past observing history', default=True)

parser.add_argument('-r','--run_extra_rounds', help='Run the bonus round', action='store_true')
parser.add_argument('-t','--timeout', help='Max time spent optimizing (sec)',type=int, default=300)
parser.add_argument('-s','--slot_size', help='The slot size (minutes)', type=int, default=10)
parser.add_argument('-w','--run_weather_loss', help='If False, do not simulate weather', action='store_false')
parser.add_argument('-g','--show_gurobi', help='Turn on Gurobi console print', action='store_false')
parser.add_argument('-b','--run_backups', help='Turn on plot outputs', action='store_true')
args = parser.parse_args()

# files for semester solver
request_sheet = args.folder + "inputs/Requests.csv"
allocated_nights = args.folder + "inputs/2024B_Binary_Schedule.txt"
past_database = args.folder + "inputs/queryJumpDatabase.csv"
twilight_times = args.folder + "inputs/2024B_twilight_times.csv"
access_map = args.folder + "inputs/2024B_AccessMaps_" + str(args.slot_size) + "minSlots.txt"
special_map = args.folder + 'inputs/2024B_specialMaps_' + str(args.slot_size) + 'minSlots.txt'
nonqueue_map =  args.folder + 'inputs/2024B_NonQueueMap'  + str(args.slot_size) + '.txt'

# files for plotting
folder_forecasts = args.folder + '/data/first_forecasts/'
folder_cadences = args.folder + '/outputs/' + str(args.schedule_dates[0]) + '/cadences/'
zero_out_file = 'nofilename.txt'
turn_on_off_file = args.folder + "inputs/2024B_turnOnOffDates.csv"
starmap_template_filename = args.folder + "inputs/2024B_cadenceTemplateFile.csv"
nonqueue_list = args.folder + "inputs/NonQueue.csv"
future_forecast = args.folder + "outputs/" + str(args.schedule_dates[0]) + \
                      '/raw_combined_semester_schedule_Round2.txt'

# files for ttp
backup_file = path2modules + "data/bright_backups_frame.csv"
backup_observability_file = path2modules + "data/bright_backup_observability.csv"
nightly_start_stop_times = args.folder + "inputs/2024B_NightlyStartStopTimes.csv"

if args.run_scheduler:
    ss.run_kpfcc(args.schedule_dates,
                              request_sheet,
                              allocated_nights,
                              access_map,
                              twilight_times,
                              args.folder + 'outputs/' + str(args.schedule_dates[0]) + '/',
                              args.slot_size,
                              args.run_extra_rounds,
                              past_database,
                              starmap_template_filename,
                              turn_on_off_file,
                              nonqueue_map,
                              special_map,
                              zero_out_file,
                              args.run_weather_loss,
                              args.exclude_history,
                              args.show_gurobi,
                              args.run_plots,
                              args.timeout)

if args.run_plots:
    data = ptf.DataHandler(args.schedule_dates[0], request_sheet, past_database, future_forecast,
                           twilight_times, nonqueue_list, folder_forecasts, folder_cadences,
                           args.folder + 'reports/', args.slot_size)
    ptf.run_plot_suite(data, args.folder, build_starmaps=True)

if args.run_ttp:
    ttp.run_ttp(future_forecast, nightly_start_stop_times, request_sheet,
                                    args.schedule_dates[0], args.folder)
