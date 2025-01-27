import sys
import argparse
import os

import numpy as np
import pandas as pd

# The -3 cuts of the "bin/" of the path to this current file
#path2modules = os.path.dirname(os.path.abspath(__file__))[:-3]
# sys.path.append(path2modules  + "/kpfcc/")
import kpfcc.solve_semester as ss
import kpfcc.helper_functions as hf
import kpfcc.plotting_functions as ptf
import kpfcc.processing_functions as pf
import kpfcc.ttp_functions as ttpf
from kpfcc import DATADIR

parser = argparse.ArgumentParser(description='Generate schedules with KPF-CC v2')
# Required parameters
parser.add_argument('-d','--schedule_dates',action='append',help='Date(s) to be scheduled; \
                            string in format YYYY-MM-DD. Use a -d flag between each date.')
parser.add_argument('-f','--folder', help='Folder to save all outputs',
                            default=os.environ["KPFCC_SAVE_PATH"])
# Optional parameters - High Level
parser.add_argument('-a','--run_scheduler', help='Turn off the autoscheduler', action='store_false')
parser.add_argument('-p','--run_plots', help='Turn off the plotting', action='store_false')
parser.add_argument('-r','--run_extra_rounds', help='Run the bonus round', action='store_true')
parser.add_argument('-ttp','--run_ttp', help='Turn off the TTP.', action='store_false')
parser.add_argument('-w','--run_weather_loss', help='If False, do not simulate weather', action='store_false')
# Optional parameters - Use Default 99% of the time
parser.add_argument('-t','--timeout', help='Max time spent optimizing (sec)',type=int, default=300)
parser.add_argument('-s','--slot_size', help='The slot size (minutes)', type=int, default=5)
parser.add_argument('-g','--show_gurobi', help='Turn on Gurobi console print', action='store_false')
parser.add_argument('-b','--run_backups', help='Turn on plot outputs', action='store_true')
args = parser.parse_args()

# Files for semester solver
request_sheet = os.path.join(args.folder, "inputs/Requests.csv")
allocated_nights = os.path.join(args.folder, "inputs/2024B_Binary_Schedule.txt")
past_database = os.path.join(args.folder, "inputs/queryJumpDatabase.csv")
twilight_times = os.path.join(args.folder, "inputs/2024B_twilight_times.csv")
access_map = os.path.join(args.folder, "inputs/2024B_AccessMaps_" + str(args.slot_size) + "minSlots.txt")
special_map = os.path.join(args.folder, "inputs/2024B_specialMaps_" + str(args.slot_size) + "minSlots.txt")
nonqueue_map =  os.path.join(args.folder, "inputs/2024B_NonQueueMap"  + str(args.slot_size) + ".txt")
zero_out_file = os.path.join(args.folder, "inputs/zero_out.csv")

# Files for plotting
folder_forecasts = os.path.join(args.folder, "/data/first_forecasts/")
folder_cadences = os.path.join(args.folder, "/outputs/" + str(args.schedule_dates[0]) + "/cadences/")
turn_on_off_file = os.path.join(args.folder, "inputs/2024B_turnOnOffDates.csv")
starmap_template_filename = os.path.join(args.folder, "inputs/2024B_cadenceTemplateFile.csv")
nonqueue_list = os.path.join(args.folder, "inputs/NonQueue.csv")
future_forecast = os.path.join(args.folder, "outputs/" + str(args.schedule_dates[0]) + \
                      "/raw_combined_semester_schedule_Round2.txt")

# Files for the TTP
nightly_start_stop_times = os.path.join(args.folder, "inputs/2024B_NightlyStartStopTimes.csv")
backup_file = os.path.join(DATADIR,"bright_backups_frame.csv")
backup_observability_file = os.path.join(DATADIR,"bright_backup_observability.csv")

if args.run_scheduler:
    ss.run_kpfcc(args.schedule_dates,
                              request_sheet,
                              allocated_nights,
                              access_map,
                              twilight_times,
                              args.folder + "outputs/" + str(args.schedule_dates[0]) + "/",
                              args.slot_size,
                              args.run_extra_rounds,
                              past_database,
                              starmap_template_filename,
                              turn_on_off_file,
                              nonqueue_map,
                              special_map,
                              zero_out_file,
                              args.run_weather_loss,
                              args.show_gurobi,
                              args.run_plots,
                              args.timeout)

if args.run_plots:
    data = ptf.DataHandler(args.schedule_dates[0], request_sheet, past_database, future_forecast,
                           twilight_times, nonqueue_list, folder_forecasts, folder_cadences,
                           args.folder + "reports/", args.slot_size)
    ptf.run_plot_suite(data, args.folder)

if args.run_ttp:
    ttpf.run_ttp(future_forecast, nightly_start_stop_times, request_sheet,
                                    args.schedule_dates[0], args.folder)
