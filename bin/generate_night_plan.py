import sys
import argparse
import os

import numpy as np
import pandas as pd

# The -3 cuts of the "bin/" of the path to this current file
path2modules = os.path.dirname(os.path.abspath(__file__))[:-3]
sys.path.append(path2modules)
import kpfcc.admin_functions as af
import kpfcc.solve_semester as ss
import kpfcc.helper_functions as hf
import kpfcc.plotting_functions as ptf
import kpfcc.processing_functions as pf
import kpfcc.ttp_functions as ttpf
import kpfcc.star_tracker as st
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
parser.add_argument('-r','--run_round_two', help='Run the bonus round', action='store_true')
parser.add_argument('-mB','--max_bonus', help='For round 2, the maximum bonus observations as a percentage of total request', type=int, default=0.5)
parser.add_argument('-ttp','--run_ttp', help='Turn off the TTP.', action='store_false')
parser.add_argument('-w','--run_weather_loss', help='If False, do not simulate weather', action='store_false')

# Optimal Allocation Flags - optional
parser.add_argument('-opt','--run_as_optimal_allocation', help='If True, run in optimal allocation mode', action='store_true')
parser.add_argument('-aes','--run_with_aesthetic', help='If True, run optimal allocation with aesthetic constraints', action='store_false')
parser.add_argument('-mQ','--max_quarters', help='The maximum quarters that can be allocated in optimal allocation', type=int, default=1)
parser.add_argument('-mN','--max_nights', help='The maximum unique nights on sky that can be allocated in optimal allocation', type=int, default=1)
parser.add_argument('-aS','--allow_single_quarters', help='If True, in optimal allocation we can allocate single quarter nights.', action='store_false')
parser.add_argument('-mxC','--max_consecutive', help='In optimal allocation set max consecutive on sky nights allowed.', type=int, default=6)
parser.add_argument('-mnC','--min_consecutive', help='In optimal allocation set min consecutive off sky nights allowed.', type=int, default=10)
parser.add_argument('-mxB','--max_baseline', help='In optimal allocation set minimum days from start/end.', type=int, default=5)

# Optional parameters - Use Default 99% of the time
parser.add_argument('-smg','--solve_max_gap', help='Max gap between solution and ideal (%)',type=int, default=0.05)
parser.add_argument('-stl','--solve_time_limit', help='Max time spent optimizing (sec)',type=int, default=300)
parser.add_argument('-s','--slot_size', help='The slot size (minutes)', type=int, default=5)
parser.add_argument('-g','--show_gurobi', help='Turn on Gurobi console print', action='store_false')
parser.add_argument('-b','--run_backups', help='Turn on plot outputs', action='store_true')
parser.add_argument('-stmp','--build_starmaps', help='Turn on plot outputs', action='store_true')
args = parser.parse_args()

manager = af.data_admin(
                        # basic parameters
                        args.folder,
                        str(args.schedule_dates[0]),
                        args.slot_size,
                        # files for running scheduler
                        # important to name files appropriately
                        os.path.join(args.folder, "inputs/requests.csv"),
                        os.path.join(args.folder, "inputs/twilight_times.csv"),
                        os.path.join(args.folder, "inputs/queryJumpDatabase.csv"),
                        os.path.join(args.folder, "inputs/allocation_schedule.txt"),
                        os.path.join(args.folder, "inputs/AccessMaps_" + str(args.slot_size) + "minSlots.txt"),
                        os.path.join(args.folder, "inputs/specialMaps_" + str(args.slot_size) + "minSlots.txt"),
                        os.path.join(args.folder, "inputs/zero_out.csv"),
                        os.path.join(args.folder, "inputs/NonQueueMap"  + str(args.slot_size) + ".txt"),
                        os.path.join(args.folder, "inputs/NonQueue.csv"),
                        args.run_weather_loss,
                        # parameters for running optimal allocation
                        args.run_as_optimal_allocation,
                        args.run_with_aesthetic,
                        args.max_quarters,
                        args.max_nights,
                        args.allow_single_quarters,
                        args.max_consecutive,
                        args.min_consecutive,
                        args.max_baseline,
                        os.path.join(args.folder, "inputs/whiteout_dates.txt"),
                        os.path.join(args.folder, "inputs/blackout_dates.txt"),
                        # flags for running gurobi
                        args.show_gurobi,
                        args.run_plots,
                        args.solve_time_limit,
                        args.solve_max_gap,
                        args.run_round_two,
                        args.max_bonus,
                        # files for plotting
                        os.path.join(args.folder, "/data/first_forecasts/"),
                        os.path.join(args.folder, "/outputs/" + str(args.schedule_dates[0]) + "/cadences/"),
                        os.path.join(args.folder, "inputs/turnOnOffDates.csv"),
                        os.path.join(args.folder, "inputs/cadenceTemplateFile.csv"),
                        os.path.join(args.folder, "outputs/" + str(args.schedule_dates[0]) + "/raw_combined_semester_schedule_Round2.txt"),
                        args.build_starmaps,
                        # files for ttp
                        os.path.join(args.folder, "inputs/NightlyStartStopTimes.csv"),
                        os.path.join(DATADIR,"bright_backups_frame.csv"),
                        os.path.join(DATADIR,"bright_backup_observability.csv")
                        )

manager.run_admin()

if args.run_scheduler:
    ss.run_kpfcc(manager)

if args.run_plots:
    star_tracker = st.StarTracker(manager)
    ptf.write_cadence_plot_files(manager)
    ptf.run_plot_suite(star_tracker, manager)

if args.run_ttp:
    ttpf.run_ttp(manager)
