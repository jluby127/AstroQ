import sys
import argparse
import os

import numpy as np
import pandas as pd

# The -3 cuts of the "bin/" of the path to this current file
path2modules = os.path.dirname(os.path.abspath(__file__))[:-3]
sys.path.append(path2modules)
import kpfcc.management as mn
import kpfcc.schedule as sc
import kpfcc.plot as pl
import kpfcc.io as io
# import kpfcc.onsky as sk
import kpfcc.tracking as tk
from kpfcc import DATADIR

parser = argparse.ArgumentParser(description='Generate schedules with KPF-CC')
parser.add_argument('-d','--today', help="Today's date, in format YYYY-MM-DD.")
parser.add_argument('-c','--config', help="Path and filename of to config file.")
args = parser.parse_args()

manager = mn.data_admin(args.config, args.today)
manager.run_admin()

if manager.run_scheduler:
    sc.run_kpfcc(manager)

if manager.run_plots:
    forecast_frame = pd.read_csv(manager.future_forecast)
    manager.combined_semester_schedule_stars = forecast_frame.values
    star_tracker = tk.StarTracker(manager)
#    io.report_allocation_stats(manager)
    pl.write_cadence_plot_files(manager)
    pl.run_plot_suite(star_tracker, manager)

if manager.run_ttp:
    sk.run_ttp(manager)
