"""Standalone TTP runner (upstream bin/runTTP.py).

Preserved as a standalone script alongside the vendored astroq.ttp
package so that one-off optimal target-ordering runs (no semester or
night plan required) remain possible. Note: this script tracks the
upstream API; some plotting helper signatures may have drifted relative
to AstroQ's main pipeline.

Usage:
    python -m astroq.ttp.bin.runTTP -f path/to/inputs.txt
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta
import os
import imageio

from astroq.ttp import formatting, telescope, plotting, model

# You must use the -f flag and point to the inputs file
import argparse
parser = argparse.ArgumentParser(description='Generate optimal slew paths with the TTP.')
parser.add_argument('-f','--folder', help='Path to the inputs file', default=None)
args = parser.parse_args()

print("Prepare schedule for the TTP.")
inputs = formatting.readInputs(args.folder)
if not os.path.isdir(inputs['SavePath']):
    os.mkdir(inputs['SavePath'])

current_day = str(inputs['StartNight'][:10])
tel = telescope.create_tel(inputs['Telescope'])
startObs = Time(inputs['StartNight'], format='isot')
endObs = Time(inputs['EndNight'], format='isot')
total_time = np.round((endObs.jd-startObs.jd)*24,3)
print("Time in Night for Observations: " + str(total_time) + " hours.")

targlist = formatting.theTTP(inputs['StarListFile'])
print("Solving optimal slew path.")
solution = model.TTPModel(startObs, endObs, targlist, tel, inputs['SavePath'])
solutionDict = solution.schedule

print("Solution found.")
orderedList = []
for i in range(len(solutionDict['Starname'])):
    orderedList.append(solutionDict['Starname'][i])

print("Writing star list.")
plotting.writeStarList(solution.plotly, startObs, current_day, outputdir=inputs['SavePath'])
print("Plotting the 2D path.")
plotting.plot_path_2D(solution, outputdir=inputs['SavePath'])
print("Generating the Night Plan plots")
plotting.nightPlan(solution.plotly, current_day, outputdir=inputs['SavePath'])
print("Producing animated slewpath.")
plotting.animate_telescope(solution, startObs, endObs, outputdir=inputs['SavePath'], animationStep=120)
print()
print()
print("Done. Clear skies!")
