# AstroQ
This software represents Version 2.0 of the KPF-CC autoscheduler, known as AstroQ.

This package contains the code for solving two related problems:
1. Optimal Semester Scheduler -- Determine which targets should be observed on which nights
2. Optimal Instrument Scheduler -- Determine the best distribution of nights for the KPF-CC queue in a semester

# Documentation

For information on the installation and use of the ``AstroQ`` code, please visit the [documentation page](https://astroq.readthedocs.io).

# Notes on HIRES implementation

These are notes for the HIRES implementation starting in 2026A. You'll need the following extra environment variables to be set

CC_CONDA_ENV = <name of conda environment> # needed to test and develop different versions of the cod
HIRES_PROGRAM_SHEET_URLS_CSV = <csv file containing urls to google sheet of hires requests>
CC_OUTPUT_PATH = usually the same directory as where the repo is stored.

# More Info
More information on the KPF-CC program and the algorithm can be found in these places:
- The original autoscheduler paper: https://ui.adsabs.harvard.edu/abs/2024AJ....167..122H/abstract
- The original TTP paper (please cite this if you obtain observations through KPF-CC):  https://ui.adsabs.harvard.edu/abs/2024AJ....167...33H/abstract
- The updated autoscheduler paper (describes this new algorithm, please cite this if you obtain observations through KPF-CC):
- A FAQ page: https://docs.google.com/document/d/1-TfR6lNEtHO3muw_2Sc7l9Klveb3yovMUu28Jeq21t0/edit?usp=sharing
