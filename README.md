## optimalAllocation
Software for solving two related problems:
1. Determine which targets should be observed on which nights (semester-solver)
2. Determine the best distribution of nights for the KPF-CC queue in a semester, given target request info

This software represents Version 2 of the KPF-CC algorithm, now with more granular scheduling of targets into individual slots within the semester. This new framework allows additional features as well.

#Installation instructions
A number of packages are required to run this software. They can all be installed by running:
pip install -r requirements.txt
while within the main directory. It is recommended that you set up a new conda environment for running this software.

# Run instructions:

python runKPFCCv2.py -d YYYY-MM-DD

This will run the semester-solver from the supplied date (in given format) onwards to the end of the semester. Additionally, you can specify parameters via a series of flags:

-d -- specifies the date to start the semester-schedule solver. When running Optimal Allocation, only use the first day of the semester here. If you use another date, the code automatically resets the date to the first day of the semester. Note that any previous days in the semester will not be scheduled. There is no default and not specifying a date will throw an error that quits the program.

-f -- name of the folder you wish to save all the output files and figures (include path here too). The default can be hard-coded.

-l -- the maximum time, in seconds, to allow the solver to find the optimal solution. The default is 300s.

-s -- the slot size, in minutes, to use when solving the schedule. The default is 10 min.

-b -- the path and filename to the EnforcedNO list of nights/quarters. These are nights/quarters that cannot be chosen to be allocated. This is only relevant for when solving Optimal Allocation. The default is none.

-w -- the path and filename to the EnforcedYES list of nights/quarters. These are nights/quarters that must be chosen to be allocated. This is only relevant for when solving Optimal Allocation. The default is none.

-g -- Turn off the printing out of Gurobi's solver information. Default is on.

-o -- Turn on solving Optimal Allocation. Default is off.

-r -- Turn on solving Round 2 of the semester solver. Default is off.

-p -- Turn off printing out of plots and reports. Default is on.
