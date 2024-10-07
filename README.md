# KPF-CC Autoscheduler
This software represents Version 2.0 of the KPF-CC autoscheduler. The algorithm is now more granular, working on the timescale of "slots" and trying to schedule targets into specific slots within the semester. This new framework allows additional features as well.

This package contains the code for solving two related problems:
1. Optimal Semester Scheduler -- Determine which targets should be observed on which nights
2. Optimal Instrument Scheduler -- Determine the best distribution of nights for the KPF-CC queue in a semester

# Installation instructions

<!-- Use `conda` or `mamba` to install a dedicated envirnoment for KPF-CC

```
mamba env create -n kpf-cc -f environment.yml
``` -->

We highly recommend (though optional) setting everything up in a new Conda environment first. To do so, run this command and then be sure to activate the environment:
```
conda create -n kpfauto python=3.9
conda activate kpfauto
```

To install `optimalAllocation`, clone this repository:
```
$ git clone https://github.com/jluby127/optimalAllocation.git
```
Next cd in to optimalAllocation directory and you can install all relevant packages at once using pip:
```
pip install -r requirements.txt
```
This will set up the environment with the _**exception**_ of Gurobipy. Very important _**do not**_ pip install Gurobipy at this step. Before installing Gurobipy, we must first install Gurobi from their website and acquire a license. See below.


### Installing Gurobi
`ttpsolver` relies on Gurobi for solving large matrix equations efficiently. Follow these steps to install and set up Gurobi:

1. **Create an an Account** on Gurobi's [registration site](https://portal.gurobi.com/iam/register/). Select that you are an "Academic", type in your home institution, and submit the form via "Access Now". You will receive an email to complete the registration.
2. **Download Gurobi** for your OS from [this download page](https://www.gurobi.com/downloads/gurobi-software/). Run the installer file following their instructions.
3. **Request an Academic License** from your [user portal](https://portal.gurobi.com/iam/licenses/request/) *while connected to a university network*. You want the 'Named-User Academic License' which has a one year lifetime. At the end of the year, you can obtain a new license easily within your account (and for free) so long as you have maintained your academic status.
4. **Retrieve the License** by running the command from the popup window in a shell. It should look like:
```
grbgetkey 253e22f3...
```
5. **Install Gurobipy** using either pip or conda *after* completing the previous steps. The reason we must wait to run this line until after obtaining a license is that, without a specific license, Gurobipy will give you a trial license. This trial license limits the size of the models you can solve and expires quickly. We need a full, true license in order to solve this model, so you must obtain a 'Name-User Academic License' (see above) before installing Gurobipy.
```
pip install gurobipy
```
Or:
```
conda install -c gurobi gurobi
```



# Run instructions:
As there are two main problems this package solves, each has its own run file. Together they share a few flags, and then each also has its own individual flag options.



### Required Input Files
This package expects certain files in certain formats to be operational. Below is a list. The files should be saved in a directory together and then the path and filenames can be hard-coded into the generateScript.py and runOptimalAllocation.py files. These files contain information about the target stars, the requests, the telescope allocation schedule, etc. These files will need to change every semester and therefore the hard-coded paths should be updated each semester.



### Both solvers
Request Sheet - a csv containing the information on each target's coordinates, stellar parameters, and observational strategy. Specific column names are required.

Instrument Schedule - A txt file representing a python dictionary where each row is formatted as "key:value". The keys are calendar dates in format YYYY-MM-DD and the values are a string of format "X X X X" where X is either a 1 or a 0 to indicate if that quarter (1st to 4th running left to right) is allocated or not (respectively) to the queue.

Twilight Times - a csv containing pre-computed 12 degree twilight times for both evening and following morning for each calendar day of the semester. Specific column names are required.

Accessibility Maps - a pickle file of a python dictionary where keys are the target names and values are the 1D array of the precomputed accessibility map for that target.

Time Sensitive Observations - a csv of the time sensitive observations, their target names and start/stop times in UTC format

Turn On/Off Dates - a csv where each targets first and last calendar day that it is observable is recorded. For plotting purposes later. Specific column names required.

Starmap Template - a csv template of the file which will be used to produce a cadence plot of a specific target. For plotting purposes later. Specific column names required.

### Optimal Semester Schedule
Past Database - a csv copy of the database which stores all previous true observations. Specific column names are required.

### Optimal Instrument Schedule
Enforced Yes - a csv containing the calendar dates and quarters that the solver must select to be allocated to the queue. Specific column names are required.

Enforced No - a csv containing the calendar dates and quarters that the solver is forbidden from selecting to be allocated to the queue. Specific column names are required.

### Shared flags
-d -- specifies the date to start the solver from. Note that any previous days in the semester will not be scheduled but all days up until the end of the semester will be solved for. Note that for Optimal Instrument Schedule, the date specified should be the first day of the semester in question. There is no default and not specifying a date will throw an error that quits the program.

-f -- name of the folder you wish to save all the output files and figures (include path here too). It is helpful to hard-code the default, then this flag is not needed each time.

-l -- the maximum time, in seconds, to allow the solver to find the optimal solution. The default is 300s. For the Optimal Semester Schedule this is generally enough time. For Optimal Instrument Schedule you should consider upping this to at least 900s.

-s -- the slot size, in minutes, to use when solving the schedule. The default is 5 minutes. Depending on your computer's specs, 5 minute slots may not be feasible to solve. See below for guidelines on computational costs.

-g -- Turn off the printing out of Gurobi's solver information. Default is on.

-p -- Turn off printing out of plots and reports. Default is on.

### Optimal Semester Schedule only flags
To solve the Optimal Semester Schedule from the supplied date (in given format) onwards to the end of the semester and produce a script for tonight's observations, run: python generateScript.py -d YYYY-MM-DD

-r -- Turn on solving Round 2 of the semester solver. Default is off.

### Optimal Instrument Schedule only flags
To solve the Optimal Instrument Schedule, run: python runOptimalAllocation.py -d YYYY-MM-DD

# Computational Costs
Solving these models is a heavy computational problem and some standard out of the box computers may not be able to solve them. The actual runtime and computational power needed varies highly on the size of the model you are trying to solve, which in turn has many input variables including: the number of nights in a semester, the size of the slots in minutes, and the number of individual requests you are attempting to schedule.

Given a standard Keck Observatory semester of 182 nights, with 5 minute slots, and an average request load ~100, it is recommended that your machine have at least 64GB of RAM and 10+ cores. For reference, the KPF-CC program uses an Apple Mac Studio with 60 CPU cores and 128 GB of RAM and an M2 chip; it takes XX seconds to solve the model. The computational power required is significantly lessened for 10 minute slot sizes, although this produces a less realistic and more inefficient model as very short exposure time targets must occupy an entire 10 minutes of real time. Similarly, it is recommended that requests to observe a star only once in the semester be excluded from the model and instead be scheduled semi-manually. These requests are highly flexible, creating a vast symmetry to the problem space, which makes the model more difficult to solve. Plus they add an entire row to the solution matrix which increases runtime.

We have designed the algorithm to only solve the semester from the given input date onward to the end of the semester. This means that as the semester progresses, the size of the model is continually shrinking as there are fewer days in the semester to solve for.

# Solving the Traveling Telescope Problem (TTP)
Once the Optimal Semester Schedule is solved, we then wish to take the list of stars that have been chosen to be observed on a given night and further determine the optimal order within the night to observe these stars, effectively minimizing the slew time between them. To do this, we employ the Traveling Telescope Problem (TTP) software.

For more installation info and documentation on the TTP, see its Github Repo at: https://github.com/lukehandley/ttp. Be sure to update your environment variables to include the path to the TTP. The TTP must be installed in order for optimalAllocation to run fully and properly.

# More Info
More information on the KPF-CC program and the algorithm can be found in these places:
- The original autoscheduler paper: https://ui.adsabs.harvard.edu/abs/2024AJ....167..122H/abstract
- The original TTP paper (please cite this if you obtain observations through KPF-CC):  https://ui.adsabs.harvard.edu/abs/2024AJ....167...33H/abstract
- The updated autoscheduler paper (describes this new algorithm, please cite this if you obtain observations through KPF-CC):
- A FAQ page: https://docs.google.com/document/d/1-TfR6lNEtHO3muw_2Sc7l9Klveb3yovMUu28Jeq21t0/edit?usp=sharing
