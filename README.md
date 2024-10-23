# KPF-CC Autoscheduler
This software represents Version 2.0 of the KPF-CC autoscheduler. The algorithm is now more granular, working on the timescale of "slots" and trying to schedule targets into specific slots within the semester. This new framework allows additional features as well.

This package contains the code for solving two related problems:
1. Optimal Semester Scheduler -- Determine which targets should be observed on which nights
2. Optimal Instrument Scheduler -- Determine the best distribution of nights for the KPF-CC queue in a semester

# Installation instructions

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

# Run instructions

To test installation and practice the workflow of running the autoscheduler, use the supplied example. From your terminal, navigate to the /bin directory and run the following command:
```
python generateScript.py -d 2024-08-02 -f ../example/
```

# # Flags
-d specifies the date for which to produce a script, format YYYY-MM-DD (No Default)
-f specifies the folder containing the inputs and outputs directory (can be defaulted to an environment variable: key name "KPFCC_SAVE_PATH")
-s specifies the slot size, in minutes (Default: 10)
-r specifies to run the "bonus" round (Default: False)
-t specifies the maximum time to solve the model (Default: 600s)
-g specifies to print the gorubi output to terminal (Default: True)
-p specifies to plot the results (Default: True)

# Solving the Traveling Telescope Problem (TTP)
Once the Optimal Semester Schedule is solved, we then wish to take the list of stars that have been chosen to be observed on a given night and further determine the optimal order within the night to observe these stars, effectively minimizing the slew time between them. To do this, we must employ the Traveling Telescope Problem (TTP) software.

For installation info and documentation on the TTP, see its Github Repo at: https://github.com/lukehandley/ttp. Be sure to update your environment variables to include the path to the TTP. The TTP must be installed in order for optimalAllocation to run fully and properly.

# More Info
More information on the KPF-CC program and the algorithm can be found in these places:
- The original autoscheduler paper: https://ui.adsabs.harvard.edu/abs/2024AJ....167..122H/abstract
- The original TTP paper (please cite this if you obtain observations through KPF-CC):  https://ui.adsabs.harvard.edu/abs/2024AJ....167...33H/abstract
- The updated autoscheduler paper (describes this new algorithm, please cite this if you obtain observations through KPF-CC):
- A FAQ page: https://docs.google.com/document/d/1-TfR6lNEtHO3muw_2Sc7l9Klveb3yovMUu28Jeq21t0/edit?usp=sharing
