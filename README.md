# KPF-CC Autoscheduler
This software represents Version 2.0 of the KPF-CC autoscheduler. The algorithm is now more granular, working on the timescale of "slots" and trying to schedule targets into specific slots within the semester. This new framework allows additional features as well.

This package contains the code for solving two related problems:
1. Optimal Semester Scheduler -- Determine which targets should be observed on which nights
2. Optimal Instrument Scheduler -- Determine the best distribution of nights for the KPF-CC queue in a semester

# Installation instructions

The installation of this package is not straight-forward. Please read carefully and follow closely.

We _*highly recommend*_ (though optional) setting up within a new Conda environment first. To do so, run this command and then be sure to activate the environment:
```
conda create -n kpfcc python=3.9
conda activate kpfcc
```

### Installing Gurobi
Before doing anything else, you must obtain a license for Gurobi.

`kpfcc` relies on Gurobi for solving large matrix equations efficiently. Follow these steps to install and set up Gurobi:

1. **Create an an Account** on Gurobi's [registration site](https://portal.gurobi.com/iam/register/). Select that you are an "Academic", type in your home institution, and submit the form via "Access Now". You will receive an email to complete the registration.
2. **Download Gurobi** for your OS from [this download page](https://www.gurobi.com/downloads/gurobi-software/). Run the installer file following their instructions.
3. **Request an Academic License** from your [user portal](https://portal.gurobi.com/iam/licenses/request/) *while connected to a university network*. You want the 'Named-User Academic License' which has a one year lifetime. At the end of the year, you can obtain a new license easily within your account (and for free) so long as you have maintained your academic status.
4. **Retrieve the License** by running the command from the popup window in a shell. It should look like:
```
grbgetkey 253e22f3...
```

### Install optimalAllocation

Only once you have obtained a Gurobi license can clone this repository:
```
git clone https://github.com/jluby127/optimalAllocation.git
```

And then install via:
```
pip install .
```

This will set up the environment and all dependency packages. Once again, it is very important that you _**do not**_ run the pip installer until you have obtained a Gurobi license. I cannot stress this enough.

### Install ttpsolver
One of the optimalAllocation dependencies is a package we co-developed, the TTP Solver (Traveling Telescope Problem). While the optimalAllocation autoscheduler decides which requests should be observed on a given night, the TTP solves for the optimal slew path to observe all those targets within the night. When pip installing the optimalAllocation package, the ttpsolver will be cloned from its repo as well. You must navigate to its local path and similarly pip install the package. For more information on installation and documentation of the TTP, see its Github Repo at: https://github.com/lukehandley/ttp. Since you have obtained a Gurobi license as part of the installation instructions above, you may skip that section of the TTP's installation instructions. Be sure to update your environment variables to include the path to the TTP.

# Run instructions

To test installation and practice the workflow of running the autoscheduler, use the supplied example. From your terminal, run the following command:
```
python <LOCAL_DIR>/bin/generate_night_plan.py -d 2024-08-02 -f <LOCAL_DIR>/examples/
```

## Flags

### Required

-f specifies the folder containing the inputs and outputs directory (can be defaulted to an environment variable: key name "KPFCC_SAVE_PATH")

-d specifies the date for which to produce a script, format YYYY-MM-DD (No Default)

### Advanced

-a specifies to run the autoscheduler (Default: True)

-p specifies to run the plotting/reporting suite (Default: True)

-ttp specifies to turn off the ttp solver (Default: False)

-r specifies to run the "bonus" round (Default: False)

-w specifies to turn off the weather loss simulations. Now no nights will be lost to weather. (Default: False)

### Optional

-s specifies the slot size, in minutes (Default: 5)

-b specifies to turn on the back up bright star list script. Now we will produce an additional night plan of only bright (V < 8) stars. (Default: False)

-t specifies the maximum time to solve the model (Default: 300s)

-g specifies to print the gorubi output to terminal (Default: True)


# More Info
More information on the KPF-CC program and the algorithm can be found in these places:
- The original autoscheduler paper: https://ui.adsabs.harvard.edu/abs/2024AJ....167..122H/abstract
- The original TTP paper (please cite this if you obtain observations through KPF-CC):  https://ui.adsabs.harvard.edu/abs/2024AJ....167...33H/abstract
- The updated autoscheduler paper (describes this new algorithm, please cite this if you obtain observations through KPF-CC):
- A FAQ page: https://docs.google.com/document/d/1-TfR6lNEtHO3muw_2Sc7l9Klveb3yovMUu28Jeq21t0/edit?usp=sharing
