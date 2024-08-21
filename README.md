# ttpsolver
Astronomical observation planning tool that solves the Traveling Telescope Problem (TTP), a time dependent version of the Traveling Salesman Problem.

## Algorithm
`ttpsolver` helps you make an observing plan for your targets within a specified time frame. It simplifies the process by allowing you to input your observatory details, observing time constraints, and a target table. The tool plans observations to maximize the number of completed exposures and optimizes the slew path between targets.

`ttpsolver` employs a Mixed-Integer Linear Programming formulation to optimize scheduling in a time-dependent manner, similar to the Traveling Salesman Problem. It determines both the order of targets and the timing of observations to respect celestial object rise and set times. The algorithm is detailed in the Astronomical Journal [here](https://iopscience.iop.org/article/10.3847/1538-3881/ad0dfb). Please cite this work when using `ttpsolver` for planning scientific observations.

## Installation
For Conda users, start by setting up a new environment in python 3.9:
```
conda create -n kpfauto python=3.9
```

To install `ttpsolver`, clone this repository:
```
$ git clone https://github.com/lukehandley/ttp.git
```
Next you can install all relevant packages at once with pip using:
```
pip install -r requirements.txt
```
This will set up the environment with the _**exception**_ of Gurobipy. Very important _**do not**_ pip install Gurobipy at this step.


### Installing Gurobi
`ttpsolver` relies on Gurobi for solving large matrix equations efficiently. Follow these steps to install and set up Gurobi:

1. **Create an an Account** on Gurobi's [registration site](https://portal.gurobi.com/iam/register/). Select that you are an "Academic", type in your home institution, and submit the form via "Access Now". You will receive an email to complete the registration.
2. **Download Gurobi** for your OS from [this download page](https://www.gurobi.com/downloads/gurobi-software/).
3. **Request an Academic License** from your [user portal](https://portal.gurobi.com/iam/licenses/request/) *while connected to a university network*. You want the 'Named-User Academic License' which has a one year lifetime.
4. **Retrieve the License** by running the command from the popup window in a shell. It should look like:
```
grbgetkey 253e22f3...
```
5. **Install Gurobipy** using either pip or conda *after* completing the previous steps.
```
pip install gurobipy
```
Or:
```
conda install -c gurobi gurobi
```

## Tutorials
You can find a tutorial for using the package in the `/tutorial` directory. In the tutorial, we attempt to schedule 85 targets (the survey list for the TESS-Keck Survey, see [here](https://ui.adsabs.harvard.edu/abs/2022AJ....163..297C/abstract), as this is a large list of stars spread out over RA hours). Not all stars will be observable on the specified night but the `ttpsolver` handles this easily.

To run `ttpsolver`, we require an inputs file where the user specifies information about the logistics of the night (start and stop time to observations and the telescope to be used), as well as a file containing the table of targets and their observational strategy. Each target requires RA/Dec in hour angle/degrees, and an exposure time in seconds. Next each target's strategy is described by the number of exposures to take at each visit, the number of visits to make to the target and the minimum separation in time (hours) between visits to the same target. Lastly, each target requires a priority ranking from 1 to 10, with 10 being the highest priority.

Feel free to play with this tutorial in changing the exposure times, the observational strategies, the priorities, and/or the time window to compute the optimal path (both in terms of calendar date and hour window within the night).

`ttpsolver` returns a few reports and plots. First, it outputs a .txt file where the optimal slew path is laid out as each each target is listed next to the starting time of the exposure. Next, it outputs a .txt file with statistics about the solution. Then it generates 3 plots: a Night Plan (both as an interactive .html and as a static png), a 2D slew path representation, and an animation portraying the slew path over the course of the night.

We built `ttpsolver` with two telescopes pre-made: Keck1 and WIYN. However, we designed the telescope class to be easily replicated for any other custom observatory. See instructions within the telescope.py file.  
