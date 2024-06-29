# ttpsolver
Astronomical observation planning tool using the Traveling Telescope Problem (TTP) formulation. 

## Purpose
`ttpsolver` helps you plan observations for over 100 targets within a specified time frame. It simplifies the process by allowing you to input your observatory details, observing time constraints, and a target table. The tool plans observations to maximize successful exposures and optimizes the path between targets.

## Installation
To install `ttpsolver`, run:

*Update this line with package installation*

This will set up the environment with the _**exception**_ of Gurobipy.

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

## Algorithm

`ttpsolver` employs a Mixed-Integer Linear Programming formulation to optimize scheduling in a time-dependent manner, similar to the Traveling Salesman Problem. It determines both the order of targets and the timing of observations to respect celestial object rise and set times. The algorithm is detailed in the Astronomical Journal [here](https://iopscience.iop.org/article/10.3847/1538-3881/ad0dfb). Please cite this work when using ttpsolver for planning scientific observations.

## Tutorials
You can find a simple tutorial for using the package in the `/tutorial` directory.


