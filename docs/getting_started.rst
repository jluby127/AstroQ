.. _getting_started:

Getting Started
===============

.. _installation:

Installation
++++++++++++

    
Installing Gurobi
------------------

``AstroQ`` relies on the Gurobi optimization software to efficiently solve large matrix equations.
Follow these steps to install and set up Gurobi:

    1. **Create an Account** on Gurobi's `registration site <https://portal.gurobi.com/iam/register/>`_. Select that you are an "Academic," type in your home institution, and submit the form via "Access Now." You will receive an email to complete the registration.
    
    2. **Download Gurobi** for your OS from `this download page <https://www.gurobi.com/downloads/gurobi-software/>`_. Follow the instructions to run the installer file.
    
    3. **Request an Academic License** from your `user portal <https://portal.gurobi.com/iam/login/?target=https%3A%2F%2Fportal.gurobi.com%2Fiam%2Flicenses%2Frequest%2F>`_ *while connected to a university network*. You want the 'Named-User Academic License,' which has a one-year lifetime. At the end of the year, you can obtain a new license easily within your account (and for free) as long as you have maintained your academic status.
    
    4. **Retrieve the License** by running the command given in the popup window in a shell. It should look like:
        .. code-block:: bash

            $ grbgetkey 253e22f3...
        
    This command will download a file called ``gurobi.lic`` to your machine.
    
    5. **Define the GRB_LICENSE_FILE environment variable** by adding the following to your ``.zshrc`` file:
        .. code-block:: bash

            export GRB_LICENSE_FILE=/path/to/gurobi.lic
            
    6. **Optionally define some additional environment variables** while you're at it:
    
        .. code-block:: bash

            export KPFCC_JUMP_USERNAME=your_username
            export KPFCC_JUMP_PASSWORD=your_pass
            export KPFCC_SAVE_PATH=path/to/AstrQ/results


Installing AstroQ
------------------

We recommend installing ``AstroQ`` within a new 
`Anaconda <https://www.anaconda.com/download>`_ environment. 
To do so, clone the `AstroQ Github repository <https://github.com/jluby127/AstroQ>`_ 
and run the following commands in the top-level directory:

    .. code-block:: bash

        $ conda env create --name astroq --file environment.yml
        $ conda activate astroq

Then install the package:
        
    .. code-block:: bash

        $ pip install -e .
        
        
Installation Test
------------------
        
Run ``pytest`` from your terminal to execute a suite of tests verifying that the installation was successful:

    .. code-block:: bash

        $ pytest -v

This may take a few minutes on the first run. The desired result is that all tests pass successfully.
                

Example Usage
+++++++++++++

**AstroQ requires four files to run. Additional files are optional.**

* **config.ini** - Contains configuration information for the AstroQ run. See examples in the ``examples/`` directory.

* **allocation.csv** - Contains information about the nights and times when the telescope is available for observation. It must contain appropriate column header
    - "start"
    - "end"
    Times are in format "YYYY-MM-DD HH:MM"

* **request.csv** - Contains information about the targets and their observational strategies. It must contain appropriate column headers:
    - "unique_id" - A unique identifier for the target.
    - "starname" - The human-readable name of the target.
    - "program_code" - The program code of the target. 
    - "ra" - The right ascension of the target.
    - "dec" - The declination of the target.
    - "exptime" - The exposure time of the target.
    - "n_exp" - The number of exposures of the target.
    - "n_intra_max" - The desired number of visits within a night.
    - "n_intra_min" - The acceptable number of visits within a night.
    - "tau_intra" - The minimum time (hours)between visits of the target within a night.
    - "n_inter_max" - The maximum number of unique nights to observe the target.
    - "tau_inter" - The minimum time (days) between unique night observations of the target.
    More columns may be added as needed for your queue. 

* **past.csv** - Contains information about the past history of observations. It must contain appropriate column headers:
    - "id" - the unique identifier for the target.
    - "target" - the human-readable name of the target.
    - "semid" - the program code of the target.
    - "exposure_start_time" - the start time of the recorded exposure. 
    - "exposure_time" - the recorded duration of the exposure.

* **OPTIONAL: custom.csv** - Contains information about the specific time windows when targets may be observed. It must contain appropriate column headers:
    - "unique_id" - the unique identifier for the target.
    - "starname" - the human-readable name of the target.
    - "start" - the start time of the time window.
    - "stop" - the stop time of the time window.
    Times are in format "YYYY-MM-DD HH:MM"

* **OPTIONAL: programs.csv** - Contains information the awarded time to each program. Useful for plotting later. It must contain appropriate column headers:
    - "program" - the program code.
    - "hours" - the hours of awarded time to the program.
    - "nights" - the number of nights awarded to the program.

**The kpfcc_prep subcommand compiles and prepares all necessary files for the KPFCC program specifically.** Use of this command requires token access to the KPFCC database. If you are using AstroQ for a different facility, then you must prepare these files through your own means. Feel free to write your own "prep" command to automate the process.

It is strongly recommended that you create the following folder structure for your AstroQ project:
    <upstream path>/<semester>/<date>/

    where upstream path is the path to the top-level directory of your AstroQ project.
    semester is the semester ID of your project.
    date is the date of your project.

    For example, if your upstream path is ``/Desktop/``, your semester is ``2025B``, your date is ``2025-08-01``, and your band is ``1``, then your folder structure should be:
    /Desktop/2025B/2025-08-01/

    In this folder structure, you should place (or automate the creation of) the necessary files above to run AstroQ. We recommend a new folder for every run, so as to create a history of schedules. You may further create additional weather bands within the date level. Each band should get its own copy of all files.

Test your installation further by running some of the commands below. We will use ``AstroQ``'s command line interface to mock up a toy observing schedule. After performing the steps in the :ref:`installation` section, you should have access to the ``astroq`` binary at the command line. If not, ensure that you have run ``pip install -e .`` in the top-level directory of the repository.

Let's see which subcommands are available using ``astroq --help``:

    .. code-block:: bash
    
    
        $ usage: astroq [-h] [-V] {bench,plot,kpfcc_prep,webapp,plan-semester,plan-night,compare} ...

            AstroQ: Optimized observation scheduling

            optional arguments:
            -h, --help            show this help message and exit
            -V, --version         Print version number and exit.

            subcommands:
            {bench,plot,kpfcc_prep,webapp,plan-semester,plan-night,compare}
                

The AstroQ command-line interface provides the following subcommands:

* **bench** - Runs the benchmark test. Here you can compare your computer's performance with the AstroQ benchmark from our paper, Lubin et al. 2025.
* **kpfcc_prep** - Compile and prepare all necessary files for the KPFCC program. Note that if you are building a queue for a different facility, you will have to build your own "prep" command.
* **webapp** - Launch web app to view interactive plots for a given solution of AstroQ.
* **plan-semester** - Solve for the optimal semester-long schedule, determining what stars to observe on what nights. This is the heart of AstroQ.
* **plan-night** - Solve for the optimal slew path using the TTP package.
* **compare** - Compare two AstroQ solutions.

**Below are more detailed explanations of these commands and their outputs.**

To create your mock observing schedule, let's run the hello world example. Here we have supplied the necessary files. We will create a schedule for the first night when we have the telescope, 2024-08-02, as specified in the config file. Then we run:

    .. code-block:: bash
    
        $ astroq plan-semester -cf examples/hello_world/config_hello_world.ini
        
Let's take a look at the outputs produced:

    .. code-block:: bash
    
        $ ls -ltr examples/hello_world/outputs
        
            -rw-r--r--@ 1 jack  staff    3153 Oct 29 12:15 semester_plan.csv
            -rw-r--r--@ 1 jack  staff  230542 Oct 29 12:15 serialized_outputs_dense_v1.csv
            -rw-r--r--@ 1 jack  staff  253640 Oct 29 12:15 serialized_outputs_dense_v2.csv
            -rw-r--r--@ 1 jack  staff     346 Oct 29 12:15 runReport.txt
            -rw-r--r--@ 1 jack  staff     873 Oct 29 12:15 request_selected.csv
            -rw-r--r--@ 1 jack  staff  244512 Oct 29 12:15 semester_planner.h5
      
- ``semester_plan.csv`` contains the day/slot/name tuples of the scheduled observations. Grouped by target name. See example:

    .. code-block:: csv
    
        r,d,s,name
        e1,4,62,HIP1532
        e1,5,66,HIP1532
        e1,8,65,HIP1532
        e1,12,61,HIP1532
        e1,13,82,HIP1532
        e1,61,77,HIP1532
        e1,63,76,HIP1532
        e1,64,77,HIP1532
        e1,86,67,HIP1532
        e1,87,67,HIP1532
        e2,4,36,TOI-1670
        e2,5,35,TOI-1670
        e2,8,41,TOI-1670

- ``serialized_outputs_dense_v1.csv`` contains the same information, but now all slots, even those not scheduled to have an observation are included, and it is ordered by time. See example (the first 14 slots of the semester are not scheduled, as they are during day time):

    .. code-block:: csv

        d,s,r,name
        0,0,,
        0,1,,
        0,2,,
        0,3,,
        0,4,,
        0,5,,
        0,6,,
        0,7,,
        0,8,,
        0,9,,
        0,10,,
        0,11,,
        0,12,,
        0,13,,
        0,14,,
- ``serialized_outputs_dense_v2.png`` is identical to ``serialized_outputs_dense_v1.csv``, but now slots that cannot be filled (due to day time or allocation, are denoted with an "X") See example:

    .. code-block:: csv
        d,s,r,name
        0,0,X,
        0,1,X,
        0,2,X,
        0,3,X,
        0,4,X,
        0,5,X,
        0,6,X,
        0,7,X,
        0,8,X,
        0,9,X,
        0,10,X,
        0,11,X,
        0,12,X,
        
    
- `runReport.txt`: contains some basic statistics about the fullness of the schedule. See example:
    .. image:: plots/runreport.png

- `request_selected.csv`: contains a copy of the request.csv file but only for the targets that were selected to be observed tonight.:
        
- `semester_planner.h5`: is a serialized and compressed version of the splan object use to schedule. This will be used later in the webapp/plotting routines.

Now that we have the stars to be observed each night of the semester, let's determine the optimal ordering of the stars selected for tonight, mimizing slew times, to generate a night plan:

    .. code-block:: bash
    
        $ astroq plan-night -cf examples/hello_world/config_hello_world.ini 
        

Here are the new files in ``examples/hello_world/outputs/``:

    .. code-block:: bash
        
        -rw-r--r--@ 1 jack  staff     447 Oct 29 13:36 ttp_prepared.csv
        -rw-r--r--@ 1 jack  staff     433 Oct 29 13:36 TTPstatistics.txt
        -rw-r--r--@ 1 jack  staff     158 Oct 29 13:36 ObserveOrder_2018-08-05.txt
        -rw-r--r--@ 1 jack  staff    1157 Oct 29 13:36 script_2018-08-05_nominal.txt
        -rw-r--r--@ 1 jack  staff   16080 Oct 29 13:36 night_planner.h5

- ``ttp_prepared.csv`` is the input to the TTP solver. It contains the target information in the format required by the TTP solver. See example and the TTP repository for more info:
    .. code-block:: txt

        Starname,RA,Dec,Exposure Time,Exposures Per Visit,Visits In Night,Intra_Night_Cadence,Priority,First Available,Last Available
        e1,4.77317661843625,-9.964852409,360,2,1,0,10,2018-08-05 10:20,2018-08-05 15:00
        e2,259.0173367,72.16115935,1200,1,1,0,10,2018-08-05 05:50,2018-08-05 08:50
        e7,285.679422455377,50.2413060048164,50,1,4,1,10,2018-08-05 05:50,2018-08-05 12:20
        e9,348.320729001503,57.1683566176719,60,3,1,0,10,2018-08-05 08:40,2018-08-05 15:00

- ``TTPstatistics.txt`` contains some basic statistics about the TTP solution. See example:

    .. code-block:: txt
        Stats for TTP Solution
        ------------------------------------
            Model ran for 0.05 seconds
            Observations Requested: 7
            Observations Scheduled: 7
        Maximum Observations Bound: 1
        ------------------------------------
        Observing Duration (min): 560.00
        Time Spent Exposing (min): 40.57
            Time Spent Idle (min): 510.51
        Time Spent Slewing (min): 8.92
        Minimum Slew Bound (min): 4.17
        ------------------------------------

- ``ObserveOrder_2018-08-05.txt`` the timestamps at which each star is set to be observed. See example:

    .. code-block:: txt
        unique_id,Target,StartExposure
        e2,TOI-1670,05:47
        e7,Kepler-10,06:07
        e7,Kepler-10,07:07
        e7,Kepler-10,08:47
        e9,219134,08:49
        e1,HIP1532,11:41
        e7,Kepler-10,11:59

- ``script_2018-08-05_nominal.txt`` is a backwards compatible file in the style of the old HIRES queue night plans.

- ``night_planner.ht``: a serialized and compressed version of the night planner object use to schedule. This will be used later in the webapp/plotting routines.

Now we can launch the webapp to view the schedule. This will launch on a local server at http://localhost:50001. When it launches, append "/2018B/2018-08-05/band1/admin" to the URL to view the admin page. Check out and familiarize yourself with the various figures and reports.

    .. code-block:: bash
    
        $ astroq webapp -up examples/hello_world/




        
        
        
