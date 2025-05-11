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

        $ pytest

This may take a few minutes on the first run. The desired result is that all tests pass successfully.
        
Temporary note: you can also run ``coverage`` to see how much code the tests touch:

    .. code-block:: bash

        $ coverage run -m pytest
        

Example Usage
+++++++++++++

Test your installation further by running some of the commands below. We will use ``AstroQ``'s command line interface to mock up a toy observing schedule. After performing the steps in the :ref:`installation` section, you should have access to the ``astroq`` binary at the command line. If not, ensure that you have run ``pip install -e .`` in the top-level directory of the repository.

Let's see which subcommands are available using ``astroq --help``:

    .. code-block:: bash

        $ astroq --help
          usage: astroq [-h] [-V] {bench,plot,schedule,kpfcc} ...

          AstroQ: Optimized observation scheduling

          optional arguments:
              -h, --help            show this help message and exit
              -V, --version         Print version number and exit.

          subcommands:
              {bench,plot,schedule,kpfcc}


To create your mock observing schedule, first run ``kpfcc build``:

    .. code-block:: bash
    
        $ astroq kpfcc build -cf examples/hello_world/config_hello_world.ini
        
Let's take a look at the outputs produced:

    .. code-block:: bash
    
        $ ls -ltr examples/hello_world/outputs/2024-08-01
        
        -rw-r--r--@ 1 staff       0 May 11 08:54 runReport.txt
        -rw-r--r--  1 staff   35685 May 11 08:54 weather_loss_visualization.png
        -rw-r--r--@ 1 staff    4904 May 11 08:54 Weather_Simulation_Results.csv
        -rw-r--r--  1 staff  729964 May 11 08:54 request_set.json
      
- ``request_set.json`` contains the PI-provided target observation criteria used to generate the observing schedule
- ``Weather_Simulation_Results.csv`` contains the results of a simulation to estimate the number of nights lost due to poor weather.
- ``weather_loss_visualization.png`` provides a visual display of the information in ``Weather_Simulation_Results.csv``. Here is an example of what it might look like:

    .. image:: plots/weather_loss_visualization.png
- `runReport.txt`: <DESCRIPTION HERE>
        
Next, use ``request_set.json`` to generate an observing schedule:

    .. code-block:: bash
        $ astroq schedule -cf examples/hello_world/config_hello_world.ini -rf examples/hello_world/outputs/2024-08-01/request_set.json
        

Here are the new files in ``examples/hello_world/outputs/2024-08-01``:

    .. code-block:: bash
        
        -rw-r--r--  1  staff   32652 May 11 09:23 raw_combined_semester_schedule_Round1.txt
        -rw-r--r--  1  staff   31447 May 11 09:23 raw_combined_semester_schedule_available.txt
        -rw-r--r--  1  staff   32652 May 11 09:23 raw_combined_semester_schedule_Round2.txt
        -rw-r--r--  1  staff       0 May 11 09:23 Round2_Requests.txt
        -rw-r--r--  1  staff  206459 May 11 09:23 serialized_outputs_dense.csv

- ``raw_combined_semester_schedule_Round1.txt``, ``raw_combined_semester_schedule_Round2.txt``, and ``raw_combined_semester_schedule_available.txt`` give visual representations of the scheduling output. In these files, each line represents one night of the semester and each column represents a single *10-minute* slot. Slots populated with any combination of an ``X``, ``*``, or ``W`` are unavailable for observation, indicating respectively that the telescope is allocated to another group, the sky is not sufficiently dark, or that the weather simulation has predicted that slot will be weathered out. Slots allocated to a specific target are populated with that target's name. Here are a few example rows from ``raw_combined_semester_schedule_Round2.txt``:

    .. image:: plots/schedule_Round2_screenshot.png

    - Each night begins and ends with an interval of un-allocated bright time, corresponding to twilight.
    - ``Line 70`` shows a first-half-only night, during which three observations of `Kepler-10` are scheduled.
    - No time is allocated to our program on the night corresponding to ``Line 71``.
    - The night of ``Line 79`` is entirely weathered out.
    
- ``Round2_Requests.txt``: <EXPLANATION HERE>
- ``serialized_outputs_dense.csv`` presents the information in the ``raw_combined`` files (which one specifically?) in a machine-readable format. There is one row for every slot of every night in the semester, and six columns:
    - ``d``: The day of the semester, beginning at 0
    - ``s``: The slot of a given night, beginning at 0
    - ``r``: The target requested during a given slot (empty if none)
    - ``isNight``: Whether it is night (1 if yes, 0 if no)
    - ``isAlloc``: Whether the telescope is allocated to your program (1 if yes, 0 if no)
    - ``isClear``: Whether the weather is conducive to observing (1 if yes, 0 if no)
    
    .. image:: plots/outputs_dense_screenshot.png
    
    - The image above shows 11 slots during the 21st night of the semester. `HAT-P-11` is scheduled during slots 53-XX (update with dense version that shows all filled slots), followed by `219134` in slots 62-XX. The sky is dark, the telescope is allocated to our program, and the weather is good during all of these slots.      
        
        
        
