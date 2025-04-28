.. _getting_started:

Getting Started
===============

.._installation:


Installation
++++++++++++


Conda Environment
------------------

We recommend installing ``AstroQ`` within a new 
`Anaconda <https://www.anaconda.com/download>`_ environment. 
To do so, clone the `AstroQ Github repository <https://github.com/jluby127/AstroQ>`_ 
and run the following commands in the top-level directory:

    .. code-block:: bash

        $ conda create --name astroq --file environment.yml
        $ conda activate astroq
    
Installing Gurobi
------------------

``AstroQ`` relies on the Gurobi optimization software to efficiently solve large matrix equations.
Follow these steps to install and set up Gurobi:

    1. **Create an Account** on Gurobi's `registration site <https://portal.gurobi.com/iam/register/>`_. Select that you are an "Academic," type in your home institution, and submit the form via "Access Now." You will receive an email to complete the registration.
    
    2. **Download Gurobi** for your OS from `this download page <https://www.gurobi.com/downloads/gurobi-software/>`_. Follow the instructions to run the installer file.