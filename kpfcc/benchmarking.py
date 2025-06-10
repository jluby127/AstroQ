"""
Module for preparing the benchmark tests

"""
import json
import os
import sys
import numpy as np
import matplotlib.pyplot as pt
import pandas as pd
import math
import random
from configparser import ConfigParser
from argparse import Namespace

import kpfcc.request as rq
import kpfcc.driver as dr

# In the paper, we used random seed = 24.
np.random.seed(24)

def getDec(maxDec=90, minDec=-20):
    '''
    Randomly draw a declination from cosine i distribution between two values.
    The default min/max declination values are chosen based on favorable viewing from Hawaii.
    '''
    mincosdec = np.cos((90+maxDec)*(np.pi/180.))
    maxcosdec = np.cos((90+minDec)*(np.pi/180.))
    cosdec = np.random.uniform(mincosdec, maxcosdec)
    dec = (np.arccos(cosdec)*(180./np.pi))-90
    return dec

def stars_in_program(program, total_hours):
    '''
    Determine how many stars should be in a program based on that program's observing strategy and
    it's awarded time.
    '''
    e = program[1]
    n = program[2]
    v = program[3]
    total_time = (e*n*v)/3600
    n_stars = int(np.round(total_hours/total_time,0))
    return n_stars

def build_toy_model_from_paper(ns, hours_per_program = 80, plot = False):
    """
    Generate a synthetic request set based on the paper's toy model.
    Returns a pandas DataFrame with the request information.
    """
    # Define programs with their characteristics
    # [tau_inter, t_exp [seconds], n_inter_max, n_intra_max]
    program1 = [1, 300, 40, 1]  # APF-50
    program2 = [3, 600, 20, 1]  # TKS
    program3 = [10, 1200, 10, 1]  # bi-weekly
    program4 = [1, 300, 8, 5]  # Intra
    program5 = [1, 300, 40, 1]  # Constrained
    dynamic_exptime = ns*300
    program6 = [0, dynamic_exptime, 1, 1]  # Singles
    all_programs = [program1, program2, program3, program4, program5, program6]

    # Calculate number of stars per program
    stars_per_program = []
    for p in range(len(all_programs)):
        n_stars = stars_in_program(all_programs[p], hours_per_program)
        stars_per_program.append(n_stars)
        all_programs[p].append(n_stars)

   # Compute stats for the paper's table
    stars_per_program = []
    total_stars = 0
    prog_numb = []
    prog_nobs = []
    prog_inter = []
    prog_visits = []
    prog_intra = []
    prog_nexp = []
    prog_exptime = []
    prog_nstars = []
    prog_nexpos = []
    prog_nslots = []
    prog_award = []
    for p in range(len(all_programs)):
        nights = all_programs[p][2]
        inter = all_programs[p][0]
        viz = all_programs[p][3]
        exptime = all_programs[p][1]

        prog_numb.append(p)
        prog_nobs.append(nights)
        prog_inter.append(inter)
        prog_visits.append(viz)
        if viz == 1:
            prog_intra.append(0)
        else:
            prog_intra.append(1)
        prog_nexp.append(1)
        prog_exptime.append(exptime)

        n_stars = stars_in_program(all_programs[p], hours_per_program)
        prog_nstars.append(n_stars)
        prog_nexpos.append(n_stars*viz*nights)
        prog_nslots.append(n_stars*viz*nights*int(exptime/300))
        prog_award.append(np.round((n_stars*viz*nights*exptime)/3600,1))
        all_programs[p].append(n_stars)
        total_stars += n_stars

        # Metadata about toy model, currently not returned
    prog_info = pd.DataFrame({"Program #":prog_numb,
                            "# Nights":prog_nobs,
                            "Inter Cadence":prog_inter,
                            "# Visits":prog_visits,
                            "Intra Cadence":prog_intra,
                            "# Exposures":prog_nexp,
                            "Exp Time":prog_exptime,
                            "# Stars":prog_nstars,
                            "Total Exposures":prog_nexpos,
                            "Total Slots":prog_nslots,
                            "Award":prog_award,
                            })
    # Generate request data
    starname = []
    program_code = []
    RA = []
    Dec = []
    exposure_times = []
    internight_cadences = []
    intranight_cadences = []
    visits_per_night_max = []
    visits_per_night_min = []
    unique_nights = []
    exposures_per_visit = []
    n_intra_max = []
    tau_inter = []
    priority = []

    star_idx = 0
    for p, program in enumerate(all_programs):
        n_stars = program[4]  # Number of stars for this program
        for s in range(n_stars):
            starname.append(f"Star{star_idx:04d}")
            program_code.append(f"Program{p+1}")

            # Generate RA/Dec following original logic
            if p == 4:  # Constrained to Kepler field
                tmpra = np.random.uniform(19*15, 19.66*15)  # RA between 285-295 degrees
                tmpdec = np.random.uniform(40, 50)          # Dec between 40-50 degrees
            else:
                # only allow RAs that are somewhat favorable to a B semester, exclude hour anges between 12 and 18
                exclude_start = 12*15
                exclude_end = 18*15
                tmpra = np.random.uniform(0, exclude_start) if np.random.random() < (exclude_start / (360 - (exclude_end - exclude_start))) else np.random.uniform(exclude_end, 360)
                tmpdec = getDec()  # Uses cosine distribution
            RA.append(tmpra)
            Dec.append(tmpdec)

            exposure_times.append(program[1])
            internight_cadences.append(program[0])
            intranight_cadences.append(0 if program[3] == 1 else 1)
            visits_per_night_max.append(program[3])
            if program[3] == 5:
                visits_per_night_min.append(program[3]-2) # make min visits equal to 3 when max visits is 5
            else:
                visits_per_night_min.append(program[3]) # make min visits equal to max visists
            unique_nights.append(program[2])
            exposures_per_visit.append(1)
            n_intra_max.append(program[3])
            tau_inter.append(program[0])
            priority.append(np.random.randint(1, 5))
            star_idx += 1

    # Create DataFrame
    requests_data = {
        'starname': starname,
        'program_code': program_code,
        'ra': RA,
        'dec': Dec,
        'exptime': exposure_times,
        'n_exp': exposures_per_visit,
        'n_intra_max': visits_per_night_max,
        'n_intra_min': visits_per_night_min,
        'n_inter_max': unique_nights,
        'tau_inter': internight_cadences,
        'tau_intra': intranight_cadences
    }
    requests_data = pd.DataFrame(requests_data)

    if plot:
        for i in range(len(all_programs), 0, -1):
            filt = requests_data[requests_data['program'] == 'Program' + str(i)]
            if i < 5 and i != 4:
                c = 'b'
            elif i == 4:
                c = 'r'
            else:
                c = 'g'
            pt.plot(filt['ra'], filt['dec'], 'o', color=c)
        pt.xlim(0,360)
        pt.ylim(-40,90)
        pt.show()

    print("The toy model is defined! Happy benchmarking.")
    return requests_data
