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
from configparser import ConfigParser
from argparse import Namespace

import kpfcc.request as rq
import kpfcc.driver as dr

# In the paper, we used random seed = 24.
np.random.seed(24)
def getDec(maxDec=75, minDec=-30):
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

def firstN_Requests(nstar, request_set, request_file):
    '''
    Cut down the large RequestSet's number of requests into a smaller number.
    Useful for benchmarking astroq's performance as a function of number of requests.
    Note: Setting nstar < 250 will remove requests that are part of the standard toy model. This
          set of requests is designed to mimic a real semester and is paired to match the total "awarded" time
          in the Real2024B and Random allocation maps.
    '''

    unique_stars = list(set(request_set.observability.id))
    unique_stars.sort()
    keep_stars = unique_stars[:nstar]

    mask1 = request_set.observability['id'].isin(keep_stars)
    mask2 = request_set.strategy['id'].isin(keep_stars)

    request_set.observability = request_set.observability[mask1]
    request_set.observability.reset_index(inplace=True, drop=True)
    request_set.strategy = request_set.strategy[mask2]
    request_set.strategy.reset_index(inplace=True, drop=True)

    request_frame = pd.read_csv(request_file)
    request_frame.sort_values(by='starname', inplace=True)
    request_frame.reset_index(inplace=True, drop=True)
    request_frame = request_frame[:nstar]
    request_frame.to_csv(request_file[:-8] + ".csv", index=False)

    return request_set

def set_nSlots_singles(nslot, request_set, start_row=250):
    '''
    Change the number of slots required to complete one request of the extra single shot requests.
    Useful for benchmarking astroq's performance as a function of # of Slots to Complete Request

    star_row is set to 250 because this is the first row of the extra single shot requests. We don't
    want to apply this logic to any of the standard toy model requests.
    '''
    request_set.strategy.loc[request_set.strategy.iloc[start_row:].index, 't_visit'] = nslot
    return request_set

def build_toy_model_from_paper(hours_per_program = 10, plot = False, savepath = "", shortcut=0):
    """
    Generate a synthetic request set based on the paper's toy model.
    Returns a pandas DataFrame with the request information.
    """
    # Define programs with their characteristics
    # order: cadences, exptime, nobs, visits
    program0 = [1, 300, 40, 1]  # APF-50
    program1 = [5, 600, 20, 1]  # TKS
    program2 = [15, 1200, 10, 1]  # bi-weekly
    program3 = [1, 300, 8, 5]  # Intra
    program4 = [1, 300, 40, 1]  # Constrained
    program5 = [1, 3600, 1, 1]  # Singles
    program6 = [1, 300, 1, 1]  # Singles to serve as extra
    all_programs = [program0, program1, program2, program3, program4, program5, program6]

    # Calculate number of stars per program
    stars_per_program = []
    for p in range(len(all_programs)):
        if p <= 5:
            n_stars = stars_in_program(all_programs[p], hours_per_program)
        else:
            n_stars = 1350
        stars_per_program.append(n_stars)
        all_programs[p].append(n_stars)

    #prog_info = pd.DataFrame(all_programs)
    #savepath = "./examples/bench/inputs/toy_model_program_info.csv"
    #test = pd.read_csv(savepath, index_col=0)
    #savepath = "./examples/bench/inputs/Requests.csv"
    #test2 = pd.read_csv(savepath, index_col=0)


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
            program_code.append(f"Program{p}")
            
            # Generate RA/Dec following original logic
            if p == 4:  # Constrained to Kepler field
                tmpra = np.random.uniform(18*15, 20*15)  # RA between 270-300 degrees
                tmpdec = np.random.uniform(40, 50)       # Dec between 40-50 degrees
            else:
                tmpra = np.random.uniform(18*15, ((20*15)+180))  # RA between 270-450 degrees
                tmpdec = getDec()  # Uses cosine distribution
            if tmpra > 360:
                tmpra -= 360
            RA.append(tmpra)
            Dec.append(tmpdec)
            
            exposure_times.append(program[1])
            internight_cadences.append(program[0])
            intranight_cadences.append(0 if program[3] == 1 else 1)
            visits_per_night_max.append(program[3])
            visits_per_night_min.append(program[3])
            unique_nights.append(program[2])
            exposures_per_visit.append(1)
            n_intra_max.append(program[3])
            tau_inter.append(program[0])
            priority.append(np.random.randint(1, 5))
            star_idx += 1

    # Create DataFrame
    requests_data = {
        'starname': starname,
        'program': program_code,
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
    return requests_data


