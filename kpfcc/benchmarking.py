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

def build_toy_model_from_paper(hours_per_program = 100, plot = False, shortcut=0):
    # order: cadences, exptime, nobs, visits
    program0 = [1, 300, 40, 1] # APF-50
    program1 = [5, 600, 20, 1] # TKS
    program2 = [15, 1200, 10, 1] # bi-weekly
    program3 = [1, 300, 8, 5] # Intra
    program4 = [1, 300, 40, 1] # Constrained
    program5 = [1, 3600, 1, 1] # Singles
    program6 = [1, 300, 1, 1] # Singles to serve as extra; one slot as placeholder, later will edit.
    all_programs = [program0, program1, program2, program3, program4, program5, program6]

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

        # First six programs get a number of stars to fill their allocation and observing strategy
        # Program 6 is the extra requests
        if p <= 5:
            n_stars = stars_in_program(all_programs[p], hours_per_program)
        else:
            n_stars = 1350
        prog_nstars.append(n_stars)
        prog_nexpos.append(n_stars*viz*nights)
        prog_nslots.append(n_stars*viz*nights*int(exptime/300))
        prog_award.append(np.round((n_stars*viz*nights*exptime)/3600,1))
        all_programs[p].append(n_stars)
        total_stars += n_stars

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

    # Randomly generate the RA/Decs for these stars
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
    timecounter = 0
    slotcounter = 0
    starcounter = 0
    for p in range(len(all_programs)):
        for s in range(all_programs[p][4]):
            starname.append("Star" + f"{starcounter:04d}")
            starcounter += 1
            program_code.append("Program" + str(p))
            # Program 4 is the constrained to roughly the Kepler field
            if p==4:
                tmpra = np.random.uniform(18*15, 20*15)
                tmpdec = np.random.uniform(40, 50)
            else:
                tmpra = np.random.uniform(18*15, ((20*15)+180))
                tmpdec = getDec()
            if tmpra > 360:
                tmpra -= 360
            RA.append(tmpra)
            Dec.append(tmpdec)
            exposure_times.append(all_programs[p][1])
            exposures_per_visit.append(1)
            visits_per_night_max.append(all_programs[p][3])
            if all_programs[p][3] == 5:
                visits_per_night_min.append(3)
            else:
                visits_per_night_min.append(1)
            unique_nights.append(all_programs[p][2])
            internight_cadences.append(all_programs[p][0])
            if all_programs[p][3] == 1:
                intranight_cadences.append(0)
            else:
                intranight_cadences.append(1)
            timecounter += ((all_programs[p][1]*all_programs[p][3]*all_programs[p][2])/3600)
            slotcounter += ((all_programs[p][1]*all_programs[p][3]*all_programs[p][2])/300)

    toy_requests = pd.DataFrame({'starname':starname, 'program':program_code, 'ra':RA, 'dec':Dec,
                                'exptime':exposure_times,
                                'n_exp':exposures_per_visit,
                                'n_intra_max':visits_per_night_max,
                                'n_intra_min':visits_per_night_min,
                                'n_inter_max':unique_nights,
                                'tau_inter':internight_cadences,
                                'tau_intra':intranight_cadences
                                })

    if plot:
        # Plot out the stars on the sky
        for i in range(len(all_programs), 0, -1):
            filt = toy_requests[toy_requests['program'] == 'Program' + str(i)]
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

    if shortcut > 0:
        toy_requests = toy_requests[:shortcut]
    print("The toy model is defined! Happy benchmarking.")
    return toy_requests