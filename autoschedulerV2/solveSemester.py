import numpy as np
import matplotlib.pyplot as pt
import pandas as pd
import sys
import math
import time
import pickle
import copy
from collections import defaultdict
import os
import warnings
warnings.filterwarnings('ignore')

from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astroplan as apl
import astropy.units as u

import gurobipy as gp
from gurobipy import GRB

# KPF-CC specific files
dirpath = '/Users/jack/Documents/Github/optimalAllocation/autoschedulerV2/'
sys.path.append(dirpath)
import helperFunctions as hf
import twilightFunctions as tw
import reportingFunctions as rf
import processingFunctions as pf
import mappingFunctions as mf


def runKPFCCv2(current_day,
               requestsFile,
               allocationFile,
               accessibilitiesFile,
               twilightFile,
               outputDirectory,
               STEP,
               runRound2,
               pastObservationsFile,
               semesterTemplateFile,
               turnOffOnFile,
               nonqueueMap,
               specialMaps,
               zeroOutFile,
               gurobi_output = True,
               plot_results = True,
               SolveTime = 300):

    """runKPFCCv2
    Args:
        - current_day (str) = the calendar date of the night to produce a script. Sets the "first" day of the semester from which to compute the semester schedule solution from this day forward. Format: YYYY-MM-DD.
        - requestsFile (str) = the path and file name to the CSV with all the PI requests. Confirm that column names are correct.
        - allocationFile (str) = the path and file name to the binary map of allocated nights.
        - accessibilitiesFile (str) = the path and file name to the pickle file containing a dictionary of target names and associated pre-computed 1D accessibility maps of length equal to nSlotsInSemester.
        - twilightFile (str) = the path and file name to the CSV with precomputed twilight times.
        - outputDirectory (str) = the path to the directory where all outputs of this function should be saved. It is recommended that the path be outside the git repo.
        - STEP (int) = the time, in minutes, for a single slot.
        - runRound2 (boolean) = when True, run the bonus round. When False, do not run the bonus round.
        - pastObservationsFile (str) = the path and file name of the CSV containing information on all previous observations in the semester. If equal to nofilename.csv, then we are ignoring prior observations.
        - semesterTemplateFile (str) = the path and file name of the CSV containing the visual template of the semester. For plotly plotting purposes only.
        - turnOffOnFile (str) = the path and file name of the CSV containing the pre-computed first and last day of accessiblity for each target. For plotly plotting purposes only.
        - nonqueueMap (str) = the path and file name of the CSV containining a grid of nNightsInSemester by nSlotsInNight cells where the slots set aside for non-queue RM observations are filled in with the target name.
        - gurobi_output (boolean) = a flag to turn off or on the feature of Gurobi printing to the terminal as it solves the model.
        - plot_results (boolean) = a flag to turn off or on the plotting outputs.
        - SolveTime (int) = the maximum time, in seconds, to allow Gurobi to solve the model.
    Returns:
        None
    """

    # I suggest your output directory be something so that it doesn't autosave
    # to the same directory as the run files and crowds up the GitHub repo.
    check = os.path.isdir(outputDirectory)
    if not check:
        os.makedirs(outputDirectory)
    # this represents the day of the semester we are scheduling
    # for now this parameter only sets the number of days left in the semester,
    # and therefore the number of slots in the semester
    # and therefore the size of the Yns, Wns, etc variables.
    current_day = current_day[0] # Multiple days inputted is for the TTP. Here we only are concerned with the first day of the sequence.

    start_all = time.time()

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set up logistics parameters
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    semester_start_date, semester_end_date, semesterLength, semesterYear, semesterLetter = hf.getSemesterInfo(current_day)
    semesterLength = 89 # note: special for now! delete later
    all_dates_dict = hf.buildDayDateDictionary(semester_start_date, semesterLength)
    dates_in_semester = list(all_dates_dict.keys())
    nNightsInSemester = hf.currentDayTracker(current_day, all_dates_dict)# + 1 #note: delete the plus 1!!

    nQuartersInNight = 4
    nHoursInNight = 14
    nSlotsInQuarter = int(((nHoursInNight*60)/nQuartersInNight)/STEP)
    nSlotsInNight = nSlotsInQuarter*nQuartersInNight
    nSlotsInSemester = nSlotsInNight*nNightsInSemester
    semester_grid = np.arange(0,nNightsInSemester,1)
    quarters_grid = np.arange(0,nQuartersInNight,1)
    semesterSlots_grid = np.arange(0,nSlotsInSemester,1)
    nightSlots_grid = np.arange(0,nSlotsInNight,1)
    startingSlot = all_dates_dict[current_day]*nSlotsInNight
    startingNight =  all_dates_dict[current_day]


    nSlotsInSemester -= nSlotsInNight # note: delete this line later!

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Read in files and prep targets
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    print("Reading in and prepping files")

    # Read in and prep target info
    all_targets_frame = pd.read_csv(requestsFile)

    # Build dictionary where keys are target names and
    # values are how many slots are needed to complete one exposure for that target
    alwaysRoundUp = False
    slotsNeededDict = {}
    for n,row in all_targets_frame.iterrows():
        name = row['Starname']
        # schedule multi-shots and multi-visits as if a single, long exposure. When n_shots and n_visits are both 1, this reduces down to just the stated true exposure time.
        singlevisit = row['Nominal_ExpTime']*row['N_Observations_per_Visit'] + 45*(row['N_Observations_per_Visit'] - 1)
        exptime = singlevisit*row['N_Visits_per_Night']
        if alwaysRoundUp:
            slotsneededperExposure = math.ceil(exptime/(STEP*60.))
        else:
            if exptime > STEP*60.:
                slotsneededperExposure = hf.roundSlots(exptime, STEP*60.)
            else:
                slotsneededperExposure = 1
        slotsNeededDict[name] = slotsneededperExposure

    # Sub-divide target list into only those that are more than 1 unique night of observations.
    # Single shot observations may scheduled in a future Round 1.5, TBD.
    cadenced_targets = all_targets_frame[(all_targets_frame['N_Unique_Nights_Per_Semester'] > 1)]
    cadenced_targets.reset_index(inplace=True)
    ntargets = len(cadenced_targets)

    # Read in twilight times
    twilight_frame = pd.read_csv(twilightFile, parse_dates=True)

    # AvailableSlotsInGivenNight is a 1D matrix of length nights n
    # This is not a gorubi variable, but a regular python variable
    # Each element will hold an integer which represents the true number of STEP minute slots
    # in a quarter in a given night, after accounting for non-observable times due to day/twilight.
    AvailableSlotsInGivenNight = []
    for date in dates_in_semester:
        edge = tw.determineTwilightEdges(date, twilight_frame, STEP, verbose=False)
        AvailableSlotsInGivenNight.append(edge)

    # Build the twilight times maps
    twilightMap_all = np.array(mf.buildTwilightMap(AvailableSlotsInGivenNight, nSlotsInNight, invert=False))

    twilightMap_toDate = twilightMap_all.copy()
    twilightMap_toDate = twilightMap_toDate[all_dates_dict[current_day]:]

    twilightMap_all_flat = twilightMap_all.copy()
    twilightMap_all_flat = twilightMap_all_flat.flatten()

    twilightMap_toDate_flat = twilightMap_toDate.copy()
    twilightMap_toDate_flat = twilightMap_toDate_flat.flatten()

    # Read in serialized pre-computed accessibility maps and deserialize it
    rewriteFlag = False
    accessmaps_precompute = mf.readAccessibilityMapsDict(accessibilitiesFile)
    #for name in all_targets_frame['Starname']:
    for n,row in all_targets_frame.iterrows():
        name = row['Starname']
        # check that this target has a precomputed accessibility map, if not, make one and add it to the pickle file
        try:
            tmpRead = accessmaps_precompute[name]
        except:
            print(name + " not found in precomputed accessibilty maps. Running now.")
            newmap = mf.singleTargetAccessible(name, row['RA'], row['Dec'], semester_start_date, semesterLength-1, STEP) #note the -1 is just for the 2024B shortened semester, delete later
            accessmaps_precompute[name] = np.array(newmap).flatten()
            rewriteFlag = True
    if rewriteFlag:
        mf.writeAccessibilityMapsDict(accessmaps_precompute, '/Users/jack/Desktop/tmpaccessmap.txt') #accessibilitiesFile)

    # read in the customized maps for unique targets
    if specialMaps != 'nofilename.txt':
        uniqueTargetMaps = mf.readAccessibilityMapsDict(specialMaps)
    else:
        uniqueTargetMaps = {}

    # list of stars to be purposefully excluded from tonight's script
    # for generating more P1 stars as gapFillers
    if zeroOutFile != 'nofilename.txt':
        # zeroOut_names = np.loadtxt(zeroOutFile, delimiter=' ', dtype=str)
        zeroOut = pd.read_csv(zeroOutFile)
        zeroOut_names = list(zeroOut['Target'])
    else:
        zeroOut_names = []


    print("Solving the semester schedule problem.")

    # Read in allocation info
    # Then a simple conversion from human-readable to computer-readable
    # I wanted to keep the data in the Binary Schedule for easy manual editing if needed
    allocation_schedule_load = np.loadtxt(allocationFile, dtype=str)
    allocation_schedule_full = []
    for a in range(all_dates_dict[current_day], len(allocation_schedule_load)):# - all_dates_dict[current_day]): #maybe revert back without the minus
        convert = list(map(int, allocation_schedule_load[a][2:]))
        allocation_schedule_full.append(convert)

    # Repeat but include the entire semester's worth of dates
    allocation_schedule_full_semester = []
    for b in range(len(allocation_schedule_load)):
        convert = list(map(int, allocation_schedule_load[b][2:]))
        allocation_schedule_full_semester.append(convert)

    # determine the first and last allocated day of the semester, for COF plotting purposes only
    firstday = 0
    lastday = len(allocation_schedule_full_semester)
    for c in range(len(allocation_schedule_full_semester)):
        if np.sum(allocation_schedule_full_semester[c]) > 0:
            firstday = c
            break
    for e in range(len(allocation_schedule_full_semester)-1, 1, -1):
        if np.sum(allocation_schedule_full_semester[e]) > 0:
            lastday = e
            break
    startends = [semester_start_date, dates_in_semester[firstday], dates_in_semester[lastday], semester_end_date]

    # Randomly sample out 30% of future allocated quarters to simulate weather loss
    print("Sampling out weather losses")

    protectNonQueueNights = False
    if protectNonQueueNights and nonqueueMap != 'nofilename.csv':
        nonqueuemap_slots_strs = np.loadtxt(nonqueueMap, delimiter=',', dtype=str)
        nonqueuemap_slots_ints = []
        for i in range(len(nonqueuemap_slots_strs)):
            holder = []
            for j in range(len(nonqueuemap_slots_strs[i])):
                if nonqueuemap_slots_strs[i][j] == '':
                    holder.append(1)
                else:
                    holder.append(0)
            nonqueuemap_slots_ints.append(holder)
        summs = np.sum(nonqueuemap_slots_ints, axis=1)
        nonqueueNightMask = summs < len(tmp[0])
        protected = np.where(nonqueueNightMask)[0]
        protectedQuarters = []
        for k in range(len(protected)):
            protectedQuarters.append(protected[k]*4)
            protectedQuarters.append(protected[k]*4 + 1)
            protectedQuarters.append(protected[k]*4 + 2)
            protectedQuarters.append(protected[k]*4 + 3)
    else:
        protectedQuarters = []

    allocation_schedule_long = np.array(allocation_schedule_full).flatten()
    allindices = [i for i in range(len(allocation_schedule_long)) if allocation_schedule_long[i] == 1 and i not in protectedQuarters]
    weatherlosses = np.random.choice(allindices[4:], int(0.3*len(allindices)), replace=False)
    for w in weatherlosses:
       allocation_schedule_long[w] = 0
    allocation_schedule = np.reshape(allocation_schedule_long, (nNightsInSemester, 4))
    weatherDiff = allocation_schedule_full - allocation_schedule
    weatherDiff_1D = weatherDiff.flatten()

    # Build the allocation map
    allocation_map, allocation_map_NS, weathered_map = mf.buildAllocationMap(allocation_schedule, weatherDiff, AvailableSlotsInGivenNight[startingNight:], nSlotsInNight)
    allocation_map_fullsemester, allocation_map_NS_fullsemester, weathered_map_ignorethis = mf.buildAllocationMap(allocation_schedule_full_semester, np.zeros(np.shape(allocation_schedule_full_semester)), AvailableSlotsInGivenNight, nSlotsInNight)

    # create the intersection of the two
    alloAndTwiMap = allocation_map&twilightMap_toDate_flat

    # output the maps for testing
    # np.savetxt('/Users/jack/Desktop/availSlots.txt', AvailableSlotsInGivenNight, delimiter=',', fmt='%1i')
    # np.savetxt('/Users/jack/Desktop/alloMaps.txt', allocation_map_NS, delimiter=',', fmt='%1i')
    # np.savetxt('/Users/jack/Desktop/alloMaps_all.txt', allocation_map_NS_fullsemester, delimiter=',', fmt='%1i')
    # np.savetxt('/Users/jack/Desktop/twiMaps.txt', twilightMap_all, delimiter=',', fmt='%1i')
    # np.savetxt('/Users/jack/Desktop/intersection.txt', np.reshape(alloAndTwiMap, (nNightsInSemester, nSlotsInNight)), delimiter=',', fmt='%1i')
    # print("OUTPUTTED MAPS")

    # Pull the database of past observations and process.
    # For each target, determine the most recent date of observations and the number of unique days observed.
    # Also process the past to build the starmap for each target.
    pastObs_Info = {}
    if pastObservationsFile != 'nofilename.csv':
        pf.getKPFAllObservations(pastObservationsFile)
        database = pd.read_csv(pastObservationsFile)
        for i in range(len(all_targets_frame['Starname'])):
            starmask = database['star_id'] == all_targets_frame['Starname'][i]
            star_past_obs = database[starmask]
            star_past_obs.sort_values(by='utctime', inplace=True)
            star_past_obs.reset_index(inplace=True)

            total_past_observations = int(len(star_past_obs)/(all_targets_frame['N_Visits_per_Night'][i]*all_targets_frame['N_Observations_per_Visit'][i]))
            # print(all_targets_frame['Starname'][i], len(star_past_obs), total_past_observations, all_targets_frame['N_Visits_per_Night'][i], all_targets_frame['N_Observations_per_Visit'][i])

            #total_open_shutter_time = np.sum(star_past_obs['Nominal_ExpTime'])
            star_past_obs, unique_hstdates_observed, quarterObserved = pf.getUniqueNights(star_past_obs, twilight_frame)
            if len(unique_hstdates_observed) > 0:
                mostRecentObservationDate = unique_hstdates_observed[-1]
            else:
                mostRecentObservationDate = '0000-00-00'
            Nobs_on_date = pf.getNobs_on_Night(star_past_obs, unique_hstdates_observed)

            # within the pastObs_Info dictionary, each target's data is always in the given order:
            # element 0 = the calendar date of the most recent observation (HST date)
            # element 1 = the list of calendar dates with at least one observation
            # element 2 = the list of the quarter where the observation took place on the corresponding unique night from element 1's list (if multi visits, then this is for the first visit)
            # element 3 = the list of the number of observations that took place on the corresponding unique night from element 1's list
            pastObs_Info[all_targets_frame['Starname'][i]] = [mostRecentObservationDate, unique_hstdates_observed, quarterObserved, Nobs_on_date]


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set Gorubi model variables
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    print("Building variables")
    m = gp.Model('Semester_Scheduler')

    # Yns (ExposureStartSlot) is a 2D matrix of targets n and slots s
    # Slot s for target n will be 1 to indicate starting an exposure
    Yns = m.addVars(all_targets_frame['Starname'], semesterSlots_grid, vtype = GRB.BINARY, name = 'ExposureStart_Slot')

    # Wnd (ExposureStartNight) is a 2D matrix of targets t and nights d
    # Night d for target n will be 1 to indicate this target gets an exposure on this night
    # useful for tracking cadence and building the DeltaDays variable
    Wnd = m.addVars(all_targets_frame['Starname'], semesterSlots_grid, vtype = GRB.BINARY, name = 'ExposureStart_Night')

    # Thn (Theta_n) is the "shortfall" variable
    Thn = m.addVars(all_targets_frame['Starname'], name='Theta_n')


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set Gorubi constraints
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------

    print("Constraint: relating Yns to Wns")
    # Relate Yns to Wns
    # For target n, only one slot within the night can be set to 1 (one exposure per night)
    # For target n, if any slot s within night d is 1, then Wnd must be 1
    for d in range(nNightsInSemester):
        start = d*nSlotsInNight
        end = start + nSlotsInNight
        for name in all_targets_frame['Starname']:
            m.addConstr((gp.quicksum(Yns[name,s] for s in range(start, end)) <= 1), 'oneObsPerNight_' + str(name) + "_" + str(d) + "d")
            for s in range(start, end):
                m.addConstr(Yns[name, s] <=  Wnd[name, d], "related_Yns_Wnd_" + str(name) + "_" + str(d) + "_" + str(s))

    print("Constraint: no other exposures can start in multi-slot exposures")
    # Enforce no other exposures start within t slots if Yns[name, slot] is 1
    # Velibor's new logic
    for n,row in all_targets_frame.iterrows():
        name = row['Starname']
        slotsneededperExposure = slotsNeededDict[name]
        if (slotsneededperExposure > 1):
            for s in range(nSlotsInSemester):
                m.addConstr((slotsneededperExposure-1)*Yns[name,s] <= ( (slotsneededperExposure - 1) - gp.quicksum(Yns[name2,s+t] for name2 in all_targets_frame['Starname'] for t in range(1, min(slotsneededperExposure, nSlotsInSemester - s)))), 'dont_Start_Exposure_' + str(name) + "_" + str(s) + "s")

    print("Constraint: only 1 exposure per slot")
    # Enforce only one exposure per slot
    for s in range(nSlotsInSemester):
        m.addConstr(gp.quicksum(Yns[name,s] for name in all_targets_frame['Starname']) <= 1, 'ensure_singleOb_perSlot_' + str(s) + "s")

    print("Constraint: exposure can't start if it won't complete in the night/quarter.")
    # Example: a full night allocation where the 2nd quarter of the night is "weathered". In this scenario,
    # This constraint must be applied at both the end of the 1st quarter and the end of the 4th quarter.
    # But this constraint is not applied to the 3rd quarter: an exposure can span the boundary of 3 to 4 quarter.
    # Enforce that an exposure cannot start if it will not complete within the same night/quarter.
    for d in range(nNightsInSemester):
        for s in range(nSlotsInNight - 1): # -1 so that we don't hit the end of the loop. Plus the last and second to last slot should never be allocated (always set to 0) so no need to test.
            if allocation_map_NS[d][s] == 1 and allocation_map_NS[d][s+1] == 0:
                for name in all_targets_frame['Starname']:
                    if slotsNeededDict[name] > 1:
                        for e in range(0, slotsNeededDict[name]-1): # the -1 is because a target can be started if it just barely fits
                            m.addConstr(Yns[name,d*nSlotsInNight + s - e] == 0, 'dont_start_near_endof_night_' + name + "_" + str(s) + 's_' + str(e) + 'e')

    print("Constraint: inter-night cadence of future observations.")
    # Constrain the future inter-night cadence
    counter = 0
    for d in range(nNightsInSemester):
        for t,row in all_targets_frame.iterrows():
            name = row['Starname']
            internightcadence = int(row['Inter_Night_Cadence'])
            for dlta in range(1, internightcadence):
                try:
                    m.addConstr(Wnd[name,d] <= 1 - Wnd[name,d+dlta], 'enforce_internightCadence_' + str(name) + "_" + str(d) + "d_" + str(dlta) + "dlta")
                except:
                    # enter this except when the inter-night cadence brings us beyond the end of the semester (i think)
                    counter += 1
                    continue

    print("Constraint: multi-visit targets can only be attempted under certain conditions.")
    # Constrain the nights when a multi-visit target can be scheduled
    alloAndTwiMap_d = copy.deepcopy(alloAndTwiMap)
    alloAndTwiMap_d = np.reshape(alloAndTwiMap_d, (nNightsInSemester, nSlotsInNight))
    for t,row in all_targets_frame.iterrows():
        name = row['Starname']
        visits = int(row['N_Visits_per_Night'])
        if visits > 1:
            cadence = int(row['Intra_Night_Cadence'])
            # this equation is a political decision and can be modified. It states that for each visit, after the intra-night cadence time has elapsed,
            # we require a 90 minute window within which to allow for scheduling the next visit. We then assume that this next visit is scheduled
            # at the very end of this 90 minute window, which then restarts the clock for any future visits to require the same 90 min grace period after
            # the intra-night cadence time is satisfied.
            minTimeRequired = (visits - 1)*(cadence + 1.5) # hours
            minimumSlotsRequired = math.ceil(minTimeRequired/(STEP*60.))
            all_access_target = np.array(accessmaps_precompute[name])
            # newexptime = (visits-1)*200 + row['N_Observations_per_Visit']*45 + row['N_Observations_per_Visit']*row['Nominal_ExpTime']
            newexptime = visits * (row['N_Observations_per_Visit']*45 + row['N_Observations_per_Visit']*row['Nominal_ExpTime'])
            new_nslots = math.ceil(newexptime/(STEP*60.))
            for d in range(nNightsInSemester):
                possibleOpenSlots = np.sum(alloAndTwiMap_d[d]&all_access_target[d])
                if possibleOpenSlots < minimumSlotsRequired:
                    m.addConstr(Wnd[name,d] == 0, 'cannotMultiVisit_' + str(name) + "_" + str(d) + "d")


    print("Constraint: build Theta variable")
    # Construct the Theta_n variable
    extra = 5
    for t,row in all_targets_frame.iterrows():
        name = row['Starname']
        Nobs_Unique = row['N_Unique_Nights_Per_Semester']
        if pastObs_Info == {}:
            past_Unique = 0
        else:
            past_Unique = len(pastObs_Info[name][1])
        m.addConstr(Thn[name] >= 0, 'ensureGreaterThanZero_' + str(name))
        m.addConstr(Thn[name] >= ((Nobs_Unique - past_Unique) - gp.quicksum(Yns[name, s] for s in range(nSlotsInSemester))), 'ensureGreaterThanNobsShortfall_' + str(name))
        # add an upper limit to how many extra obs can be scheduled for a single target
        m.addConstr((gp.quicksum(Wnd[name,d] for d in range(nNightsInSemester)) <= (Nobs_Unique - past_Unique) + extra), 'max_Nobs_Unique_Nights_' + str(name))

    print("Constraint: enforce twilight times")
    print("Constraint: only observe if accessible")
    print("Constraint: enforce no observations when not allocated.")
    print("Constraint: enforce custom maps.")
    uniqueTargetMap_names = list(uniqueTargetMaps.keys())
    for name in all_targets_frame['Starname']:
        # all_access = np.array(accessmaps_precompute[name]).flatten()
        all_access = accessmaps_precompute[name]
        startslot = (all_dates_dict[current_day])*nSlotsInNight # plus 1 to account for python indexing?
        access = all_access[startslot:]

        alloAndTwiMap_short = alloAndTwiMap[:-nSlotsInNight] # note: temporary, delete this later.
        #alloAndTwiMap_short = alloAndTwiMap[startslot:]

        # fullmap = alloAndTwiMap&access[:len(alloAndTwiMap)]
        # fullmap = alloAndTwiMap&access
        if name in uniqueTargetMap_names:
            # customMap = uniqueTargetMaps[name][startslot:-nSlotsInNight]
            customMap = uniqueTargetMaps[name][startslot:]
        else:
            customMap = np.array([1]*nSlotsInSemester)

        # enforce that certain targets not be allowed to be scheduled tonight
        # the idea is this creates a second script from which we can pull more P1 targets
        # to fill gaps in the schedule.
        # This is because the TTP is too good at saving us time!
        zeroMap = np.array([1]*nSlotsInSemester)
        if name in zeroOut_names:
            zeroMap[:nSlotsInNight] = np.array([0]*nSlotsInNight)

        #print(name, np.shape(alloAndTwiMap), np.shape(alloAndTwiMap_short), np.shape(access), np.shape(customMap), np.shape(zeroMap))
        fullmap = alloAndTwiMap_short&access&customMap&zeroMap
        # fullmap = alloAndTwiMap&access&customMap&zeroMap
        for s in range(nSlotsInSemester):
            m.addConstr(Yns[name,s] <= fullmap[s], 'enforceMaps_' + str(name) + "_" + str(s) + "s")

    print("Constraint: first observation of new schedule can't violate past inter-night cadence.")
    for t,row in all_targets_frame.iterrows():
        name = row['Starname']
        internightcadence = int(row['Inter_Night_Cadence'])
        if pastObs_Info == {}:
            dateLastObserved = '0000-00-00'
        else:
            dateLastObserved = pastObs_Info[name][0]
        if dateLastObserved != '0000-00-00':
            dateLastObserved_number = all_dates_dict[dateLastObserved]
            today_number = all_dates_dict[current_day]
            diff = today_number - dateLastObserved_number
            if diff < internightcadence:
                doublediff = internightcadence - diff
                for e in range(doublediff):
                    m.addConstr(Wnd[name,e] == 0, 'enforce_internightCadence_fromStart_' + str(name) + "_" + str(e) + "e")

    if nonqueueMap != 'nofilename.csv':
        print("Constraint: certain slots on allocated nights must be zero to accommodate Non-Queue observations.")
        nonqueuemap_slots_strs = np.loadtxt(nonqueueMap, delimiter=',', dtype=str)
        nonqueuemap_slots_ints = []
        for i in range(len(nonqueuemap_slots_strs)):
            holder = []
            for j in range(len(nonqueuemap_slots_strs[i])):
                if nonqueuemap_slots_strs[i][j] == '':
                    holder.append(1)
                else:
                    holder.append(0)
            nonqueuemap_slots_ints.append(holder)
        nonqueuemap_slots_ints = np.array(nonqueuemap_slots_ints).flatten()
        for s in range(nSlotsInSemester):
            nonqueueslot = int(nonqueuemap_slots_ints[s + startingSlot])
            for name in all_targets_frame['Starname']:
                m.addConstr(Yns[name,s] <= nonqueueslot, 'enforce_NonQueueSlots_' + str(name) + "_" + str(s) + "s")
    else:
        print("No non-queue observations are scheduled.")

    endtime2 = time.time()
    print("Total Time to build constraints: ", np.round(endtime2-start_all,3))

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set Gorubi Objective and Run
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    print("Setting objective.")

    # Minimize Theta_n
    m.setObjective(gp.quicksum(Thn[name] for name in all_targets_frame['Starname']), GRB.MINIMIZE)

    m.params.TimeLimit = SolveTime
    m.Params.OutputFlag = gurobi_output
    m.params.MIPGap = 0.05 # can stop at 5% gap or better to prevent it from spending lots of time on marginally better solutions
    m.params.Presolve = 2 # more aggressive presolve gives better solution in shorter time
    m.update()
    m.optimize()

    if m.Status == GRB.INFEASIBLE:
        print('Model remains infeasible. Searching for invalid constraints')
        search = m.computeIIS()
        print("Printing bad constraints:")
        for c in m.getConstrs():
            if c.IISConstr:
                print('%s' % c.ConstrName)
        for c in m.getGenConstrs():
            if c.IISGenConstr:
                print('%s' % c.GenConstrName)
    else:
        print("")
        print("")
        print("Round 1 Model Solved.")

    endtime3 = time.time()
    print("Total Time to finish solver: ", np.round(endtime3-start_all,3))


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Retrieve data from solution
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    file = open(outputDirectory + "runReport.txt", "w")
    file.close()

    print("Building human readable schedule.")
    combined_semester_schedule = hf.buildHumanReadableSchedule(Yns, twilightMap_toDate, all_targets_frame, nNightsInSemester, nSlotsInNight, AvailableSlotsInGivenNight, nSlotsInQuarter, all_dates_dict, current_day, allocation_map_NS, weathered_map, slotsNeededDict, nonqueueMap)

    round = 'Round 1'
    rf.buildFullnessReport(allocation_map, twilightMap_all, combined_semester_schedule, nSlotsInQuarter, nSlotsInSemester, all_targets_frame, outputDirectory, STEP, round)
    np.savetxt(outputDirectory + 'raw_combined_semester_schedule_Round1.txt', combined_semester_schedule, delimiter=',', fmt="%s")

    print("Writing Report.")
    filename = open(outputDirectory + "runReport.txt", "a")
    theta_n_var = []
    counter = 0
    for v in Thn.values():
        varname = v.VarName
        varval = v.X
        counter += varval
    print("Sum of Theta: " + str(counter))
    filename.write("Sum of Theta: " + str(counter) + "\n")

    if plot_results:
        print("Plotting results.")

        all_starmaps = {}
        # Only generate the COF for results of Round 1.
        for i in range(len(all_targets_frame)):
            if pastObs_Info != {}:
                starmap = rf.buildObservedMap_past(pastObs_Info[all_targets_frame['Starname'][i]][1], pastObs_Info[all_targets_frame['Starname'][i]][2], pastObs_Info[all_targets_frame['Starname'][i]][3], semesterTemplateFile, weatherDiff_1D)
            else:
                starmap = pd.read_csv(semesterTemplateFile)
                starmap_columns = starmap.columns
                if 'Observed' not in starmap_columns:
                    starmap['Observed'] = [False]*len(starmap)
                if 'N_obs' not in starmap_columns:
                    starmap['N_obs'] = [0]*len(starmap)
                unique_hstdates_observed = []
            starmap_updated = rf.buildObservedMap_future(all_targets_frame['Starname'][i], slotsNeededDict[all_targets_frame['Starname'][i]], combined_semester_schedule, starmap, all_dates_dict[current_day], slotsNeededDict, np.array(allocation_schedule_full_semester), weatherDiff_1D, nSlotsInNight)
            #if optimalAllocation:
            #    pd.to_csv(outputDirectory + "/FirstForecasts/Forecast_" + str(all_targets_frame['Starname'][i]) + "_semester.csv", index=False)
            all_starmaps[all_targets_frame['Starname'][i]] = starmap_updated
            rf.writeCadencePlotFile(all_targets_frame['Starname'][i], i, starmap, turnOffOnFile, all_targets_frame, outputDirectory, unique_hstdates_observed, current_day)

        rf.buildAllocationPicture(allocation_schedule_full, nNightsInSemester, nQuartersInNight, startingNight, all_dates_dict, outputDirectory)
        sum_schedule = np.sum(allocation_schedule_full_semester, axis=1)
        allocated_days = np.nonzero(sum_schedule)
        firstDayAllocated = dates_in_semester[allocated_days[0][0]]
        lastDayAllocated = dates_in_semester[allocated_days[0][-1]]
        notable_dates = [semester_start_date, firstDayAllocated, lastDayAllocated, semester_end_date]
        rf.buildCOF(outputDirectory, current_day, startends, all_targets_frame, all_dates_dict, all_starmaps, notable_dates, False)

    endtime4 = time.time()
    print("Total Time to complete Round 1: " + str(np.round(endtime4-start_all,3)))
    filename.write("Total Time to complete Round 1: " + str(np.round(endtime4-start_all,3)) + "\n")
    filename.write("\n")
    filename.write("\n")
    filename.close()
    print("Round 1 complete.")


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Initiate Round 2 Scheduling
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------

    if runRound2:
        print("Beginning Round 2 Scheduling.")
        # Extract theta values:
        first_stage_objval = m.objval
        eps = 5 # allow a tolerance

        m.params.TimeLimit = SolveTime
        m.Params.OutputFlag = gurobi_output
        m.params.MIPGap = 0.05 # stop at 5% gap, this is good enough and leaves some room for including standard stars
        m.addConstr(gp.quicksum(Thn[name] for name in all_targets_frame['Starname']) <= first_stage_objval + eps)
        m.setObjective(gp.quicksum(slotsNeededDict[name]*Yns[name,s] for name in all_targets_frame['Starname'] for s in range(nSlotsInSemester)), GRB.MAXIMIZE)
        m.update()
        m.optimize()

        if m.Status == GRB.INFEASIBLE:
            print('Model remains infeasible. Searching for invalid constraints')
            search = m.computeIIS()
            print("Printing bad constraints:")
            for c in m.getConstrs():
                if c.IISConstr:
                    print('%s' % c.ConstrName)
            for c in m.getGenConstrs():
                if c.IISGenConstr:
                    print('%s' % c.GenConstrName)
        else:
            print("")
            print("")
            print("Round 2 Model Solved.")

        endtime5 = time.time()
        print("Total Time to finish round 2 solver: ", np.round(endtime5-start_all,3))

        combined_semester_schedule = hf.buildHumanReadableSchedule(Yns, twilightMap_toDate, all_targets_frame, nNightsInSemester, nSlotsInNight, AvailableSlotsInGivenNight, nSlotsInQuarter, all_dates_dict, current_day, allocation_map_NS, weathered_map, slotsNeededDict, nonqueueMap)

        round = 'Round 2'
        rf.buildFullnessReport(allocation_map, twilightMap_all, combined_semester_schedule, nSlotsInQuarter, nSlotsInSemester, all_targets_frame, outputDirectory, STEP, round)
        np.savetxt(outputDirectory + 'raw_combined_semester_schedule_Round2.txt', combined_semester_schedule, delimiter=',', fmt="%s")


        scheduleR1 = np.loadtxt(outputDirectory + 'raw_combined_semester_schedule_Round1.txt', delimiter=',', dtype=str)
        scheduleR2 = np.loadtxt(outputDirectory + 'raw_combined_semester_schedule_Round2.txt', delimiter=',', dtype=str)
        gapStars = hf.getGapFillerTargets(scheduleR1, scheduleR2, all_dates_dict[current_day])
        np.savetxt(outputDirectory + 'gapFillerTargets.txt', gapStars, delimiter=',', fmt="%s")

        filename = open(outputDirectory + "runReport.txt", "a")
        theta_n_var = []
        counter = 0
        for v in Thn.values():
            varname = v.VarName
            varval = v.X
            counter += varval
        filename.write("Sum of Theta: " + str(counter) + "\n")

        endtime6 = time.time()
        filename.write("Total Time to complete Round 2: " + str(np.round(endtime5-start_all,3)) + "\n")
        filename.write("\n")
        filename.write("\n")
        filename.close()
        print("Round 2 complete.")
    else:
        print("Not running Round 2. Duplicating Raw Schedule as dummy file.")
        np.savetxt(outputDirectory + 'raw_combined_semester_schedule_Round2.txt', combined_semester_schedule, delimiter=',', fmt="%s")
        gapStars = []
        np.savetxt(outputDirectory + 'gapFillerTargets.txt', gapStars, delimiter=',', fmt="%s")

    print("done")
