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
               slotsize,
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
        - slotsize (int) = the time, in minutes, for a single slot.
        - runRound2 (boolean) = when True, run the bonus round. When False, do not run the bonus round.
        - pastObservationsFile (str) = the path and file name of the CSV containing information on all previous observations in the semester. If file does not exist, then we are ignoring prior observations.
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
    # and therefore the size of the Yrs, Wrd, etc variables.
    current_day = current_day[0] # Multiple days inputted is for the TTP. Here we only are concerned with the first day of the sequence.

    start_the_clock = time.time()

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set up logistics parameters
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    semester_start_date, semester_end_date, semesterLength, semesterYear, semesterLetter = hf.getSemesterInfo(current_day)
    all_dates_dict = hf.buildDayDateDictionary(semester_start_date, semesterLength)
    all_dates_array = list(all_dates_dict.keys())
    nNightsInSemester = hf.currentDayTracker(current_day, all_dates_dict)
    print("There are " + str(nNightsInSemester) + " calendar nights remaining in this semester.")

    nQuartersInNight = 4
    nHoursInNight = 14
    nSlotsInQuarter = int(((nHoursInNight*60)/nQuartersInNight)/slotsize)
    nSlotsInNight = nSlotsInQuarter*nQuartersInNight
    nSlotsInSemester = nSlotsInNight*nNightsInSemester
    semester_grid = np.arange(0,nNightsInSemester,1)
    quarters_grid = np.arange(0,nQuartersInNight,1)
    semesterSlots_grid = np.arange(0,nSlotsInSemester,1)
    nightSlots_grid = np.arange(0,nSlotsInNight,1)
    # the slot and night from which represents this date in the semester onward
    startingSlot = all_dates_dict[current_day]*nSlotsInNight
    startingNight =  all_dates_dict[current_day]

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Read in files and prep targets
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    print("Reading in and prepping files")

    # Read in and prep target info
    requests_frame = pd.read_csv(requestsFile)

    # Build dictionary where keys are target names and
    # values are how many slots are needed to complete one exposure for that target
    slotsNeededForExposure_dict = {}
    for n,row in requests_frame.iterrows():
        name = row['Starname']
        # schedule multi-shots and multi-visits as if a single, long exposure.
        # When n_shots and n_visits are both 1, this reduces down to just the stated true exposure time.
        singlevisit = row['Nominal Exposure Time [s]']*row['# of Exposures per Visit'] + 45*(row['# Visits per Night'] - 1)
        exptime = singlevisit*row['# Visits per Night']
        slotsNeededForExposure_dict[name] = hf.slotsRequired(exptime, slotsize)

    # Sub-divide target list into only those that are more than 1 unique night of observations.
    # Currently this is not used anywhere else. May be deleted. - Jack Oct. 29, 2024
    cadenced_targets = requests_frame[(requests_frame['# of Nights Per Semester'] > 1)]
    cadenced_targets.reset_index(inplace=True)
    ntargets = len(cadenced_targets)

    # Read in twilight times
    twilight_frame = pd.read_csv(twilightFile, parse_dates=True)

    # AvailableSlotsInGivenNight is a 1D matrix of length nights n
    # This is not a gorubi variable, but a regular python variable
    # Each element will hold an integer which represents the true number of slotsize minute slots
    # in a quarter in a given night, after accounting for non-observable times due to day/twilight.
    AvailableSlotsInGivenNight = []
    for date in all_dates_array:
        slotsTonight = tw.determineTwilightEdges(date, twilight_frame, slotsize, verbose=False)
        AvailableSlotsInGivenNight.append(slotsTonight)

    # Build the twilight times maps
    twilightMap_all = np.array(mf.buildTwilightMap(AvailableSlotsInGivenNight, nSlotsInNight, invert=False))
    twilightMap_toDate = twilightMap_all[all_dates_dict[current_day]:]
    twilightMap_toDate_flat = twilightMap_toDate.copy().flatten()

    # Read in serialized pre-computed accessibility maps and deserialize it
    rewriteFlag = False
    accessmaps_precompute = mf.readAccessibilityMapsDict(accessibilitiesFile)
    for n,row in requests_frame.iterrows():
        name = row['Starname']
        # check that this target has a precomputed accessibility map, if not, make one and add it to the pickle file
        try:
            tmpRead = accessmaps_precompute[name]
        except:
            print(name + " not found in precomputed accessibilty maps. Running now.")
            newmap = mf.singleTargetAccessible(name, row['RA'], row['Dec'], semester_start_date, semesterLength-1, slotsize) #note the -1 is to account for python indexing
            accessmaps_precompute[name] = np.array(newmap).flatten()
            rewriteFlag = True
    if rewriteFlag:
        # overwrite with the updated file
        mf.writeAccessibilityMapsDict(accessmaps_precompute, accessibilitiesFile)

    # read in the customized maps for unique targets
    if os.path.exists(specialMaps):
        customRequestMaps = mf.readAccessibilityMapsDict(specialMaps)
    else:
        customRequestMaps = {}

    # List of requests to be purposefully excluded from tonight's script
    if os.path.exists(zeroOutFile):
        zeroOut = pd.read_csv(zeroOutFile)
        zeroOut_names = list(zeroOut['Target'])
    else:
        zeroOut_names = []

    # Read in allocation info
    # Conversion from human to computer-readable
    allocation_raw = np.loadtxt(allocationFile, dtype=str)
    allocation_toDate = []
    allocation_all = []
    for a in range(semesterLength):
        convert = list(map(int, allocation_raw[a][2:]))
        allocation_all.append(convert)
        if a >= all_dates_dict[current_day]:
            allocation_toDate.append(convert)

    # Randomly sample out 30% of future allocated quarters to simulate weather loss
    print("Sampling out weather losses")
    protectNonQueueNights = False
    if protectNonQueueNights and os.path.exists(nonqueueMap):
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
    allocation_toDate_1D = np.array(allocation_toDate).flatten()
    allindices = [i for i in range(len(allocation_toDate_1D)) if allocation_toDate_1D[i] == 1 and i not in protectedQuarters]
    weatherlosses = np.random.choice(allindices[4:], int(0.3*len(allindices)), replace=False)
    for w in weatherlosses:
       allocation_toDate_1D[w] = 0

    allocation_postLoss = np.reshape(allocation_toDate_1D, (nNightsInSemester, 4))
    weatherDiff = allocation_toDate - allocation_postLoss
    weatherDiff_1D = weatherDiff.flatten()

    # Build the allocation map
    allocation_map_1D, allocation_map_2D, weathered_map = mf.buildAllocationMap(allocation_postLoss, weatherDiff, AvailableSlotsInGivenNight[startingNight:], nSlotsInNight)
    # allocation_map_fullsemester, allocation_map_NS_fullsemester, weathered_map_ignorethis = mf.buildAllocationMap(allocation_all, np.zeros(np.shape(allocation_all)), AvailableSlotsInGivenNight, nSlotsInNight)

    # create the intersection of the two: this is the queue's "observability_1D"
    observability_1D = allocation_map_1D&twilightMap_toDate_flat
    observability_2D = copy.deepcopy(observability_1D)
    observability_2D = np.reshape(observability_2D, (nNightsInSemester, nSlotsInNight))

    # Pull the database of past observations and process.
    # For each target, determine the most recent date of observations and the number of unique days observed.
    # Also process the past to build the starmap for each target.
    pastObs_Info = {}
    pf.getKPFAllObservations(pastObservationsFile)
    if os.path.exists(pastObservationsFile):
        database = pd.read_csv(pastObservationsFile)
        for i in range(len(requests_frame['Starname'])):
            starmask = database['star_id'] == requests_frame['Starname'][i]
            star_past_obs = database[starmask]
            star_past_obs.sort_values(by='utctime', inplace=True)
            star_past_obs.reset_index(inplace=True)
            total_past_observations = int(len(star_past_obs)/(requests_frame['# Visits per Night'][i]*requests_frame['# of Exposures per Visit'][i]))
            star_past_obs, unique_hstdates_observed, quarterObserved = pf.getUniqueNights(star_past_obs, twilight_frame)
            if len(unique_hstdates_observed) > 0:
                mostRecentObservationDate = unique_hstdates_observed[-1]
            else:
                mostRecentObservationDate = '0000-00-00'
            Nobs_on_date = pf.getNobs_on_Night(star_past_obs, unique_hstdates_observed)
            # within the pastObs_Info dictionary, each target's data is always in the given order:
            # element 0 = the calendar date of the most recent observation (HST date)
            # element 1 = the list of calendar dates with at least one observation
            # element 2 = the list of the quarters where the observation took place on the corresponding unique night from element 1's list (if multi visits, then this is for the first visit)
            # element 3 = the list of the number of observations that took place on the corresponding unique night from element 1's list
            pastObs_Info[requests_frame['Starname'][i]] = [mostRecentObservationDate, unique_hstdates_observed, quarterObserved, Nobs_on_date]


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set Gorubi model variables
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    print("Building variables")
    m = gp.Model('Semester_Scheduler')

    # Yrs is a 2D matrix of requests r and slots s
    # Slot s for request r will be 1 to indicate starting an exposure for that request in that slot
    Yrs = m.addVars(requests_frame['Starname'], semesterSlots_grid, vtype = GRB.BINARY, name = 'Requests_Slots')

    # Wrd is a 2D matrix of requests r and nights d
    # Night d for request r will be 1 to indicate this requests gets an exposure on this night
    # Useful for tracking cadence and building the DeltaDays variable
    Wrd = m.addVars(requests_frame['Starname'], semesterSlots_grid, vtype = GRB.BINARY, name = 'Requests_Nights')

    # theta is the "shortfall" variable
    theta = m.addVars(requests_frame['Starname'], name='Shortfall')


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set Gorubi constraints
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------

    print("Constraint: relating Yrs to Wrd")
    # Relate Yrs to Wrd
    # For request r, only one slot within the night can be set to 1 (one exposure per night)
    # For request r, if any slot s within night d is 1, then Wrd must be 1
    for d in range(nNightsInSemester):
        start = d*nSlotsInNight
        end = start + nSlotsInNight
        for name in requests_frame['Starname']:
            m.addConstr((gp.quicksum(Yrs[name,s] for s in range(start, end)) <= 1), 'oneObservationPerNight_' + str(name) + "_" + str(d) + "d")
            for s in range(start, end):
                m.addConstr(Yrs[name, s] <=  Wrd[name, d], "relate_Yrs_Wrd_" + str(name) + "_" + str(d) + "_" + str(s))

    print("Constraint: no other exposures can start in multi-slot exposures")
    # Enforce no other exposures start within so many slots if Yrs[name, slot] is 1
    for n,row in requests_frame.iterrows():
        name = row['Starname']
        print(name)
        slotsNeeded = slotsNeededForExposure_dict[name]
        if (slotsNeeded > 1):
            for s in range(nSlotsInSemester):
                m.addConstr((slotsNeeded-1)*Yrs[name,s] <= ((slotsNeeded - 1) - gp.quicksum(Yrs[name2,s+t] for name2 in requests_frame['Starname'] for t in range(1, min(slotsNeeded, nSlotsInSemester - s)))), 'UsedByMultiSlotExposure_' + str(name) + "_" + str(s) + "s")

    print("Constraint: only 1 exposure per slot")
    # Enforce only one request can be in a slot
    for s in range(nSlotsInSemester):
        m.addConstr(gp.quicksum(Yrs[name,s] for name in requests_frame['Starname']) <= 1, 'ensureSingleObservationPerSlot_' + str(s) + "s")

    print("Constraint: exposure can't start if it won't complete in the night/quarter.")
    # Example: a full night allocation where the 2nd quarter of the night is "weathered". In this scenario,
    # This constraint must be applied at both the end of the 1st quarter and the end of the 4th quarter.
    # But this constraint is not applied to the 3rd quarter: an exposure can span the boundary of 3 to 4 quarter.
    # Enforce that an exposure cannot start if it will not complete within the same night/quarter.
    for d in range(nNightsInSemester):
        for s in range(nSlotsInNight - 1): # -1 so that we don't hit the end of the loop. Plus the last and second to last slot should never be allocated (always set to 0) so no need to test.
            if allocation_map_2D[d][s] == 1 and allocation_map_2D[d][s+1] == 0:
                for name in requests_frame['Starname']:
                    slotsNeeded = slotsNeededForExposure_dict[name]
                    if slotsNeeded > 1:
                        for e in range(0, slotsNeeded-1): # the -1 is because a target can be started if it just barely fits
                            try:
                                m.addConstr(Yrs[name,d*nSlotsInNight + s - e] == 0, 'mustCompleteWithinNight_' + name + "_" + str(s) + 's_' + str(e) + 'e')
                            except KeyError:
                                # hit the end of the array
                                continue
                            except:
                                print("Non-Key Error. Manually check: ", name, s, e)

    print("Constraint: inter-night cadence with respect to most recent observation.")
    if pastObs_Info != {}:
        for t,row in requests_frame.iterrows():
            name = row['Starname']
            dateLastObserved = pastObs_Info[name][0]
            if dateLastObserved != '0000-00-00':
                dateLastObserved_number = all_dates_dict[dateLastObserved]
                today_number = all_dates_dict[current_day]
                diff = today_number - dateLastObserved_number
                if diff < int(row['Minimum Inter-Night Cadence']):
                    blockUpcomingDays = int(row['Minimum Inter-Night Cadence']) - diff
                    for e in range(blockUpcomingDays):
                        m.addConstr(Wrd[name,e] == 0, 'enforcePastInterNightCadence_' + str(name) + "_" + str(e) + "e")

    print("Constraint: inter-night cadence of future observations.")
    # Constrain the future inter-night cadence
    for d in range(nNightsInSemester):
        for t,row in requests_frame.iterrows():
            name = row['Starname']
            for dlta in range(1, int(row['Minimum Inter-Night Cadence'])):
                try:
                    m.addConstr(Wrd[name,d] <= 1 - Wrd[name,d+dlta], 'enforceFutureInterNightCadence_' + str(name) + "_" + str(d) + "d_" + str(dlta) + "dlta")
                except:
                    # enter this except when the inter-night cadence brings us beyond the end of the semester
                    continue

    print("Constraint: multi-visit targets can only be attempted under certain conditions.")
    # Constrain the nights when a multi-visit target can be scheduled
    for t,row in requests_frame.iterrows():
        name = row['Starname']
        if int(row['# Visits per Night']) > 1:
            # this equation is a political decision and can be modified. It states that for each visit, after the intra-night cadence time has elapsed,
            # we require a 90 minute window within which to allow for scheduling the next visit. We then assume that this next visit is scheduled
            # at the very end of this 90 minute window, which then restarts the clock for any future visits to require the same 90 min grace period after
            # the intra-night cadence time is satisfied.
            minTimeRequired = ((int(row['# Visits per Night']) - 1)*(int(row['Minimum Intra-Night Cadence']) + 1.5))*3600 #convert hours to seconds
            minimumSlotsRequired = hf.slotsRequired(minTimeRequired, slotsize)
            accessibilty_r = np.array(accessmaps_precompute[name])
            for d in range(nNightsInSemester):
                possibleOpenSlots = np.sum(observability_2D[d]&accessibilty_r[d])
                if possibleOpenSlots < minimumSlotsRequired:
                    m.addConstr(Wrd[name,d] == 0, 'cannotMultiVisitTonight_' + str(name) + "_" + str(d) + "d")

    print("Constraint: build Theta variable")
    # Construct the Theta_n variable
    B = 5
    for t,row in requests_frame.iterrows():
        name = row['Starname']
        if pastObs_Info == {}:
            past_Unique = 0
        else:
            past_Unique = len(pastObs_Info[name][1])
        m.addConstr(theta[name] >= 0, 'ensureGreaterThanZero_' + str(name))
        m.addConstr(theta[name] >= ((row['# of Nights Per Semester'] - past_Unique) - gp.quicksum(Yrs[name, s] for s in range(nSlotsInSemester))), 'ensureGreaterThanNobsShortfall_' + str(name))
        m.addConstr((gp.quicksum(Wrd[name,d] for d in range(nNightsInSemester)) <= (row['# of Nights Per Semester'] - past_Unique) + B), 'maximumUniqueNightsPerTarget_' + str(name))

    print("Constraint: enforce twilight times")
    print("Constraint: only observe if accessible")
    print("Constraint: enforce no observations when not allocated.")
    print("Constraint: enforce custom maps.")
    uniqueTargetMap_names = list(customRequestMaps.keys())
    for name in requests_frame['Starname']:
        accessibility_r = accessmaps_precompute[name]
        startslot = (all_dates_dict[current_day])*nSlotsInNight
        access = accessibility_r[startslot:]
        if name in uniqueTargetMap_names:
            customMap = customRequestMaps[name][startslot:]
        else:
            customMap = np.array([1]*nSlotsInSemester)
        # enforce that certain targets not be allowed to be scheduled tonight
        # the idea is this creates a second script from which we can pull more P1 targets
        # to fill gaps in the schedule. This is because the TTP is too good at saving us time!
        zeroMap = np.array([1]*nSlotsInSemester)
        if name in zeroOut_names:
            zeroMap[:nSlotsInNight] = np.array([0]*nSlotsInNight)
        #print(name, np.shape(observability_1D), np.shape(alloAndTwiMap_short), np.shape(access), np.shape(customMap), np.shape(zeroMap))
        fullmap = observability_1D&access&customMap&zeroMap
        for s in range(nSlotsInSemester):
            m.addConstr(Yrs[name,s] <= fullmap[s], 'enforceIntersectionOfMaps_' + str(name) + "_" + str(s) + "s")

    if os.path.exists(nonqueueMap):
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
            for name in requests_frame['Starname']:
                m.addConstr(Yrs[name,s] <= nonqueueslot, 'enforceNonQueueSlots_' + str(name) + "_" + str(s) + "s")
    else:
        print("No non-queue observations are scheduled.")

    endtime2 = time.time()
    print("Total Time to build constraints: ", np.round(endtime2-start_the_clock,3))

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Set Gorubi Objective and Run
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    print("Setting objective.")

    # Minimize Theta_n
    m.setObjective(gp.quicksum(theta[name] for name in requests_frame['Starname']), GRB.MINIMIZE)

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
    print("Total Time to finish solver: ", np.round(endtime3-start_the_clock,3))


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #               Retrieve data from solution
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    file = open(outputDirectory + "runReport.txt", "w")
    file.close()

    print("Building human readable schedule.")
    combined_semester_schedule = hf.buildHumanReadableSchedule(Yrs, twilightMap_toDate, requests_frame, nNightsInSemester, nSlotsInNight, AvailableSlotsInGivenNight, nSlotsInQuarter, all_dates_dict, current_day, allocation_map_2D, weathered_map, slotsNeededForExposure_dict, nonqueueMap)

    round = 'Round 1'
    rf.buildFullnessReport(allocation_map_1D, twilightMap_all, combined_semester_schedule, nSlotsInQuarter, nSlotsInSemester, requests_frame, outputDirectory, slotsize, round)
    np.savetxt(outputDirectory + 'raw_combined_semester_schedule_Round1.txt', combined_semester_schedule, delimiter=',', fmt="%s")

    print("Writing Report.")
    filename = open(outputDirectory + "runReport.txt", "a")
    theta_n_var = []
    counter = 0
    for v in theta.values():
        varname = v.VarName
        varval = v.X
        counter += varval
    print("Sum of Theta: " + str(counter))
    filename.write("Sum of Theta: " + str(counter) + "\n")

    if plot_results:
        # Only generate the plots for results of Round 1.
        print("Plotting results.")
        all_starmaps = {}
        for i in range(len(requests_frame)):
            if pastObs_Info != {}:
                starmap = rf.buildObservedMap_past(pastObs_Info[requests_frame['Starname'][i]][1], pastObs_Info[requests_frame['Starname'][i]][2], pastObs_Info[requests_frame['Starname'][i]][3], semesterTemplateFile, weatherDiff_1D)
            else:
                starmap = pd.read_csv(semesterTemplateFile)
                starmap_columns = starmap.columns
                if 'Observed' not in starmap_columns:
                    starmap['Observed'] = [False]*len(starmap)
                if 'N_obs' not in starmap_columns:
                    starmap['N_obs'] = [0]*len(starmap)
                unique_hstdates_observed = []
            starmap_updated = rf.buildObservedMap_future(requests_frame['Starname'][i], slotsNeededForExposure_dict[requests_frame['Starname'][i]], combined_semester_schedule, starmap, all_dates_dict[current_day], slotsNeededForExposure_dict, np.array(allocation_all), weatherDiff_1D, nSlotsInNight)
            all_starmaps[requests_frame['Starname'][i]] = starmap_updated
            rf.writeCadencePlotFile(requests_frame['Starname'][i], i, starmap, turnOffOnFile, requests_frame, outputDirectory, unique_hstdates_observed, current_day)

        sum_schedule = np.sum(allocation_all, axis=1)
        allocated_days = np.nonzero(sum_schedule)
        firstDayAllocated = all_dates_array[allocated_days[0][0]]
        lastDayAllocated = all_dates_array[allocated_days[0][-1]]
        notable_dates = [semester_start_date, firstDayAllocated, lastDayAllocated, semester_end_date]
        rf.buildCOF(outputDirectory, current_day, notable_dates, requests_frame, all_dates_dict, all_starmaps, False)

    endtime4 = time.time()
    print("Total Time to complete Round 1: " + str(np.round(endtime4-start_the_clock,3)))
    filename.write("Total Time to complete Round 1: " + str(np.round(endtime4-start_the_clock,3)) + "\n")
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
        m.addConstr(gp.quicksum(theta[name] for name in requests_frame['Starname']) <= first_stage_objval + eps)
        m.setObjective(gp.quicksum(slotsNeededForExposure_dict[name]*Yrs[name,s] for name in requests_frame['Starname'] for s in range(nSlotsInSemester)), GRB.MAXIMIZE)
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
        print("Total Time to finish round 2 solver: ", np.round(endtime5-start_the_clock,3))

        combined_semester_schedule = hf.buildHumanReadableSchedule(Yrs, twilightMap_toDate, requests_frame, nNightsInSemester, nSlotsInNight, AvailableSlotsInGivenNight, nSlotsInQuarter, all_dates_dict, current_day, allocation_map_2D, weathered_map, slotsNeededForExposure_dict, nonqueueMap)

        round = 'Round 2'
        rf.buildFullnessReport(allocation_map_1D, twilightMap_all, combined_semester_schedule, nSlotsInQuarter, nSlotsInSemester, requests_frame, outputDirectory, slotsize, round)
        np.savetxt(outputDirectory + 'raw_combined_semester_schedule_Round2.txt', combined_semester_schedule, delimiter=',', fmt="%s")

        scheduleR1 = np.loadtxt(outputDirectory + 'raw_combined_semester_schedule_Round1.txt', delimiter=',', dtype=str)
        scheduleR2 = np.loadtxt(outputDirectory + 'raw_combined_semester_schedule_Round2.txt', delimiter=',', dtype=str)
        R2_requests = hf.getGapFillerTargets(scheduleR1, scheduleR2, all_dates_dict[current_day])
        np.savetxt(outputDirectory + 'Round2_Requests.txt', R2_requests, delimiter=',', fmt="%s")

        filename = open(outputDirectory + "runReport.txt", "a")
        theta_n_var = []
        counter = 0
        for v in theta.values():
            varname = v.VarName
            varval = v.X
            counter += varval
        filename.write("Sum of Theta: " + str(counter) + "\n")

        endtime6 = time.time()
        filename.write("Total Time to complete Round 2: " + str(np.round(endtime5-start_the_clock,3)) + "\n")
        filename.write("\n")
        filename.write("\n")
        filename.close()
        print("Round 2 complete.")
    else:
        print("Not running Round 2. Duplicating Raw Schedule as dummy file.")
        np.savetxt(outputDirectory + 'raw_combined_semester_schedule_Round2.txt', combined_semester_schedule, delimiter=',', fmt="%s")
        R2_requests = []
        np.savetxt(outputDirectory + 'Round2_Requests.txt', R2_requests, delimiter=',', fmt="%s")

    print("All done, clear skies!")
