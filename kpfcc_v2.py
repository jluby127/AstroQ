import numpy as np
import matplotlib.pyplot as pt
import pandas as pd
import sys
import math
import time
import pickle
from collections import defaultdict

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

# set a few parameters and flags
# these should be the only things that need changing as we test CPU usage
runOptimalAllocation = False
if runOptimalAllocation:
    SolveTime = 1200
    runRound2 = False
else:
    SolveTime = 300
    runRound2 = False # false for now
# The slot size, in minutes. Also determines the access_map used.
STEP = 10

# set your output directory here
# I suggest specifying something so that it doesn't autosave
# to the same directory as the run files and crowds up the GitHub repo.
outputdir = "/Users/jack/Desktop/"
# this represents the day of the semester we are scheduling
# for now this parameter only sets the number of days left in the semester,
# and therefore the number of slots in the semester
# and therefore the size of the Yns, Wns, etc variables.
# For testing CPU usage, always set this to '2024-02-01' as this
# was the first day of the 2024A semester.
current_day = '2024-07-01'
semesterLetter = 'A'
semesterYear = current_day[:4]


start_all = time.time()

# -------------------------------------------------------------------
# -------------------------------------------------------------------
#               Set up logistics parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------
semester_start_date, semesterLength = hf.getSemesterInfo(semesterYear, semesterLetter)
all_dates_dict = hf.buildDayDateDictionary(semester_start_date, semesterLength)
dates_in_semester = list(all_dates_dict.keys())
nNightsInSemester = hf.currentDayTracker(current_day, all_dates_dict)

accessmapfile = "kpfcc_v2_files/2024A_AccessibilityMaps__" + str(STEP) + "minSlot_14Hr.pkl"

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

# -------------------------------------------------------------------
# -------------------------------------------------------------------
#               Read in files and prep targets
# -------------------------------------------------------------------
# -------------------------------------------------------------------
print("Reading in and prepping files")

# Read in and prep target info
all_targets_frame = pd.read_csv("kpfcc_v2_files/requests.csv")

# build dictionary where keys are target names and
# values are how many slots are needed to complete one exposure for that target
slotsNeededDict = {}
for n,row in all_targets_frame.iterrows():
    name = row['Starname']
    exptime = row['Exposure_Time']
    slotsneededperExposure = math.ceil(exptime/(STEP*60.))
    slotsNeededDict[name] = slotsneededperExposure
    #if slotsneededperExposure > 1:
    #    print(name, slotsneededperExposure)

# sub-divide target list into only those that are more than 1 unique night of observations
# later, single shot observations are scheduled in Round 3
cadenced_targets = all_targets_frame[(all_targets_frame['N_Unique_Nights_Per_Semester'] > 1)]
cadenced_targets.reset_index(inplace=True)
ntargets = len(cadenced_targets)

# Read in twilight times
twilight_frame = pd.read_csv("kpfcc_v2_files/twilight_times.csv", parse_dates=True, index_col=0)

# AvailableSlotsInGivenNight is a 1D matrix of length nights n
# This is not a gorubi variable, but a regular python variable
# Each element will hold an integer which represents the true number of STEP minute slots
# in a quarter in a given night, after accounting for non-observable times due to day/twilight.
AvailableSlotsInGivenNight = []
for date in dates_in_semester:
    edge = tw.determineTwilightEdges(date, twilight_frame, verbose=False)
    AvailableSlotsInGivenNight.append(edge)

# set different constraints depending on if we are running the
# optimal allocation algorithm or a regular semester schedule solver
if runOptimalAllocation == False:
    print("Solve the regular semester schedule.")
    # Read in allocation info
    # Then a simple conversion from human-readable to computer-readable
    # I wanted to keep the data in the Binary Schedule for easy manual editing if needed
    allocation_schedule_load = np.loadtxt("kpfcc_v2_files/2024A_Binary_Schedule.txt", dtype=str)
    allocation_schedule_full = []
    for a in range(all_dates_dict[current_day], len(allocation_schedule_load)):
        convert = list(map(int, allocation_schedule_load[a][2:]))
        allocation_schedule_full.append(convert)

    print("Sampling out weather losses")
    # Randomly sample out 30% of allocated quarters to simulate weather loss
    allocation_schedule_long = np.array(allocation_schedule_full).flatten()
    allindices = [i for i in range(len(allocation_schedule_long)) if allocation_schedule_long[i] == 1]
    weatherlosses = np.random.choice(allindices, int(0.3*len(allindices)), replace=False)
    for w in weatherlosses:
        allocation_schedule_long[w] = 0
    allocation_schedule = np.reshape(allocation_schedule_long, (nNightsInSemester, 4))
    weatherDiff = allocation_schedule_full - allocation_schedule
    allocation_map, allocation_map_NS, weathered_map = hf.buildNonAllocatedMap(allocation_schedule, weatherDiff, AvailableSlotsInGivenNight, nSlotsInSemester, nNightsInSemester, nQuartersInNight, nSlotsInQuarter, nSlotsInNight)


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

if runOptimalAllocation == True:
    # Anq is a 2D matrix of N_nights_in_semester by N_quarters_in_night
    # element will be 1 if that night/quarter is allocated to KPF and 0 otherwise
    Anq = m.addVars(semester_grid, quarters_grid, vtype = GRB.BINARY, name = 'Allocation_Map')

    # Un is a 1D matrix of N_nights_in_semester
    # element will be 1 if at least one quarter in that night is allocated
    Un = m.addVars(semester_grid, vtype = GRB.BINARY, name = 'UniqueAllocation_Map')
    m.update()


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
            m.addConstr(Yns[name, s] <= Wnd[name, d], "related_Yns_Wnd_" + str(name) + "_" + str(d) + "_" + str(s))

print("Constraint: no other exposures can start in multi-slot exposures")
# Enforce no other exposures start within t slots if Yns[name, slot] is 1
# Velibor's new logic
for n,row in all_targets_frame.iterrows():
    name = row['Starname']
    slotsneededperExposure = slotsNeededDict[name]
    if (slotsneededperExposure > 1):
        for s in range(nSlotsInSemester):
            m.addConstr((slotsneededperExposure-1)*Yns[name,s] <= ( (slotsneededperExposure - 1) - gp.quicksum(Yns[name2,s+t] for name2 in all_targets_frame['Starname'] for t in range(1, min(slotsneededperExposure, nSlotsInSemester - s)))), 'dont_Start_Exposure_' + str(name) + "_" + str(s) + "s")

print("Constraint: only 1 exposures per slot")
# Enforce only one exposure per slot
for s in range(nSlotsInSemester):
    m.addConstr(gp.quicksum(Yns[name,s] for name in all_targets_frame['Starname']) <= 1, 'ensure_singleOb_perSlot_' + str(s) + "s")

if runOptimalAllocation == False:
    print("Constraint: enforce allocation map")
    # Enforce zeros when we are not allocated the telescope
    # Enforce zeros when target cannot be completed before an allocation ends
    for s in range(nSlotsInSemester):
        for k,row in all_targets_frame.iterrows():
            name = row['Starname']
            exptime = row['Exposure_Time']
            slotsneededperExposure = math.ceil(exptime/(STEP*60.))

            m.addConstr(Yns[name,s] <= allocation_map[s], 'allocationMatch_' + name + "_" + str(s) + "s")

            if slotsneededperExposure > 1:
                if np.sum(allocation_map[s:s+slotsneededperExposure]) != slotsneededperExposure:
                    m.addConstr(Yns[name,s] == 0, 'cantComplete' + name + "_" + str(s) + "s")

print("Constraint: enforce twilight times")
# Enforce twilight times within each quarter
for n in range(nNightsInSemester):
    for q in range(nQuartersInNight):
        start = n*nSlotsInQuarter*nQuartersInNight + q*nSlotsInQuarter
        end  = start + nSlotsInQuarter
        m.addConstr(gp.quicksum(Yns[name,s] for s in range(start,end) for name in all_targets_frame['Starname']) <= int(AvailableSlotsInGivenNight[n]/nQuartersInNight), "ensure_max_slots_used_in_quarter_" + str(n) + "n_" + str(q) + "q")

print("Constraint: exposure can't start if it won't complete in the night")
# Enforce that an exposure cannot start if it will not complete within the same night
for d in range(nNightsInSemester):
    end_night_slot = d*nSlotsInNight + nSlotsInNight - int(AvailableSlotsInGivenNight[d]/nQuartersInNight) - 1 # the -1 is to account for python indexing
    for t,row in all_targets_frame.iterrows():
        name = row['Starname']
        exptime = row['Exposure_Time']
        slotsneededperExposure = math.ceil(exptime/(STEP*60.))
        if slotsneededperExposure > 1:
            for e in range(0, slotsneededperExposure-1): # the -1 is so because a target can be started if it just barely fits
                m.addConstr(Yns[name,end_night_slot-e] == 0, 'dont_start_near_endof_night_' + name + "_" + str(end_night_slot) + 's_' + str(e) + 'e')

print("Constraint: inter-night cadence")
# Constrain the inter-night cadence
counter = 0
for d in range(nNightsInSemester):
    for t,row in all_targets_frame.iterrows():
        name = row['Starname']
        internightcadence = int(row['Inter_Night_Cadence'])
        for dlta in range(1, internightcadence):
            try:
                m.addConstr(Wnd[name,d] <= 1 - Wnd[name,d+dlta], 'enforce_internightCadence_' + str(name) + "_" + str(d) + "_" + str(dlta) + "dlta")
            except:
                # enter this except when the internightcadence brings us beyond the end of the semester (i think)
                counter += 1
                continue

print("Constraint: only observe if accessible")
# Ensure that a target cannot be scheduled at a time when it is not accessible
# Read in serialized pre-computed accessibility maps and deserialize it
with open(accessmapfile, 'rb') as file:
    accessmaps_precompute = pickle.load(file)
for name in all_targets_frame['Starname']:
    if name[:5] == 'Stand':
        target_access_map = np.array([1]*len(np.array(accessmaps_precompute['KOI-8107']).flatten()))
    else:
        target_access_map = np.array(accessmaps_precompute[name]).flatten()
    for s in range(startingSlot, nSlotsInSemester):
        m.addConstr((Yns[name,s] <= target_access_map[s]), 'constr_accessibility_' + str(name) + "_" + str(s))

print("Constraint: build Theta variable")
# Construct the Theta_n variable
extra = 10
for t,row in all_targets_frame.iterrows():
    name = row['Starname']
    Nobs_Unique = row['N_Unique_Nights_Per_Semester']
    m.addConstr(Thn[name] >= 0, 'ensureGreaterThanZero_' + str(name))
    # normalized theta value
    # m.addConstr(Thn[name] >= (Nobs_Unique - gp.quicksum(Yns[name, s] for s in range(nSlotsInSemester)))/Nobs_Unique, 'ensureGreaterThanNobsShortfall_' + str(name))
    # unnormalized theta value
    m.addConstr(Thn[name] >= (Nobs_Unique - gp.quicksum(Yns[name, s] for s in range(nSlotsInSemester))), 'ensureGreaterThanNobsShortfall_' + str(name))
    # add an upper limit to how many extra obs can be scheduled for a single target
    m.addConstr((gp.quicksum(Wnd[name,d] for d in range(nNightsInSemester)) <= Nobs_Unique + extra), 'max_Nobs_Unique_Nights_' + str(name))

if runOptimalAllocation == True:
    print("Solve the optimal allocation.")

    # Constraints on the way the allocation map can be filled
    maxQuarters = 150
    maxNights = 60
    minQuarterSelection = 5

    print("Constraint: setting max number of quarters allocated.")
    # No more than a maximum number of quarters can be allocated
    m.addConstr(gp.quicksum(Anq[d,q] for d in range(nNightsInSemester) for q in range(nQuartersInNight)) <= maxQuarters, "maximumQuartersAllocated")

    print("Constraint: relating allocation map and unique night allocation map.")
    # relate unique_allocation and allocation
    # if any one of the q in allocation[given date, q] is 1, then unique_allocation[given date] must be 1, zero otherwise
    for d in range(nNightsInSemester):
        for q in range(nQuartersInNight):
            m.addConstr(Un[d] >= Anq[d,q], "relatedUnique_andNonUnique_lowerbound_" + str(d) + "d_" + str(q) + "q")
        m.addConstr(Un[d] <= gp.quicksum(Anq[d,q] for q in range(nQuartersInNight)), "relatedUnique_andNonUnique_upperbound_" + str(d) + "d")

    print("Constraint: setting max number of unique nights allocated.")
    # No more than a maximum number of unique nights can be allocated
    m.addConstr(gp.quicksum(Un[d] for d in range(nNightsInSemester)) <= maxNights, "maximumNightsAllocated")

    print("Constraint: setting min number each quarter to be allocated.")
    # Minimum number of each quarter must be allocated
    m.addConstr(gp.quicksum(Anq[d,0] for d in range(nNightsInSemester)) >= minQuarterSelection, "minQuarterSelection_0q")
    m.addConstr(gp.quicksum(Anq[d,1] for d in range(nNightsInSemester)) >= minQuarterSelection, "minQuarterSelection_1q")
    m.addConstr(gp.quicksum(Anq[d,2] for d in range(nNightsInSemester)) >= minQuarterSelection, "minQuarterSelection_2q")
    m.addConstr(gp.quicksum(Anq[d,3] for d in range(nNightsInSemester)) >= minQuarterSelection, "minQuarterSelection_3q")

    print("Constraint: forbid certain patterns of quarter night allocations within night.")
    # Disallow certain patterns of quarters selected within same night
    for d in range(nNightsInSemester):
        # Cannot have 1st and 3rd quarter allocated without also allocating 2nd quarter (no gap), regardless of if 4th quarter is allocated or not
        m.addConstr(Anq[d,0] + (Un[d]-Anq[d,1]) + Anq[d,2] <= 2*Un[d], "NoGap2_" + str(d) + "d")
        # Cannot have 2nd and 4th quarter allocated without also allocating 3rd quarter (no gap), regardless of if 1st quarter is allocated or not
        m.addConstr(Anq[d,1] + (Un[d]-Anq[d,2]) + Anq[d,3] <= 2*Un[d], "NoGap3_" + str(d) + "d")
        # Cannot have only 2nd and 3rd quarters allocated (no middle half)
        m.addConstr((Un[d]-Anq[d,0]) + Anq[d,1] + Anq[d,2] + (Un[d]-Anq[d,3]) <= 3*Un[d], "NoMiddleHalf_" + str(d) + "d")
        # Cannot have only 1st and 4th quarters allocated (no end-cap half)
        m.addConstr(Anq[d,0] + (Un[d]-Anq[d,1]) + (Un[d]-Anq[d,2]) + Anq[d,3] <= 3*Un[d], "NoEndCapHalf_" + str(d) + "d")

    # Aesthetic constraints:
    # -----------------
    aesthetic = True
    if aesthetic:
        print("Yes include aesthetic constraints.")

        print("Constraint: setting max number of consecutive unique nights allocated.")
        # Don't allocate more than X consecutive nights
        consecMax = 6
        for d in range(nNightsInSemester - consecMax):
            m.addConstr(gp.quicksum(Un[d + t] for t in range(consecMax)) <= consecMax - 1, "consecutiveNightsMax_" + str(d) + "d")

        print("Constraint: setting min number of consecutive unique nights not allocated.")
        # Enforce at least one day allocated every X days (no large gaps)
        maxGap = 10
        for d in range(nNightsInSemester - maxGap):
            m.addConstr(gp.quicksum(Un[d + t] for t in range(maxGap)) >= 2, "noLargeGaps_" + str(d) + "d")

        print("Constraint: maximize the baseline of unique nights allocated.")
        # Enforce a night to be allocated within the first X nights and the last X nights of the semester (max baseline)
        maxBase = 5
        m.addConstr(gp.quicksum(Un[0 + t] for t in range(maxBase)) >= 1, "maxBase_early")
        m.addConstr(gp.quicksum(Un[nNightsInSemester - maxBase + t] for t in range(maxBase)) >= 1, "maxBase_late")

    else:
        print("No aesthetic constraints.")




endtime2 = time.time()
print("Total Time to build constraints: ", round(endtime2-start_all,3))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
#               Set Gorubi Objective and Run
# -------------------------------------------------------------------
# -------------------------------------------------------------------
print("Setting objective.")

# Minimize Theta_n
m.setObjective(gp.quicksum(Thn[name] for name in all_targets_frame['Starname']), GRB.MINIMIZE)

m.params.TimeLimit = SolveTime
m.params.MIPGap = 0.01 # can stop at 1% gap or better
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
print("Total Time to finish solver: ", round(endtime3-start_all,3))


# -------------------------------------------------------------------
# -------------------------------------------------------------------
#               Retrieve data from solution
# -------------------------------------------------------------------
# -------------------------------------------------------------------
file = open(outputdir + "runReport.txt", "w")
file.close()

twilightMap = hf.buildTwilightMap(AvailableSlotsInGivenNight, nSlotsInQuarter)

if runOptimalAllocation == True:
    allocation_schedule_1d = []
    for v in Anq.values():
        if np.round(v.X,0) == 1:
            allocation_schedule_1d.append(1)
        else:
            allocation_schedule_1d.append(0)
    allocation_schedule = np.reshape(allocation_schedule_1d, (nNightsInSemester, nQuartersInNight))
    allocation_schedule_full = allocation_schedule
    holder = np.zeros(np.shape(allocation_schedule))
    allocation_map, allocation_map_NS, weathered_map = hf.buildNonAllocatedMap(allocation_schedule, holder, AvailableSlotsInGivenNight, nSlotsInSemester, nNightsInSemester, nQuartersInNight, nSlotsInQuarter, nSlotsInNight)
    combined_semester_schedule = hf.buildHumanReadableSchedule(Yns, twilightMap, all_targets_frame, nNightsInSemester, nSlotsInNight, AvailableSlotsInGivenNight, nSlotsInQuarter, all_dates_dict, current_day, allocation_map_NS, weathered_map, slotsNeededDict)
else:
    allocation_schedule_full
    combined_semester_schedule = hf.buildHumanReadableSchedule(Yns, twilightMap, all_targets_frame, nNightsInSemester, nSlotsInNight, AvailableSlotsInGivenNight, nSlotsInQuarter, all_dates_dict, current_day, allocation_map_NS, weathered_map, slotsNeededDict)

round = 'Round 1'
rf.buildFullnessReport(allocation_schedule, twilightMap, combined_semester_schedule, nSlotsInQuarter, nSlotsInSemester, all_targets_frame, outputdir, STEP, round)

filename = open(outputdir + "runReport.txt", "a")
theta_n_var = []
counter = 0
for v in Thn.values():
    varname = v.VarName
    varval = v.X
    #print(varname, varval)
    counter += varval
filename.write("Sum of Theta: " + str(counter) + "\n")

# Only generate the COF for results of Round 1.
rf.buildCOF(outputdir, current_day, all_targets_frame, all_dates_dict, combined_semester_schedule, dates_in_semester)
rf.buildAllocationPicture(allocation_schedule_full, nNightsInSemester, nQuartersInNight, startingNight, outputdir)

endtime4 = time.time()
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

    m.params.TimeLimit = 900
    m.params.MIPGap = 0.25 # stop at 5% gap, this is good enough and leaves some room for standards
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
        print("Round 1 Model Solved.")

    endtime5 = time.time()
    print("Total Time to finish round 2 solver: ", np.round(endtime5-start_all,3))

    combined_semester_schedule = hf.buildHumanReadableSchedule(Yns, twilightMap, all_targets_frame, nNightsInSemester, nSlotsInNight, AvailableSlotsInGivenNight, nSlotsInQuarter, all_dates_dict, current_day, allocation_map_NS, weathered_map, slotsNeededDict)

    round = 'Round 2'
    rf.buildFullnessReport(allocation_schedule, twilightMap, combined_semester_schedule, nSlotsInQuarter, nSlotsInSemester, all_targets_frame, outputdir, STEP, round)

    filename = open(outputdir + "runReport.txt", "a")
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

print("done")
