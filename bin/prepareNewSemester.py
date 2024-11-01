import numpy as np
import matplotlib.pyplot as pt
import pandas as pd
import astropy as apy
import astroplan as apl
from astropy.time import Time
from astropy.time import TimeDelta
import os
# import pickle
import json
import warnings
warnings.filterwarnings('ignore')
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__))[:-3] + "kpfcc/")
import mappingFunctions as mf

import argparse
parser = argparse.ArgumentParser(description='Generate schedules with KPF-CC v2')
parser.add_argument('-f','--folder', help='Folder to save generated scripts and plots', default=os.environ["KPFCC_SAVE_PATH"])
args = parser.parse_args()

# This file produces all the meta-data files needed to run the autoscheduler.
# This file should be run as soon as the allocation schedule is announced by the observatory.

# When working with the observatory scheduler to produce the allocation schedule, if you want to test a specific
# allocation schedule (as opposed to solving the Optimal Instrument Schedule where the schedule is itself the data product),
# then you will need to run this file but use the flag to turn off the allocation

# To prepare for a new semester, 3 files are needed:
# --- The allocation schedule that Keck Observatory produces each semester, see here and filter by KPF: https://www2.keck.hawaii.edu/observing/keckSchedule/queryForm.php
# --- The set of PI requests for cadenced queue observations, downloaded from the Request Submission Webform
# --- The set of PI requests for time-sensitive non-queue observations, downloaded from the Request Submission Webform which include start/stop times for the event windows

requests = pd.read_csv(args.folder + 'inputs/Requests.csv')
nonqueues = pd.read_csv(args.folder + 'inputs/NonQueue.csv', comment='#')
allocation = mf.reformatKeckAllocationData(args.folder + 'inputs/AllocationSchedule.csv')


# Generate dictionary between calendar day and day of semester
# -----------------------------------------------------------------------------------------
print("Preparing meta data. ")
# These are the only values that should ever be manually changed.
# And even then, really only change the semester and start date.
semester = '2024B'
start_date = '2024-08-01'
end_date = '2025-01-31'
slotStartTimestamp = '17:30:00' #HST = 03:30 UTC
slotEndTimestamp = '07:30:00' #HST = 17:30 UTC
nQuartersInNight = 4
nHoursInNight = 14
stepsizes = [10] # list of slot sizes, in minutes, to compute accessibility maps and non-queue maps.

# Build a dictionary relating calendar day to day of semester.
all_dates = []
date_formal = Time(start_date,format='iso',scale='utc')
date = str(date_formal)[:10]
all_dates.append(date)
nNightsInSemester = 1
while True:
    date_formal += TimeDelta(1,format='jd')
    date = str(date_formal)[:10]
    all_dates.append(date)
    nNightsInSemester += 1
    if date==end_date:
        break
    assert nNightsInSemester < 1000

# Generate the binary map for allocations this semester
# -----------------------------------------------------------------------------------------
print("Generating binary map of allocated nights/quarters.")
allocationMap = []
allocationMap_ints = []
uniqueDays = 0
for j in range(len(all_dates)):
    datemask = allocation['Date'] == all_dates[j]
    oneNight = allocation[datemask]
    if np.sum(datemask) == 0:
        # for night that is not allocated
        map1 = "0 0 0 0"
        map2 = [0, 0, 0, 0]
    elif np.sum(datemask) == 1:
        # for night where only one program is allocated (regardless of length of allocation)
        oneNight.reset_index(inplace=True)
        map1 = mf.translator(oneNight['Start'][0], oneNight['Stop'][0])
        map2 = [int(map1[0]), int(map1[2]), int(map1[4]), int(map1[6])]
        uniqueDays += 1
    elif np.sum(datemask) >= 1:
        # for night where multiple programs are allocated (regardless of their lengths)
        oneNight.reset_index(inplace=True)
        last = len(oneNight)
        map1 = mf.translator(oneNight['Start'][0], oneNight['Stop'][last-1])
        map2 = [int(map1[0]), int(map1[2]), int(map1[4]), int(map1[6])]
        uniqueDays += 1
    else:
        print("We have a problem, error code 5.")
        map1 = "0 0 0 0"
        map2 = [0, 0, 0, 0]
    allocationMap.append(map1)
    allocationMap_ints.append(map2)
print("Total number of quarters allocated: ", np.sum(allocationMap_ints))
print("Total unique nights allocated: ", uniqueDays)

# Write the binary allocation map to file
filename = args.folder + "inputs/" + semester + '_Binary_Schedule.txt'
file = open(filename, 'w')
for a in range(len(allocationMap)):
    line = all_dates[a] + " : " + str(allocationMap[a])
    file.write(str(line) + "\n")
file.close()


# Produce and write the start and stop times of each night to file.
# For the TTP.
# -----------------------------------------------------------------------------------------
print("Generate the nightly start/stop times for observing.")
listdates = list(allocation['Date'])
processed_dates = []
starts = []
stops = []
for t in range(len(allocation)):
    date = allocation['Date'][t]
    if date in processed_dates:
        continue
    start = allocation['Time'][t][:5]
    datecount = listdates.count(date)
    if datecount > 1:
        stop = allocation['Time'][t+datecount-1][8:13]
    else:
        stop = allocation['Time'][t][8:13]
    processed_dates.append(date)
    starts.append(start)
    stops.append(stop)
allocation_frame = pd.DataFrame({'Date':processed_dates, 'Start':starts, 'Stop':stops})
allocation_frame.to_csv(args.folder + "inputs/" + semester + '_NightlyStartStopTimes.csv', index=False)


# Create the template file for the cadence plots including true weather map
# -----------------------------------------------------------------------------------------
print("Generate the cadence plot template file.")
dateslist = []
quarterlist = []
for d in range(len(all_dates)):
    for q in range(4):
        dateslist.append(all_dates[d])
        quarterlist.append(0.5+q)
falselist = [False]*len(dateslist)
template_frame = pd.DataFrame({'Date':dateslist, 'Quarter':quarterlist,'Allocated':falselist,'Weathered':falselist})
template_frame.to_csv(args.folder + "inputs/" + semester + "_cadenceTemplateFile.csv", index=False)


# Create the csv file of twilight times for the semester
# -----------------------------------------------------------------------------------------
print("Generate twilight times for each day of semester.")
keck = apl.Observer.at_site('W. M. Keck Observatory')
twilight_frame = pd.DataFrame({'time_utc':all_dates})
eighteen_deg_evening = []
twelve_deg_evening = []
six_deg_evening = []
eighteen_deg_morning = []
twelve_deg_morning = []
six_deg_morning = []
for day in twilight_frame['time_utc']:
    as_day = Time(day,format='iso',scale='utc')
    eighteen_deg_evening.append(keck.twilight_evening_astronomical(as_day,which='next'))
    twelve_deg_evening.append(keck.twilight_evening_nautical(as_day,which='next'))
    six_deg_evening.append(keck.twilight_evening_civil(as_day,which='next'))
    eighteen_deg_morning.append(keck.twilight_morning_astronomical(as_day,which='next'))
    twelve_deg_morning.append(keck.twilight_morning_nautical(as_day,which='next'))
    six_deg_morning.append(keck.twilight_morning_civil(as_day,which='next'))
twilight_frame['18_evening'] = eighteen_deg_evening
twilight_frame['12_evening'] = twelve_deg_evening
twilight_frame['6_evening'] = six_deg_evening
twilight_frame['18_morning'] = eighteen_deg_morning
twilight_frame['12_morning'] = twelve_deg_morning
twilight_frame['6_morning'] = six_deg_morning
twilight_frame.to_csv(args.folder + "inputs/" + semester + '_twilight_times.csv', index=False)


# Generate the accessibility maps for each target.
# Repeat for different slot sizes.
# Accessibility is True for a given slot if target position at that time passes two tests:
#       - target's altitude is above the telescope pointing limits (see: https://www2.keck.hawaii.edu/inst/common/TelLimits.html)
#       - target's RA/Dec are >30deg from Moon's RA/Dec that night
# Note: no check for day/twilight/night time is made. This is handled in the twilight map.

# Also generate the turn on/off dates for each quarter for each star.
# This is the first and last date of the semester where the target is accessible in a given quarter.
#      For plotting purposes later. It is easist to precompute these values.
# -----------------------------------------------------------------------------------------
print("Generate accessibilty map and target turn on/off date for each target.")
Starname = []
ondate_q1 = []
offdate_q1 = []
ondate_q2 = []
offdate_q2 = []
ondate_q3 = []
offdate_q3 = []
ondate_q4 = []
offdate_q4 = []
for s in range(len(stepsizes)):
    all_maps = {}
    print("Computing for " + str(stepsizes[s]) + " minute slots.")
    for n,row in requests.iterrows():
        print("Calculating access map for: ", row['Starname'])
        access_map, turnonoff = mf.singleTargetAccessible(row['Starname'], row['RA'], row['Dec'], start_date, nNightsInSemester, stepsizes[s], turnonoff=True)

        flat_access_map = np.array(access_map).flatten()
        #all_maps[row['Starname']] = access_map
        # because you can't json serialize an array/list or integers
        # must create a string "list". Later decode it back to a real array
        stringmap = '['
        for e in range(len(flat_access_map)):
            stringmap += str(flat_access_map[e]) + ','
        stringmap += ']'
        all_maps[row['Starname']] = stringmap
        if s == 0:
            Starname.append(row['Starname'])
            ondate_q1.append(all_dates[turnonoff[0][0]])
            offdate_q1.append(all_dates[turnonoff[0][1]])
            ondate_q2.append(all_dates[turnonoff[1][0]])
            offdate_q2.append(all_dates[turnonoff[1][1]])
            ondate_q3.append(all_dates[turnonoff[2][0]])
            offdate_q3.append(all_dates[turnonoff[2][1]])
            ondate_q4.append(all_dates[turnonoff[3][0]])
            offdate_q4.append(all_dates[turnonoff[3][1]])

    # Serialize the dictionary and write it to a file
    with open(args.folder + "inputs/" + semester + "_AccessMaps_" + str(stepsizes[s]) + "minSlots.txt", 'w') as convert_file:
         convert_file.write(json.dumps(all_maps))

    # Create and write out the Turn On/Off dataframe
    turnOnOff = pd.DataFrame({'Starname':Starname,
                              'Q1_on_date':ondate_q1, 'Q1_off_date':offdate_q1,
                              'Q2_on_date':ondate_q2, 'Q2_off_date':offdate_q2,
                              'Q3_on_date':ondate_q3, 'Q3_off_date':offdate_q3,
                              'Q4_on_date':ondate_q4, 'Q4_off_date':offdate_q4,})
    turnOnOff.to_csv(args.folder + "inputs/" + semester + "_turnOnOffDates.csv", index=False)


# Generate the array of slots which must be occupied by non-queue observations
# -----------------------------------------------------------------------------------------
print("Reserving slots for time-sensitive multi-hour observations.")
startRounds = []
endRounds = []
startSlots = []
endSlots = []
durations_m = []
durations_h = []
daynumbers = []
for s in range(len(stepsizes)):
    print("Computing for " + str(stepsizes[s]) + " minute slots.")
    nSlotsInQuarter = int(((nHoursInNight*60)/nQuartersInNight)/stepsizes[s])
    nSlotsInNight = nSlotsInQuarter*nQuartersInNight

    slot2time = []
    time_formal = Time(start_date + 'T' + slotStartTimestamp, format='isot',scale='utc')
    timestr = str(time_formal)[11:16]
    slot2time.append(timestr)
    for i in range(nSlotsInNight):
        time_formal += TimeDelta(6.9444444*10**(-4)*stepsizes[s],format='jd') # 6.9444444*10**(-4) days = 1 minute
        timestr = str(time_formal)[11:16]
        slot2time.append(timestr)

    fullsemester_NonQueueMap_str = np.array([['Supercalifragilisticexpialidocious'] * nSlotsInNight] * nNightsInSemester, dtype=str)
    # For reasons I don't understand, the maximum length of string you can put into a slot is equal to the length of the
    # string that initializes the slot, so I'm putting the longest word i can think of here.

    onedaymap_str_all = []
    for i in range(len(nonqueues)):
        starttime = mf.convertUTC2HST(nonqueues['Start'][i][11:16])
        endtime = mf.convertUTC2HST(nonqueues['Stop'][i][11:16])

        startRound = mf.roundToSlot(starttime, stepsizes[s])
        startRounds.append(startRound)
        endRound = mf.roundToSlot(endtime, stepsizes[s])
        endRounds.append(endRound)

        startSlot = slot2time.index(startRound)
        startSlots.append(startSlot)
        endSlot = slot2time.index(endRound)
        endSlots.append(endSlot)
        dur = (endSlot - startSlot)*stepsizes[s]
        durations_m.append(dur)
        durations_h.append(round(dur/60,2))

        # build map
        date = nonqueues['Start'][i][:10]
        dayNumberOfSemester = all_dates.index(date)
        if dayNumberOfSemester == 0:
            print("Error. UTC to HST conversion error. Window begins in previous semester!")
        daynumbers.append(dayNumberOfSemester-1) # -1 to account for UTC to HST time change
        onedaymap = np.ones(nSlotsInNight, dtype=int)
        for j in range(startSlot, endSlot):
            onedaymap[j] = 0

        onedaymap_str = []
        for t in range(len(onedaymap)):
            if onedaymap[t] == 0:
                onedaymap_str.append('RM___' + nonqueues['Starname'][i])
            else:
                onedaymap_str.append('')
        onedaymap_str_all.append(onedaymap_str)

    # post process to remove the supercalis
    counter = 0
    for x in range(len(fullsemester_NonQueueMap_str)):
        if x in list(daynumbers):
            # continue
            fullsemester_NonQueueMap_str[x] = onedaymap_str_all[counter]
            counter += 1

        else:
            holder = ['']*nSlotsInNight
            fullsemester_NonQueueMap_str[x] = holder

    # Write out the file
    np.savetxt(args.folder + "inputs/" + semester + '_NonQueueMap' + str(stepsizes[s]) + '.txt', fullsemester_NonQueueMap_str, delimiter=',', fmt="%s")

print("All files written to directory: " + args.folder + 'inputs/ complete.')
print("done.")
