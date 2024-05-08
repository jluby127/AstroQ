import numpy as np
import matplotlib.pyplot as pt
import pandas as pd
import sys
import math
import time
import pickle
import os
import sys
from collections import defaultdict
from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astroplan as apl
import astropy.units as u
import plotly.graph_objects as go
import plotly.express as px
sys.path.append("/Users/jack/Desktop/")
sys.path.append("/Users/jack/Documents/Github/optimalAllocation/")
import helperFunctions as hf
import twilightFunctions as tw


def buildFullnessReport(allocation_schedule, twilightMap, combined_semester_schedule, nSlotsInQuarter, nSlotsInSemester, all_targets_frame, outputdir, STEP, round):
    file = open(outputdir + "runReport.txt", "a")
    file.write("Stats for " + str(round) + "\n")
    file.write("------------------------------------------------------" + "\n")
    # Compute how full the semester is
    nights_on_sky = []
    twilight_slots_on_sky_nights = 0
    for b in range(len(allocation_schedule)):
        if np.sum(allocation_schedule[b]) > 0:
            nights_on_sky.append(1)
            twilight_slots_on_sky_nights += ((np.sum(twilightMap[b])/4)*np.sum(allocation_schedule[b]))
        else:
            nights_on_sky.append(0)

    listnames = list(all_targets_frame['Starname'])
    slot_used_counter = 0
    for c in range(len(combined_semester_schedule)):
        for d in range(len(combined_semester_schedule[c])):
            if combined_semester_schedule[c][d] in listnames:
            #if combined_semester_schedule[c][d] != '' and combined_semester_schedule[c][d] != 'X' and combined_semester_schedule[c][d] != 'X*' and combined_semester_schedule[c][d] != '*' and combined_semester_schedule[n][s] != 'XW' and combined_semester_schedule[n][s] != 'X*W':
                slot_used_counter += 1

    total_available_slots = np.sum(np.array(allocation_schedule).flatten())*nSlotsInQuarter - twilight_slots_on_sky_nights
    file.write("N slots in semester:" + str(nSlotsInSemester) + "\n")
    file.write("N available slots: " + str(int(total_available_slots)) + "\n")
    file.write("N slots scheduled: " + str(slot_used_counter) + "\n")

    totalslotsrequested = 0
    for i in range(len(all_targets_frame)):
        totalslotsrequested += all_targets_frame['N_Unique_Nights_Per_Semester'][i]*math.ceil(all_targets_frame['Exposure_Time'][i]/(STEP*60.))
    file.write("N slots requested (total): " + str(totalslotsrequested) + "\n")
    percentage = np.round((slot_used_counter*100)/total_available_slots,3) # round((slot_used_counter*100)/total_available_slots,3)
    file.write("Percent full: " + str(percentage) + "%." + "\n")
    file.close()

    ff = open(outputdir + "semester_schedule.txt", "w")
    for ind in range(len(combined_semester_schedule)):
        ff.write("This is day " + str(ind) + "\n")
        ff.write("Q1: " + str(combined_semester_schedule[ind][nSlotsInQuarter*0:nSlotsInQuarter*1]) + "\n")
        ff.write("Q2: " + str(combined_semester_schedule[ind][nSlotsInQuarter*1:nSlotsInQuarter*2]) + "\n")
        ff.write("Q3: " + str(combined_semester_schedule[ind][nSlotsInQuarter*2:nSlotsInQuarter*3]) + "\n")
        ff.write("Q4: " + str(combined_semester_schedule[ind][nSlotsInQuarter*3:nSlotsInQuarter*4]) + "\n")
        ff.write("--------------------------------------------------------------------------------------" + "\n")
    ff.close()

def buildCOF(outputdir, current_day, all_targets_frame, all_dates_dict, starmap_maps, allocation_map_NS_fullsemester):
    dates_in_semester = list(all_dates_dict.keys())
    x = []
    y = []
    prog = []
    totobs = []
    commentsfile = open(outputdir + "ProgramData.csv", 'w')

    for program in all_targets_frame['Program_Code'].unique():

        programMask = all_targets_frame['Program_Code'] == program
        programDict = all_targets_frame[programMask]
        programDict.reset_index(inplace=True)
        if program == 'S001':
            tot_obs = len(programDict)*extra
        else:
            tot_obs = 0
            for t in range(len(programDict)):
                tot_obs += programDict['N_Unique_Nights_Per_Semester'][t]*programDict['N_Visits_per_Night'][t]*programDict['N_Observations_per_Visit'][t]

        runningObsList = [0.]*len(starmap_maps[all_targets_frame['Starname'][0]])
        for day in range(len(runningObsList)):
            for targ in programDict['Starname']:
                runningObsList[day] += starmap_maps[targ]['N_obs'][day]

        newrunning = 0.
        #uniqueAllocatedDays = list(starmap_maps[all_targets_frame['Starname'][0]]['Dates'].unique())
        for e in range(len(runningObsList)):
            x.append(starmap_maps[all_targets_frame['Starname'][0]]['Date'][e])
            newrunning += runningObsList[e]
            y.append(round((newrunning/tot_obs)*100,2))
            prog.append(program)
            totobs.append(tot_obs)

        commentsfile.write('#' + str(program) + '_trueComplete:' + str(round(y[-1],2)) + '\n')

    # # now run through one more time with dummy date to build a 1-to-1 line
    # newrunning = 0
    # for e in range(len(dates_in_semester)):
    #     x.append(dates_in_semester[e])
    #     newrunning += 1
    #     y.append(round((newrunning/len(dates_in_semester))*100,2))
    #     prog.append('Even Burn Rate Line')
    #     totobs.append(len(dates_in_semester))

    programdata = pd.DataFrame({"Program":prog, "Date":x, "Percent Complete (Observations)":y, "Total Obs Requested":totobs})
    programdata.to_csv(commentsfile, index=False)

    # fig = px.line(programdata, x="Date", y="Percent Complete (Observations)", hover_data=['Total Obs Requested'],
    #             color='Program',title='Cumulative Observation Function - N_Obs')

    originalforecast = pd.read_csv('/Users/jack/Documents/KPF_CC/PlotlyTesting/feb2024B_forecast.csv', comment='#')
    Program = ['Even Burn Rate', 'Even Burn Rate', 'Even Burn Rate', 'Even Burn Rate']
    Date = ['2024-02-01', '2024-02-24', '2024-07-29', '2024-07-31']
    PercentComplete = [0, 0, 100, 100]
    TotalObsRequested = [100, 100, 100, 100]
    burnrateline = pd.DataFrame({"Program":Program,"Date":Date,"Percent Complete (Observations)":PercentComplete,"Total Obs Requested":TotalObsRequested})

    fig = px.line(burnrateline, x="Date", y="Percent Complete (Observations)", hover_data=['Total Obs Requested'],
                    color='Program',title='Cumulative Observation Function - N_Obs')

    # there are 30 colors here, hopefully this is enough for a given semester.
    # colors = ['red', 'blue', 'green', 'yellow', 'cyan', 'purple', 'darkorange', 'brown',
    #           'firebrick', 'cornflowerblue', 'lime', 'gold', 'magenta', 'gray',
    #            'pink', 'royalblue', 'forestgreen', 'khaki', 'paleturquoise', 'blueviolet', 'moccasin',
    #           'salmon', 'steelblue', 'yellowgreen', 'darkkhaki', 'aquamarine', 'orchid', 'sandybrown', 'peru', 'gray'
    #          ]
    # light colors don't look good. Only use same core colors and just repeat as needed
    colors = ['red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            'red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            'red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            'red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            'red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            'red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            ]

    progs = originalforecast['Program'].unique()
    for i in range(len(progs)):
        progmask = programdata['Program'] == progs[i]
        progmask_original = originalforecast['Program'] == progs[i]

        fig.add_trace(go.Scatter(x=originalforecast['Date'][progmask_original], y=originalforecast['Percent Complete (Observations)'][progmask_original], name=progs[i],
                             line=dict(color=colors[i], width=2, dash='dash')))

        fig.add_trace(go.Scatter(x=programdata['Date'][progmask], y=programdata['Percent Complete (Observations)'][progmask], name=progs[i],
                             line=dict(color=colors[i], width=2)))

    fig.add_vrect(
            x0=current_day,
            x1=current_day,
            annotation_text="Today",
            line_dash="dash",
            fillcolor=None,
            line_width=2,
            line_color='black',
            annotation_position="bottom left"
        )
    fig.write_html(outputdir + "/COF_Nobs_" + str(current_day) + ".html")



def buildAllocationPicture(allocation_schedule, nNightsInSemester, nQuartersInNight, startingNight, all_dates_dict, outputdir):
    dateslist = list(all_dates_dict.keys())
    ff = open(outputdir + "AllocationPicture_stats.txt", "w")

    fig = pt.figure(figsize=(12,5))
    q1s = 0
    q2s = 0
    q3s = 0
    q4s = 0
    count0 = 0
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0
    for j in range(len(allocation_schedule)):
        if allocation_schedule[j][0] == 1.:
            q1s += 1
            pt.axvline(startingNight + j, ymin=0., ymax=0.25, color='b')
            date_info = dateslist[startingNight + j] + " - q0"
            ff.write(date_info + "\n")

        if allocation_schedule[j][1] == 1.:
            q2s += 1
            pt.axvline(startingNight + j, ymin=0.25, ymax=0.5, color='b')
            date_info = dateslist[startingNight + j] + " - q1"
            ff.write(date_info + "\n")

        if allocation_schedule[j][2] == 1.:
            q3s += 1
            pt.axvline(startingNight + j, ymin=0.5, ymax=0.75, color='b')
            date_info = dateslist[startingNight + j] + " - q2"
            ff.write(date_info + "\n")

        if allocation_schedule[j][3] == 1.:
            q4s += 1
            pt.axvline(startingNight + j, ymin=0.75, ymax=1.0, color='b')
            date_info = dateslist[startingNight + j] + " - q3"
            ff.write(date_info + "\n")


        allocated_quarters = np.sum(allocation_schedule[j])
        if allocated_quarters == 0:
            count0 += 1
        if allocated_quarters == 1:
            count1 += 1
        if allocated_quarters == 2:
            count2 += 1
        if allocated_quarters == 3:
            count3 += 1
        if allocated_quarters == 4:
            count4 += 1

    size=15
    pt.xlabel('Day')
    pt.ylabel('Quarter')
    pt.xlim(-5, 188)
    pt.yticks(range(4), ['Q1', 'Q2', 'Q3', 'Q4'])
    pt.ylim(0,4)
    pt.axhline(1, color='k', linestyle='-')
    pt.axhline(2, color='k', linestyle='-')
    pt.axhline(3, color='k', linestyle='-')
    pt.axhline(4, color='k', linestyle='-')
    pt.savefig(outputdir + "AllocationPicture.png", dpi=300, bbox_inches='tight', facecolor='w')

    ff.write("\n")
    ff.write("\n")
    ff.write("\n")
    ff.write("\n")
    ff.write("There are " + str(q1s) + " first quarters." + "\n")
    ff.write("There are " + str(q2s) + " second quarters." + "\n")
    ff.write("There are " + str(q3s) + " third quarters." + "\n")
    ff.write("There are " + str(q4s) + " fourth quarters." + "\n")
    ff.write("\n")
    ff.write("There are " + str(count0) + " no quarter nights." + "\n")
    ff.write("There are " + str(count1) + " 1 quarter nights."+ "\n")
    ff.write("There are " + str(count2) + " 2 quarter nights."+ "\n")
    ff.write("There are " + str(count3) + " 3 quarter nights."+ "\n")
    ff.write("There are " + str(count4) + " 4 quarter nights."+ "\n")
    ff.write("\n")
    total_quarters = count1 + 2*count2 + 3*count3 + 4*count4
    total_nights = count1 + count2 + count3 + count4
    ff.write("Total quarters allocated: " + str(total_quarters) + "\n")
    ff.write("Total unique nights allocated: " + str(total_nights) + "\n")
    ff.close()


def buildBinaryAllocationMap(outputdir, allocation_schedule):
    # Build the allocation map for the auto-scheduler
    # Example: run this to create a map of the results of optimal allocaiton
    # that then can be used to run the semester scheduler algorithm
    filename = outputdir + 'testing_Binary_Schedule.txt'
    file = open(filename, 'w')
    for a in range(len(allocation_schedule)):
        line = all_dates[a] + " : " + str(allocation_schedule[a])[1:-1]
        file.write(str(line) + "\n")
    file.close()

def buildObservedMap_past(unique_hstdates_observed, quarterObserved, Nobs_on_date, starmap_template_filename):
    starmap = pd.read_csv(starmap_template_filename)
    observed = [False]*len(starmap)
    N_observed = [0]*len(starmap)
    for i in range(len(unique_hstdates_observed)):
        ind = list(starmap['Date']).index(unique_hstdates_observed[i])
        observed[ind + int(quarterObserved[i]-0.5)] = True
        N_observed[ind + int(quarterObserved[i]-0.5)] = Nobs_on_date[i]
    starmap['Observed'] = observed
    starmap['N_obs'] = N_observed
    return starmap

def quarterDeterminer(value, nSlotsPerNight):
    if value <= int(nSlotsPerNight)/4:
        quart = 0
    elif value <= 2*int(nSlotsPerNight)/4 and value > int(nSlotsPerNight)/4:
        quart = 1
    elif value <= 3*int(nSlotsPerNight)/4 and value > 2*int(nSlotsPerNight)/4:
        quart = 2
    elif value <= 4*int(nSlotsPerNight)/4 and value > 3*int(nSlotsPerNight)/4:
        quart = 3
    else:
        quart = 0
        print("Houston, we've had a problem. No valid quarter.")
    return quart

def buildObservedMap_future(targetname, slotsPerExposure, combined_semester_schedule, starmap, current_day_number):
    observed = [False]*len(starmap)
    N_observed = [0]*len(starmap)
    for i in range(len(combined_semester_schedule)):
        if targetname in combined_semester_schedule[i]:
            ind = list(combined_semester_schedule[i]).index(targetname)
            quart = quarterDeterminer(ind, 84)
            starmap['Observed'][i*4 + quart] = True
            starmap['N_obs'][i*4 + quart] = list(combined_semester_schedule[i]).count(targetname) #int(combined_semester_schedule[i].count(targetname)/slotsPerExposure)
    return starmap

def writeCadencePlotFile(targetname, target_counter, starmap, turnFile, all_targets_frame, outputdir, unique_hstdates_observed, current_day):
    turnOnOffs = pd.read_csv(turnFile)
    request_id = all_targets_frame.index[all_targets_frame['Starname']==str(targetname)][0]
    request_name = all_targets_frame.loc[request_id,'Starname']
    program_code = all_targets_frame.loc[request_id,'Program_Code']
    print("Plotting cadence for star [" + str(request_name) + "] in program [" + str(program_code) + "]...target #" + str(target_counter) + " of " + str(len(all_targets_frame)) + ".")

    n_obs_desired = all_targets_frame.loc[request_id,'Total_Observations_Requested']
    n_obs_taken = len(unique_hstdates_observed)
    n_obs_scheduled = np.sum(starmap['N_obs'] - n_obs_taken)
    cadence = all_targets_frame.loc[request_id,'Inter_Night_Cadence']
    turnIDnumber = turnOnOffs.index[turnOnOffs['Starname']==str(targetname)][0]
    turns = [[turnOnOffs['Q1_on_date'][turnIDnumber], turnOnOffs['Q1_off_date'][turnIDnumber]], [turnOnOffs['Q2_on_date'][turnIDnumber], turnOnOffs['Q2_off_date'][turnIDnumber]], [turnOnOffs['Q3_on_date'][turnIDnumber], turnOnOffs['Q3_off_date'][turnIDnumber]], [turnOnOffs['Q4_on_date'][turnIDnumber], turnOnOffs['Q4_off_date'][turnIDnumber]]]
    savepath = outputdir + "/Cadences/" + str(program_code) + "/"
    os.makedirs(savepath, exist_ok = True)
    commentsfile = open(savepath + targetname + "_Cadence_Interactive.csv", 'w')
    commentsfile.write('#starname:' + str(targetname) + '\n')
    commentsfile.write('#programcode:' + str(program_code) + '\n')
    commentsfile.write('#Nobs_scheduled:' + str(n_obs_scheduled) + '\n')
    commentsfile.write('#Nobs_desired:' + str(n_obs_desired) + '\n')
    commentsfile.write('#Nobs_taken:' + str(n_obs_taken) + '\n')
    commentsfile.write('#cadence:' + str(cadence) + '\n')
    commentsfile.write('#q1_start:' + str(turns[0][1]) + '\n')
    commentsfile.write('#q1_end:' + str(turns[0][0]) + '\n')
    commentsfile.write('#q2_start:' + str(turns[1][1]) + '\n')
    commentsfile.write('#q2_end:' + str(turns[1][0]) + '\n')
    commentsfile.write('#q3_start:' + str(turns[2][1]) + '\n')
    commentsfile.write('#q3_end:' + str(turns[2][0]) + '\n')
    commentsfile.write('#q4_start:' + str(turns[3][1]) + '\n')
    commentsfile.write('#q4_end:' + str(turns[3][0]) + '\n')
    commentsfile.write('#current_day:' + str(current_day) + '\n')
    starmap = starmap[starmap.Allocated == True]
    starmap.to_csv(commentsfile, index=False)
    commentsfile.close()
