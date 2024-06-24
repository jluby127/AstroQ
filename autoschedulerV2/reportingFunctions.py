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

    listnames = list(all_targets_frame['Starname'])
    unavailable = 0
    unused = 0
    used = 0
    for b in range(len(combined_semester_schedule)):
        before = ''
        for c in range(len(combined_semester_schedule[b])):
            if "X" in combined_semester_schedule[b][c] or "*" in combined_semester_schedule[b][c]:
                unavailable += 1
            if combined_semester_schedule[b][c] == "":
                unused += 1
            if combined_semester_schedule[b][c] in listnames:
                used += 1
    available = unused + used
    available2 = np.sum(allocation_schedule.flatten())
    file.write("N slots in semester:" + str(nSlotsInSemester) + "\n")
    file.write("N available slots: " + str(available) + "\n")
    file.write("N available slots (recalc):" + str(available2) + "\n")
    file.write("N slots scheduled: " + str(used) + "\n")
    file.write("N slots left empty: " + str(unused) + "\n")

    totalslotsrequested = 0
    for i in range(len(all_targets_frame)):
        totalslotsrequested += all_targets_frame['N_Unique_Nights_Per_Semester'][i]*math.ceil(all_targets_frame['Nominal_ExpTime'][i]/(STEP*60.))
    file.write("N slots requested (total): " + str(totalslotsrequested) + "\n")
    percentage = np.round((used*100)/available,3) # round((slot_used_counter*100)/total_available_slots,3)
    percentage2 = np.round((used*100)/available2,3) # round((slot_used_counter*100)/total_available_slots,3)
    file.write("Percent full: " + str(percentage) + "%." + "\n")
    file.write("Percent full (recalc): " + str(percentage2) + "%." + "\n")

    file.close()

    # ff = open(outputdir + "semester_schedule.txt", "w")
    # for ind in range(len(combined_semester_schedule)):
    #     ff.write("This is day " + str(ind) + "\n")
    #     ff.write("Q1: " + str(combined_semester_schedule[ind][nSlotsInQuarter*0:nSlotsInQuarter*1]) + "\n")
    #     ff.write("Q2: " + str(combined_semester_schedule[ind][nSlotsInQuarter*1:nSlotsInQuarter*2]) + "\n")
    #     ff.write("Q3: " + str(combined_semester_schedule[ind][nSlotsInQuarter*2:nSlotsInQuarter*3]) + "\n")
    #     ff.write("Q4: " + str(combined_semester_schedule[ind][nSlotsInQuarter*3:nSlotsInQuarter*4]) + "\n")
    #     ff.write("--------------------------------------------------------------------------------------" + "\n")
    # ff.close()

def buildCOF(outputdir, current_day, all_targets_frame, all_dates_dict, starmap_maps, notable_dates, compare_original):
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
            tot_obs = len(programDict)#*extra
        else:
            tot_obs = 0
            for t in range(len(programDict)):
                #tot_obs += programDict['N_Unique_Nights_Per_Semester'][t]*programDict['N_Visits_per_Night'][t]*programDict['N_Observations_per_Visit'][t]
                tot_obs += programDict['N_Unique_Nights_Per_Semester'][t]

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

    # notable_dates is a 4 element list containing the following information in the following order
    # element 1: start day of semester
    # element 2: first day in semester that is allocated
    # element 3: last day in semester that is allocated
    # element 4: last day of the semester
    BurnProg = ['Even Burn Rate']*4
    PercentComplete = [0, 0, 100, 100]
    TotalObsRequested = [100, 100, 100, 100]
    #notable_dates = ['2024-08-01', '2024-08-01', '2025-01-31', '2025-01-31']
    notable_dates = ['2024-02-01', '2024-02-24', '2024-07-28', '2024-07-31']
    for b in range(len(notable_dates)):
        prog.append(BurnProg[b])
        x.append(notable_dates[b])
        y.append(PercentComplete[b])
        totobs.append(TotalObsRequested[b])

    programdata = pd.DataFrame({"Program":prog, "Date":x, "Percent Complete (Observations)":y, "Total Obs Requested":totobs})
    programdata.to_csv(commentsfile, index=False)

    # light colors don't look good. Only use same core colors and just repeat as needed
    colors = ['red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            'red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            'red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            'red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            'red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            'red', 'blue', 'green', 'purple', 'darkorange', 'brown', 'gold', 'forestgreen', 'firebrick', 'royalblue',
            ]

    if compare_original:
        burnrate = pd.DataFrame({"Program":BurnProg, "Date":notable_dates, "Percent Complete (Observations)":PercentComplete, "Total Obs Requested":TotalObsRequested})
        fig = px.line(burnrate, x="Date", y="Percent Complete (Observations)", hover_data=['Total Obs Requested'],
                color='Program',title='Cumulative Observation Function - N_Obs')

        originalforecast = pd.read_csv('/Users/jack/Documents/KPF_CC/PlotlyTesting/feb2024B_forecast.csv', comment='#')
        progs = originalforecast['Program'].unique()
        for i in range(len(progs)):
            progmask = programdata['Program'] == progs[i]
            progmask_original = originalforecast['Program'] == progs[i]
            fig.add_trace(go.Scatter(x=originalforecast['Date'][progmask_original], y=originalforecast['Percent Complete (Observations)'][progmask_original], name=progs[i],
                                 line=dict(color=colors[i], width=2, dash='dash')))
            fig.add_trace(go.Scatter(x=programdata['Date'][progmask], y=programdata['Percent Complete (Observations)'][progmask], name=progs[i],
                                 line=dict(color=colors[i], width=2)))
    else:
        fig = px.line(programdata, x="Date", y="Percent Complete (Observations)", hover_data=['Total Obs Requested'],
                    color='Program',title='Cumulative Observation Function - N_Obs')

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
    # fig.write_html(outputdir + "/COF_Nobs_" + str(current_day) + ".html")
    fig.write_html(outputdir + "/COF_Plot_V-.html")



def buildAllocationPicture(allocation_schedule, nNightsInSemester, nQuartersInNight, startingNight, all_dates_dict, outputdir):
    dateslist = list(all_dates_dict.keys())
    ff = open(outputdir + "AllocationPicture_stats.txt", "w")
    fff = open(outputdir + "AllocationMap.txt", "w")
    fff.write('Date,Q1,Q2,Q3,Q4\n')

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
        holder = [0, 0, 0, 0]

        if allocation_schedule[j][0] == 1.:
            q1s += 1
            pt.axvline(startingNight + j, ymin=0., ymax=0.25, color='b')
            date_info = dateslist[startingNight + j] + " - q0"
            ff.write(date_info + "\n")
            holder[0] = 1

        if allocation_schedule[j][1] == 1.:
            q2s += 1
            pt.axvline(startingNight + j, ymin=0.25, ymax=0.5, color='b')
            date_info = dateslist[startingNight + j] + " - q1"
            ff.write(date_info + "\n")
            holder[1] = 1

        if allocation_schedule[j][2] == 1.:
            q3s += 1
            pt.axvline(startingNight + j, ymin=0.5, ymax=0.75, color='b')
            date_info = dateslist[startingNight + j] + " - q2"
            ff.write(date_info + "\n")
            holder[2] = 1

        if allocation_schedule[j][3] == 1.:
            q4s += 1
            pt.axvline(startingNight + j, ymin=0.75, ymax=1.0, color='b')
            date_info = dateslist[startingNight + j] + " - q3"
            ff.write(date_info + "\n")
            holder[3] = 1


        line = str(dateslist[startingNight + j]) + "," + str(holder[0]) + "," + str(holder[1]) + "," + str(holder[2]) + "," + str(holder[3]) + "\n"
        fff.write(line)

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
    fff.close()


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

def buildObservedMap_past(unique_hstdates_observed, quarterObserved, Nobs_on_date, starmap_template_filename, weathermap):
    starmap = pd.read_csv(starmap_template_filename)
    observed = [False]*len(starmap)
    N_observed = [0]*len(starmap)
    for i in range(len(unique_hstdates_observed)):
        ind = list(starmap['Date']).index(unique_hstdates_observed[i])
        observed[ind + int(quarterObserved[i]-0.5)] = True
        N_observed[ind + int(quarterObserved[i]-0.5)] = 1#Nobs_on_date[i]
    starmap['Observed'] = observed
    starmap['N_obs'] = N_observed
    return starmap

def quarterDeterminer(value, nSlotsInNight):
    if value <= int(nSlotsInNight*(1/4.)):
        quart = 0
    elif value <= int(nSlotsInNight*(2/4.)) and value > int(nSlotsInNight*(1/4.)):
        quart = 1
    elif value <= int(nSlotsInNight*(3/4.)) and value > int(nSlotsInNight*(2/4.)):
        quart = 2
    elif value <= nSlotsInNight and value > int(nSlotsInNight*(3/4.)):
        quart = 3
    elif value <= int(nSlotsInNight*(5/4.)) and value > nSlotsInNight:
        quart = 3
    else:
        quart = 0
        print("Houston, we've had a problem. No valid quarter: ", value, nSlotsInNight)
    return quart

def buildObservedMap_future(targetname, slotsPerExposure, combined_semester_schedule, starmap, current_day_number, slotsNeededDict, allocationMap, weatherDiff_1D, nSlotsInNight):
    observed = [False]*len(starmap)
    N_observed = [0]*len(starmap)
    for i in range(len(combined_semester_schedule)):
        if targetname in combined_semester_schedule[i]:
            exptimeslots = slotsNeededDict[targetname]
            ind = list(combined_semester_schedule[i]).index(targetname)
            quart = quarterDeterminer(ind, nSlotsInNight)
            starmap['Observed'][i*4 + quart] = True
            starmap['N_obs'][i*4 + quart] = list(combined_semester_schedule[i]).count(targetname)/exptimeslots #int(combined_semester_schedule[i].count(targetname)/slotsPerExposure)
    allocationMap_flat = allocationMap.flatten()
    alloBool = []
    for a in range(len(allocationMap_flat)):
        # if targetname == '10700':
        #     print(a, allocationMap_flat[a], bool(allocationMap_flat[a]))
        alloBool.append(bool(allocationMap_flat[a]))
    starmap['Allocated'] = alloBool
    for w in range(len(weatherDiff_1D)):
        if weatherDiff_1D[w] == 1:
            starmap['Weathered'][w] = True
    return starmap

def writeCadencePlotFile(targetname, target_counter, starmap, turnFile, all_targets_frame, outputdir, unique_hstdates_observed, current_day):

    turnOnOffs = pd.read_csv(turnFile)
    request_id = all_targets_frame.index[all_targets_frame['Starname']==str(targetname)][0]
    request_name = all_targets_frame.loc[request_id,'Starname']
    program_code = all_targets_frame.loc[request_id,'Program_Code']
    #print("Plotting cadence for star [" + str(request_name) + "] in program [" + str(program_code) + "]...target #" + str(target_counter) + " of " + str(len(all_targets_frame)) + ".")

    # n_obs_desired = all_targets_frame.loc[request_id,'Total_Observations_Requested']
    n_obs_desired = all_targets_frame.loc[request_id,'N_Unique_Nights_Per_Semester']
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


def makeTemplateFile(all_dates_dict, filename):
    # filename = "/Users/jack/Desktop/template_star_observed.csv"
    date = []
    quart = []
    for i in range(len(all_dates_dict)):
        for j in range(4):
            date.append(all_dates_dict[i])
            quart.append(j + 0.5)
    starmap = pd.DataFrame({"Date":date, "Quarter":quart})
    starmap['Weathered'] = [False]*len(date)
    starmap['Allocated'] = [False]*len(date)

    starmap.to_csv(filename, index=False)
    print("Template file saved to ", filename)
