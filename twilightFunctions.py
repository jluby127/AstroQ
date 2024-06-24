import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
import astroplan as apl


def createTwilightTimesFile(semester_start_date, semesterLength, outputdir):
    #Initialize the astroplan observer object
    keck = apl.Observer.at_site('W. M. Keck Observatory')

    #note from jack 04/17/2024 -- computing the end date is untested.....might need fixing!
    startdate = Time(semester_start_date, format='iso')
    startdate.jd
    semester_end_date = startdate + TimeDelta(semesterLength,format='jd')
    semester_end_date.iso

    #Create a range of dates to calculate (some before and some after semester in this case)
    twilight_frame = pd.date_range(semester_start_date,str(semester_end_date)[:10]).to_frame()
    twilight_frame = twilight_frame.rename(columns={0:'time_utc'})
    eighteen_deg_evening = []
    twelve_deg_evening = []
    six_deg_evening = []
    eighteen_deg_morning = []
    twelve_deg_morning = []
    six_deg_morning = []
    for day in twilight_frame.index.strftime(date_format='%Y-%m-%d').tolist():
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

    #Create simple csv with date and twilight evening/morning times
    twilight_frame.to_csv(outputdir + 'twilight_times.csv')


def generateTwilightTimes(startDate, endDate):
    # This function computes the evening and morning twilight times for each day in the semester, in units of BJD
    # Note: both input variables should be strings, of the form "2024-02-01"

    #Initialize the astroplan observer object
    keck = apl.Observer.at_site('W. M. Keck Observatory')

    #Create a range of dates to calculate (some before and some after semester in this case)
    twilight_frame = pd.date_range(startDate, endDate).to_frame()
    twilight_frame = twilight_frame.rename(columns={0:'time_utc'})
    eighteen_deg_evening = []
    twelve_deg_evening = []
    six_deg_evening = []
    eighteen_deg_morning = []
    twelve_deg_morning = []
    six_deg_morning = []
    for day in twilight_frame.index.strftime(date_format='%Y-%m-%d').tolist():
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

    #Create simple csv with date and twilight evening/morning times
    twilight_frame.to_csv('data/twilight_times.csv')
    return twilight_frame



def getAvailableSlotsPerNight(startDate, nNightsInSemester):
    # This function generates a 1D list of length equal to the number of nights in the semester where
    # each element is the integer number of slots that are available to the scheduler in that night

    # startDate is a string with format: '2024-02-01'
    # nNightsInSemester is an integer

    AvailableSlotsInGivenNight = []
    date_formal = Time(startDate,format='iso',scale='utc')
    date = str(date_formal)[:10]
    for d in range(nNightsInSemester):
        edge = determineTwilightEdges(date, twilight_frame, verbose=False)
        AvailableSlotsInGivenNight.append(edge)
        date_formal += TimeDelta(1,format='jd')
        date = str(date_formal)[:10]
    return AvailableSlotsInGivenNight



def determineTwilightEdges(date, twilight_frame, STEP, verbose=False):
    # Using the pre-computed twilight times, determine how many slots within each night are available to the scheduler.
    # This is because each "night" in this framework is 14 hours long, designed to be longer than the longest night of the year
    # so that we always have the same fixed matrix size, and we can always "remove" slots to match the true, known length of a given night
    # by constraining that certain slots in a single night are not available.
    #
    # Returns an integer representing how many slots within that night are available to the scheduler for placing targets on given night
    # Note that the slot size is given in units of minutes by the STEP variable

    # find date in frame
    ind = twilight_frame.index[twilight_frame['time_utc'] == date].tolist()

    # get 12deg twilight evening time
    # This is to ensure that the map of X min intervals of twilight times matches the same X min intervals in the Accessibilty map
    start_time_bjd = twilight_frame.loc[ind,'12_evening'].values[0]
    start_time_bjd2 = Time(start_time_bjd,format='jd')

    start = date + "T03:30:00"
    startime = Time(start)
    totaltime = 14. # hours

    slots = int((totaltime*60)/STEP)
    binedges = []
    binedges.append(str(startime))
    for i in range(slots):
        startime += TimeDelta(STEP/(24*60),format='jd')
        binedges.append(str(startime))

    for i in range(len(binedges)):
        if start_time_bjd2 < binedges[i]:
            edge = i
            break

    availableWindowsInNight = len(binedges) - 2*edge

    if verbose:
        print('total windows in night: ', availableWindowsInNight)
        print('windows in quarter: ', int(availableWindowsInNight/4))
        print("minutes in quarter: ", int(availableWindowsInNight/4)*STEP)
        print("hours in quarter: ", round(((availableWindowsInNight/4)*STEP)/60,3))
        print("hours in night: ", 4*round(((availableWindowsInNight/4)*STEP)/60,3))

    # the max number of windows to assign in each of the four slots for a given night
    return availableWindowsInNight



def reorderAccessibility(AccessibilityMap, AvailableSlotsInGivenNight, nSlotsInNight):
    # This function re-orders the slots within a night.
    # In reality, all slots that are constrained as unobservable due to twilight occur at the beginning and end of the night
    # However, we want to redistribute those twilight slots to have an even number within each of the 4 quarters of the night
    # The idea being that each of the 4 quarters must have the same number of available slots, ie the same length of time.
    #
    # This function returns the same information encoded in the given AccessibilityMap, but with the slots redistributed evenly within the 4 quarters of each night.

    reorderedMap = []
    for n in range(len(AccessibilityMap)-1): #this minus one is simply because my pregenerated AccessibilityMap(s) have an extra day by accident. delete this later

        edge = int((nSlotsInNight - int(AvailableSlotsInGivenNight[n]))/2)
        flag1 = (AvailableSlotsInGivenNight[n]/4)%1 != 0
        subset = AccessibilityMap[n][edge:-edge]
        nightwindowsPer_quarter = len(subset)/4
        flag2 = nightwindowsPer_quarter%1 != 0
        nightwindowsPer_quarter = int(nightwindowsPer_quarter)

        # lets do this manually for simplicity
        q1 = []
        for w in range(nightwindowsPer_quarter*0, nightwindowsPer_quarter*1):
            q1.append(subset[w])
        q2 = []
        for w in range(nightwindowsPer_quarter*1, nightwindowsPer_quarter*2):
            q2.append(subset[w])
        q3 = []
        for w in range(nightwindowsPer_quarter*2, nightwindowsPer_quarter*3):
            q3.append(subset[w])
        q4 = []
        for w in range(nightwindowsPer_quarter*3, nightwindowsPer_quarter*4):
            q4.append(subset[w])

        comblen = len(q1)+len(q2)+len(q3)+len(q4)
        bufflen = int((nSlotsInNight-comblen)/4)
        buffer0 = [0]*bufflen
        q1.extend(buffer0)
        q2.extend(buffer0)
        q3.extend(buffer0)
        q4.extend(buffer0)

        combine = q1 + q2 + q3 + q4
        reorderedMap.extend(combine)

    #reorderedMap = np.array(reorderedMap).flatten()
    return reorderedMap


def reorderAccessibilityOneDay(OneDayMap, AvailableSlotsInTheNight, nSlotsInNight):
    # This function re-orders the slots within a night.
    # In reality, all slots that are constrained as unobservable due to twilight occur at the beginning and end of the night
    # However, we want to redistribute those twilight slots to have an even number within each of the 4 quarters of the night
    # The idea being that each of the 4 quarters must have the same number of available slots, ie the same length of time.
    #
    # This function returns the same information encoded in the given AccessibilityMap, but with the slots redistributed evenly within the 4 quarters of each night.
    #
    # This version of the function operates on a single day at a time, specifically for the non-queue observations
    # because it fills slots that the NonQueue Obs YES needs as 0's and the slots that it DOESN'T need as 1's.

    edge = int((nSlotsInNight - int(AvailableSlotsInTheNight))/2)
    flag1 = (AvailableSlotsInTheNight/4)%1 != 0
    subset = OneDayMap[edge-1:-edge]
    nightwindowsPer_quarter = len(subset)/4
    flag2 = nightwindowsPer_quarter%1 != 0
    nightwindowsPer_quarter = int(nightwindowsPer_quarter)

    # lets do this manually for simplicity
    q1 = []
    for w in range(nightwindowsPer_quarter*0, nightwindowsPer_quarter*1):
        q1.append(subset[w])
    q2 = []
    for w in range(nightwindowsPer_quarter*1, nightwindowsPer_quarter*2):
        q2.append(subset[w])
    q3 = []
    for w in range(nightwindowsPer_quarter*2, nightwindowsPer_quarter*3):
        q3.append(subset[w])
    q4 = []
    for w in range(nightwindowsPer_quarter*3, nightwindowsPer_quarter*4):
        q4.append(subset[w])

    comblen = len(q1)+len(q2)+len(q3)+len(q4)
    bufflen = int((nSlotsInNight-comblen)/4)
    buffer0 = [1]*bufflen
    q1.extend(buffer0)
    q2.extend(buffer0)
    q3.extend(buffer0)
    q4.extend(buffer0)

    combine = q1 + q2 + q3 + q4

    return np.array(combine)
