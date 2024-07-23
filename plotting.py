import numpy as np
import pandas as pd
import astropy as ap
import astroplan as apl
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.time import TimeDelta
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import os
import star

def writeStarList(orderData, starttime, current_day, outputdir=''):
    """Write a text file with the optimized schedule.

    Args:
        orderData (dict): Formatted and ordered schedule from the
            model object
        starttime (astropy time object): Beginning of observing interval
        current_day (string): Date of observations, to be included in the
            filename
        outputdir (str): Folder in which to write the text output
    """

    f = open(os.path.join(outputdir,"ObserveOrder_" + current_day + ".txt"), "w")
    f.write("Target,StartExposure\n")
    for i in range(len(orderData['Starname'])):
        adjusted_timestamp = TimeDelta(orderData['Start Exposure'][i]*60,format='sec') + starttime
        formatted_timestamp = str(adjusted_timestamp)[11:16]
        f.write(orderData['Starname'][i] + "," + formatted_timestamp + "\n")
    f.close()


def nightPlan(orderData, current_day, outputdir='plots'):
    """Create an interactive plot which illustrates the solution.

    Args:
        orderData (dict): Formatted and ordered schedule from the
            model object
        current_day (string): Date of observations, to be included in the
            filename
        outputdir (str): Folder in which to save the html link to the plot
    """

    fig = px.scatter(orderData, x='Minutes the from Start of the Night', y="Starname", hover_data=['First Available', 'Last Available', 'Exposure Time (min)', "N_shots", "Total Exp Time (min)"] ,title='Night Plan') #color='Program'

    new_already_processed = []
    for i in range(len(orderData['Starname'])):
        if orderData['Starname'][i] not in new_already_processed:
            counter1 = list(orderData['Starname']).count(orderData['Starname'][i])
            for j in range(counter1):
                if i == 0 and j == 0:
                    fig.add_shape(type="rect", x0=orderData['Start Exposure'][i+j], x1=orderData['Start Exposure'][i+j] + orderData["Total Exp Time (min)"][i+j], y0=i-0.5, y1=i+0.5, fillcolor='red', showlegend=True, name='Expose')
                else:
                    fig.add_shape(type="rect", x0=orderData['Start Exposure'][i+j], x1=orderData['Start Exposure'][i+j] + orderData["Total Exp Time (min)"][i+j], y0=i-0.5, y1=i+0.5, fillcolor='red')
            if i == 0:
                fig.add_shape(type="rect", x0=orderData['First Available'][i], x1=orderData['Last Available'][i], y0=i-0.5, y1=i+0.5, fillcolor='lime', opacity=0.3, showlegend=True, name='Accessible')
            else:
                fig.add_shape(type="rect", x0=orderData['First Available'][i], x1=orderData['Last Available'][i], y0=i-0.5, y1=i+0.5, fillcolor='lime', opacity=0.3, showlegend=False)

            new_already_processed.append(orderData['Starname'][i])

    fig.write_html(outputdir + "/NightPlan_" + str(current_day) + ".html")
    # fig.show()


# Old imported code from kpfautoscheduler repo
def plot_path_2D(model,outputdir='plots'):
    """Plot the telescope path (alt/az) over time

    Args:
        model (object): A solved TTP model object
        outputdir (str): Folder in which to save the plot
    """

    names = [s for s in model.schedule['Starname']]
    times = model.times
    az_path = model.az_path
    alt_path = model.alt_path

    # convert obs_time from JD to UTC
    obs_time = [t.jd for t in times]

    start = Time(obs_time[0], format='jd')
    end = Time(obs_time[-1], format='jd')
    step = TimeDelta(30*60.,format='sec')
    ttemp = np.arange(start.jd, end.jd, step.jd)
    new_times = Time(ttemp,format='jd')
    new_obs_time = []
    for t in range(len(ttemp)):
        new_tt = str(new_times[t].isot)[11:16]
        new_obs_time.append(new_tt)

    fig, axs = plt.subplots(2,sharex=True,sharey=False,figsize = (20,8))
    fig.patch.set_alpha(1)
    axs[0].plot(obs_time,az_path,color = 'indigo')
    axs[0].set_xticks(ttemp, new_obs_time)
    axs[0].vlines(obs_time,0,360,linestyle = 'dashed', alpha = .5, color = 'gray')

    # Mark the given wrap angle
    wrap = model.observatory.wrapLimitAngle
    if wrap:
        axs[0].hlines(wrap,obs_time[0],obs_time[-1],linestyle='dashed', alpha=0.5, color='red')
        axs[0].text(obs_time[-1],wrap,f'Wrap={wrap}',color='red', fontsize=8, va='center')
    axs[0].set_yticks([0,120,240,360],[0,120,240,360])
    ax2 = axs[0].twiny()
    ax2.set_xlim(axs[0].get_xlim())

    topticks = []
    index = 0
    while index < len(obs_time):
        val = (obs_time[index+1]+obs_time[index])/2
        topticks.append(val)
        index+=2

    # Overlay the target names on the top axis
    ax2.set_xticks(topticks)
    ax2.set_xticklabels(names,rotation=45)
    axs[1].plot(obs_time,alt_path,color = 'seagreen')
    axs[1].vlines(obs_time,0,90,linestyle = 'dashed', alpha = .5, color = 'gray')
    #axs[2].plot(obs_time,airmass_path,color = 'firebrick')
    #axs[2].vlines(obs_time,1,3.5,linestyle = 'dashed', alpha = .5, color = 'gray')

    bottomticks = []
    bottomticks.append(obs_time[0])
    for i in range(1,4):
        val = obs_time[0] + i*(obs_time[-1]-obs_time[0])/4
        bottomticks.append(val)
    bottomticks.append(obs_time[-1])

    i = 0
    while i < len(az_path):
        axs[0].fill_betweenx([0,360],obs_time[i],obs_time[i+1],color = 'orange',alpha = .25)
        axs[1].fill_betweenx([0,90],obs_time[i],obs_time[i+1],color = 'orange',alpha = .25)
        #axs[2].fill_betweenx([1,3.5],obs_time[i],obs_time[i+1],color = 'orange',alpha = .25)
        i += 2

    axs[0].set(ylabel='Azimuth Angle (Deg)')
    axs[1].set(ylabel='Elevation Angle (Deg)')
    #axs[2].set(ylabel='Airmass')
    axs[1].set(xlabel='Observation Time (JD)')
    plt.title('Telescope Path Over Time')
    if outputdir != '':
        plt.savefig(os.path.join(outputdir,'Telescope_Path.png'))
    else:
        plt.show()
    plt.close()

def plot_slew_histogram(model,bins=30,outputdir='plots'):
    """Make histograms of the estimated/real slews in the schedule
    
    Args:
        model (object): A solved TTP model object
        bins (int): The number of bins to display in the histogram
        outputdir (str): Folder in which to save the plot
    """

    est_slews = model.estimated_slews
    real_slews = model.real_slews
    plt.style.use('ggplot')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6),sharey=True)

    # Left plot will have the worst case estimate of the slews which were used in the optimization
    ax1.hist(est_slews, bins=bins, color='skyblue', edgecolor='black', alpha=0.7)
    median1 = np.median(est_slews)
    ax1.axvline(median1, color='red', linestyle='dashed', linewidth=2)
    ax1.text(median1, ax1.get_ylim()[1] * 0.9, f'Median: {median1:.2f}', color='red', fontsize=12, ha='center', va='center', backgroundcolor='white')
    ax1.set_title(r'Worst Case Slew from $\tau$', fontsize=16, fontweight='bold')
    ax1.set_xlabel('Slew Duration (Minutes)', fontsize=14)
    ax1.set_ylabel('Number', fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Right plot will have the real slews (hopefully lower than before) when evaluated at the exact time of the slew
    ax2.hist(real_slews, bins=bins, color='skyblue', edgecolor='black', alpha=0.7)
    median2 = np.median(real_slews)
    ax2.axvline(median2, color='red', linestyle='dashed', linewidth=2)
    ax2.text(median2, ax2.get_ylim()[1] * 0.9, f'Median: {median2:.2f}', color='red', fontsize=12, ha='center', va='center', backgroundcolor='white')
    ax2.set_title(r'Real Slew at Time $t_{ijm}$', fontsize=16, fontweight='bold')
    ax2.set_xlabel('Slew Duration (Minutes)', fontsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    filename = os.path.join(outputdir,'SlewStats.png')
    plt.savefig(filename,dpi=200)