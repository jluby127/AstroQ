import numpy as np
import pandas as pd
import astropy as ap
import astroplan as apl
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.time import TimeDelta
import plotly
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import os
import ttp.star as star
import imageio
import base64


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

    # reverse the order so that the plot flows from top to bottom with time
    orderData = pd.DataFrame.from_dict(orderData)
    orderData = orderData.iloc[::-1]
    orderData.reset_index(inplace=True)

    # Each priority gets a different color. Make sure that each priority is actually included here or the plot will break. Recall bigger numbers are higher priorities.
    colordict = {'10':'red',
                 '9':'tomato',
                 '8':'darkorange',
                 '9':'sandybrown',
                 '7':'gold',
                 '6':'olive',
                 '5':'green',
                 '4':'cyan',
                 '3':'darkviolet',
                 '2':'magenta',
                '1':'blue'}

    # build the outline of the plot, add dummy points that are not displyed within the x/y limits so as to fill in the legend
    fig = px.scatter(orderData, x='Minutes the from Start of the Night', y="Starname", hover_data=['First Available', 'Last Available', 'Exposure Time (min)', "N_shots", "Total Exp Time (min)"] ,title='Night Plan', width=800, height=1000) #color='Program'
    fig.add_shape(type="rect", x0=-100, x1=-80, y0=-0.5, y1=0.5, fillcolor='red', showlegend=True, name='Expose P1')
    fig.add_shape(type="rect", x0=-100, x1=-80, y0=-0.5, y1=0.5, fillcolor='blue', showlegend=True, name='Expose P3')
    fig.add_shape(type="rect", x0=-100, x1=-80, y0=-0.5, y1=0.5, fillcolor='lime', opacity=0.3, showlegend=True, name='Accessible')

    new_already_processed = []
    ifixer = 0 # for multi-visit targets, it throws off the one row per target plotting...this fixes it
    for i in range(len(orderData['Starname'])):
        if orderData['Starname'][i] not in new_already_processed:
            counter1 = list(orderData['Starname']).count(orderData['Starname'][i])
            # find all the times in the night when the star is being visited
            indices = [k for k in range(len(orderData['Starname'])) if orderData['Starname'][k] == orderData['Starname'][i]]
            for j in range(len(indices)):
                fig.add_shape(type="rect", x0=orderData['Start Exposure'][indices[j]], x1=orderData['Start Exposure'][indices[j]] + orderData["Total Exp Time (min)"][indices[j]], y0=i+ifixer-0.5, y1=i+ifixer+0.5, fillcolor=colordict[str(orderData['Priority'][indices[j]])])
                if j == 0:
                    # only do this once, otherwise the green bar gets discolored compared to other rows
                    fig.add_shape(type="rect", x0=orderData['First Available'][indices[j]], x1=orderData['Last Available'][indices[j]], y0=i+ifixer-0.5, y1=i+ifixer+0.5, fillcolor='lime', opacity=0.3, showlegend=False)
            new_already_processed.append(orderData['Starname'][i])
        else:
            # if we already did this star, it is a multi-visit star and we need to adjust the row counter for plotting purposes
            ifixer -= 1

    fig.update_layout(xaxis_range=[0,orderData['Start Exposure'][0] + orderData["Total Exp Time (min)"][0]])

    if outputdir != '':
        # Save as both a .html intereactive plot and a low-res png plot
        fig.write_html(outputdir + "NightPlan_" + str(current_day) + ".html")
        fig.write_image(outputdir + 'NightPlan_' + str(current_day) + ".png")
    return fig

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
    # else:
    #     plt.show()
    plt.close()
    return fig

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


def animate_telescope(model, outputdir, animationStep=120):
    '''
    Produce the animation slew path GIF

    model (object) - the model object returned out from the TTP's model.py
    outputdir (str) - the path to save the animation
    animationStep (int) - the time, in seconds, between animation still frames. Default to 120s.
    '''

    tel_ims_dir = os.path.join(outputdir,'tel_ims')
    if not os.path.isdir(tel_ims_dir):
        os.mkdir(tel_ims_dir)

    # Gather info on the slew path from the TTP solution
    names = [s for s in model.schedule['Starname']]
    times = model.times

    # set the timestamps at which to generate a frame for the animation
    t = np.arange(model.nightstarts.jd, model.nightends.jd, TimeDelta(animationStep,format='sec').jd)
    t = Time(t,format='jd')

    # we need the astropy target objects within the TTP's Star object and we need them in the order to be observed
    list_targets = []
    for n in range(len(model.schedule['Starname'])):
        for s in model.stars:
            if s.name == model.schedule['Starname'][n]:
                list_targets.append(s.target)

    # compute the alt/az of each star at each animation still
    AZ = model.observatory.observer.altaz(t, list_targets, grid_times_targets=True)
    alt = []
    az = []
    for m in range(len(AZ)):
        tmpalt = []
        tmpaz = []
        for n in range(len(AZ[m])):
            tmpalt.append(np.round(AZ[m,n].az.rad,2))
            tmpaz.append(90-np.round(AZ[m,n].alt.deg,2))
        alt.append(tmpalt)
        az.append(tmpaz)

    # Call the custom code to orient targets in the correct path
    # Compute the telescope's pointing at each frame
    stamps = [0]*len(t)
    slewPath = createTelSlewPath(stamps, model.schedule['Time'], list_targets)

    AZ1 = model.observatory.observer.altaz(t, slewPath, grid_times_targets=False)
    tel_az = np.round(AZ1.az.rad,2)
    tel_zen = 90 - np.round(AZ1.alt.deg,2)

    # Build the animation frame by frame
    theta = np.arange(model.observatory.deckAzLim1/180, model.observatory.deckAzLim2/180, 1./180)*np.pi
    allaround = np.arange(0., 360/180, 1./180)*np.pi
    plotlowlim = 80
    # theta = np.arange(5.3/180, 146.2/180, 1./180)*np.pi
    filenames = []
    for i in range(len(t)):
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.set_ylim(0,plotlowlim)
        ax.set_title(t.isot[i])
        ax.set_yticklabels([])
        # ax.fill_between(theta,56.7,80,color = 'red',alpha=.7)
        ax.fill_between(theta,90-model.observatory.deckAltLim,plotlowlim,color = 'red',alpha=.7)
        ax.fill_between(allaround,90-model.observatory.vigLim,plotlowlim,color = 'red',alpha=.7)
        ax.fill_between(allaround,90-model.observatory.zenLim,0,color = 'red',alpha=.7)
        ax.set_theta_zero_location('N')

        # if arrival time less than current frame time, star "was observed"
        wasObserved =  np.array(model.schedule['Time']) <= float(t[i].jd)
        observed_list = np.array(model.schedule['Starname'])[wasObserved]

        # place stars on the plot
        for j in range(len(model.schedule['Starname'])):
            if names[j] in observed_list:
                ax.scatter(alt[j][i],az[j][i],color='orange',marker='*')
            else:
                ax.scatter(alt[j][i],az[j][i],color='white',marker='*')

        # place lines connecting stars on plot
        ax.plot(tel_az[:i],tel_zen[:i],color='orange')

        ax.set_facecolor('black')

        # create file name and append it to a list
        filename = f'{i}.png'
        filenames.append(filename)

        # save frame
        plt.savefig(os.path.join(tel_ims_dir,filename),dpi=100)
        plt.close()

    # build gif
    with imageio.get_writer(os.path.join(outputdir,'Observing_Animation.gif'), mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(os.path.join(tel_ims_dir,filename))
            writer.append_data(image)

    with open(os.path.join(outputdir,'Observing_Animation.gif'), 'rb') as f:
        gif_data = f.read()
    gif_base64 = base64.b64encode(gif_data).decode('utf-8')
    html = f'<img src="data:image/gif;base64,{gif_base64}" alt="Observing Animation"/>'

    # Remove files
    for filename in set(filenames):
        os.remove(os.path.join(tel_ims_dir,filename))
    try:
        os.remove(tel_ims_dir)
    except:
        print('Cannot remove redundant tel_ims directory due to file permissions')
    print("Animation plot complete. Still frame files deleted.")
    return html


def createTelSlewPath(stamps, changes, pointings, animationStep=120):
    '''
    Correctly assign each frame of the animation to the telescope pointing at that time

    stamps (list of zeros) - the list where each element represents a frame of the animation. We manipulate and return this at the end.
    changes (list) - the times at which the telescope pointing changes (in order of the slew path)
    poitings (list) - the astropy target objects of for the stars to be observed, in order of the slew path
    animationStep (int) - the time, in seconds, between frames

    return
        stamps - now a list where element holds the pointing of the telescope (aka the star object) at that frame

    '''
    # determine how many minutes each frame of the animation represents
    minPerStep = int(animationStep/60)
    mins = int(60/minPerStep)

    # offset the timestamps of the observations to a zero point
    changes = (changes - changes[0])*24*mins
    for c in range(len(changes)):
        changes[c] = int(changes[c])

    # determine telescope pointing at each frame
    for i in range(len(changes)-1):
        for j in range(len(stamps)):
            if j >= changes[i] and j < changes[i+1]:
                stamps[j] = pointings[i]

    # Add edge cases of the first and last telescope pointing
    k = 0
    while stamps[k] == 0:
        stamps[k] = pointings[0]
        k += 1
    l = len(stamps)-1
    while stamps[l] == 0:
        stamps[l] = pointings[-1]
        l -= 1

    return stamps
