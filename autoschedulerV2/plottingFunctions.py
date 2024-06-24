from collections import defaultdict
import matplotlib.pyplot as plt
import imageio
import os
import numpy as np
import astropy as apy
import astroplan as apl
from astropy.time import Time
from astropy.time import TimeDelta
import astropy.units as u
import pandas as pd
import math

def plot_path_2D(obs_time,az_path,alt_path,names,targ_list,outputdir,current_day):

    # convert obs_time from JD to UTC
    start = Time(obs_time[0], format='jd')
    end = Time(obs_time[-1], format='jd')
    step = TimeDelta(30*60.,format='sec')
    ttemp = np.arange(start.jd, end.jd, step.jd)
    new_times = Time(ttemp,format='jd')
    new_obs_time = []
    for t in range(len(ttemp)):
        # print(new_times[t].isot)
        new_tt = str(new_times[t].isot)[11:16]
        new_obs_time.append(new_tt)

    fig, axs = plt.subplots(2,sharex=True,sharey=False,figsize = (14,8))
    fig.patch.set_alpha(1)
    axs[0].plot(obs_time,az_path,color = 'indigo')
    axs[0].set_xticks(ttemp, new_obs_time)
    axs[0].vlines(obs_time,0,360,linestyle = 'dashed', alpha = .5, color = 'gray')
    axs[0].set_yticks([0,120,240,360],[0,120,240,360])
    ax2 = axs[0].twiny()
    ax2.set_xlim(axs[0].get_xlim())

    topticks = []
    index = 0
    while index < len(obs_time):
        val = (obs_time[index+1]+obs_time[index])/2
        topticks.append(val)
        index+=2

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
    axs[1].set(xlabel='Observation Time (UTC)')
    plt.title('Telescope Path {}'.format(current_day))
    plt.savefig(os.path.join(outputdir,'Telescope_Path'))
    plt.close()

def animate_telescope(time_strings,total_azimuth_list,total_zenith_list,tel_az,tel_zen,observed_at_time,plotpath):

    theta = np.arange(5.3/180, 146.2/180, 1./180)*np.pi
    total_azimuth_list = np.array(total_azimuth_list)
    total_zenith_list = np.array(total_zenith_list)
    tel_ims_dir = os.path.join(plotpath,'tel_ims')
    if not os.path.isdir(tel_ims_dir):
        os.mkdir(tel_ims_dir)

    filenames = []
    for i in range(len(time_strings)):
        if i % 60 == 0:
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.set_ylim(0,70)
            ax.set_title(time_strings[i])
            ax.set_yticklabels([])
            ax.fill_between(theta,56.7,70,color = 'red',alpha=.7)
            #ax.fill_between(np.arange(0,2,1/180)*np.pi,50,70,color= 'red',alpha=.4)
            ax.set_theta_zero_location('N')
            observed_list = observed_at_time[:i]
            for j in set(observed_list):
                ax.scatter(total_azimuth_list[i][j],total_zenith_list[i][j],color='orange',marker='*')
            for j in set(observed_at_time):
                if j not in set(observed_list):
                    ax.scatter(total_azimuth_list[i][j],total_zenith_list[i][j],color='white',marker='*')
            ax.plot(tel_az[:i],tel_zen[:i],color='orange')
            ax.set_facecolor('black')

            # create file name and append it to a list
            filename = f'{i}.png'
            filenames.append(filename)

            # save frame
            plt.savefig(os.path.join(tel_ims_dir,filename),dpi=100)
            plt.close()

    # build gif
    with imageio.get_writer(os.path.join(plotpath,'Observing_Animation.gif'), mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(os.path.join(tel_ims_dir,filename))
            writer.append_data(image)

    # Remove files
    for filename in set(filenames):
        os.remove(os.path.join(tel_ims_dir,filename))
    try:
        os.remove(tel_ims_dir)
    except:
        print('Cannot remove redundant tel_ims directory due to file permissions')
