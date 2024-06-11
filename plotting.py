import numpy as np
import pandas as pd
import astropy as ap
import astroplan as apl
import astropy.units as u
from astropy.coordinates import SkyCoord
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import os

import star

def nightPlan(orderData):

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

    # fig.write_html(outputdir + "/CDF_" + str(current_day) + ".html")
    fig.show()


# Old imported code from kpfautoscheduler repo
def plot_path_2D(model,outputdir=None):
    names = [s.name for s in model.schedule['Star']]
    times = model.times
    az_path = model.az_path
    alt_path = model.alt_path

    obs_time = [t.jd for t in times]
    fig, axs = plt.subplots(2,sharex=True,sharey=False,figsize = (20,8))
    fig.patch.set_alpha(1)
    axs[0].plot(obs_time,az_path,color = 'indigo')
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
    axs[1].set(xlabel='Observation Time (JD)')
    plt.title('Telescope Path Over Time')
    if outputdir:
        plt.savefig(os.path.join(outputdir,'Telescope_Path'))
    else:
        plt.show()
    plt.close()
