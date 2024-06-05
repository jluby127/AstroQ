import numpy as np
import pandas as pd
import astropy as ap
import astroplan as apl
import astropy.units as u
from astropy.coordinates import SkyCoord
import plotly.graph_objects as go
import plotly.express as px

import star

def getOrderedList(listOfStars, order):

    already_processed = []
    starname = []
    starearliest = []
    starlatest = []
    starendObs = []
    starexptime = []
    midexposure = []
    totalexptimes = []
    nExposures = []

    for i in range(len(order)):
        for j in range(len(listOfStars)):
            if order[i] == listOfStars[j].name and order[i] not in already_processed:
                counter = already_processed.count(order[i])
                for k in range(len(listOfStars[j].tend)):
                    starname.append(order[i])
                    starearliest.append(listOfStars[j].te - np.round(listOfStars[j].exptime/60,1))
                    starlatest.append(listOfStars[j].tl)
                    nExposures.append(listOfStars[j].shots)
                    totalexptime = (np.round(listOfStars[j].exptime/60,1)*listOfStars[j].shots + ((45/60)*(listOfStars[j].shots-1)))
                    totalexptimes.append(totalexptime)
                    starendObs.append(listOfStars[j].tend[k] - totalexptime)
                    midexposure.append(listOfStars[j].tend[k] - totalexptime/2)
                    starexptime.append(np.round(listOfStars[j].exptime/60,1))
                already_processed.append(order[i])

    orderData = pd.DataFrame({'Starname':starname, 'First Available From Start of Night (min)':starearliest, 'Last Available From Start of Night (min)':starlatest, 'Start Exposure From Start of Night (min)':starendObs, 'Minutes From Start Of Night':midexposure, 'Exposure Time (min)':starexptime, "Total Exp Time (min)":totalexptimes, "N_Shots":nExposures})

    fig = px.scatter(orderData, x="Minutes From Start Of Night", y="Starname", hover_data=['First Available From Start of Night (min)', 'Last Available From Start of Night (min)', 'Exposure Time (min)', "N_Shots"] ,title='Night Plan') #color='Program'

    new_already_processed = []
    for i in range(len(order)):
        if order[i] not in new_already_processed:
            counter1 = list(orderData['Starname']).count(orderData['Starname'][i])
            for j in range(counter1):
                if i == 0 and j == 0:
                    fig.add_shape(type="rect", x0=orderData['Start Exposure From Start of Night (min)'][i+j], x1=orderData['Start Exposure From Start of Night (min)'][i+j] + orderData["Total Exp Time (min)"][i+j], y0=i-0.5, y1=i+0.5, fillcolor='red', showlegend=True, name='Expose')
                else:
                    fig.add_shape(type="rect", x0=orderData['Start Exposure From Start of Night (min)'][i+j], x1=orderData['Start Exposure From Start of Night (min)'][i+j] + orderData["Total Exp Time (min)"][i+j], y0=i-0.5, y1=i+0.5, fillcolor='red')
            if i == 0:
                fig.add_shape(type="rect", x0=orderData['First Available From Start of Night (min)'][i], x1=orderData['Last Available From Start of Night (min)'][i], y0=i-0.5, y1=i+0.5, fillcolor='lime', opacity=0.3, showlegend=True, name='Accessible')
            else:
                fig.add_shape(type="rect", x0=orderData['First Available From Start of Night (min)'][i], x1=orderData['Last Available From Start of Night (min)'][i], y0=i-0.5, y1=i+0.5, fillcolor='lime', opacity=0.3, showlegend=False)

            new_already_processed.append(orderData['Starname'][i])

    # fig.write_html(outputdir + "/CDF_" + str(current_day) + ".html")
    fig.show()
