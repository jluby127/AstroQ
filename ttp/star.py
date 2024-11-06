import numpy as np
import pandas as pd
import astropy as ap
import astroplan as apl
import astropy.units as u
from astropy.coordinates import SkyCoord

# from astroplan import download_IERS_A
# download_IERS_A()

class star(object):
    """Object which describes a star and the requested observation scheme

    Args:
        row (pandas series): A row from the input .csv file as a dataframe row
        index (int): The row identifier of the star in the input .csv file

    Attributes:
        name (string): Given name of the star
        ra (float): Right ascension of the star in degrees
        dec (float): Declination of the star in degrees
        target (astroplan FixedTarget): Object representation of the star in
            the astroplan module
        exptime (float): Duration of each requested exposure of the star
            in minutes
        shots (int): How many exposures must be collected each time the
            star is visited by the telescope
        expwithreadout (float): Total duration (minutes) of a visit of the star,
            after considering the number of exposures and the readout delay
            between each
        visits (int): Number of unique times the telescope should return to the
            star during the night
        intra_night_cadence (float): Minimum requested time (hours) between any
            two unique visits to the star
        index (int): The row identifier of the star in the input .csv file

        NOTE: The following are given as minutes from the beginning of the
            observing interval, and are only determined after solving a model

        te (int): Earliest time  which the star can finish an observation
        tl (int): Latest time which the star can finish an observation
        tend (list): The time(s) at which a unique visit of the star is
            scheduled to finish
        orderInNight (list): The place(s) in the schedule in which the star
            will be visited
    """

    def __init__(self, row, index):

        # must hard code in data types because, for whatever reason, if ALL characters in ALL starnames are numbers, then these values all get read in as floats 
        self.name = str(row['Starname'])
        # ignore these notes for now!!!!
        # NOTE: RA must be in following format: XXhYYmZZ.Zs
        # NOTE: Dec must be in following format: +XXdYYmZZ.Zs
        self.ra = float(row['RA'])
        self.dec = float(row['Dec'])
        coords = SkyCoord(self.ra*u.deg, self.dec*u.deg, frame='icrs')
        self.priority = int(row['Priority'])
        self.target = apl.FixedTarget(name=self.name, coord=coords)
        self.exptime = int(row['Exposure Time']/60) # Seconds to minutes
        self.shots = int(row['Exposures Per Visit'])
        self.visits = int(row['Visits In Night'])
        self.expwithreadout = int(self.exptime*self.shots + (45/60)*(self.shots-1))
        self.intra_night_cadence = int(row['Intra_Night_Cadence']) # Hours
        self.index = index
        print(str(self.index) + ". " + str(self.name))

        # 4 additional attributes will be set after the model is solved
        self.te = None
        self.tl = None
        self.tend = []
        self.orderInNight = []

    def nObs(self):
        """Calculates the total number of exposures of the star object"""

        return self.shots*self.visits

    def timeRequired(self):
        """Calculates the total time (minutes) requested of exposing on
        the star"""

        return nObs(self)*self.exptime

    def printStar(self):
        """Prints the relevant attributes of the star"""

        print("Starname: ", self.name)
        print('RA: ', np.round(self.target.ra.deg,3))
        print('Dec: ', np.round(self.target.dec.deg,3))
        print("Exp Time: ", self.exptime)
        print("Make " + str(self.visits) + " visit(s) and take " + str(self.shots) + " exposure(s) each.")
        print('Index: ', self.index)
        print()
