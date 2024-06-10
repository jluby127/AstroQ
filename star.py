import numpy as np
import pandas as pd
import astropy as ap
import astroplan as apl
import astropy.units as u
from astropy.coordinates import SkyCoord

# from astroplan import download_IERS_A
# download_IERS_A()

class star(object):
    """docstring for ."""

    def __init__(self, row, index):

        self.name = row['Starname']
        # ignore these notes for now!!!!
        # NOTE: RA must be in following format: XXhYYmZZ.Zs
        # NOTE: Dec must be in following format: +XXdYYmZZ.Zs
        self.ra = row['RA']
        self.dec = row['Dec']
        coords = SkyCoord(self.ra*u.deg, self.dec*u.deg, frame='icrs')
        self.target = apl.FixedTarget(name=self.name, coord=coords)
        self.exptime = row['Exposure Time']/60 # Seconds to minutes
        self.shots = row['Exposures Per Visit']
        self.visits = row['Visits In Night']
        self.expwithreadout = self.exptime*self.shots + (45/60)*(self.shots-1)
        self.intra_night_cadence = row['Intra_Night_Cadence']
        self.index = index
        print(str(self.index) + ". " + str(self.name))

        # 4 additional attributes will be set after the model is solved
        self.te = None
        self.tl = None
        self.tend = None
        self.orderInNight = None

    def nObs(self):
        return self.shots*self.visits

    def timeRequired(self):
        return nObs(self)*self.exptime

    def printStar(self):
        print("Starname: ", self.name)
        print('RA: ', np.round(self.target.ra.deg,3))
        print('Dec: ', np.round(self.target.dec.deg,3))
        print("Exp Time: ", self.exptime)
        print("Make " + str(self.visits) + " visit(s) and take " + str(self.shots) + " exposure(s) each.")
        print('Index: ', self.index)
        print()
