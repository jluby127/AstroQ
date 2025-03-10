"""
Module for generating backup bright star list for a night. Designed to be only run as a function
call from the generateScript.py script.

Example usage:
    import backup_star_functions as bsf
"""
from astropy.time import Time
from astropy.time import TimeDelta
import astropy as apy
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import astroplan as apl
import astropy.units as u

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

import kpfcc.mapping_functions as mf
