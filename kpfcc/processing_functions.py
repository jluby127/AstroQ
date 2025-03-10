"""
Module for processin the inputs and outputs of the autoscheduler to/from various sources.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import processing_functions as pf
"""
import os

from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.time import TimeDelta

import numpy as np
import pandas as pd
JUMP_USERNAME = os.environ['KPFCC_JUMP_USERNAME']
JUMP_PASSWORD = os.environ['KPFCC_JUMP_PASSWORD']
