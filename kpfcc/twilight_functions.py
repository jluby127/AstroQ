"""
Module for performing functions relating to twilight time determination.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import twilight_functions as tf
"""
import numpy as np
import pandas as pd

from astropy.time import Time
from astropy.time import TimeDelta
import astroplan as apl

import kpfcc.mapping_functions as mf
