import astropy.units as u
from astroplan import Observer, FixedTarget
import numpy as np

'''
Note to future users:
We have designed the TTP to take in a telescope object which can be custom made for your specific observatory.
Follow the template below for Keck1 to create your own named observatory.
Important features are the telescope pointing limits, the slew rate, and a wrap limit (only one limit is supported).

Increasing the number of "nSlots" increases the number of times within the night where the TTP computes the slew path lengths
between every target. Larger numbers greatly increase the computation time. We recommend nSlots no greater than 4 unless
you have a particularly powerful computer and/or a strong need.
'''

def create_tel(user_input):

    class_ = tel_map.get(user_input)

    if class_:
        tel = class_()
        return tel
    else:
        print("This observatory is not in our pre-built classes. See the telescope.py file for instructions on how to add your own custom built telescope class. Exiting, this will break code.")
        return None

# Thank you to Ryan Rubenzahl for most of this class
class Keck1(object):
    """
    Define astropy objects specific to Keck observatory

    Also define useful functions (pointing limits, etc.)
    """

    def __init__(self, nSlots=3):
        self.observer = Observer.at_site('Keck Observatory', name='Keck', timezone='US/Hawaii')
        self.readOutTime = 45.
        self.slewrate = 6./10. # degrees per second
        self.wrapLimitAngle = 270 # degrees azimuth
        self.deckAzLim1 = 5.3
        self.deckAzLim2 = 146.2
        self.deckAltLim = 33.3
        self.vigLim = 18
        self.zenLim = 84 # Note, this cannot be 0 or no observations can be scheduled.
        self.nSlots = nSlots

    def pointing_limits(self, az, unvignetted=True):
        '''
        Return the elevations that are observable from Keck I at some azimuth
        https://www2.keck.hawaii.edu/inst/common/TelLimits.html
        '''
        zenith_limit = self.zenLim # to keep up with guiding, 88.9 deg physically
        if self.deckAzLim1 < az and az < self.deckAzLim2:
            return [self.deckAltLim, zenith_limit]
        elif unvignetted:
            return [self.vigLim, zenith_limit]  # unvignetted
        else:
            # lower 18 deg vignetted by bottom shutter
            # actual lower limit set by software at some positive elevation
            return [0, zenith_limit]

    def changeSlots(self, newSlotsVal):
        self.nSlots = newSlotsVal

    def changeSlewRate(self, newRateVal):
        self.slewrate = newRateVal

    def is_up(self, alt, az, unvignetted=True):
        # yesno is a list of binaries where if the target is up based on pointing limits, the value is 1 and if it is not up, the value is 0.
        yesno = []
        for i in range(len(alt)):
            val = int((alt[i] > self.pointing_limits(az[i], unvignetted)[0]) and (alt[i] < self.pointing_limits(az[i], unvignetted)[1]))
            if val:
                yesno.append(1)
            else:
                yesno.append(0)
        return np.array(yesno)

    def is_up_highElevation(self, alt, az, elvLimit=40.0):
        # yesno is a list of binaries where if the target is above a desireable elevation. The value is 1 it true, and 0 otherwise.
        yesno = []
        for i in range(len(alt)):
            val = int(alt[i] > elvLimit)
            if val:
                yesno.append(1)
            else:
                yesno.append(0)
        return np.array(yesno)

class WIYN(object):
    """
    Define astropy objects specific to the WIYN telescope at Kitt Peak National Observatory

    Also define useful functions (pointing limits, etc.)
    """

    def __init__(self, nSlots=3):
        self.observer = Observer.at_site('KPNO', name='WIYN', timezone='US/Arizona')
        self.readOutTime = 45.
        self.slewrate = 10/6. # degrees per second
        self.wrapLimitAngle = 270 # degrees azimuth
        self.deckAzLim1 = 0
        self.deckAzLim2 = 0
        self.deckAltLim = 0
        self.vigLim = 20
        self.zenLim = 85 # Note, this cannot be 0 or no observations can be scheduled.
        self.nSlots = nSlots

    def pointing_limits(self, az, unvignetted=True):
        '''
        Return the elevations that are observable from Keck I at some azimuth
        https://www2.keck.hawaii.edu/inst/common/TelLimits.html
        '''
        zenith_limit = self.zenLim # to keep up with guiding, 88.9 deg physically
        if self.deckAzLim1 < az and az < self.deckAzLim2:
            return [self.deckAltLim, zenith_limit]
        elif unvignetted:
            return [self.vigLim, zenith_limit]  # unvignetted
        else:
            # lower 18 deg vignetted by bottom shutter
            # actual lower limit set by software at some positive elevation
            return [0, zenith_limit]

    def changeSlots(self, newSlotsVal):
        self.nSlots = newSlotsVal

    def changeSlewRate(self, newRateVal):
        self.slewrate = newRateVal

    def is_up(self, alt, az, unvignetted=True):
        # yesno is a list of binaries where if the target is up based on pointing limits, the value is 1 and if it is not up, the value is 0.
        yesno = []
        for i in range(len(alt)):
            val = int((alt[i] > self.pointing_limits(az[i], unvignetted)[0]) and (alt[i] < self.pointing_limits(az[i], unvignetted)[1]))
            if val:
                yesno.append(1)
            else:
                yesno.append(0)
        return np.array(yesno)

    def is_up_highElevation(self, alt, az, elvLimit=40.0):
        # yesno is a list of binaries where if the target is above a desireable elevation. The value is 1 it true, and 0 otherwise.
        yesno = []
        for i in range(len(alt)):
            val = int(alt[i] > elvLimit)
            if val:
                yesno.append(1)
            else:
                yesno.append(0)
        return np.array(yesno)

# If you create a new telescope class, be sure to add it here.
tel_map = {
    'Keck1': Keck1,
    'WIYN': WIYN,
}
