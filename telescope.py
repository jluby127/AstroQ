import astropy.units as u
from astroplan import Observer, FixedTarget
import numpy as np

# Thank you to Ryan Rubenzahl for most of this class
class Keck1(object):
    """
    Define astropy objects specific to Keck observatory

    Also define useful functions (Keck I pointing limits, etc.)
    """

    def __init__(self, nSlots=3):
        self.observer = Observer.at_site('Keck Observatory', name='Keck', timezone='US/Hawaii')
        self.slewrate = 10/6. # degrees per second
        self.wrapLimitAngle = 270 # degrees azimuth
        self.nSlots = nSlots

    def pointing_limits(self, az, unvignetted=True):
        '''
        Return the elevations that are observable from Keck I at some azimuth
        https://www2.keck.hawaii.edu/inst/common/TelLimits.html
        '''
        zenith_limit = 84 # to keep up with guiding, 88.9 deg physically
        if 5.3 < az and az < 146.2:
            return [33.3, zenith_limit]
        elif unvignetted:
            return [18, zenith_limit]  # unvignetted
        else:
            # lower 18 deg vignetted by bottom shutter
            # actual lower limit set by software at some positive elevation
            return [0, zenith_limit]

    def changeSlots(self, newSlotsVal):
        self.nSlots = newSlotsVal

    def changeSlewRate(self, newRateVal):
        self.slewrate = newRateVal

    def is_up(self, alt, az):
        # yesno is a list of binaries where if the target is up based on pointing limits, the value is 1 and if it is not up, the value is 0.
        yesno = []
        for i in range(len(alt)):
            val = int((alt[i] > self.pointing_limits(az[i])[0]) and (alt[i] < self.pointing_limits(az[i])[1]))
            if val:
                yesno.append(1)
            else:
                yesno.append(0)
        return np.array(yesno)
