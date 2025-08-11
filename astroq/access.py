"""
Module for computing accessibility and observability of targets, including telescope pointing,
twilight, moon constraints, and allocation maps.

This module combines functionality for:
- Basic accessibility calculations
- Comprehensive observability maps
- Allocation and scheduling constraints
"""

# Standard library imports
import logging

# Third-party imports
import astropy as apy
import astropy.units as u
import astroplan as apl
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
import os

# Local imports
import astroq.weather as wh

logs = logging.getLogger(__name__)

# Define list of observatories which are currently supported.
# To add an obseratory/telescope, add the Astroplan resolvable name to the list in generate_night_plan
# Then add the same name to the appropriate element of the locations dictionary.
# If necessary, add a location to the locations dictionary, if so, add the location to each of the
# pre_sunrise and post_sunrise dictionaries. Ensure times are 14 hours apart, at least one hour
# before the earliest sunset and one hour after the latest sunrise of the year.
locations = {"Keck Observatory":"Hawaii", "Kitt Peak National Observatory":"Arizona"}
pre_sunset = {'Hawaii':'03:30', 'Arizona':'05:00'}
post_sunrise = {'Hawaii':'17:30', 'Arizona':'14:00'}

class Access:
    """
    Access class that provides accessibility computation functionality.
    
    This class encapsulates all the parameters needed for accessibility computation
    and provides an object-oriented interface to the accessibility computation.
    """
    
    def __init__(self, semester_start_date, semester_length, slot_size, observatory, 
                 current_day, all_dates_dict, all_dates_array, n_nights_in_semester,
                 custom_file, allocation_file, past_history, today_starting_night, 
                 slots_needed_for_exposure_dict, run_weather_loss, output_directory, run_band3):
        """
        Initialize the Access object with explicit parameters.
        
        Args:
            semester_start_date: Start date of the semester
            semester_length: Total number of nights in the semester
            slot_size: Size of each time slot in minutes
            observatory: Observatory name/location
            current_day: Current day identifier
            all_dates_dict: Dictionary mapping dates to day numbers
            all_dates_array: Array of date strings for the semester
            n_nights_in_semester: Number of remaining nights in the semester
            custom_file: Path to custom times file
            allocation_file: Path to allocation file
            past_history: Past observation history
            today_starting_night: Starting night number for today
            slots_needed_for_exposure_dict: Dictionary mapping star names to required slots
            run_weather_loss: Whether to run weather loss simulation
            output_directory: Directory for output files
        """
        self.semester_start_date = semester_start_date
        self.semester_length = semester_length
        self.slot_size = slot_size
        self.observatory = observatory
        self.current_day = current_day
        self.all_dates_dict = all_dates_dict
        self.all_dates_array = all_dates_array
        self.n_nights_in_semester = n_nights_in_semester
        self.custom_file = custom_file
        self.allocation_file = allocation_file
        self.past_history = past_history
        self.today_starting_night = today_starting_night
        self.slots_needed_for_exposure_dict = slots_needed_for_exposure_dict
        self.run_weather_loss = run_weather_loss
        self.output_directory = output_directory
        self.run_band3 = run_band3
    
    def produce_ultimate_map(self, rf, running_backup_stars=False):
        """
        Combine all maps for a target to produce the final map

        Args:
            rf (dataframe): request frame
            running_backup_stars (bool): if true, then do not run the extra map of stepping back in time to account for the starting slot fitting into the night

        Returns:
            available_indices_for_request (dictionary): keys are the starnames and values are a 1D array
                                                      the indices where available_slots_for_request is 1.
        """
        # Prepatory work
        start_date = Time(self.semester_start_date,format='iso',scale='utc')
        ntargets = len(rf)
        nnights = self.semester_length # total nights in the full semester
        nslots = int((24*60)/self.slot_size) # slots in the night
        slot_size_time = TimeDelta(self.slot_size*u.min)

        # Determine observability
        coords = apy.coordinates.SkyCoord(rf.ra * u.deg, rf.dec * u.deg, frame='icrs')
        targets = apl.FixedTarget(name=rf.starname, coord=coords)
        keck = apl.Observer.at_site(self.observatory)

        # Set up time grid for one night, first night of the semester
        daily_start = Time(start_date, location=keck.location)
        daily_end = daily_start + TimeDelta(1.0, format='jd') # full day from start of first night
        t = Time(np.arange(daily_start.jd, daily_end.jd, slot_size_time.jd), format='jd',location=keck.location)
        t = t[np.argsort(t.sidereal_time('mean'))] # sort by lst

        # Compute base alt/az pattern, shape = (ntargets, nslots)
        coord0 = keck.altaz(t, targets, grid_times_targets=True)
        alt0 = coord0.alt.deg
        az0 = coord0.az.deg
        # 2D mask (n targets, n slots))
        is_altaz0 = np.ones_like(alt0, dtype=bool)
        is_altaz0 &= ~((5.3 < az0 ) & (az0 < 146.2) & (alt0 < 33.3)) # remove nasymth deck
        # remove min elevation for mid declination stars
        ismiddec = ((-30 < targets.dec.deg) & (targets.dec.deg < 75))
        fail = ismiddec[:,np.newaxis] & (alt0 < 30) # broadcast declination array
        is_altaz0 &= ~fail
        # all stars must be between 18 and 85 deg
        fail = (alt0 < 18) | (alt0 > 85)
        # fail = (alt0 < 38) | (alt0 > 85)
        is_altaz0 &= ~fail
        # computing slot midpoint for all nights in semester 2D array (slots, nights)
        slotmidpoint0 = daily_start + (np.arange(nslots) + 0.5) *  self.slot_size * u.min
        days = np.arange(nnights) * u.day
        slotmidpoint = (slotmidpoint0[np.newaxis,:] + days[:,np.newaxis])
        # 3D mask
        is_altaz = np.empty((ntargets, nnights, nslots),dtype=bool)
        # Pre-compute the sidereal times for interpolation
        x = t.sidereal_time('mean').value
        x_new = slotmidpoint.sidereal_time('mean').value
        idx = np.searchsorted(x, x_new, side='left')
        idx = np.clip(idx, 0, len(x)-1) # Handle edge cases
        is_altaz = is_altaz0[:,idx]

        is_future = np.ones((ntargets, nnights, nslots),dtype=bool)
        today_daynumber = self.all_dates_dict[self.current_day]
        is_future[:,:today_daynumber,:] = False
        
        # Compute moon accessibility
        is_moon = np.ones_like(is_altaz, dtype=bool)
        moon = apy.coordinates.get_moon(slotmidpoint[:,0] , keck.location)
        # Reshaping uses broadcasting to achieve a (ntarget, night) array
        ang_dist = apy.coordinates.angular_separation(
            targets.ra.reshape(-1,1), targets.dec.reshape(-1,1),
            moon.ra.reshape(1,-1), moon.dec.reshape(1,-1),
        ) # (ntargets)
        is_moon = is_moon & (ang_dist.to(u.deg) > 30*u.deg)[:, :, np.newaxis]

        is_custom = np.ones((ntargets, nnights, nslots), dtype=bool)
        # Handle case where custom file doesn't exist
        if os.path.exists(self.custom_file):
            custom_times_frame = pd.read_csv(self.custom_file)
            # Check if the file has any data rows (not just header)
            if len(custom_times_frame) > 0:
                starname_to_index = {name: idx for idx, name in enumerate(rf['starname'])}
                custom_times_frame['start'] = custom_times_frame['start'].apply(Time)
                custom_times_frame['stop'] = custom_times_frame['stop'].apply(Time)
                for _, row in custom_times_frame.iterrows():
                    starname = row['starname']
                    # Skip if the star is not in the current requests frame
                    if starname not in starname_to_index:
                        print(f"Warning: Star '{starname}' in custom times file not found in requests frame, skipping")
                        continue
                    mask = (slotmidpoint >= row['start']) & (slotmidpoint <= row['stop'])
                    star_ind = starname_to_index[starname]
                    current_map = is_custom[star_ind]
                    if np.all(current_map):  # If all ones, first interval: restrict with AND
                        is_custom[star_ind] = mask
                    else:  # Otherwise, union with OR
                        is_custom[star_ind] = current_map | mask
        else:
            print(f"Custom times file not found: {self.custom_file}. Using no custom constraints.")

        # TODO add in logic to remove stars that are not observable, currently code is a no-op

        # Set to False if internight cadence is violated
        is_inter = np.ones((ntargets, nnights, nslots),dtype=bool)
        for itarget in range(ntargets):
            name = rf.iloc[itarget]['starname']
            if name in self.past_history and rf.iloc[itarget]['tau_inter'] > 1:
                inight_start = self.all_dates_dict[self.past_history[name].date_last_observed] - self.today_starting_night
                inight_stop = min(inight_start + rf.iloc[itarget]['tau_inter'],nnights)
                is_inter[itarget,inight_start:inight_stop,:] = False

        allocated_times_frame = pd.read_csv(self.allocation_file)
        allocated_times_frame['start'] = allocated_times_frame['start'].apply(Time)
        allocated_times_frame['stop'] = allocated_times_frame['stop'].apply(Time)
            
        allocated_times_map = []
        allocated_mask = np.zeros((nnights, nslots), dtype=bool)
        for i in range(len(allocated_times_frame)):
            start_time = allocated_times_frame['start'].iloc[i]
            stop_time = allocated_times_frame['stop'].iloc[i]
            mask = (slotmidpoint >= start_time) & (slotmidpoint <= stop_time)
            allocated_mask |= mask
        is_alloc = allocated_mask
        is_alloc = np.ones_like(is_altaz, dtype=bool) & is_alloc[np.newaxis,:,:] # shape = (ntargets, nnights, nslots)

        # run the weather loss simulation
        is_clear = np.ones_like(is_altaz, dtype=bool)
        if self.run_weather_loss:
            loss_stats_this_semester = wh.get_loss_stats(self.all_dates_array)
            is_clear = wh.simulate_weather_losses(
                self.semester_length, 
                self.n_nights_in_semester, 
                self.slot_size, 
                loss_stats_this_semester, 
                covariance=0.14
            )
        else:
            logs.info("Pretending weather is always clear!")

        is_observable_now = np.logical_and.reduce([
            is_altaz,
            is_future,
            is_moon,
            is_custom,
            is_inter,
            is_alloc,
            is_clear,
        ])

        # the target does not violate any of the observability limits in that specific slot, but
        # it does not mean it can be started at the slot. retroactively grow mask to accomodate multishot exposures.
        # Is observable now,
        is_observable = is_observable_now.copy()
        if running_backup_stars == False:
            for itarget in range(ntargets):
                e_val = self.slots_needed_for_exposure_dict[rf.iloc[itarget]['starname']]
                if e_val == 1:
                    continue
                for shift in range(1, e_val):
                    # shifts the is_observable_now array to the left by shift
                    # for is_observable to be true, it must be true for all shifts
                    is_observable[itarget, :, :-shift] &= is_observable_now[itarget, :, shift:]

        access = {
            'is_altaz': is_altaz,
            'is_future': is_future,
            'is_moon': is_moon,
            'is_custom':is_custom,
            'is_inter': is_inter,
            'is_alloc': is_alloc,
            'is_clear': is_clear,
            'is_observable_now': is_observable_now,
            'is_observable': is_observable
        }
        access = np.rec.fromarrays(list(access.values()), names=list(access.keys()))
        return access

    def observability(self, requests_frame, access=None):
        """
        Extract available indices dictionary from the record array returned by produce_ultimate_map
        
        Args:
            requests_frame: DataFrame containing request information with starname column
            access: Optional record array from produce_ultimate_map (if None, will compute it)
            
        Returns:
            dict: Dictionary where keys are target names and values are lists of available slots per night
        """
        if access is None:
            access = self.produce_ultimate_map(requests_frame)
        ntargets, nnights, nslots = access.shape
        
        # specify indeces of 3D observability array
        itarget, inight, islot = np.mgrid[:ntargets,:nnights,:nslots]

        # define flat table to access maps
        df = pd.DataFrame(
            {'itarget':itarget.flatten(),
             'inight':inight.flatten(),
             'islot':islot.flatten()}
        )
        df['is_observable'] = access.is_observable.flatten()
        df = pd.merge(requests_frame[['starname']].reset_index(drop=True),df,left_index=True,right_on='itarget')
        namemap = {'starname':'id','inight':'d','islot':'s'}
        df = df.query('is_observable').rename(columns=namemap)[namemap.values()]
        return df
