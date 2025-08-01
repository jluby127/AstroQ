"""
Module for computing accessibility and observability of targets, including telescope pointing,
twilight, moon constraints, and allocation maps.

This module combines functionality for:
- Basic accessibility calculations
- Comprehensive observability maps
- Allocation and scheduling constraints
"""

# Standard library imports
from typing import Dict, List, Optional, Union

# Third-party imports
import astropy as apy
import astropy.units as u
import astroplan as apl
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta

# Local imports
import astroq.weather as wh
import astroq.history as hs

class Access:
    """
    Access class that provides accessibility computation functionality.
    
    This class encapsulates all the parameters needed for accessibility computation
    and provides an object-oriented interface to the accessibility computation.
    """
    
    def __init__(self, semester_planner):
        """
        Initialize the Access object from a SemesterPlanner object.
        
        Args:
            semester_planner: SemesterPlanner object containing all configuration
        """
        # Extract core parameters from semester_planner
        self.semester_start_date = semester_planner.semester_start_date
        self.semester_length = semester_planner.semester_length
        self.slot_size = semester_planner.slot_size
        self.observatory = semester_planner.observatory
        self.current_day = semester_planner.current_day
        self.custom_file = semester_planner.custom_file
        self.allocation_file = semester_planner.allocation_file
        self.past_history_file = semester_planner.past_file
        self.today_starting_night = semester_planner.today_starting_night
        self.slots_needed_for_exposure_dict = semester_planner.slots_needed_for_exposure_dict
        self.run_weather_loss = semester_planner.run_weather_loss
        
        # Extract pre-computed date attributes (no need to recompute)
        self.all_dates_dict = semester_planner.all_dates_dict
        self.all_dates_array = semester_planner.all_dates_array
        self.n_nights_in_semester = semester_planner.n_nights_in_semester
        
        # Load past history
        self.past_history = hs.process_star_history(self.past_history_file)

    
    def _setup_observing_context(self, rf):
        """
        Set up shared observing context used by all constraint methods.
        
        Args:
            rf: Request frame containing target information
        """
        # Basic dimensions and timing
        start_date = Time(self.semester_start_date, format='iso', scale='utc')
        self.ntargets = len(rf)
        self.nnights = self.semester_length
        self.nslots = int((24*60)/self.slot_size)
        self.slot_size_time = TimeDelta(self.slot_size*u.min)
        
        # Target coordinates and observer
        coords = apy.coordinates.SkyCoord(rf.ra * u.deg, rf.dec * u.deg, frame='icrs')
        self.targets = apl.FixedTarget(name=rf.starname, coord=coords)
        self.keck = apl.Observer.at_site(self.observatory)
        
        # Set up time grid for one night, first night of the semester
        daily_start = Time(start_date, location=self.keck.location)
        daily_end = daily_start + TimeDelta(1.0, format='jd')
        self.t = Time(np.arange(daily_start.jd, daily_end.jd, self.slot_size_time.jd), 
                     format='jd', location=self.keck.location)
        self.t = self.t[np.argsort(self.t.sidereal_time('mean'))]  # sort by lst
        
        # Computing slot midpoint for all nights in semester 2D array (nights, slots)
        slotmidpoint0 = daily_start + (np.arange(self.nslots) + 0.5) * self.slot_size * u.min
        days = np.arange(self.nnights) * u.day
        self.slotmidpoint = (slotmidpoint0[np.newaxis,:] + days[:,np.newaxis])

    def _compute_altaz_mask(self):
        """
        Compute altitude/azimuth telescope pointing constraints.
        
        Returns:
            np.ndarray: Boolean mask of shape (ntargets, nnights, nslots)
        """
        # Compute base alt/az pattern, shape = (ntargets, nslots)
        coord0 = self.keck.altaz(self.t, self.targets, grid_times_targets=True)
        alt0 = coord0.alt.deg
        az0 = coord0.az.deg
        
        # 2D mask (ntargets, nslots)
        is_altaz0 = np.ones_like(alt0, dtype=bool)
        
        # Remove nasmyth deck constraint
        is_altaz0 &= ~((5.3 < az0) & (az0 < 146.2) & (alt0 < 33.3))
        
        # Remove min elevation for mid declination stars
        ismiddec = ((-30 < self.targets.dec.deg) & (self.targets.dec.deg < 75))
        fail = ismiddec[:,np.newaxis] & (alt0 < 30)  # broadcast declination array
        is_altaz0 &= ~fail
        
        # All stars must be between 18 and 85 deg
        fail = (alt0 < 18) | (alt0 > 85)
        is_altaz0 &= ~fail
        
        # Expand to 3D mask using sidereal time interpolation
        is_altaz = np.empty((self.ntargets, self.nnights, self.nslots), dtype=bool)
        
        # Pre-compute the sidereal times for interpolation
        x = self.t.sidereal_time('mean').value
        x_new = self.slotmidpoint.sidereal_time('mean').value
        idx = np.searchsorted(x, x_new, side='left')
        idx = np.clip(idx, 0, len(x)-1)  # Handle edge cases
        is_altaz = is_altaz0[:,idx]
        
        return is_altaz

    def _compute_future_mask(self):
        """
        Compute future time constraints (can't observe in the past).
        
        Returns:
            np.ndarray: Boolean mask of shape (ntargets, nnights, nslots)
        """
        is_future = np.ones((self.ntargets, self.nnights, self.nslots), dtype=bool)
        today_daynumber = self.all_dates_dict[self.current_day]
        is_future[:,:today_daynumber,:] = False
        return is_future

    def _compute_moon_mask(self):
        """
        Compute moon distance constraints (30° minimum separation).
        
        Returns:
            np.ndarray: Boolean mask of shape (ntargets, nnights, nslots)
        """
        is_moon = np.ones((self.ntargets, self.nnights, self.nslots), dtype=bool)
        moon = apy.coordinates.get_moon(self.slotmidpoint[:,0], self.keck.location)
        
        # Reshaping uses broadcasting to achieve a (ntarget, night) array
        ang_dist = apy.coordinates.angular_separation(
            self.targets.ra.reshape(-1,1), self.targets.dec.reshape(-1,1),
            moon.ra.reshape(1,-1), moon.dec.reshape(1,-1),
        )
        is_moon = is_moon & (ang_dist.to(u.deg) > 30*u.deg)[:, :, np.newaxis]
        return is_moon

    def _compute_custom_times_mask(self, rf):
        """
        Compute custom time window constraints from CSV file.
        
        Args:
            rf: Request frame for target name mapping
            
        Returns:
            np.ndarray: Boolean mask of shape (ntargets, nnights, nslots)
        """
        # Start with all times allowed (True everywhere)
        is_custom = np.ones((self.ntargets, self.nnights, self.nslots), dtype=bool)
        df = pd.read_csv(self.custom_file)
        
        # Check if the file has any data rows (not just header)
        if len(df) == 0:
            return is_custom
            
        # Create mapping from starname to target index
        starname_to_index = {name: idx for idx, name in enumerate(rf['starname'])}
        
        # Convert time columns to astropy Time objects
        df['start'] = df['start'].apply(Time)
        df['stop'] = df['stop'].apply(Time)
        
        # Check for stars in custom file that are not in request frame
        custom_stars = set(df['starname'].unique())
        request_stars = set(starname_to_index.keys())
        missing_stars = custom_stars - request_stars
        
        if missing_stars:
            missing_list = sorted(list(missing_stars))
            raise ValueError(
                f"Custom times file contains {len(missing_stars)} star(s) not found in request frame: "
                f"{', '.join(missing_list)}. Please check that star names in custom.csv match those in request.csv"
            )
        
        # All stars are valid, so no filtering needed
        # But keep the structure in case we want to add filtering logic later
        
        # Process constraints by star using groupby for cleaner logic
        for starname, star_group in df.groupby('starname'):
            star_idx = starname_to_index[starname]
            
            # Initialize this star's mask as False (no time allowed initially)
            star_mask = np.zeros((self.nnights, self.nslots), dtype=bool)
            
            # Union all time windows for this star
            for _, interval in star_group.iterrows():
                interval_mask = ((self.slotmidpoint >= interval['start']) & 
                               (self.slotmidpoint <= interval['stop']))
                star_mask |= interval_mask
            
            # Apply the combined mask for this star
            is_custom[star_idx] = star_mask
        
        return is_custom

    def _compute_inter_night_cadence_mask(self, rf):
        """
        Compute inter-night cadence constraints.
        
        Args:
            rf: Request frame containing tau_inter values
            
        Returns:
            np.ndarray: Boolean mask of shape (ntargets, nnights, nslots)
        """
        is_inter = np.ones((self.ntargets, self.nnights, self.nslots), dtype=bool)
        
        for itarget in range(self.ntargets):
            name = rf.iloc[itarget]['starname']
            if name in self.past_history and rf.iloc[itarget]['tau_inter'] > 1:
                inight_start = self.all_dates_dict[self.past_history.date_last_observed] - self.today_starting_night
                inight_stop = min(inight_start + rf.iloc[itarget]['tau_inter'], self.nnights)
                is_inter[itarget, inight_start:inight_stop, :] = False
        
        return is_inter

    def _compute_allocation_mask(self):
        """
        Compute allocated time window constraints from CSV file.
        
        Returns:
            np.ndarray: Boolean mask of shape (ntargets, nnights, nslots)
        """
        allocated_times_frame = pd.read_csv(self.allocation_file)
        allocated_times_frame['start'] = allocated_times_frame['start'].apply(Time)
        allocated_times_frame['stop'] = allocated_times_frame['stop'].apply(Time)
        
        allocated_mask = np.zeros((self.nnights, self.nslots), dtype=bool)
        for i in range(len(allocated_times_frame)):
            start_time = allocated_times_frame['start'].iloc[i]
            stop_time = allocated_times_frame['stop'].iloc[i]
            mask = (self.slotmidpoint >= start_time) & (self.slotmidpoint <= stop_time)
            allocated_mask |= mask
        
        # Expand to 3D array for all targets
        is_alloc = np.ones((self.ntargets, self.nnights, self.nslots), dtype=bool) & allocated_mask[np.newaxis,:,:]
        return is_alloc

    def _compute_weather_mask(self):
        """
        Compute weather clearness mask using weather simulation.
        
        Returns:
            np.ndarray: Boolean mask of shape (ntargets, nnights, nslots)
        """
        is_clear = np.ones((self.ntargets, self.nnights, self.nslots), dtype=bool)
        
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
            print("Pretending weather is always clear!")
        
        return is_clear

    def _compute_multi_slot_exposure_mask(self, rf, is_observable_now):
        """
        Compute multi-slot exposure constraints by checking consecutive slot availability.
        
        Args:
            rf: Request frame containing starname for lookup
            is_observable_now: Base observability mask before multi-slot adjustment
            
        Returns:
            np.ndarray: Boolean mask of shape (ntargets, nnights, nslots)
        """
        is_observable = is_observable_now.copy()
        
        for itarget in range(self.ntargets):
            e_val = self.slots_needed_for_exposure_dict[rf.iloc[itarget]['starname']]
            if e_val == 1:
                continue
                
            for shift in range(1, e_val):
                # Shifts the is_observable_now array to the left by shift
                # For is_observable to be true, it must be true for all shifts
                is_observable[itarget, :, :-shift] &= is_observable_now[itarget, :, shift:]
        
        return is_observable

    def produce_ultimate_map(self, rf):
        """
        Combine all maps for a target to produce the final map.

        Args:
            rf (dataframe): request frame

        Returns:
            np.rec.array: Record array containing all constraint masks and final observability
        """
        # Set up shared observing context
        self._setup_observing_context(rf)
        
        # Compute individual constraint masks
        is_altaz = self._compute_altaz_mask()
        is_future = self._compute_future_mask()
        is_moon = self._compute_moon_mask()
        is_custom = self._compute_custom_times_mask(rf)
        is_inter = self._compute_inter_night_cadence_mask(rf)
        is_alloc = self._compute_allocation_mask()
        is_clear = self._compute_weather_mask()
        
        # Combine all constraints
        is_observable_now = np.logical_and.reduce([
            is_altaz,
            is_future,
            is_moon,
            is_custom,
            is_inter,
            is_alloc,
            is_clear
        ])
        
        # Apply multi-slot exposure constraints
        is_observable = self._compute_multi_slot_exposure_mask(rf, is_observable_now)
        
        # Package results
        access = {
            'is_altaz': is_altaz,
            'is_future': is_future,
            'is_moon': is_moon,
            'is_custom': is_custom,
            'is_inter': is_inter,
            'is_alloc': is_alloc,
            'is_clear': is_clear,
            'is_observable_now': is_observable_now,
            'is_observable': is_observable
        }
        access = np.rec.fromarrays(list(access.values()), names=list(access.keys()))
        return access

    def observability(self, requests_frame):
        """
        Extract available indices dictionary from the record array returned by produce_ultimate_map
        
        Args:
            access: Record array from produce_ultimate_map containing observability masks
            requests_frame: DataFrame containing request information with starname column
            
        Returns:
            dict: Dictionary where keys are target names and values are lists of available slots per night
        """
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