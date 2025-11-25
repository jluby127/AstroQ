"""
Module for computing the intersection of the various accessibility maps for all targets for the following constraints:
    - telescope pointing
    - telescope allocation
    - sky brightness
    - moon separation
    - internight cadence from past history
    - PI custom windows
    - simulated weather loss
    - enough time to complete the exposure tonight

The Access class is saved as an attribute of the splan object and used again in plotting.
"""

# Standard library imports
import logging

# Third-party imports
from astropy.utils.iers import conf
conf.auto_max_age = None
import astropy as apy
import astropy.units as u
import astroplan as apl
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
import os

DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),'data')

logs = logging.getLogger(__name__)

# Keck limits -- lets talk about making this a subclass
nays_az_low = 5.3
nays_az_high = 146.2
nays_alt = 33.3
tel_min = 18
tel_max = 85

class Access:
    """
    The Access class encapsulates all the parameters needed for accessibility computation
    and provides an object-oriented interface to the accessibility computation.
    """
    
    def __init__(self, 
                 semester_start_date, 
                 semester_length, 
                 n_nights_in_semester,
                 today_starting_night,  
                 current_day, 
                 all_dates_dict, 
                 all_dates_array, 
                 slot_size, 
                 slots_needed_for_exposure_dict, 
                 custom_file, 
                 allocation_file, 
                 past_history, 
                 output_directory, 
                 run_weather_loss, 
                 run_band3, 
                 observatory_string, 
                 request_frame, 
                 ):
        """
        Initialize the Access object with explicit parameters.
        
        Args:
            semester_start_date: Start date of the semester
            semester_length: Total number of nights in the semester
            n_nights_in_semester: Number of remaining nights in the semester
            today_starting_night: Starting night number for today
            current_day: Current day identifier
            all_dates_dict: Dictionary mapping dates to day numbers
            all_dates_array: Array of date strings for the semester
            slot_size: Size of each time slot in minutes
            slots_needed_for_exposure_dict: Dictionary mapping star names to required slots
            custom_file: Path to custom times file
            allocation_file: Path to allocation file
            past_history: Past observation history
            output_directory: Directory for output files
            run_weather_loss: Whether to run weather loss simulation
            run_band3: Whether to run band 3 (used for not peforming the is_observble step for the football plot)
            observatory_string: Observatory name/location string
            request_frame: DataFrame containing request information
        """
        # parameters 
        self.semester_start_date = semester_start_date
        self.semester_length = semester_length
        self.slot_size = slot_size
        self.current_day = current_day
        self.all_dates_dict = all_dates_dict
        self.all_dates_array = all_dates_array
        self.n_nights_in_semester = n_nights_in_semester
        self.today_starting_night = today_starting_night
        self.slots_needed_for_exposure_dict = slots_needed_for_exposure_dict
        self.run_weather_loss = run_weather_loss
        self.run_band3 = run_band3

        # files
        self.custom_file = custom_file
        self.allocation_file = allocation_file
        self.past_history = past_history
        self.output_directory = output_directory
        self.request_frame = request_frame

        # Prepatory work
        self.start_date = Time(self.semester_start_date,format='iso',scale='utc')
        self.ntargets = len(self.request_frame)
        self.nnights = self.semester_length # total nights in the full semester
        self.nslots = int((24*60)/self.slot_size) # slots in the night
        self.slot_size_time = TimeDelta(self.slot_size*u.min)
        self.observatory = apl.Observer.at_site(observatory_string)
        coords = apy.coordinates.SkyCoord(self.request_frame.ra * u.deg, self.request_frame.dec * u.deg, frame='icrs')
        self.targets = apl.FixedTarget(name=self.request_frame.unique_id, coord=coords)

        # Set up time grid for one night, first night of the semester
        self.daily_start = Time(self.start_date, location=self.observatory.location)
        self.daily_end = self.daily_start + TimeDelta(1.0, format='jd') # full day from start of first night
        self.timegrid = Time(np.arange(self.daily_start.jd, self.daily_end.jd, self.slot_size_time.jd), format='jd',location=self.observatory.location)
        self.timegrid = self.timegrid[np.argsort(self.timegrid.sidereal_time('mean'))] # sort by lst
    
        # computing slot midpoint for all nights in semester 2D array (slots, nights)
        self.slotmidpoints_oneday = self.daily_start + (np.arange(self.nslots) + 0.5) *  self.slot_size * u.min
        days = np.arange(self.nnights) * u.day
        self.slotmidpoints = (self.slotmidpoints_oneday[np.newaxis,:] + days[:,np.newaxis])

    def compute_altaz(self, tel_min):
        """
        Compute boolean mask of is_altaz for targets according to a minimum elevation. 
        May be superceded by a specific compute_altaz method for a specific observatory, see astroq/queue/ modules.

        Args:
            tel_min (float): the minimum elevation for the telescope

        Returns:
            is_altaz (array): boolean mask of is_altaz for targets
        """
        # Compute base alt/az pattern, shape = (ntargets, nslots)
        altazes = self.observatory.altaz(self.timegrid, self.targets, grid_times_targets=True)
        alts = altazes.alt.deg
        min_elevation = self.request_frame['minimum_elevation'].values  # Get PI-desired minimum elevation values
        min_elevation = np.maximum(min_elevation, tel_min)  # Ensure minimum elevation is at least tel_min
        
        # Pre-compute the sidereal times for interpolation
        x = self.timegrid.sidereal_time('mean').value
        x_new = self.slotmidpoints.sidereal_time('mean').value
        idx = np.searchsorted(x, x_new, side='left')
        idx = np.clip(idx, 0, len(x)-1) # Handle edge cases

        is_altaz0 = np.ones_like(alts, dtype=bool)        
        fail = (alts < tel_min)
        is_altaz0 &= ~fail
        # self.is_altaz = np.empty((self.ntargets, self.nnights, self.nslots),dtype=bool)
        self.is_altaz = is_altaz0[:,idx]

    def compute_future(self):
        """
        Compute boolean mask of is_future for all targets according to today's current_day. 

        Args:
        Returns:
            is_altaz (array): boolean mask of is_altaz for targets
        """
        self.is_future = np.ones((self.ntargets, self.nnights, self.nslots),dtype=bool)
        today_daynumber = self.all_dates_dict[self.current_day]
        self.is_future[:,:today_daynumber,:] = False
    
    def compute_moon(self):
        """
        Compute boolean mask of is_moon for all targets according to the moon's position.
        """
        self.is_moon = np.ones_like(self.is_altaz, dtype=bool)
        moon = apy.coordinates.get_moon(self.slotmidpoints[:,0] , self.observatory.location)
        # Reshaping uses broadcasting to achieve a (ntarget, night) array
        ang_dist = apy.coordinates.angular_separation(
            self.targets.ra.reshape(-1,1), self.targets.dec.reshape(-1,1),
            moon.ra.reshape(1,-1), moon.dec.reshape(1,-1),
        )
        # Use per-row minimum_moon_separation values instead of hardcoded 30 degrees
        min_moon_sep = self.request_frame['minimum_moon_separation'].values * u.deg  # Convert to degrees
        self.is_moon = self.is_moon & (ang_dist.to(u.deg) > min_moon_sep[:, np.newaxis])[:, :, np.newaxis]

    def compute_inter(self):
        """
        Compute boolean mask of is_inter for all targets according to the internight cadence.
        """
        # Set to False if internight cadence is violated
        self.is_inter = np.ones((self.ntargets, self.nnights, self.nslots),dtype=bool)
        for itarget in range(self.ntargets):
            name = self.request_frame.iloc[itarget]['unique_id']
            if name in self.past_history and self.request_frame.iloc[itarget]['tau_inter'] > 1:
                inight_start = self.all_dates_dict[self.past_history[name].date_last_observed]
                inight_stop = min(inight_start + self.request_frame.iloc[itarget]['tau_inter'],self.nnights)
                self.is_inter[itarget,inight_start:inight_stop,:] = False

    def compute_custom(self):
        """
        Compute boolean mask of is_custom for all targets according to the custom times.
        """
        self.is_custom = np.ones((self.ntargets, self.nnights, self.nslots), dtype=bool)
        # Handle case where custom file doesn't exist
        if os.path.exists(self.custom_file):
            custom_times_frame = pd.read_csv(self.custom_file)
            # Check if the file has any data rows (not just header)
            if len(custom_times_frame) > 0:
                starid_to_index = {name: idx for idx, name in enumerate(self.request_frame['unique_id'])}
                custom_times_frame['start'] = custom_times_frame['start'].apply(Time)
                custom_times_frame['stop'] = custom_times_frame['stop'].apply(Time)
                for _, row in custom_times_frame.iterrows():
                    starid = row['unique_id']
                    # Skip if the star is not in the current requests frame
                    if starid not in starid_to_index:
                        print(f"Warning: Star {row['starname']} with unique_id '{starid}' in custom times file not found in requests frame, skipping")
                        continue
                    mask = (self.slotmidpoints >= row['start']) & (self.slotmidpoints <= row['stop'])
                    star_ind = starid_to_index[starid]
                    current_map = self.is_custom[star_ind]
                    if np.all(current_map):  # If all ones, first interval: restrict with AND
                        self.is_custom[star_ind] = mask
                    else:  # Otherwise, union with OR
                        self.is_custom[star_ind] = current_map | mask
        else:
            print(f"Custom times file not found: {self.custom_file}. Using no custom constraints.")

    def compute_allocated(self):
        """
        Compute boolean mask of is_allocated for all targets according to the allocated times.
        """
        allocated_times_frame = pd.read_csv(self.allocation_file)
        allocated_times_frame['start'] = allocated_times_frame['start'].apply(Time)
        allocated_times_frame['stop'] = allocated_times_frame['stop'].apply(Time)
            
        allocated_times_map = []
        allocated_mask = np.zeros((self.nnights, self.nslots), dtype=bool)
        for i in range(len(allocated_times_frame)):
            start_time = allocated_times_frame['start'].iloc[i]
            stop_time = allocated_times_frame['stop'].iloc[i]
            mask = (self.slotmidpoints >= start_time) & (self.slotmidpoints <= stop_time)
            allocated_mask |= mask
        self.is_allocated_mask = allocated_mask
        self.is_allocated = np.ones_like(self.is_altaz, dtype=bool) & self.is_allocated_mask[np.newaxis,:,:] # shape = (ntargets, nnights, nslots)

    def compute_clear(self,weather_loss_file=None):
        """
        Compute boolean mask of is_clear for all targets according to the clear times.

        Args:
            weather_loss_file: Path to file with weather loss statistics information
        """
        self.is_clear = np.ones_like(self.is_altaz, dtype=bool)
        if self.run_weather_loss:
            if weather_loss_file is None:
                raise ValueError("weather_loss_file is required when run_weather_loss is True")
            logs.info("Running weather loss model.")
            self.get_loss_stats(weather_loss_file)
            self.is_clear = self.simulate_weather_losses(covariance=0.14)
            self.is_clear = np.tile(self.is_clear[np.newaxis, :, :], (self.ntargets, 1, 1))
        else:
            logs.info("Pretending weather is always clear!")
            self.is_clear = np.ones((self.ntargets, self.nnights, self.nslots), dtype=bool)

    def produce_ultimate_map(self, running_backup_stars=False):
        """
        Compute boolean mask of is_observable for all targets according to the ultimate map.
        """
        self.compute_altaz(33)
        self.compute_future()
        self.compute_moon()
        self.compute_custom()
        self.compute_inter()
        self.compute_allocated()
        self.compute_clear()
        
        self.is_observable_now = np.logical_and.reduce([
            self.is_altaz,
            self.is_future,
            self.is_moon,
            self.is_custom,
            self.is_inter,
            self.is_allocated,
            self.is_clear,
        ])

        # the target does not violate any of the observability limits in that specific slot, but
        # it does not mean it can be started at the slot. retroactively grow mask to accomodate multishot exposures.
        # Is observable now,
        self.is_observable = self.is_observable_now.copy()
        if running_backup_stars == False:
            for itarget in range(self.ntargets):
                e_val = self.slots_needed_for_exposure_dict[self.request_frame.iloc[itarget]['unique_id']]
                if e_val == 1:
                    continue
                for shift in range(1, e_val):
                    # shifts the is_observable_now array to the left by shift
                    # for is_observable to be true, it must be true for all shifts
                    self.is_observable[itarget, :, :-shift] &= self.is_observable_now[itarget, :, shift:]

        access = {
            'is_altaz': self.is_altaz,
            'is_future': self.is_future,
            'is_moon': self.is_moon,
            'is_custom':self.is_custom,
            'is_inter': self.is_inter,
            'is_alloc': self.is_allocated,
            'is_clear': self.is_clear,
            'is_observable_now': self.is_observable_now,
            'is_observable': self.is_observable
        }
        access_record = np.rec.fromarrays(list(access.values()), names=list(access.keys()))
        return access_record

    def observability(self, requests_frame, access=None):
        """
        Extract a dictionary of the available indices from the record array returned by produce_ultimate_map
        
        Args:
            requests_frame: DataFrame containing request information
            access: Optional record array from produce_ultimate_map (if None, this function will compute it)
            
        Returns:
            df (dict): Dictionary where keys are target names and values are lists of available slots per night
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
        df = pd.merge(requests_frame[['unique_id']].reset_index(drop=True),df,left_index=True,right_on='itarget')
        namemap = {'starid':'unique_id','inight':'d','islot':'s'}
        df = df.query('is_observable').rename(columns=namemap)[namemap.values()]
        return df

    def get_loss_stats(self, weather_loss_file):
        """
        Gather the loss probabilities for each night in the semester from the saved historical weather data.
        """
        historical_weather_data = pd.read_csv(os.path.join(DATADIR,weather_loss_file))
        loss_stats_this_semester = []
        for i, item in enumerate(self.all_dates_array):
            ind = historical_weather_data.index[historical_weather_data['Date'] == \
                self.all_dates_array[i][5:]].tolist()[0]
            loss_stats_this_semester.append(historical_weather_data['% Total Loss'][ind])
        self.loss_stats_this_semester = loss_stats_this_semester

    def simulate_weather_losses(self, covariance=0.14):
        """
        Simulate nights totally lost to weather using historical data

        Args:
            covariance (float): the added percent chance that tomorrow will be lost if today is lost

        Returns:
            is_clear (array): Trues represent clear nights, Falses represent weathered nights
        """
        previous_day_was_lost = False
        is_clear = np.ones((self.semester_length, int((24*60)/ self.slot_size)),dtype=bool)
        for i in range(len(self.loss_stats_this_semester)):
            value_to_beat = self.loss_stats_this_semester[i]
            if previous_day_was_lost:
                value_to_beat += covariance
            roll_the_dice = np.random.uniform(0.0,1.0)

            if roll_the_dice < value_to_beat:
                # the night is simulated a total loss
                is_clear[i] = np.zeros(is_clear.shape[1])  # Set all slots to False
                previous_day_was_lost = True
            else:
                previous_day_was_lost = False
        logs.info(f"Total nights simulated as weathered out: {np.sum(~np.any(is_clear, axis=1))} of {len(is_clear)} nights remaining.")
        return is_clear

def build_twilight_allocation_file(semester_planner):
    """
    Build an allocation.csv file where every night of the semester is allocated 
    from evening to morning 12-degree twilight times. 
    This is used exclusively by the football plot in the webapp.
    
    Args:
        semester_planner (SemesterPlanner): a semester planner object from splan.py
        
    Returns:
        twilight_file (str): Path to the created allocation.csv file
    """
    
    # Create the filename based on semester
    semester = semester_planner.semester_start_date[:4] + semester_planner.semester_letter
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
    twilight_file = os.path.join(data_dir, f"{semester}_twilights.csv")
    
    # Check if file already exists
    if os.path.exists(twilight_file):
        return twilight_file
    
    # Create data directory if it doesn't exist
    os.makedirs(data_dir, exist_ok=True)
    
    # Get observatory location
    observatory = apl.Observer.at_site(semester_planner.observatory)
    
    # Create allocation data
    allocation_data = []
    
    for date_str in semester_planner.all_dates_dict.keys():
        # Parse the date
        date = Time(date_str, format='iso', scale='utc')
        
        # Get 12-degree twilight times for this night
        evening_12 = observatory.twilight_evening_nautical(date, which='next')
        morning_12 = observatory.twilight_morning_nautical(date, which='next')
        
        # Add to allocation data
        allocation_data.append({
            'start': evening_12.isot,
            'stop': morning_12.isot
        })
    
    # Create DataFrame and save to CSV
    twilight_df = pd.DataFrame(allocation_data)
    twilight_df.to_csv(twilight_file, index=False)
    
    return twilight_file