"""
Module for night-level observation planning and optimization.
Uses the vendored TTP MILP solver (``astroq.ttp.model.TTPModel``) to optimize
nightly observation sequences.
"""

# Standard library imports
import os

# Third-party imports
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time, TimeDelta

# Local imports
from astroq.splan import SemesterPlanner
from astroq.ttp import model

# HDF5 schema version for night_planner.h5. Bumped when the layout changes so
# stale files fail-fast with a clear message rather than silently misloading.
NIGHT_PLANNER_H5_SCHEMA = 2

class NightPlanner(object):
    """
    The NightPlanner object is responsible for preparing, running, and outputting the TTP slew path optimization. 
    It is built from the config file and requires a semester_planner object to have been created and saved to an h5 file first.
    """
    
    def __init__(self, config_file):
        """
        Initialize the Night Planner with a config file.
        
        Args:
            config_file: Path to configuration file
        """
        
        # Parse config file directly for paths (following SemesterPlanner pattern)
        from configparser import ConfigParser
        config = ConfigParser()
        config.read(config_file)
        
        # Get workdir from global section
        workdir = str(config.get('global', 'workdir'))
        self.upstream_path = workdir
        self.semester_directory = self.upstream_path
        self.current_day = str(config.get('global', 'current_day'))
        self.output_directory = os.path.join(self.upstream_path, "outputs")
        self.reports_directory = os.path.join(self.upstream_path, "outputs")

        # Get night plan specific parameters
        self.max_solve_gap = config.getfloat('night', 'max_solve_gap')
        self.max_solve_time = config.getint('night', 'max_solve_time')
        self.show_gurobi_output = config.getboolean('night', 'show_gurobi_output')
        
        # Set up allocation file path from data section
        allocation_file_config = str(config.get('data', 'allocation_file'))
        if os.path.isabs(allocation_file_config):
            self.allocation_file = allocation_file_config
        else:
            self.allocation_file = os.path.join(self.semester_directory, allocation_file_config)
            
        # Set up backup file path
        filler_file_config = str(config.get('data', 'filler_file'))
        if os.path.isabs(filler_file_config):
            self.filler_file = filler_file_config
        else:
            self.filler_file = os.path.join(self.semester_directory, filler_file_config)

        # Set up custom file path from data section
        custom_file_config = str(config.get('data', 'custom_file'))
        if os.path.isabs(custom_file_config):
            self.custom_file = custom_file_config
        else:
            self.custom_file = os.path.join(self.semester_directory, custom_file_config)
        
        # Load SemesterPlanner from pickle file instead of creating new one
        config = ConfigParser()
        config.read(config_file)
        workdir = os.path.join(str(config.get('global', 'workdir')), "outputs")

        semester_planner_h5 = os.path.join(workdir, 'semester_planner.h5')
        self.semester_planner = SemesterPlanner.from_hdf5(semester_planner_h5)

        # The night plan inherits the queue from the semester plan to guarantee
        # they use identical telescope/instrument descriptions.
        self.queue = self.semester_planner.queue

        # Pull properties from SemesterPlanner for consistency
        self.semester_start_date = self.semester_planner.semester_start_date
        self.semester_length = self.semester_planner.semester_length
        self.all_dates_dict = self.semester_planner.all_dates_dict
        self.all_dates_array = self.semester_planner.all_dates_array
        self.today_starting_night = self.semester_planner.today_starting_night
        self.past_history = self.semester_planner.past_history
        self.slots_needed_for_exposure_dict = self.semester_planner.slots_needed_for_exposure_dict
        self.run_weather_loss = self.semester_planner.run_weather_loss

        # Compute tonight's allocation gaps: runs of unallocated slots (zeros) between allocated slots (ones)
        access_record = self.semester_planner.access_record
        tonight_index = self.today_starting_night
        slot_size = self.semester_planner.slot_size  # minutes per slot
        # is_alloc shape (ntargets, nnights, nslots); allocation is same for all targets
        tonight_allocated = access_record.is_alloc[0, tonight_index, :]  # 1D: 1=allocated, 0=not
        allocated = tonight_allocated.astype(np.int8)
        diff = np.diff(allocated)
        # Gap = run of zeros between ones (exclude leading/trailing zeros). diff==-1: 1->0 (gap start); diff==1: 0->1 (gap end)
        gap_start_slots = np.where(diff == -1)[0] + 1  # first zero slot of each potential gap
        gap_end_slots = np.where(diff == 1)[0]         # last zero slot before each 0->1 transition
        total_slots_in_night = len(allocated)
        total_allocated_slots = int(np.sum(allocated))
        total_nonallocated_slots = int(np.sum(1 - allocated))
        self.tonight_allocation_gaps = []
        for start_slot in gap_start_slots:
            # Pair with next gap_end that is >= start_slot (excludes trailing zeros)
            candidates = gap_end_slots[gap_end_slots >= start_slot]
            if len(candidates) > 0:
                end_slot = candidates[0]
                n_slots_in_gap = end_slot - start_slot + 1
                start_minutes = start_slot * slot_size
                end_minutes = (end_slot + 1) * slot_size  # end of last zero slot
                gap_length = n_slots_in_gap * slot_size
                gap_start_time = f"{int(start_minutes // 60):02d}:{int(start_minutes % 60):02d}"
                gap_stop_time = f"{int(end_minutes // 60):02d}:{int(end_minutes % 60):02d}"
                self.tonight_allocation_gaps.append({
                    'total_slots_in_night': total_slots_in_night,
                    'total_allocated_slots': total_allocated_slots,
                    'total_nonallocated_slots': total_nonallocated_slots,
                    'n_slots_in_gap': n_slots_in_gap,
                    'gap_start_slot': start_slot,
                    'gap_start_time': gap_start_time,
                    'gap_stop_slot': end_slot,
                    'gap_stop_time': gap_stop_time,
                    'gap_length': gap_length,
                })
        self.tonight_total_unallocated_slots = int(np.sum(1 - allocated))
        self.tonight_total_unallocated_minutes = self.tonight_total_unallocated_slots * slot_size
        self.tonight_gap_unallocated_slots = sum(g['n_slots_in_gap'] for g in self.tonight_allocation_gaps)
        self.tonight_gap_unallocated_minutes = self.tonight_gap_unallocated_slots * slot_size

    def run_ttp(self):
        """
        Prepare the TTP input dataframe by parsing the request_selected.csv file. Ensure data is in the correct format for TTP.
        Then run the TTP optimization to produce the solution which is then saved out as an hdf5 file.
        If no targets are selected, the function will gracefully return without running the TTP.

        Args:
            None

        Returns:
            None
        """

        observers_path = os.path.join(self.semester_directory, 'outputs/')
        if not os.path.isdir(observers_path):
            os.makedirs(observers_path)

        try:
            observation_start_time, observation_stop_time = get_nightly_times_from_allocation(self.allocation_file, self.current_day)
            total_time = np.round((observation_stop_time.jd-observation_start_time.jd)*24,3)
            print("Time in Night for Observations: " + str(total_time) + " hours.")
        except ValueError:
            print(f"No allocation times found for date {self.current_day}. Not running TTP. No night_planner.h5 file will be created.")
            return False

        selected_path = os.path.join(self.output_directory, 'request_selected.csv')
        if not os.path.exists(selected_path):
            raise FileNotFoundError(f"{selected_path} not found. Please run the scheduler first.")
        selected_df = pd.read_csv(selected_path)
        if len(selected_df) == 0:
            print(f"No targets found in {selected_path}. Not running TTP. No night_planner.pkl file will be created.")
            return

        first_available, last_available = self.get_first_last_indices(selected_df)
        selected_df['first_available'] = first_available
        selected_df['last_available'] = last_available

        # Webform data hygiene: coerce "None" strings to defaults.
        selected_df['n_intra_max'] = selected_df['n_intra_max'].replace('None', np.nan).fillna(1)
        selected_df['n_intra_min'] = selected_df['n_intra_min'].replace('None', np.nan).fillna(1)
        selected_df['tau_intra'] = selected_df['tau_intra'].replace('None', np.nan).fillna(0.0)
        selected_df['jmag'] = selected_df['jmag'].replace('None', np.nan).fillna(0.0)
        selected_df['Vmag'] = selected_df['Vmag'].replace('None', np.nan).fillna(0.0)
        selected_df['pmra'] = selected_df['pmra'].replace('None', np.nan).fillna(0.0)
        selected_df['pmdec'] = selected_df['pmdec'].replace('None', np.nan).fillna(0.0)
        selected_df['epoch'] = selected_df['epoch'].replace('None', np.nan).fillna(0.0)

        # Build TTP input frame directly in AstroQ-native vocabulary
        # (TTPModel.REQUIRED_COLUMNS). No CSV intermediate; we hand the
        # DataFrame to TTPModel directly.
        to_ttp = selected_df[[
            "unique_id", "ra", "dec", "exptime",
            "n_exp", "n_intra_max", "tau_intra",
            "first_available", "last_available",
        ]].copy()
        to_ttp["priority"] = 10  # default per-request priority

        # Synthetic "Gap N" rows block out unallocated gaps inside the night.
        # The dummy unique_id prefix "Gap " is matched verbatim by
        # `drop_gap_rows` below to scrub these from the user-facing outputs.
        if len(self.tonight_allocation_gaps) > 0:
            avg_ra = selected_df["ra"].mean()
            avg_dec = selected_df["dec"].mean()
            tonight_date = self.current_day
            gap_rows = []
            for i, gap in enumerate(self.tonight_allocation_gaps, start=1):
                first_av = f"{tonight_date} {gap['gap_start_time']}"
                last_av = f"{tonight_date} {gap['gap_stop_time']}"
                exposure_time_sec = (gap['gap_length'] - self.semester_planner.slot_size) * 60
                gap_rows.append({
                    "unique_id": f"Gap {i}",
                    "ra": avg_ra,
                    "dec": avg_dec,
                    "exptime": exposure_time_sec,
                    "n_exp": 1,
                    "n_intra_max": 1,
                    "tau_intra": 0.0,
                    "priority": 20,  # high enough to guarantee selection
                    "first_available": first_av,
                    "last_available": last_av,
                })
            to_ttp = pd.concat([to_ttp, pd.DataFrame(gap_rows)], ignore_index=True)

        tm = model.TTPModel(
            requests_frame=to_ttp,
            night_start=observation_start_time,
            night_end=observation_stop_time,
            outdir=observers_path,
            observer=self.queue.observer,
            slew_rate=self.queue.slew_rate,
            wrap_limit=self.queue.wrap_limit,
            readout_time=self.queue.readout_time,
            n_slots=self.queue.nSlots,
            inaccessible_zones=self.queue.inaccessible_zones,
        )
        tm.build_nodes()
        tm.build_slew_tensor()
        tm.build_model()
        tm.model.params.TimeLimit = self.max_solve_time
        tm.model.params.MIPGap = self.max_solve_gap
        tm.model.update()
        tm.run_model()

        del tm.model             # remove attribute so object is hdf5 compatable

        # Compute gap stats BEFORE scrubbing (for adjusted TTP statistics)
        gap_exposure_min = 0.0
        gap_count = 0
        if len(self.tonight_allocation_gaps) > 0:
            plotly_exp = tm.plotly.get('Total Exp Time (min)', [])
            for i, name in enumerate(tm.plotly.get('Starname', [])):
                if str(name).startswith('Gap '):
                    gap_count += 1
                    gap_exposure_min += float(plotly_exp[i]) if i < len(plotly_exp) else 0
            if tm.extras and tm.extras.get('Starname'):
                extras_exp = tm.extras.get('Total Exp Time (min)', [])
                for j, name in enumerate(tm.extras['Starname']):
                    if str(name).startswith('Gap '):
                        gap_count += 1
                        gap_exposure_min += float(extras_exp[j]) if j < len(extras_exp) else 0
        gap_total_min = self.tonight_gap_unallocated_minutes if len(self.tonight_allocation_gaps) > 0 else 0.0
        n_gap_targets = len(self.tonight_allocation_gaps)

        # Scrub dummy "Gap N" rows so they never appear in user-facing outputs.
        def drop_gap_rows(d):
            keep = [i for i, s in enumerate(d['Starname']) if not str(s).startswith('Gap ')]
            return {k: [v[i] for i in keep] for k, v in d.items()}
        tm.plotly = drop_gap_rows(tm.plotly)
        if tm.extras is not None:
            if isinstance(tm.extras, pd.DataFrame):
                tm.extras = tm.extras[
                    ~tm.extras['Starname'].astype(str).str.startswith('Gap ')
                ]
            elif len(tm.extras.get('Starname', [])) > 0:
                tm.extras = drop_gap_rows(tm.extras)
        if getattr(tm, 'schedule', None) is not None and isinstance(tm.schedule, dict) and 'Starname' in tm.schedule:
            tm.schedule = drop_gap_rows(tm.schedule)
        # Drop "Gap N" rows from requests_frame too; downstream consumers iterate it.
        if getattr(tm, 'requests_frame', None) is not None:
            mask = ~tm.requests_frame['unique_id'].astype(str).str.startswith('Gap ')
            tm.requests_frame = tm.requests_frame[mask].reset_index(drop=True)

        # Update TTP stats to exclude Gap observations (observing duration, exposing, idle)
        if gap_total_min > 0 or gap_exposure_min > 0:
            tm.dur = max(0, tm.dur - gap_total_min)
            tm.time_exposing = max(0, tm.time_exposing - gap_exposure_min)
            tm.time_idle = max(0, tm.dur - tm.time_exposing - tm.time_slewing)
            tm.num_scheduled = tm.num_scheduled - gap_count
            # Re-print and overwrite TTPstatistics.txt with gap-adjusted stats
            ttp_stats_path = os.path.join(observers_path, 'TTPstatistics.txt')
            with open(ttp_stats_path, 'w') as f:
                f.write("Stats for TTP Solution (Gap observations excluded)\n")
                f.write("------------------------------------\n")
                f.write(f'    Model ran for {tm.solve_time:.2f} seconds\n')
                f.write(f'     Observations Requested: {tm.N - 2 - n_gap_targets}\n')
                f.write(f'     Observations Scheduled: {tm.num_scheduled}\n')
                f.write("------------------------------------\n")
                f.write(f'   Observing Duration (min): {tm.dur:.2f}\n')
                f.write(f'  Time Spent Exposing (min): {tm.time_exposing:.2f}\n')
                f.write(f'      Time Spent Idle (min): {tm.time_idle:.2f}\n')
                f.write(f'   Time Spent Slewing (min): {tm.time_slewing:.2f}\n')
                f.write("------------------------------------\n")
            print('\n------------------------------------')
            print(' (Gap observations excluded from stats)')
            print('------------------------------------')
            print(f'     Observations Requested: {tm.N - 2 - n_gap_targets}')
            print(f'     Observations Scheduled: {tm.num_scheduled}')
            print('------------------------------------')
            print(f'   Observing Duration (min): {tm.dur:.2f}')
            print(f'  Time Spent Exposing (min): {tm.time_exposing:.2f}')
            print(f'      Time Spent Idle (min): {tm.time_idle:.2f}')
            print(f'   Time Spent Slewing (min): {tm.time_slewing:.2f}')
            print('------------------------------------')

        # add human readable starname to the solution so that it can be used in the plotting functions
        id_to_name = dict(zip(selected_df['unique_id'], selected_df['starname']))
        tm.plotly['human_starname'] = [
            id_to_name.get(uid, "NO MATCHING NAME") for uid in tm.plotly['Starname']
        ]
        self.solution = [tm]

        tm.plotly['UTC Start Time'] = [0]*len(tm.plotly['Start Exposure'])
        for i in range(len(tm.plotly['Start Exposure'])):
            tm.plotly['UTC Start Time'][i] = str(TimeDelta(tm.plotly['Start Exposure'][i]*60,format='sec') + observation_start_time)[11:16]
        numeric_columns = ['Start Exposure', 'First Available', 'Last Available', 'Minutes the from Start of the Night']
        for col in numeric_columns:
            if col in tm.plotly:
                tm.plotly[col] = np.round(np.array(tm.plotly[col]), 2).tolist()

        # Convert tm.plotly to a DataFrame for easier handling
        observe_order_file = os.path.join(observers_path, f"ObserveOrder_{self.current_day}.txt")
        plotly_df = pd.DataFrame(tm.plotly)
        use_starnames = []
        use_star_ids = []
        use_start_exposures = []
        for i in range(len(plotly_df)):
            adjusted_timestamp = TimeDelta(plotly_df['Start Exposure'].iloc[i]*60,format='sec') + observation_start_time
            use_start_exposures.append(str(adjusted_timestamp)[11:16])
            use_starnames.append(selected_df[selected_df['unique_id'] == plotly_df['Starname'].iloc[i]]['starname'].iloc[0])
            use_star_ids.append(str(plotly_df['Starname'].iloc[i]))

        extras_df = pd.DataFrame(tm.extras)
        for j in range(len(extras_df)):
            use_start_exposures.append('24:00')
            use_star_ids.append(str(extras_df['Starname'].iloc[j]))
            use_starnames.append(selected_df[selected_df['unique_id'] == extras_df['Starname'].iloc[j]]['starname'].iloc[0])
        use_frame = pd.DataFrame({'unique_id': use_star_ids, 'Target': use_starnames, 'StartExposure': use_start_exposures})
        use_frame.to_csv(observe_order_file, index=False)

        self.queue.write_starlist(selected_df, tm.plotly, observation_start_time, tm.extras,
                            [], str(self.current_day), observers_path,
                            all_active_requests=self.semester_planner.requests_frame,
                            past_history=self.past_history)
        print("The optimal path through the sky for the selected stars is found. Clear skies!")

        return True

    def get_first_last_indices(self, selected_df):
        """
        Get the first and last available time slots for each target in selected_df.
        
        Args:
            selected_df (pd.DataFrame): DataFrame containing selected targets with unique_id column
            
        Returns:
            first_available_list (list) - Lists of time strings in HH:MM format for each target's first available slot
            last_available_list (list) - Lists of time strings in HH:MM format for each target's last available slot
        """

        # Get tonight's index from the all_dates_dict
        tonight_index = self.all_dates_dict[self.current_day]
        
        # Get the access record from semester planner
        access_record = self.semester_planner.access_record
        
        # Create mapping from unique_id to target index in the access record
        # The access record was created from the original requests_frame, so we need to map
        # selected_df targets back to their indices in the original requests_frame
        target_to_index = {}
        for idx, row in self.semester_planner.requests_frame.iterrows():
            target_to_index[row['unique_id']] = idx
        
        # Initialize the new columns
        first_available = []
        last_available = []
        
        # For each target in selected_df, find first and last available slots tonight
        for _, row in selected_df.iterrows():
            target_id = row['unique_id']
            
            # Get the target's index in the access record
            if target_id in target_to_index:
                target_idx = target_to_index[target_id]
                
                # Get tonight's observability array for this target (shape: nslots)
                tonight_observable = access_record.is_observable[target_idx, tonight_index, :]
                
                # Find first and last True indices
                true_indices = np.where(tonight_observable)[0]
                
                if len(true_indices) > 0:
                    first_slot = true_indices[0]
                    last_slot = true_indices[-1]
                    
                    # Convert slot indices to time strings (assuming slot_size is in minutes)
                    first_time_minutes = first_slot * self.semester_planner.slot_size
                    last_time_minutes = last_slot * self.semester_planner.slot_size
                    
                    first_hour = first_time_minutes // 60
                    first_minute = first_time_minutes % 60
                    last_hour = last_time_minutes // 60
                    last_minute = last_time_minutes % 60
                    
                    first_available.append(f"{self.current_day} {first_hour:02d}:{first_minute:02d}")
                    last_available.append(f"{self.current_day} {last_hour:02d}:{last_minute:02d}")
                else:
                    # No available slots tonight, use dummy values so TTP doesn't break
                    last_hour = 23
                    last_minute = 59
                    first_available.append(f"{self.current_day} {last_hour:02d}:{last_minute:02d}")
                    last_available.append(f"{self.current_day} {last_hour:02d}:{last_minute:02d}")
            else:
                # Target not found in original requests_frame
                last_hour = 23
                last_minute = 59
                first_available.append(f"{self.current_day} {last_hour:02d}:{last_minute:02d}")
                last_available.append(f"{self.current_day} {last_hour:02d}:{last_minute:02d}")
        
        return first_available, last_available

    def to_hdf5(self, hdf5_path=None):
        """
        Save the NightPlanner object to an HDF5 file.
        
        Args:
            hdf5_path (str, optional): Path to save the HDF5 file. 
                                      If None, saves to output_directory/night_planner.h5
        """
        import h5py
        import json
        
        if hdf5_path is None:
            hdf5_path = os.path.join(self.output_directory, 'night_planner.h5')
        # Remove existing file if it exists
        if os.path.exists(hdf5_path):
            os.remove(hdf5_path)
        
        # Define serialization mappings
        # Format: (hdf5_key, object_path, data_type, conversion_func)
        # data_type: 'scalar', 'string', 'array', 'time', 'dict_json', 'dataframe', 'stars'
        # object_path: attribute path like 'solution.plotly' or 'self.upstream_path'
        
        # NightPlanner scalar/string attributes
        nightplanner_attrs = [
            ('upstream_path', 'self.upstream_path', 'string', None),
            ('semester_directory', 'self.semester_directory', 'string', None),
            ('current_day', 'self.current_day', 'string', None),
            ('output_directory', 'self.output_directory', 'string', None),
            ('reports_directory', 'self.reports_directory', 'string', None),
            ('max_solve_gap', 'self.max_solve_gap', 'scalar', None),
            ('max_solve_time', 'self.max_solve_time', 'scalar', None),
            ('show_gurobi_output', 'self.show_gurobi_output', 'scalar', None),
            ('allocation_file', 'self.allocation_file', 'string', None),
            ('filler_file', 'self.filler_file', 'string', None),
            ('custom_file', 'self.custom_file', 'string', None),
        ]
        
        # Solution object attributes. Schema v2 (NIGHT_PLANNER_H5_SCHEMA):
        # the legacy per-Star `solution_star_*` datasets are replaced by a
        # single `solution_requests_frame` table on disk; the in-memory
        # `coord` column is stripped before save and rebuilt on load.
        solution = self.solution[0]
        solution_attrs = [
            ('solution_plotly_json', 'solution.plotly', 'dict_json', None),
            ('solution_times_jd', 'solution.times', 'time_list', None),
            ('night_start_jd', 'solution.night_start', 'time', None),
            ('night_end_jd', 'solution.night_end', 'time', None),
            ('solution_schedule_json', 'solution.schedule', 'dict_json', None),
            ('solution_az_path', 'solution.az_path', 'array', None),
            ('solution_alt_path', 'solution.alt_path', 'array', None),
        ]
        
        # Save solution.extras first (special case - DataFrame or dict)
        extras_is_dict = isinstance(solution.extras, dict)
        if isinstance(solution.extras, pd.DataFrame):
            # Save DataFrame (even if empty)
            extras_df = solution.extras
        elif isinstance(solution.extras, dict):
            # Convert dict to DataFrame (handles empty dicts with empty lists)
            # pd.DataFrame() creates empty DataFrame with columns when all lists are empty
            extras_df = pd.DataFrame(solution.extras)
        
        if extras_df.empty:
            extras_df.to_hdf(hdf5_path, key='solution_extras', mode='a', format='fixed')
        else:
            extras_df.to_hdf(hdf5_path, key='solution_extras', mode='a', format='table')

        # Persist solution.requests_frame, stripping the derived `coord` column
        # (SkyCoord pickling inside HDF5 is fragile across astropy versions;
        # we rebuild from ra/dec on load).
        rf = solution.requests_frame.drop(columns=['coord'], errors='ignore')
        if rf.empty:
            rf.to_hdf(hdf5_path, key='solution_requests_frame', mode='a', format='fixed')
        else:
            rf.to_hdf(hdf5_path, key='solution_requests_frame', mode='a', format='table')

        with h5py.File(hdf5_path, 'a') as f:
            f.attrs['schema_version'] = NIGHT_PLANNER_H5_SCHEMA
            f.attrs['extras_was_dict'] = extras_is_dict
            
            # Save solution attributes
            for hdf5_key, obj_path, data_type, extra in solution_attrs:
                obj = solution
                for attr in obj_path.split('.')[1:]:  # Skip 'solution' part
                    obj = getattr(obj, attr)
                
                if data_type == 'dict_json':
                    # Convert dict with arrays/lists to JSON-serializable format (native Python types)
                    def _to_native(x):
                        if isinstance(x, np.ndarray):
                            return _to_native(x.tolist())
                        if isinstance(x, (list, tuple)):
                            return [_to_native(v) for v in x]
                        if isinstance(x, dict):
                            return {k: _to_native(v) for k, v in x.items()}
                        if isinstance(x, (np.integer, np.int64, np.int32)):
                            return int(x)
                        if isinstance(x, (np.floating, np.float64, np.float32)):
                            return float(x)
                        if isinstance(x, (np.bool_, bool)):
                            return bool(x)
                        return x
                    serializable = {k: _to_native(v) for k, v in obj.items()}
                    f.attrs[hdf5_key] = json.dumps(serializable)
                
                elif data_type == 'time_list':
                    # Convert list of Time objects to array of JD
                    times_jd = np.array([t.jd for t in obj])
                    f.create_dataset(hdf5_key, data=times_jd)
                
                elif data_type == 'time':
                    # Convert Time object to JD scalar
                    f.attrs[hdf5_key] = obj.jd
                
                elif data_type == 'array':
                    f.create_dataset(hdf5_key, data=np.array(obj))
            
            # Save NightPlanner attributes
            for hdf5_key, obj_path, data_type, _ in nightplanner_attrs:
                attr_name = obj_path.split('.')[-1]
                value = getattr(self, attr_name)
                f.attrs[hdf5_key] = value
            
            # Save path to semester_planner.h5 file
            semester_planner_h5_path = os.path.join(self.output_directory, 'semester_planner.h5')
            f.attrs['semester_planner_h5_path'] = semester_planner_h5_path
        
        return hdf5_path

    @classmethod
    def from_hdf5(cls, hdf5_path):
        """
        Load a NightPlanner object from an HDF5 file.
        
        Args:
            hdf5_path (str): Path to the HDF5 file
            
        Returns:
            NightPlanner: Reconstructed NightPlanner object
        """
        import h5py
        import json

        instance = cls.__new__(cls)
        
        # Define deserialization mappings (inverse of to_hdf5)
        # Format: (hdf5_key, attribute_name, data_type, conversion_func)
        nightplanner_attrs = [
            ('upstream_path', 'upstream_path', 'string', None),
            ('semester_directory', 'semester_directory', 'string', None),
            ('current_day', 'current_day', 'string', None),
            ('output_directory', 'output_directory', 'string', None),
            ('reports_directory', 'reports_directory', 'string', None),
            ('max_solve_gap', 'max_solve_gap', 'scalar', None),
            ('max_solve_time', 'max_solve_time', 'scalar', None),
            ('show_gurobi_output', 'show_gurobi_output', 'scalar', None),
            ('allocation_file', 'allocation_file', 'string', None),
            ('filler_file', 'filler_file', 'string', None),
            ('custom_file', 'custom_file', 'string', None),
        ]
        
        solution_attrs = [
            ('solution_plotly_json', 'plotly', 'dict_json', None),
            ('solution_times_jd', 'times', 'time_list', None),
            ('night_start_jd', 'night_start', 'time', None),
            ('night_end_jd', 'night_end', 'time', None),
            ('solution_schedule_json', 'schedule', 'dict_json', None),
            ('solution_az_path', 'az_path', 'array', None),
            ('solution_alt_path', 'alt_path', 'array', None),
        ]

        with h5py.File(hdf5_path, 'r') as f:
            schema = int(f.attrs.get('schema_version', 0))
            if schema != NIGHT_PLANNER_H5_SCHEMA:
                raise ValueError(
                    f"night_planner.h5 schema_version={schema} but this build "
                    f"expects {NIGHT_PLANNER_H5_SCHEMA}. Re-run plan-night to regenerate."
                )
            if 'solution_extras' not in f:
                raise AttributeError("solution.extras not found in HDF5 file")

        solution_extras_df = pd.read_hdf(hdf5_path, key='solution_extras')
        solution_requests_frame = pd.read_hdf(hdf5_path, key='solution_requests_frame')
        
        # Reconstruct solution object
        class SolutionContainer:
            pass
        
        solution = SolutionContainer()
        
        with h5py.File(hdf5_path, 'r') as f:
            # Load NightPlanner attributes
            for hdf5_key, attr_name, data_type, _ in nightplanner_attrs:
                setattr(instance, attr_name, f.attrs[hdf5_key])
            
            # Load semester_planner
            semester_planner_h5_path = f.attrs['semester_planner_h5_path']
            if not os.path.exists(semester_planner_h5_path):
                raise FileNotFoundError(f"semester_planner.h5 not found at {semester_planner_h5_path}")
            instance.semester_planner = SemesterPlanner.from_hdf5(semester_planner_h5_path)
            
            # Pull properties from SemesterPlanner
            instance.semester_start_date = instance.semester_planner.semester_start_date
            instance.semester_length = instance.semester_planner.semester_length
            instance.all_dates_dict = instance.semester_planner.all_dates_dict
            instance.all_dates_array = instance.semester_planner.all_dates_array
            instance.today_starting_night = instance.semester_planner.today_starting_night
            instance.past_history = instance.semester_planner.past_history
            instance.slots_needed_for_exposure_dict = instance.semester_planner.slots_needed_for_exposure_dict
            instance.run_weather_loss = instance.semester_planner.run_weather_loss
            instance.queue = instance.semester_planner.queue
            
            # Load solution attributes
            for hdf5_key, attr_name, data_type, extra in solution_attrs:
                if data_type == 'dict_json':
                    data = json.loads(f.attrs[hdf5_key])
                    # Convert lists back to numpy arrays
                    restored = {}
                    for key, value in data.items():
                        if isinstance(value, list):
                            restored[key] = np.array(value)
                        else:
                            restored[key] = value
                    setattr(solution, attr_name, restored)
                
                elif data_type == 'time_list':
                    times_jd = f[hdf5_key][:]
                    setattr(solution, attr_name, [Time(jd, format='jd') for jd in times_jd])
                
                elif data_type == 'time':
                    jd = f.attrs[hdf5_key]
                    setattr(solution, attr_name, Time(jd, format='jd'))
                
                elif data_type == 'array':
                    data = f[hdf5_key][:]
                    setattr(solution, attr_name, data)
            
            # Rebuild solution.requests_frame; derive the `coord` cache from
            # ra/dec to avoid serializing SkyCoord (astropy-version fragile).
            solution.requests_frame = solution_requests_frame.reset_index(drop=True).copy()
            solution.requests_frame['coord'] = list(SkyCoord(
                solution.requests_frame.ra.values * u.deg,
                solution.requests_frame.dec.values * u.deg,
                frame='icrs',
            ))

            # Re-attach queue-derived primitives that the plotter reads.
            # These are NOT serialized to h5 (the queue is the canonical source).
            solution.observer = instance.queue.observer
            solution.slew_rate = instance.queue.slew_rate
            solution.wrap_limit = instance.queue.wrap_limit
            solution.readout_time = instance.queue.readout_time
            solution.n_slots = instance.queue.nSlots
            solution.inaccessible_zones = instance.queue.inaccessible_zones
        
        # Load solution.extras (convert back to dict if needed)
        with h5py.File(hdf5_path, 'r') as f:
            extras_was_dict = f.attrs['extras_was_dict']
        
        if extras_was_dict:
            solution.extras = solution_extras_df.to_dict('list')
        else:
            solution.extras = solution_extras_df

        # Belt-and-suspenders gap scrub (run-time scrubber may have already
        # cleared these; harmless on already-clean data).
        def _drop_gap_from_dict(d):
            if not isinstance(d, dict) or 'Starname' not in d:
                return d
            keep = [i for i, s in enumerate(d['Starname']) if not str(s).startswith('Gap ')]
            return {k: ([v[i] for i in keep] if isinstance(v, (list, np.ndarray)) else v) for k, v in d.items()}
        solution.plotly = _drop_gap_from_dict(solution.plotly)
        solution.schedule = _drop_gap_from_dict(solution.schedule)
        if solution.requests_frame is not None and 'unique_id' in solution.requests_frame.columns:
            solution.requests_frame = solution.requests_frame[
                ~solution.requests_frame['unique_id'].astype(str).str.startswith('Gap ')
            ].reset_index(drop=True)
        if solution.extras is not None:
            if isinstance(solution.extras, pd.DataFrame):
                solution.extras = solution.extras[
                    ~solution.extras['Starname'].astype(str).str.startswith('Gap ')
                ]
            elif isinstance(solution.extras, dict) and solution.extras.get('Starname'):
                solution.extras = _drop_gap_from_dict(solution.extras)
        
        instance.solution = [solution]
        
        return instance

def get_nightly_times_from_allocation(allocation_file, current_day):
    """
    Extract start and stop times for a specific date from allocation.csv.
    
    Args:
        allocation_file (str): path to the allocation file
        current_day (str): the date to look for in YYYY-MM-DD format
        
    Returns:
       start_time (Time object): the start time of the allocation for the current day
       stop_time (Time object): the stop time of the allocation for the current day
    """
    allocated_times_frame = pd.read_csv(allocation_file)
    allocated_times_frame['start'] = allocated_times_frame['start'].apply(Time)
    allocated_times_frame['stop'] = allocated_times_frame['stop'].apply(Time)
    
    # Filter for the current day
    current_day_str = str(current_day)
    day_allocations = []
    for _, row in allocated_times_frame.iterrows():
        start_datetime = str(row['start'])[:10]  # Extract date part (YYYY-MM-DD)
        if start_datetime == current_day_str:
            day_allocations.append(row)
    
    if not day_allocations:
        raise ValueError(f"No allocation found for date {current_day_str}")
    
    # For multiple allocations on the same day, use the earliest start and latest stop
    start_times = [row['start'] for row in day_allocations]
    stop_times = [row['stop'] for row in day_allocations]
    
    earliest_start = min(start_times)
    latest_stop = max(stop_times)
    
    return earliest_start, latest_stop

        