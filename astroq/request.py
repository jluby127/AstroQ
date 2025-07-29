"""
Module that defines the RequestSet class. This encodes all the astronomy information into one place, and
prepares it for the Scheduler object. Useful as a checkpoint to the code.

"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import json
from configparser import ConfigParser

import astroq.access as ac
import astroq.management as mn

class RequestSet(object):
    """
    Request Set

    Specifies set of request to send to schedule. May be saved and restored from json
    """
    def __init__(self, config_file):
        """
        Initialize RequestSet from a config file.
        
        Args:
            config_file: Path to the configuration file
        """
        # Create manager and run admin to set up all necessary data
        manager = mn.data_admin(config_file)
        manager.run_admin()
        
        # Build meta data from config
        config = ConfigParser()
        config.read(config_file)
        daily_starting_time = str(config.get('other', 'daily_starting_time'))
        current_day = str(config.get('required', 'current_day'))
        slot_size = int(config.get('other', 'slot_size'))
        self.meta = {"s1_time":daily_starting_time, "d1_date":current_day, "slot_duration":slot_size}
        
        # Build strategy and observable data
        strategy, observable = self._define_indices_for_requests(manager)
        
        # enforce order
        strategy_cols = 'id t_visit n_intra_min n_intra_max tau_intra n_inter_max tau_inter'.split()
        strategy = strategy[strategy_cols]
        self.strategy = strategy
        self.observability = observable

    def _define_indices_for_requests(self, manager):
        """
        Using the dictionary of indices where each request is available, define a dataframe for which
        we will use to cut/filter/merge r,d,s tuples
        """
        # Create Access object with required parameters
        access_obj = ac.Access(
            semester_start_date=manager.semester_start_date,
            semester_length=manager.semester_length,
            slot_size=manager.slot_size,
            observatory=manager.observatory,
            current_day=manager.current_day,
            all_dates_dict=manager.all_dates_dict,
            custom_file=manager.custom_file,
            allocation_file=manager.allocation_file,
            past_history=manager.past_history,
            today_starting_night=manager.today_starting_night,
            slots_needed_for_exposure_dict=manager.slots_needed_for_exposure_dict
        )
        access = access_obj.produce_ultimate_map(manager.requests_frame)
        # Convert to list of available indices
        available_indices_for_request = access_obj.extract_available_indices_from_record(access, manager.requests_frame)

        observability_keys = []
        strategy_keys = []
        for n,row in manager.requests_frame.iterrows():
            name = row['starname']
            if name in list(available_indices_for_request.keys()):
                max_n_visits = int(row['n_intra_max'])
                min_n_visits = int(row['n_intra_min'])
                intra = int(row['tau_intra'])
                nnights = int(row['n_inter_max'])
                inter = int(row['tau_inter'])
                slots_needed = manager.slots_needed_for_exposure_dict[name]
                strategy_keys.append([name, slots_needed, min_n_visits, max_n_visits, intra, nnights, inter])
                for d in range(len(available_indices_for_request[name])):
                    for s in available_indices_for_request[name][d]:
                        observability_keys.append((name, d, s))
        strategy = pd.DataFrame(strategy_keys, columns =['id', 't_visit', 'n_intra_min', 'n_intra_max',
                                                         'tau_intra', 'n_inter_max', 'tau_inter'])
        observability = pd.DataFrame(observability_keys, columns =['id', 'd', 's'])
        return strategy, observability

def cull_from_weather(request_set, weather_loss_map):
    request_set.original_observability = request_set.observability
    request_set.observability = request_set.observability[~request_set.observability['d'].isin(weather_loss_map)]
    return request_set
