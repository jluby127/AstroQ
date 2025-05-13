"""
Module that defines the RequestSet class. This encodes all the astronomy information into one place, and
prepares it for the Scheduler object. Useful as a checkpoint to the code.

"""
import sys
import time
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import json
from configparser import ConfigParser

import kpfcc.maps as mp

class RequestSet(object):
    """
    Request Set

    Specifies set of request to send to schedule. May be saved and restored from json
    """
    def __init__(self, meta, strategy, observable):
        self.meta = meta
        # enforce order
        strategy_cols = 'id t_visit n_intra_min n_intra_max tau_intra n_inter_max tau_inter'.split()
        strategy = strategy[strategy_cols]
        self.strategy = strategy
        self.observability = observable

    def __str__(self, n=10):
        s = "# Request Set # \n"
        s += "## Meta Data ## \n"
        s += '\n'*2
        s += self.meta.to_string()
        s += '\n'*2
        s += "## Strategy ##\n"
        s += self.strategy.head(n).to_string()
        s += '\n'*2
        s += "## Observable ##\n"
        s += self.observability.head(n).to_string()
        return s

    def to_json(self, fn):
        """
        Save request set to json
        """
        data = {}
        if isinstance(self.meta, dict):
            data['meta'] = self.meta
        else:
            data['meta'] = self.meta.to_dict()
        data['strategy'] = self.strategy.to_dict()
        data['observable'] = self.observability.to_dict()

        with open(fn, "w") as f:
            json.dump(data,f,indent=4) # Pretty-printed for readability

def read_json(fn):
    """
    Read request set from json
    """
    with open(fn, "r") as f:
        data = json.load(f)

    meta = pd.Series(data['meta'])
    strategy = pd.DataFrame(data['strategy'])
    observable = pd.DataFrame(data['observable'])
    rs = RequestSet(meta, strategy, observable)
    return rs

def define_indices_for_requests(manager):
    """
    Using the dictionary of indices where each request is available, define a dataframe for which
    we will use to cut/filter/merge r,d,s tuples
    """
    # Define the tuples of request and available slot for each request.
    # This becomes the grid over which the Gurobi variables are defined.
    # Now, slots that were never possible for scheduling are not included in the model.

    available_indices_for_request = mp.produce_ultimate_map(manager)

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

def build_meta(config_file):

    config = ConfigParser()
    config.read(config_file)
    daily_starting_time = str(config.get('other', 'daily_starting_time'))
    current_day = str(config.get('required', 'current_day'))
    slot_size = int(config.get('other', 'slot_size'))

    meta = {"s1_time":daily_starting_time, "d1_date":current_day, "slot_duration":slot_size}
    return meta

def cull_from_weather(request_set, weather_loss_map):
    request_set.original_observability = request_set.observability
    request_set.observability = request_set.observability[~request_set.observability['d'].isin(weather_loss_map)]
    return request_set

def convert_slot_to_quarter(d, s, twilight_map_remaining_2D_d):
    '''
    Determine the slot numbers within the night that breaks the night into "equal" length quarters
    Take extra precaution when the total number of slots between twilight times is not easily
    divisable by 4.
    '''

    n_available_slots_in_quarter_tonight = int(np.sum(twilight_map_remaining_2D_d)/4)
    extra_slots = np.sum(twilight_map_remaining_2D_d)%4

    edge_slot = np.argmax(twilight_map_remaining_2D_d)
    split_2nd3rd = int(len(twilight_map_remaining_2D_d)/2) - 1
    split_1st2nd = int((split_2nd3rd - edge_slot)/2) + edge_slot
    split_3rd4th = int((len(twilight_map_remaining_2D_d) - edge_slot - split_2nd3rd)/2) + split_2nd3rd

    if s < split_1st2nd:
        q = 0
    elif s >= split_1st2nd and s < split_2nd3rd:
        q = 1
    elif s >= split_2nd3rd and s < split_3rd4th:
        q = 2
    elif s >= split_3rd4th:
        q = 3
    else:
        q = 100
        print("Houston, we've had a problem.")

    return q
