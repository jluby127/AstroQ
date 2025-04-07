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

import kpfcc.maps as mp

class RequestSet(object):
    """A RequestSet object, in which we define the set of possible observations.

       Args:
            manager (obj): an instance of the manager
            available_indices_for_request (dictionary): keys are the starnames and values are a 1D array
                                                      the indices where available_slots_for_request is 1.
    """
    def __init__(self, manager):

        print("Building the RequestSet.")
        self.manager = manager

        self.Aset = None
        self.Aframe = None
        self.Wset = None
        self.schedulable_requests = None
        self.first_slot_of_night_time = None
        self.first_day_of_semester = None
        self.size_of_slots = None


    def define_indices_for_requests(self):
        """
        Using the dictionary of indices where each request is available, define a dataframe for which
        we will use to cut/filter/merge r,d,s tuples
        """
        # Define the tuples of request and available slot for each request.
        # This becomes the grid over which the Gurobi variables are defined.
        # Now, slots that were never possible for scheduling are not included in the model.

        available_indices_for_request = mp.produce_ultimate_map(manager, manager.allocation_map_1D.flatten(),
                                                                manager.twilight_map_remaining_2D.flatten())

        Aset = []
        Wset = []
        Aframe_keys = []
        for n,row in self.manager.requests_frame.iterrows():
            name = row['Starname']
            if name in list(available_indices_for_request.keys()):
                max_n_visits = int(row['Desired Visits per Night'])
                min_n_visits = int(row['Accepted Visits per Night'])
                intra = int(row['Minimum Intra-Night Cadence'])
                inter = int(row['Minimum Inter-Night Cadence'])
                slots_needed = self.manager.slots_needed_for_exposure_dict[name]
                for d in range(len(available_indices_for_request[name])):
                    Wset.append((name, d))
                    for s in available_indices_for_request[name][d]:
                        Aset.append((name, d, s))
                        Aframe_keys.append([name, d, s, slots_needed, max_n_visits, min_n_visits, intra, inter])

        Aframe = pd.DataFrame(Aframe_keys, columns =['r', 'd', 's', 'e', 'maxv', 'minv', 'tra', 'ter'])
        schedulable_requests = list(Aframe['r'].unique())
        for name in list(self.manager.requests_frame['Starname']):
            if name not in schedulable_requests:
                print("WARNING: Target " + name + " has no valid day/slot pairs and therefore is effectively removed from the model.")
        # duplicate columns for easy indexing later
        Aframe['rr'] = Aframe['r']
        Aframe['dd'] = Aframe['d']
        Aframe['ss'] = Aframe['s']

        self.Aframe = Aframe
        self.Aset = Aset
        self.Wset = Wset
        self.schedulable_requests = schedulable_requests

    def read_from_json(self, filename):
    """Define a RequestSet object's attributes from reading in a previously saved json file.

    Args:
        - filename (str) = the path and file name to the request_set's json file.
    Returns:
        None
    """
        with open(filename, 'r') as file:
            data = json.load(file)

        self.Aframe = pd.DataFrame(data['cadence'])

        Aset_frame = pd.DataFrame(data['access'])
        self.Aset = [tuple(x) for x in Aset_frame[['id', 'd', 's']].values]

        unique_d_for_each_id = Aset.groupby('id')['d'].unique()
        self.schedulable_requests = list(unique_d_for_each_id.keys())
        Wset = []
        for starname in schedulable_requests:
            Wset.append([(starname, idx) for idx in list(unique_d_for_each_id[starname])])
        self.Wset = [item for sublist in Wset for item in sublist]

        self.first_slot_of_night_time = data['meta']['s1_time']
        self.first_day_of_semester = data['meta']['d1_date']
        self.slot_duration = data['meta']['slot_duration']

    def write_to_json(self, filename):
    """Define a RequestSet object's attributes from reading in a previously saved json file.

    Args:
        - filename (str) = the path and file name to the request_set's json file.
    Returns:
        None
    """

    # note that as it currently stands:
    #       - all info for "meta" is in the manager (need to add in first slot time to config file)
    #       - self.Aframe = "strategy" is a pandas dataframe
    #       - self.Aset = "observable" = is a list of tuples [(starname1, d1, s1), (starname2, d2, s2),...]
    # would need to write these out as dictionaries appropriately

    # Erik's code to write out goes here
    return


# Example usage:
rs = RequestSet(manager)

#option 1: produce "strategy" and "observable"
rs.define_indices_for_requests()
rs.write_to_json("path_to_file_here.json")

# option 2: read in from previously produced file
rs.read_from_json("path_to_file_here.json")

# Either way, call Scheduler
scheduler = sc.Scheduler(rs)
