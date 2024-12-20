import pandas as pd
import numpy as np
# new set_functions.py with pandas commands

def get_Dr(df, r):
    """
    Quick return the set of days d where request r has at least one available slot

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
        - r (str) = the request ID
    """
    return df['d'][df.r==r].unique()

def get_Rds(df, d, s):
    """
    Quick return the requests r that are observable in slot s of day d

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
        - d (int) = the day of the semester
        - s (int) = the slot in the night
    """
    return list(df['r'][(df.d==d) & (df.s==s)])#.sort()

def get_Sr(df, r):
    """
    Quick return the set of d,s tuples where requests r is observable

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
        - r (str) = the request ID
    """
    return df[['d', 's']][df['r'] == r].values

def get_Srd(df, r, d):
    """
    Quick return the set of slots in day d where request r is observable

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
        - r (str) = the request ID
        - d (int) = the day of the semester
    """
    return df['s'][(df.r==r) & (df.d==d)]

def get_Tds(df, n_nights_in_semester, n_slots_in_night):
    """
    Quick return the set of slots/day pairs where at least one request is observable

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
    """
    yes = []
    no = []
    for d in range(n_nights_in_semester):
        for s in range(n_slots_in_night):
            val = get_Rds(df, d, s)
            if len(val) > 0:
                yes.append([d, s])
            else:
                no.append([d, s])
    return yes, no

def get_Ir(df):
    """
    Quick return the set of requests that want multiple visits per night

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
    """
    return df['r'][df.i>1].unique()

def get_Jr(df):
    """
    Quick return the set of requests that want only one visit per night

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
    """
    return list(df['r'][df.i==1].unique())

def get_multi_slot_requests(df):
    """
    Quick return the set of requests that want only one visit per night

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
    """
    return df['r'][df.e>1].unique()#.sort()

def get_single_slot_requests(df):
    """
    Quick return the set of requests that want only one visit per night

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
    """
    return df['r'][df.e==1].unique()#.sort()

def get_Rds_for_multislot(df, d, s, explen):
    """
    Quick return all the tuple indices for range (r, d, s) to (r, d, s + e) for all r.
    Used in Contraint 2: Reserve slots for multi-slot exposure

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
        - d (int) = the day of the semester
        - s (int) = the slot in the night
        - explen (int) = the exposure length in units of slots
    """
    all_indices = []
    for e in range(1, explen):
        one_slot = df[['r', 'd', 's']][(df.d==d) & (df.s==s + e)].values#.index
        for o in one_slot:
            all_indices.append(o)
    all_indices = sorted(all_indices, key=lambda x: x[0])
    return all_indices
