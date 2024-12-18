import pandas as pd

def get_Dr(group_frame, r):
    """
    Quick return the set of days d where request r has at least one available slot

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
        - r (str) = the request ID
    """
    return list(group_frame[group_frame.r==r].groupby('d').first().index)

def get_Rds(group_frame, d, s):
    """
    Quick return the requests r that are observable in slot s of day d

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
        - d (int) = the day of the semester
        - s (int) = the slot in the night
    """
    return list(group_frame[(group_frame.d==d) & (group_frame.s==s)].r)

def get_Sr(group_frame, r):
    """
    Quick return the set of d,s tuples where requests r is observable

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
        - r (str) = the request ID
    """
    return list(zip(list(group_frame[(group_frame.r==r)].d), list(group_frame[(group_frame.r==r)].s)))

def get_Srd(group_frame, r, d):
    """
    Quick return the set of slots in day d where request r is observable

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
        - r (str) = the request ID
        - d (int) = the day of the semester
    """
    return list(group_frame[(group_frame.r==r) & (group_frame.d==d)].s)

def get_Tds(group_frame, n_nights_in_semester, n_slots_in_night):
    """
    Quick return the set of slots/day pairs where at least one request is observable

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
    """
    yes = []
    no = []
    for d in range(n_nights_in_semester):
        for s in range(n_slots_in_night):
            val = get_Rds(group_frame, d, s)
            if len(val) > 0:
                yes.append([d, s])
            else:
                no.append([d, s])
    return yes, no

def get_Ir(group_frame):
    """
    Quick return the set of requests that want multiple visits per night

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
    """
    return list(group_frame[group_frame.i>1].groupby('r').first().index)

def get_Jr(group_frame):
    """
    Quick return the set of requests that want only one visit per night

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
    """
    return list(group_frame[group_frame.i==1].groupby('r').first().index)

def get_multi_slot_requests(group_frame):
    """
    Quick return the set of requests that want only one visit per night

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
    """
    return list(group_frame[group_frame.e>1].groupby('r').first().index)

def get_single_slot_requests(group_frame):
    """
    Quick return the set of requests that want only one visit per night

    Args:
        - group_frame (obj) = a dataframe groupped by 'rr','dd','ss'
    """
    return list(group_frame[group_frame.e==1].groupby('r').first().index)

def get_Rds_for_multislot(group_frame, d, s, explen):
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
        one_slot = list(group_frame[(group_frame.d==d) & (group_frame.s==s + e)].index)
        for o in one_slot:
            all_indices.append(o)
    return all_indices
