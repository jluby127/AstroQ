import os
import json
import pandas as pd
import numpy as np

from configparser import ConfigParser


import astroq.splan as splan
import astroq.management as mn
import astroq.benchmarking as bn
import astroq.blocks as ob
import astroq.plot as pl
import astroq.io as io
import astroq.nplan as nplan
import astroq.history as hs
import astroq.webapp as app
import astroq.weather as wh
import astroq.dynamic as dn
import astroq.access as ac
from datetime import datetime
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.ERROR)

def bench(args):
    cf = args.config_file
    number_slots = args.number_slots
    thin = args.thin
    print(f'bench function: config_file is {cf}')
    print(f'bench function: number_slots is {number_slots}')
    print(f'bench function: thin is {thin}')

    # Load the requests frame and thin it
    config = ConfigParser()
    config.read(cf)
    upstream_path = eval(config.get('required', 'folder'), {"os": os})
    semester_directory = upstream_path
    requests_frame = pd.read_csv(os.path.join(semester_directory, "inputs/requests.csv"))
    original_size = len(requests_frame)
    requests_frame = requests_frame.iloc[::thin]
    new_size = len(requests_frame)
    print(f'Request frame thinned: {original_size} -> {new_size} rows (factor of {thin})')
    
    requests_frame.to_csv(os.path.join(semester_directory, "inputs/Requests.csv"))
    
    # Run the schedule directly from config file
    schedule = splan.SemesterPlanner(cf)
    schedule.run_model()
    return

def kpfcc(args):
    """
    Main KPFCC command line interface.
    """
    if args.kpfcc_subcommand is None:
        print("run astroq kpfcc --help for helpfile")
        return

    return

def kpfcc_prep(args):
    cf = args.config_file
    print(f'kpfcc_prep function: config_file is {cf}')
    config = ConfigParser()
    config.read(cf)

    # Get workdir from global section
    workdir = str(config.get('global', 'workdir'))
    savepath = workdir
    
    semester = str(config.get('global', 'semester'))
    start_date = str(config.get('global', 'semester_start_day'))
    end_date = str(config.get('global', 'semester_end_day'))
    start = datetime.strptime(start_date, "%Y-%m-%d")
    end = datetime.strptime(end_date, "%Y-%m-%d")
    n_days = (end - start).days

    # First capture the allocation info
    allo_source = args.allo_source
    allocation_file = str(config.get('data', 'allocation_file'))
    if allo_source == 'db':
        print(f'Pulling allocation information from database')
        awarded_programs = ob.pull_allocation_info(start_date, n_days, 'KPF-CC', savepath+allocation_file)
        awarded_programs = [semester + "_" + val for val in awarded_programs if val != 'U268'] #temporary because PI made a mistake. 
    else:
        print(f'Using allocation information from file: {allo_source}')
        awarded_programs = ac.format_keck_allocation_info(allo_source, savepath+allocation_file)
        awarded_programs = [semester + "_" + val for val in awarded_programs if val != 'U268'] #temporary because PI made a mistake. 

    # Next get the request sheet
    request_file = str(config.get('data', 'request_file'))
    OBs = ob.pull_OBs(semester)
    good_obs, bad_obs_values, bad_obs_hasFields, bad_obs_count_by_semid, bad_field_histogram = ob.get_request_sheet(OBs, awarded_programs, savepath + request_file)

    # Next get the past history 
    past_source = args.past_source
    past_file = str(config.get('data', 'past_file'))
    if past_source == 'db':
        print(f'Pulling past history information from database')
        raw_history = hs.pull_OB_histories(semester)
        obhist = hs.write_OB_histories_to_csv(raw_history, savepath + past_file)
    else:
        print(f'Using past history information from file: {past_source}')
        obhist = hs.write_OB_histories_to_csv_JUMP(good_obs, past_source, savepath + past_file)

    # This is where the custom times info pull will go
    custom_file = str(config.get('data', 'custom_file'))

    return

def kpfcc_data(args):
    pull_file = args.pull_file
    savepath = args.database_file
    print(f'kpfcc_data function: pull_file is {pull_file}')
    print(f'kpfcc_data function: saving to {savepath}')

    with open(pull_file, "r") as f:
        data = json.load(f)
    semester = data["semester"]
    awarded_programs = data["awarded_programs"]
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    OBs = ob.pull_OBs(semester)
    good_obs, bad_obs_values, bad_obs_hasFields, bad_obs_count_by_semid, bad_field_histogram = ob.get_request_sheet(OBs, awarded_programs, savepath + "/Requests.csv")

    send_emails_with = []
    for i in range(len(bad_obs_values)):
        if bad_obs_values['metadata.semid'][i] in awarded_programs:
            send_emails_with.append(ob.inspect_row(bad_obs_hasFields, bad_obs_values, i))
    '''
    this is where code to automatically send emails will go. Not implemented yet.
    '''
    return

def kpfcc_webapp(args):
    """
    Launch web app to view interactive
    plots.
    """
    config_file = args.config_file
    app.launch_app(config_file)
    return

def kpfcc_plan_semester(args):
    """
    Plan a semester's worth of observations using optimization.
    Builds the request set on-the-fly from input data, then runs optimization.
    """
    cf = args.config_file
    print(f'kpfcc_plan_semester function: config_file is {cf}')

    # Run the semester planner directly from config file
    semester_planner = splan.SemesterPlanner(cf)
    semester_planner.run_model()
    return


def plot_pkl(args):
    cf = args.config_file
    print(f'plot function: config_file is {cf}')
    pl.build_semester_webapp_pkl(cf)
    return

def plot_static(args):
    cf = args.config_file
    print(f'plot function: config_file is {cf}')
    pl.build_static_plots(cf)
    return

def ttp(args):
    cf = args.config_file
    print(f'ttp function: config_file is {cf}')
    
    # Use the new NightPlanner class for object-oriented night planning
    night_planner = nplan.NightPlanner(cf)
    night_planner.run_ttp()
    return


def get_dynamics(args):
    cf = args.config_file
    print(f'get_dynamics function: config_file is {cf}')

    
    # Use NightPlanner for reports_directory
    night_planner = nplan.NightPlanner(cf)
    data_astroq = pl.read_star_objects(night_planner.reports_directory + "star_objects.pkl")
    # all_stars_list = [star_obj for star_obj_list in data_astroq[0].values() for star_obj in star_obj_list]
    
    # Create SemesterPlanner for get_cof, get_birdseye, and get_requests_frame
    semester_planner = splan.SemesterPlanner(cf)
    table_reqframe_html = dn.get_requests_frame(semester_planner, filter_condition=None)
    fig_cof_html = dn.get_cof(semester_planner, list(data_astroq[1].values()))
    fig_birdseye_html = dn.get_birdseye(semester_planner, data_astroq[2], list(data_astroq[1].values()))

    ttp_path = os.path.join(night_planner.reports_directory , "ttp_data.pkl")
    if os.path.exists(ttp_path):
        data_ttp = pl.read_star_objects(ttp_path)
        script_table_html = dn.get_script_plan(cf, data_ttp)
        ladder_html = dn.get_ladder(data_ttp)
        slew_animation_html = dn.get_slew_animation(data_ttp, animationStep=120)
        slew_path_html = dn.plot_path_2D_interactive(data_ttp)

    # TODO: need to fix this.
    #dn.get_tau_inter_line(list(data_astroq[0].values())[0])

def requests_vs_schedule(args):
    cf = args.config_file
    sf = args.schedule_file
    print(f'requests_vs_schedule function: config_file is {cf}')
    print(f'requests_vs_schedule function: schedule_file is {sf}')
    # Create semester planner to get strategy data
    semester_planner = splan.SemesterPlanner(cf)
    semester_planner.run_model()
    req = semester_planner.strategy
    sch = pd.read_csv(sf)
    sch = sch.sort_values(by=['d', 's']).reset_index(drop=True) # Re-order into the real schedule
    # First, ensure no repeated day/slot pairs (does allow missing pairs)
    no_duplicate_slot_err = ("'No duplicate slot' condition violated: "
                             "At least one pair of rows corresponds to "
                             "the same day and slot.")
    assert sch.groupby(['d','s']).size().max()<=1, no_duplicate_slot_err
    for star in req.id:
        star_request = req.query(f"id=='{star}'")
        star_schedule = sch.query(f"r=='{star}'") # Only the slots with the star listed

        # A star might not be scheduled at all. This does not violate constraints, but should be noted.
        if len(star_schedule)==0:
            print(f"{star} not scheduled" )
            continue
    # 1) t_visit: No stars scheduled during another star's slot
        t_visit = star_request.t_visit.values[0] # Number of slots needed to complete observation
        star_inds = star_schedule.index
        day_slot = sch[['d', 's']]
        # Check the number of slots between consecutive obs. If they're on the same day, demand a minimum separation
        if star_inds.max() == day_slot.index.max(): # Special case: if this star includes the last obs in the whole schedule
            star_inds = star_inds[:-1] # Exclude the very last observation to avoid index err. That obs can't be overlapped by a later target anyway
        day_slot_diffs = day_slot.iloc[star_inds+1].reset_index() - day_slot.iloc[star_inds].reset_index()
        if len(day_slot_diffs.query('d==0'))==0: # d==0 when next obs is on the same day. If the target is always the last obs of the night, pass
            pass
        else:
            closest_slot_separation = day_slot_diffs.query('d==0').s.min()
            t_visit_err = ("t_visit violated: "
                           "Two stars are scheduled too close together: "
                          f"{star} requires {t_visit} slots but another "
                          f"star is scheduled after only {closest_slot_separation}.")
            assert closest_slot_separation >= t_visit, t_visit_err

    # 2) n_inter_max: Total number of nights a target is scheduled in the semester is less than n_inter_max
        n_inter_max = star_request['n_inter_max'].values[0]
        n_inter_sch = len(set(star_schedule.d)) # All unique nights with scheduled obs
        # Now make sure the number of visits is less than the limit
        n_inter_max_err = ("n_inter_max violated: "
                          f"{star} is scheduled too many times in the semester "
                          f"(scheduled: {n_inter_sch} obs; required: {n_inter_max} obs)")
        assert n_inter_sch <= n_inter_max, n_inter_max_err

    # 3) n_intra_min, n_intra_max: N obs per day is between n_intra_min and n_intra_max
        # t_visit, the number of slots required to complete a single observation (aka visit)
        t_visit = req[req.id==star].t_visit.values
        # Upper/lower limits on N obs per day
        n_intra_min, n_intra_max = star_request[['n_intra_min', 'n_intra_max']].values[0]
        # Scheduled min/max number of obs per day
        n_intra_groupby = star_schedule.groupby(['d']).size() # The numerator gives the sum of all starting slots in which the target is observed in a day.
        n_intra_min_sch, n_intra_max_sch = n_intra_groupby.min(), n_intra_groupby.max()
        # Ensure the target is never scheduled too few/many times in one night
        n_intra_min_err = ("n_intra_min violated: "
                          f"{star} is scheduled too few times in one night "
                          f"(scheduled: {n_intra_min_sch} obs; required: {n_intra_min} obs)")
        assert n_intra_min <= n_intra_min_sch, n_intra_min_err

        n_intra_max_err = ("n_intra_max violated: "
                          f"{star} is scheduled too many times in one night "
                          f"(scheduled: {n_intra_max_sch} obs; required: {n_intra_max} obs)")
        assert n_intra_max_sch <= n_intra_max, n_intra_max_err

    # 4) tau_inter: There must be at least tau_inter nights between successive nights during which a target is observed
        tau_inter = star_request['tau_inter'].values[0] # min num of nights before another obs
        if tau_inter > 0: # only run this test if the intention is to schedule more than once
            unique_days = np.sort(np.array(list(set(star_schedule.d))))
            min_day_gaps = np.min(unique_days[1:] - unique_days[:-1])
            if n_inter_max <= 1: # If only 1 obs per semester, no risk of spacing obs too closely
                pass
            else:
                # Require that all gaps are greater than the min gap
                tau_inter_err = ("tau_inter violated: "
                                f"two obs of {star} are not spaced by enough days "
                                f"(scheduled: {min_day_gaps} days; required: {tau_inter} days)")
                assert min_day_gaps >= tau_inter, tau_inter_err

    # 5) tau_intra: There must be at least tau_intra slots between successive observations of a target in a single night
        slot_duration = semester_planner.slot_size # Slot duration in minutes
        slots_per_hour = 60/slot_duration
        tau_intra_hrs = star_request['tau_intra'].values[0] # min num of hours before another obs
        tau_intra_slots = tau_intra_hrs * slots_per_hour
        min_slot_diffs = star_schedule.groupby('d').s.diff().min() # Group by day, then find successive differences between slot numbers in the same day. Differences are not computed between the last slot of one night and the first slot of the next night (those values are NaN). The differences must all be AT LEAST tau_intra.
        if n_intra_max <= 1: # If only 1 obs per night, no risk of spacing obs too closely
            pass
        else:
            tau_intra_err = ("tau_intra_violated: "
                            f"two obs of {star} are not spaced by enough slots "
                            f"(scheduled: {min_slot_diffs} slots; required: {tau_intra_slots} slots)")
            assert min_slot_diffs >= tau_intra_slots, tau_intra_err

