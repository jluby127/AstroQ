import argparse
import os
import sys
import json
import pandas as pd
import numpy as np
import math
from configparser import ConfigParser
from argparse import Namespace

import astroq.splan as splan
import astroq.request as rq
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

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.ERROR)
# log.debug("Debug message")
# log.info("Info message")
# log.warning("Warning message here")
# log.error("Error message")
# log.critical("Critical message")

def bench(args):
    nS = args.number_slots
    cf = args.config_file
    thin = getattr(args, 'thin', 1)  # Default to 1 if thin not provided
    print(f'bench function: config_file is {cf}')
    print(f'bench function: n_Slots for single visits is {nS}')
    print(f'bench function: thinning factor is {thin}')

    # Initialize manager and compute request set on the fly
    # This is a hacky workaround. run_admin needs this file to exist. This can
    # lead to race conditions if benchmarking is run in parallel.
    config = ConfigParser()
    config.read(cf)
    upstream_path = eval(config.get('required', 'folder'), {"os": os})
    semester_directory = upstream_path
    requests_frame = bn.build_toy_model_from_paper(nS)
    
    # Apply thinning if specified
    if thin > 1:
        original_size = len(requests_frame)
        requests_frame = requests_frame.iloc[::thin].reset_index(drop=True)
        new_size = len(requests_frame)
        print(f'Request frame thinned: {original_size} -> {new_size} rows (factor of {thin})')
    
    requests_frame.to_csv(os.path.join(semester_directory, "inputs/Requests.csv"))
    manager = mn.data_admin(cf)
    manager.run_admin()

    # Build observability maps and request set
    strategy, observable = rq.define_indices_for_requests(manager)
    meta = rq.build_meta(cf)
    request_set = rq.RequestSet(meta, strategy, observable)
    # print out the last rows of strategy to ensure the size of the model looks right
    print("Last 5 rows of Request Set here: ")
    print(request_set.strategy[-5:])
    request_set.to_json(os.path.join(manager.output_directory, "request_set.json"))

    # Run the schedule
    schedule = splan.SemesterPlanner(request_set, cf)
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

def kpfcc_build(args):
    cf = args.config_file
    print(f'kpfcc_build function: config_file is {cf}')

    manager = mn.data_admin(cf)
    manager.run_admin()
    strategy, observable = rq.define_indices_for_requests(manager)
    meta = rq.build_meta(cf)
    request_set = rq.RequestSet(meta, strategy, observable)
    request_set.to_json(manager.output_directory + "request_set.json")
    return

def kpfcc_prep(args):
    cf = args.config_file
    print(f'kpfcc_prep function: config_file is {cf}')
    mn.prepare_new_semester(cf)
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

    # Build request set on the fly (no longer need separate build step)
    manager = mn.data_admin(cf)
    manager.run_admin()
    
    # Build observability maps and request set
    strategy, observable = rq.define_indices_for_requests(manager)
    meta = rq.build_meta(cf)
    request_set = rq.RequestSet(meta, strategy, observable)
    
    # Run the semester planner
    semester_planner = splan.SemesterPlanner(request_set, cf)
    semester_planner.run_model()
    return

def schedule(args):
    rf = args.request_file
    cf = args.config_file
    print(f'schedule function: request_file is {rf}')
    print(f'schedule function: config_file is {cf}')
    request_set = rq.read_json(rf)
    schedule = splan.SemesterPlanner(request_set, cf)
    schedule.run_model()
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
    manager = mn.data_admin(cf)
    manager.run_admin()
    
    # Use the new NightPlanner class for object-oriented night planning
    night_planner = nplan.NightPlanner(manager)
    night_planner.run_ttp()
    # Removed: night_planner.produce_bright_backups()
    return

def get_history(args):
    cf = args.config_file
    print(f'get_history function: config_file is {cf}')
    manager = mn.data_admin(cf)
    manager.run_admin()
    database_info_dict = hs.build_past_history(manager.past_database_file, manager.requests_frame, manager.twilight_frame)
    return

def get_dynamics(args):
    cf = args.config_file
    print(f'get_dynamics function: config_file is {cf}')
    manager = mn.data_admin(cf)
    manager.run_admin()
    data_astroq = pl.read_star_objects(manager.reports_directory + "admin/" + manager.current_day + "/star_objects.pkl")
    # all_stars_list = [star_obj for star_obj_list in data_astroq[0].values() for star_obj in star_obj_list]
    table_reqframe_html = dn.get_requests_frame(manager, filter_condition=None)
    fig_cof_html = dn.get_cof(manager, list(data_astroq[1].values()))
    fig_birdseye_html = dn.get_birdseye(manager, data_astroq[2], list(data_astroq[1].values()))
    # dn.interactive_sky_with_static_heatmap(manager, 'U001') # This is broken, need to rethink what it means to have a twilight map.

    ttp_path = os.path.join(manager.reports_directory,"observer",manager.current_day,"ttp_data.pkl")
    if os.path.exists(ttp_path):
        data_ttp = pl.read_star_objects(ttp_path)
        script_table_html = dn.get_script_plan(manager, data_ttp)
        ladder_html = dn.get_ladder(manager, data_ttp)
        slew_animation_html = dn.get_slew_animation(manager, data_ttp, animationStep=120)
        slew_path_html = dn.plot_path_2D_interactive(manager, data_ttp)
    dn.get_tau_inter_line(manager, list(data_astroq[0].values())[0])

def requests_vs_schedule(args):
    rf = args.request_file
    sf = args.schedule_file
    print(f'kpfcc_data function: request_set_file is {rf}')
    print(f'kpfcc_data function: serialized_output_file to {sf}')

    req_object = rq.read_json(rf)
    req = req_object.strategy
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
        slot_duration = req_object.meta.slot_duration # Slot duration in minutes
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
#
# def make_simulated_history(args):
#
#     cf = args.config_file
#     print(f'make_simulated_history function: config_file is {cf}')
#     config = ConfigParser()
#     config.read(cf)
#     upstream_path = eval(config.get('required', 'folder'), {"os": os})
#     semester_directory = upstream_path
#     requests_frame = bn.build_toy_model_from_paper(12)
#     requests_frame.to_csv(os.path.join(semester_directory, "inputs/Requests.csv"))
#     manager = mn.data_admin(cf)
#     manager.run_admin()
#
#     # Run a weather-loss model that will serve as truth for all runs of the semester
#     loss_stats_remaining = wh.get_loss_stats(manager)
#     loss_stats_remaining = [x + 0.2 for x in loss_stats_remaining]
#     allocation_remaining_post_weather_loss, weather_diff_remaining, weather_diff_remaining_1D, \
#         days_lost = wh.simulate_weather_losses(manager.allocation_remaining, loss_stats_remaining, \
#         covariance=0.14, run_weather_loss=True, plot=False, outputdir=manager.output_directory)
#
#     reshaped = weather_diff_remaining_1D.reshape((184, 4))
#     weathered_days = np.where(np.any(reshaped == 1, axis=1))[0].tolist()
#     print("# of on-sky days lost to weather: ", len(weathered_days))
#
#     # Run from start of semester with no weather loss to get a first forecast
#     # this is equivalent to running the bench, so let's use the function
#     myargs = Namespace(config_file='bench/config_benchmark.ini', number_slots=12, number_requests=1)
#     bench(myargs)
#     pl.run_plot_suite(cf)
#     # Save the original forecast
#     forecast = pd.read_csv(manager.output_directory + 'raw_combined_semester_schedule_Round2.txt')
#     forecast = np.array(forecast.values)
#
#     test_dates=["2018-08-15", "2018-09-01", "2018-09-15", "2018-10-01", "2018-10-15", "2018-11-01", "2018-11-15", "2018-12-01", "2018-12-15","2019-01-01","2019-01-15"]
#     for i in range(len(test_dates)):
#         # Construct the simulated past history, including the lost observations due to weather
#         all_dfs = []
#         day_index = manager.all_dates_dict[test_dates[i]]
#         for j in range(day_index):
#             if j not in weathered_days:
#                 stars_for_night = [x for x in forecast[j] if x not in ("X", "*X", "", "Past") and not (isinstance(x, float) and math.isnan(x))]
#                 stars_for_night = list(set(stars_for_night))
#                 if stars_for_night != []:
#                     today = manager.all_dates_array[j]
#                     stamps = [today + ' 12:00:00.000000+00:00']*len(stars_for_night) #time doesn't matter, date does.
#                     today_frame = pd.DataFrame({'star_id':stars_for_night, 'utctime':stamps})
#                     all_dfs.append(today_frame)
#         #update the day, then run admin to create new folders
#         line_number = 7
#         new_date = test_dates[i]
#         new_content = f"current_day = {new_date}"
#         with open(cf, "r") as f:
#             lines = f.readlines()
#         lines[line_number - 1] = new_content + "\n"
#         with open(cf, "w") as f:
#             f.writelines(lines)
#         manager = mn.data_admin(cf)
#         manager.run_admin()
#
#         # write out the simulated past history twice, the inputs version will be over-written, the outputs version will not
#         final_df = pd.concat(all_dfs, ignore_index=True)
#         final_df.to_csv(manager.upstream_path + "inputs/queryJumpDatabase.csv", index=False)
#         final_df.to_csv(manager.output_directory + "queryJumpDatabase.csv", index=False)
#
#         # rerun the build and scheduler from the new date forward
#         strategy, observable = rq.define_indices_for_requests(manager)
#         meta = rq.build_meta(cf)
#         request_set = rq.RequestSet(meta, strategy, observable)
#         # print out the last rows of strategy to ensure the size of the model looks right
#         request_set.to_json(os.path.join(manager.output_directory, "request_set.json"))
#         # Run the schedule
#         schedule = splan.SemesterPlanner(request_set, cf)
#         schedule.run_model()
#         # Run the plot suite
#         pl.run_plot_suite(cf)
#     return
