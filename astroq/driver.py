"""
Module for executing AstroQ functions based on command line interface inputs.
"""

# Standard library imports
import logging
import os
from datetime import datetime
from configparser import ConfigParser

# Third-party imports
import numpy as np
import pandas as pd
import astroplan as apl
import plotly.io as pio
from astropy.time import Time, TimeDelta
from io import BytesIO
import imageio.v3 as iio
import base64

# Local imports
import astroq.access as ac
import astroq.benchmarking as bn
import astroq.prep.kpfcc as kpfcc
import astroq.history as hs
import astroq.io as io
import astroq.nplan as nplan
import astroq.plot as pl
import astroq.splan as splan
import astroq.webapp as app

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)  # Lower level to capture more messages

def bench(args):
    """
    Benchmark the AstroQ pipeline using a toy model.
    
    Args:
        args (argparse.Namespace): the command line arguments with flags:
            -cf (str): the path to the config file.
            -ns (int): the number of slots needed to complete each of the single shot requests in program 6.
            -thin (int): the factor to thin the request frame by (for very fast testing).
    Returns:
        None
    """
    cf = args.config_file
    number_slots = args.number_slots
    thin = args.thin
    print(f'bench function: config_file is {cf}')
    print(f'bench function: number_slots is {number_slots}')
    print(f'bench function: thin is {thin}')

    # Load the requests frame and thin it
    config = ConfigParser()
    config.read(cf)
    semester_directory = config.get('global', 'workdir')
    requests_frame = bn.build_toy_model_from_paper(number_slots)    
    # Ensure starname column is always interpreted as strings
    if 'starname' in requests_frame.columns:
        requests_frame['starname'] = requests_frame['starname'].astype(str)
    original_size = len(requests_frame)
    requests_frame = requests_frame.iloc[::thin]
    new_size = len(requests_frame)
    print(f'Request frame thinned: {original_size} -> {new_size} rows (factor of {thin})')
    requests_frame.to_csv(os.path.join(semester_directory, "request.csv"))

    # Run the schedule directly from config file
    schedule = splan.SemesterPlanner(cf, run_band3=False)
    schedule.run_model()
    return

def kpfcc_prep(args):
    """
    Prepare the KPF-CC program for a new semester. This function is specific to the KPF-CC program and the observatory's infrastructure.
    If you are adapting AstroQ for a new observatory, you will need to write your own module to connect to a new prep <your observatory> command.
    
    Args:
        args (argparse.Namespace): the command line arguments with flags:
            -cf (str): the path to the config file.
            -band_number (int): the band number to filter the request.csv by.
            -is_full_band (bool): whether this is a full-band that should update allocation.csv.
            -allo_source (str): the source of the allocation information, either 'db' or a file path.
            -past_source (str): the source of the past history information, either 'db' or a file path.
            -request_source (str): the source of the request information, either 'db' or a file path.
            -filler_programs (str): the semester ID for the filler program. Ex. 2025B_E473.
    
    Returns:
        None
    """
    cf = args.config_file
    print(f'kpfcc_prep function: config_file is {cf}')
    config = ConfigParser()
    config.read(cf)
    band_number = args.band_number
    is_full_band = args.is_full_band

    # Get workdir from global section
    workdir = str(config.get('global', 'workdir'))
    savepath = workdir
    semester = str(config.get('global', 'semester'))
    start_date = str(config.get('global', 'semester_start_day'))
    end_date = str(config.get('global', 'semester_end_day'))
    current_date = str(config.get('global', 'current_day'))
    start = datetime.strptime(start_date, "%Y-%m-%d")
    end = datetime.strptime(end_date, "%Y-%m-%d")
    n_days = (end - start).days

    # CAPTURE ALLOCATION INFORMATION AND PROCESS
    # --------------------------------------------
    # --------------------------------------------
    allo_source = args.allo_source
    allocation_file = str(config.get('data', 'allocation_file'))
    # pull the allocation 
    if allo_source == 'db':
        print(f'Pulling allocation information from database')
        allocation_frame, hours_by_program, nights_by_program = kpfcc.pull_allocation_info(start_date, n_days, 'KPF-CC')
        awarded_programs = [semester + "_" + val for val in list(hours_by_program.keys())] 
        programmatics = pd.DataFrame({'program': awarded_programs, 'hours': list(hours_by_program.values()), 'nights': list(nights_by_program.values())})
        programmatics.to_csv(os.path.join(savepath, 'programs.csv'), index=False)
    else:
        print(f'Using allocation information from Keck Observatory Instrument Plan (KOIP) file: {allo_source}')
        allocation_frame, hours_by_program, nights_by_program = kpfcc.format_keck_allocation_info(allo_source)
        awarded_programs = [semester + "_" + val for val in list(hours_by_program.keys())]
        programmatics = pd.DataFrame({'program': awarded_programs, 'hours': list(hours_by_program.values()), 'nights': list(nights_by_program.values())})
        programmatics.to_csv(os.path.join(savepath, 'programs.csv'), index=False)
    # else:
    #     print(f'Using allocation information from file: {allo_source}')
    #     # Validate that the file has the correct columns
    #     expected_columns = ['start', 'stop']
    #     if os.path.exists(allo_source):
    #         df = pd.read_csv(allo_source, nrows=0)
    #         actual_columns = set(df.columns)
    #         expected_set = set(expected_columns)
    #         missing_columns = expected_set - actual_columns
    #         if missing_columns:
    #             logging.warning(f"Allocation file '{allo_source}' is missing required columns: {missing_columns}")
    #             print("Note: if inputing your own formatted allocation file, you must also provide your own programs.csv file.")
    #         else:
    #             print(f"Allocation file columns validated: all required columns present")
    #             awarded_programs = df['program'].unique().tolist()
    #     else:
    #         logging.warning(f"Allocation file '{allo_source}' does not exist")

    allocation_frame['comment'] = [''] * len(allocation_frame)
    # Update allocation times for tonight if this is a full-band
    if is_full_band:
        print("Updating allocation.csv for full-band")
        allocation_frame = kpfcc.update_allocation_file(allocation_frame, current_date)
    allocation_frame.sort_values(by='start', inplace=True)
    allocation_frame.to_csv(os.path.join(savepath, allocation_file), index=False)

    # CAPTURE REQUEST INFORMATION AND PROCESS
    # --------------------------------------------
    # --------------------------------------------
    if args.request_source == 'db':
        # Add filler programs if specified
        fillers = args.filler_programs
        # temporarily comment out this block for 2025B.
        if fillers is not None:
            print(f'Adding filler program to awarded_programs: {fillers}')
            awarded_programs.append(fillers)
        # Pull the request sheet
        request_file = str(config.get('data', 'request_file'))
        OBs = kpfcc.pull_OBs(semester)
        good_obs, bad_obs_values, bad_obs_hasFields, bad_obs_count_by_semid, bad_field_histogram = kpfcc.get_request_sheet(OBs, awarded_programs, os.path.join(savepath, request_file))
        # Filter the request sheet by weather band
        filtered_good_obs = kpfcc.filter_request_csv(good_obs, band_number)
        # If no exposure meter threshold set, then OB can only be part of band 1
        if band_number != 1:
            filtered_good_obs = filtered_good_obs[filtered_good_obs['exp_meter_threshold'] != -1.0]
        filtered_good_obs.reset_index(drop=True, inplace=True)
        # Compute nominal exposure times and increase exposure times for different bands
        slowdown_factors = {1: 1.0, 2: 2.0, 3: 4.0}
        slow = slowdown_factors[band_number]
        new_exptimes = kpfcc.recompute_exposure_times(filtered_good_obs, slow)
        filtered_good_obs['exptime'] = new_exptimes
        filtered_good_obs.to_csv(os.path.join(savepath, request_file), index=False)
    
        # CAPTURE CUSTOM INFORMATION AND PROCESS
        # --------------------------------------------
        # --------------------------------------------
        custom_file = str(config.get('data', 'custom_file'))
        custom_frame = kpfcc.format_custom_csv(OBs)
        custom_frame.to_csv(os.path.join(savepath, custom_file), index=False)

        # CAPTURE FILLER REQUEST INFORMATION AND PROCESS
        # --------------------------------------------
        # --------------------------------------------
        # Now get the bright backup stars information from the filler program
        filler_file = str(config.get('data', 'filler_file'))
        if fillers is not None:
            print(f'Generating filler.csv from program: {fillers}')
            good_obs_backup, bad_obs_values_backup, bad_obs_hasFields_backup, bad_obs_count_by_semid_backup, bad_field_histogram_backup = kpfcc.get_request_sheet(OBs, [fillers], os.path.join(savepath, filler_file))
        else:
            print(f'No fillers specified, creating blank filler.csv file.')
            good_obs_backup = pd.DataFrame(columns=good_obs.columns)
        filtered_good_obs_backup = kpfcc.filter_request_csv(good_obs_backup, band_number)
        filtered_good_obs_backup.to_csv(os.path.join(savepath, filler_file), index=False)

        # CAPTURE EMAIL INFORMATION AND PROCESS
        # --------------------------------------------
        # --------------------------------------------
        send_emails_with = []
        for i in range(len(bad_obs_values)):
            if bad_obs_values['metadata.semid'][i] in awarded_programs:
                send_emails_with.append(kpfcc.inspect_row(bad_obs_hasFields, bad_obs_values, i))
        '''
        this is where code to automatically send emails will go. Not implemented yet.
        '''
    else:
        print(f'User specified request source: {args.request_source}')
        print("No action taken on request.csv")
        print("User must also supply a custom.csv file and a filler.csv file.")


    # CAPTURE PAST HISTORY INFORMATION AND PROCESS
    # --------------------------------------------
    # --------------------------------------------
    past_source = args.past_source
    past_file = str(config.get('data', 'past_file'))
    if past_source == 'db':
        print(f'Pulling past history information from database')
        raw_history = kpfcc.pull_OB_histories(semester)
        obhist = hs.write_OB_histories_to_csv(raw_history)
        obhist.to_csv(os.path.join(savepath, past_file), index=False)
    else:
        print(f'Using past history information from file: {past_source}')
        # Validate that the file has the correct columns
        expected_columns = ['id', 'target', 'semid', 'timestamp', 'exposure_start_time', 'exposure_time', 'observer']
        if os.path.exists(past_source):
            df = pd.read_csv(past_source, nrows=0)
            actual_columns = set(df.columns)
            expected_set = set(expected_columns)
            missing_columns = expected_set - actual_columns
            if missing_columns:
                logging.warning(f"Past history file '{past_source}' is missing required columns: {missing_columns}")
            else:
                print(f"Past history file columns validated: all required columns present")
        else:
            logging.warning(f"Past history file '{past_source}' does not exist")

    return

def kpfcc_webapp(args):
    """
    Launch web app to view interactive plots.
    
    Args:
        args (argparse.Namespace): the command line arguments with flags:
            -uptree_path (str): the path to the uptree directory below which the folder structure is <semester_code>/<date>/<band>/.
    
    Returns:
        None
    """
    uptree_path = args.uptree_path
    app.launch_app(uptree_path)
    return

def plan_semester(args):
    """
    Run the core optimization algorithm for determining what stars to observe on what nights.
    
    Args:
        args (argparse.Namespace): the command line arguments with flags:
            -cf (str): the path to the config file.
            -run_band3 (bool): whether to run the band 3 filler program.
    
    Returns:
        None
    """
    cf = args.config_file
    print(f'plan_semester function: config_file is {cf}')
    b3 = args.run_band3
    print(f'plan_semester function: b3 is {b3}')
    semester_planner = splan.SemesterPlanner(cf, b3)
    semester_planner.run_model()
    return

def plan_night(args):
    """
    Run the slew path optimization using the TTP package for a given night's selected targets.
    
    Args:
        args (argparse.Namespace): the command line arguments with flags:
            -cf (str): the path to the config file.
    
    Returns:
        None
    """
    cf = args.config_file
    print(f'plan_night function: config_file is {cf}')
    night_planner = nplan.NightPlanner(cf)
    night_planner.run_ttp()
    night_planner.to_hdf5()
    return

def plot(args):
    """
    Generate html and png files of the standard AstroQ output plots. Mirrors those produced by the webapp. 
    
    Args:
        args (argparse.Namespace): the command line arguments with flags:
            -cf (str): the path to the config file.
    
    Returns:
        None
    """

    cf = args.config_file
    print(f'plot function: using config file from {cf}')
    config = ConfigParser()
    config.read(cf)
    semester_directory = config.get('global', 'workdir')

    if os.path.exists(os.path.join(semester_directory, 'outputs', 'semester_planner.h5')):
        semester_planner = splan.SemesterPlanner.from_hdf5(os.path.join(semester_directory, 'outputs', 'semester_planner.h5'))
        saveout = os.path.join(semester_planner.output_directory, "saved_plots")
        os.makedirs(saveout, exist_ok = True)
        
        data_astroq = pl.process_stars(semester_planner)
        all_stars_from_all_programs = np.concatenate(list(data_astroq[0].values()))

        # build the plots
        request_df = pl.get_request_frame(semester_planner, all_stars_from_all_programs)
        request_table_html = pl.dataframe_to_html(request_df)
    
        fig_cof = pl.get_cof(semester_planner, list(data_astroq[1].values()))
        fig_birdseye = pl.get_birdseye(semester_planner, data_astroq[2], list(data_astroq[1].values()))
        fig_football = pl.get_football(semester_planner, all_stars_from_all_programs, use_program_colors=True)
        fig_tau_inter_line = pl.get_tau_inter_line(semester_planner, all_stars_from_all_programs, use_program_colors=True)

        # write the html versions 
        fig_cof_html = pio.to_html(fig_cof, full_html=True, include_plotlyjs='cdn')
        fig_birdseye_html = pio.to_html(fig_birdseye, full_html=True, include_plotlyjs='cdn')
        fig_football_html = pio.to_html(fig_football, full_html=True, include_plotlyjs='cdn')
        fig_tau_inter_line_html = pio.to_html(fig_tau_inter_line, full_html=True, include_plotlyjs='cdn')

        # write out the html files 
        with open(os.path.join(saveout, "request_table.html"), "w") as f:
            f.write(request_table_html)
        with open(os.path.join(saveout, "all_programs_COF.html"), "w") as f:
            f.write(fig_cof_html)
        with open(os.path.join(saveout, "all_programs_birdseye.html"), "w") as f:
            f.write(fig_birdseye_html)
        with open(os.path.join(saveout, "all_programs_football.html"), "w") as f:
            f.write(fig_football_html)
        with open(os.path.join(saveout, "all_programs_tau_inter_line.html"), "w") as f:
            f.write(fig_tau_inter_line_html)
    else:
        print(f'No semester_planner.h5 found in {semester_directory}/outputs/. No plots will be generated.')

    night_planner_h5 = os.path.join(semester_directory, 'outputs', 'night_planner.h5')
    if os.path.exists(night_planner_h5):
        night_planner = nplan.NightPlanner.from_hdf5(night_planner_h5)
        data_ttp = night_planner.solution
        
        # Get the night start time from allocation file (this is "Minute 0")
        from astroq.nplan import get_nightly_times_from_allocation
        night_start_time, _ = get_nightly_times_from_allocation(
            night_planner.allocation_file, 
            night_planner.current_day
        )

        # build the plots
        script_table_df = pl.get_script_plan(night_planner)
        timepie_fig = pl.get_timepie(semester_planner, all_stars_from_all_programs, use_program_colors=False)
        ladder_fig = pl.get_ladder(data_ttp, night_start_time)
        slew_animation_fig = pl.get_slew_animation_plotly(data_ttp, os.path.join(semester_directory, "request.csv"), animationStep=120)
        slew_path_fig = pl.plot_path_2D_interactive(data_ttp, night_start_time=night_start_time)

        # write the html versions 
        script_table_html = pl.dataframe_to_html(script_table_df)
        timepie_html = pio.to_html(timepie_fig, full_html=True, include_plotlyjs='cdn')
        ladder_html = pio.to_html(ladder_fig, full_html=True, include_plotlyjs='cdn')
        slew_path_html = pio.to_html(slew_path_fig, full_html=True, include_plotlyjs='cdn')
        slew_animation_html = pio.to_html(slew_animation_fig, full_html=True, include_plotlyjs='cdn')

        # write out the html files 
        with open(os.path.join(saveout, "script_table.html"), "w") as f:
            f.write(script_table_html)
        with open(os.path.join(saveout, "timepie_plot.html"), "w") as f:
            f.write(timepie_html)
        with open(os.path.join(saveout, "ladder_plot.html"), "w") as f:
            f.write(ladder_html)
        with open(os.path.join(saveout, "slew_animation_plot.html"), "w") as f:
            f.write(slew_animation_html)
        with open(os.path.join(saveout, "slew_path_plot.html"), "w") as f:
            f.write(slew_path_html)
    else:
        print(f'No night_planner.pkl found in {semester_directory}/outputs/. No plots will be generated.')
    return

def requests_vs_schedule(args):
    """
    Compare the request.csv file to the schedule.csv file to ensure that the schedule is valid. This is a sanity check to ensure that the schedule is not violating any of the constraints.
    
    Args:
        args (argparse.Namespace): the command line arguments with flags:
            -cf (str): the path to the config file.
            -schedule_file (str): the path to the schedule file.
    
    Returns:
        None
    """
    cf = args.config_file
    sf = args.schedule_file
    print(f'requests_vs_schedule function: config_file is {cf}')
    print(f'requests_vs_schedule function: schedule_file is {sf}')
    # Create semester planner to get strategy data
    semester_planner = splan.SemesterPlanner(cf, run_band3=False)
    semester_planner.run_model()
    req = semester_planner.strategy
    sch = pd.read_csv(sf)
    sch = sch.sort_values(by=['d', 's']).reset_index(drop=True) # Re-order into the real schedule
    # First, ensure no repeated day/slot pairs (does allow missing pairs)
    no_duplicate_slot_err = ("'No duplicate slot' condition violated: "
                             "At least one pair of rows corresponds to "
                             "the same day and slot.")
    assert sch.groupby(['d','s']).size().max()<=1, no_duplicate_slot_err
    for star in req.unique_id:
        star_request = req.query(f"unique_id=='{star}'")
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
        t_visit = req[req.unique_id==star].t_visit.values
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
            if len(unique_days) <= 1: # If only 1 obs or 1 day, no risk of spacing obs too closely
                pass
            else:
                min_day_gaps = np.min(unique_days[1:] - unique_days[:-1])
                # Require that all gaps are greater than the min gap
                tau_inter_err = ("tau_inter violated: "
                                f"two obs of {star} are not spaced by enough days "
                                f"(scheduled: {min_day_gaps} days; required: {tau_inter} days)")
                assert min_day_gaps >= tau_inter, tau_inter_err

    # 5) tau_intra: There must be at least tau_intra slots between successive observations of a target in a single night
        slot_duration = semester_planner.slot_size # Slot duration in minutes
        slots_per_hour = 60/slot_duration
        tau_intra_slots = star_request['tau_intra'].values[0] # recall that the tau_intra is already in units of slots
        min_slot_diffs = star_schedule.groupby('d').s.diff().min() # Group by day, then find successive differences between slot numbers in the same day. Differences are not computed between the last slot of one night and the first slot of the next night (those values are NaN). The differences must all be AT LEAST tau_intra.
        if n_intra_max <= 1: # If only 1 obs per night, no risk of spacing obs too closely
            pass
        else:
            tau_intra_err = ("tau_intra_violated: "
                            f"two obs of {star} are not spaced by enough slots "
                            f"(scheduled: {min_slot_diffs} slots; required: {tau_intra_slots} slots)")
            assert min_slot_diffs >= tau_intra_slots, tau_intra_err

