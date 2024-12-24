import sys
import argparse
import os

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.time import TimeDelta

# the -3 cuts of the "bin/" of the path to this current file
path2modules = os.path.dirname(os.path.abspath(__file__))[:-3]
sys.path.append(path2modules  + "/kpfcc/")
import solve_semester as ss
import helper_functions as hf
import plotting_functions as ptf
import processing_functions as pf

sys.path.append(os.environ["TTP_PATH"] + "ttp/")
import formatting
import telescope
import plotting
import model

named_colors = ['blue', 'red', 'green', 'gold', 'maroon', 'gray', 'orange', 'magenta', 'purple']

parser = argparse.ArgumentParser(description='Generate schedules with KPF-CC v2')
parser.add_argument('-d','--schedule_dates',action='append',help='Date(s) to be scheduled; \
                            string in format YYYY-MM-DD. Use a -d flag between each date.')
parser.add_argument('-f','--folder', help='Folder to save all outputs',
                            default=os.environ["KPFCC_SAVE_PATH"])
parser.add_argument('-r','--run_extra_rounds', help='Run the bonus round', action='store_true')
parser.add_argument('-t','--timeout', help='Max time spent optimizing (sec)',type=int, default=300)
parser.add_argument('-s','--slot_size', help='The slot size (minutes)', type=int, default=10)
parser.add_argument('-w','--dont_run_weather_loss', help='If True, do not simulate weather losses', default=True)
parser.add_argument('-ttp','--run_ttp', help='If True, run the TTP.', action='store_false')
parser.add_argument('-g','--show_gurobi', help='Turn on Gurobi console print', action='store_false')
parser.add_argument('-p','--show_plots', help='Turn on plot outputs', action='store_false')
parser.add_argument('-b','--run_backups', help='Turn on plot outputs', action='store_true')

args = parser.parse_args()

request_sheet = args.folder + "inputs/Requests.csv"
allocated_nights = args.folder + "inputs/2024B_Binary_Schedule.txt"
past_database = args.folder + "nopasthistoryplease.txt" #"inputs/queryJumpDatabase.csv"
twilight_times = args.folder + "inputs/2024B_twilight_times.csv"
access_map = args.folder + "inputs/2024B_AccessMaps_" + str(args.slot_size) + "minSlots.txt"
turn_on_off_file = args.folder + "inputs/2024B_turnOnOffDates.csv"
starmap_template_filename = args.folder + "inputs/2024B_cadenceTemplateFile.csv"
nonqueue_map =  args.folder + 'inputs/2024B_NonQueueMap'  + str(args.slot_size) + '.txt'
special_map = args.folder + 'inputs/2024B_specialMaps_' + str(args.slot_size) + 'minSlots.txt'
nightly_start_stop_times = pd.read_csv(args.folder + "inputs/2024B_NightlyStartStopTimes.csv")
zero_out_file = 'nofilename.txt'
folder_forecasts = args.folder + '/data/first_forecasts/'
folder_cadences = args.folder + '/outputs/' + str(args.schedule_dates[0]) + '/cadences/'
backup_file = path2modules + "data/bright_backups_frame.csv"
backup_observability_file = path2modules + "data/bright_backup_observability.csv"

ss.run_kpfcc(args.schedule_dates,
                          request_sheet,
                          allocated_nights,
                          access_map,
                          twilight_times,
                          args.folder + 'outputs/' + str(args.schedule_dates[0]) + '/',
                          args.slot_size,
                          args.run_extra_rounds,
                          past_database,
                          starmap_template_filename,
                          turn_on_off_file,
                          nonqueue_map,
                          special_map,
                          zero_out_file,
                          args.dont_run_weather_loss,
                          args.show_gurobi,
                          args.show_plots,
                          args.timeout)

print("Semester is scheduled. Complete.")
if args.show_plots:
    future_forecast = args.folder + "outputs/" + str(args.schedule_dates[0]) + \
                      '/raw_combined_semester_schedule_Round2.txt'
    nonqueue_list = args.folder + "inputs/NonQueue.csv"

    data = ptf.DataHandler(args.schedule_dates[0], request_sheet, past_database, future_forecast,
                           twilight_times, nonqueue_list, folder_forecasts, folder_cadences,
                           args.folder + 'reports/', args.slot_size)

    all_program_cofs = []
    all_program_info = {}
    all_program_first_forecast = {}

    all_star_maps_all_programs = []
    all_star_colors_all_programs = []
    all_star_nObs_all_programs = []
    all_stars_in_program = []

    all_stars = []
    all_completions = []
    for p, item in enumerate(data.programs):
        program_COF, program_total_obs_request, program_first_forecast, stars_in_program, stars_completion = \
            ptf.generate_single_program_plot_suite(data, data.programs[p])
        all_program_cofs.append(program_COF)
        all_stars.append(stars_in_program)
        all_completions.append(stars_completion)

        all_program_info[data.programs[p]] = program_total_obs_request
        all_program_first_forecast[data.programs[p]] = program_first_forecast

        all_star_maps, all_star_cols, all_star_nObs, stars_in_program = \
            ptf.create_single_program_birds_eye(data, data.programs[p])
        all_star_maps_all_programs.append(all_star_maps)

        program_color = named_colors[p%len(named_colors)]
        all_star_colors_all_programs += [program_color]*len(all_star_maps)
        all_star_nObs_all_programs.append(all_star_nObs)

        all_stars_in_program.append(stars_in_program)

        print("Building Birds Eye View")
        ptf.generate_birds_eye(data, all_star_maps, all_star_cols, all_star_nObs,
                             data.programs[p], stars_in_program)

    flat_all_star_maps_all_programs = np.concatenate(all_star_maps_all_programs)
    flat_all_star_nObs_all_programs = np.concatenate(all_star_nObs_all_programs)
    flat_all_stars_in_program = np.concatenate(all_stars_in_program)

    ptf.generate_admin_view_plot_suite(data, all_program_cofs, all_program_info,
                                       all_program_first_forecast, flat_all_star_maps_all_programs,
                                       all_star_colors_all_programs,
                                       flat_all_star_nObs_all_programs, flat_all_stars_in_program)
    all_stars = np.array(all_stars).flatten()
    all_completions = np.array(all_completions).flatten()
    complete_frame = pd.DataFrame({"Starname":all_stars, 'CompletetionPercent':all_completions})
    complete_frame.to_csv(args.folder + 'reports/admin/' + str(args.schedule_dates[0]) + '/cofs/final_completion_rates.csv', index=False)

if args.run_ttp:
    savepath = args.folder + 'reports/observer/' + str(args.schedule_dates[0]) + '/'
    check1 = os.path.isdir(savepath)
    if not check1:
        os.makedirs(savepath)
    output_path = args.folder + 'outputs/' + str(args.schedule_dates[0]) + "/"
    check2 = os.path.isdir(output_path)
    if not check2:
        os.makedirs(output_path)

    print("Prepare schedule for the TTP.")
    observatory = telescope.Keck1()
    semester_start_date, semester_end_date, semester_length, semester_year, semester_letter = \
            hf.get_semester_info(args.schedule_dates[0])
    all_dates_dict = hf.build_date_dictionary(semester_start_date, semester_length)
    the_schedule = np.loadtxt(output_path + 'raw_combined_semester_schedule_Round2.txt',
                                delimiter=',', dtype=str)

    for n in range(len(args.schedule_dates)):
        day_in_semester = all_dates_dict[args.schedule_dates[n]]
        idx = nightly_start_stop_times[nightly_start_stop_times['Date'] == str(args.schedule_dates[n])].index[0]
        night_start_time = nightly_start_stop_times['Start'][idx]
        night_stop_time = nightly_start_stop_times['Stop'][idx]
        observation_start_time = Time(str(args.schedule_dates[n]) + "T" + str(night_start_time),
            format='isot')
        observation_stop_time = Time(str(args.schedule_dates[n]) + "T" + str(night_stop_time),
            format='isot')
        total_time = np.round((observation_stop_time.jd-observation_start_time.jd)*24,3)
        print("Time in Night for Observations: " + str(total_time) + " hours.")

        round_two_requests = np.loadtxt(output_path + 'Round2_Requests.txt', dtype=str)
        request_frame = pd.read_csv(request_sheet)
        send_to_ttp = pf.prepare_for_ttp(request_frame, the_schedule[day_in_semester],
                                            round_two_requests)

        filename = output_path + 'Selected_' + str(args.schedule_dates[n]) + ".txt"
        if args.run_backups and os.path.exists(backup_file) and os.path.exists(backup_observability_file):
            print("Generating filler bright stars and backup weather script")
            check = os.path.isdir(savepath + 'Backups/')
            if not check:
                os.makedirs(savepath + 'Backups/')
            import backup_star_functions as bsf
            backup_starlist = pd.read_csv(backup_file)
            backup_observability_frame = pd.read_csv(backup_observability_file)
            backups = bsf.get_stars_for_tonight(backup_starlist, backup_observability_frame,
                    args.schedule_dates[0][5:], minimum_up_time=4.0)
            backups.drop(columns='index', inplace=True)
            fill_time_hours = 0.
            if fill_time_hours > 0:
                subsetBackups = bsf.get_times_worth(backups, fill_time_hours)
                send_to_ttp_all = pd.concat([send_to_ttp, subsetBackups], ignore_index=True)
            else:
                send_to_ttp_all = send_to_ttp
            send_to_ttp_all.to_csv(filename, index=False)
            filename_backups = savepath + '/Backups_' + str(args.schedule_dates[n]) + ".txt"
            backups_full = bsf.get_times_worth(backups, total_time+1)
            backups_full.to_csv(filename_backups, index=False)
        else:
            send_to_ttp.to_csv(filename, index=False)

        target_list = formatting.theTTP(filename)

        solution = model.TTPModel(observation_start_time, observation_stop_time, target_list,
                                    observatory, savepath)

        plotting.writeStarList(solution.plotly, observation_start_time, args.schedule_dates[n],
                            outputdir=savepath)
        plotting.plot_path_2D(solution,outputdir=savepath)
        plotting.nightPlan(solution.plotly, args.schedule_dates[n], outputdir=savepath)

        obs_and_times = pd.read_csv(savepath + 'ObserveOrder_' + str(args.schedule_dates[0]) + ".txt")

        pf.write_starlist(request_frame, solution.plotly, observation_start_time, solution.extras,
                            round_two_requests, str(args.schedule_dates[0]), savepath)

        if args.run_backups and os.path.exists(backup_file) and os.path.exists(backup_observability_file):
            print('Generating independent backup weather script.')
            targlist_backups = formatting.theTTP(filename_backups)

            solution_b = model.TTPModel(observation_start_time, observation_stop_time,
                                        targlist_backups, observatory, savepath + "Backups/",
                                        runtime=600, optgap=0.05)

            plotting.writeStarList(solution_b.plotly, observation_start_time, args.schedule_dates[n],
                                        outputdir=savepath + "Backups/")
            plotting.plot_path_2D(solution_b,outputdir=savepath+ "Backups/")
            plotting.nightPlan(solution_b.plotly, args.schedule_dates[n],
                                        outputdir=savepath + "Backups/")

            obs_and_times_b = pd.read_csv(savepath + 'Backups/ObserveOrder_' + \
                                str(args.schedule_dates[0]) + ".txt")
            # pf.write_starlist(backup_starlist,
            #                     obs_and_times_b,
            #                     solution_b.extras,
            #                     [],
            #                     'Backups',
            #                     str(args.schedule_dates[0]),
            #                     savepath + "Backups/")

            pf.write_starlist(backup_starlist, solution_b.plotly, observation_start_time,
                                solution_b.extras, [], str(args.schedule_dates[0]),
                                savepath + "Backups/")

    print("Semester is scheduled and script is generated. Complete.")
