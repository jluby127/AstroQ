"""
Module for plotting.
Designed to be only run as a function call from the generateScript.py script.

Example usage:
    import plotting_functions as ptf
"""
import os
import warnings
warnings.filterwarnings('ignore')
import plotly.graph_objects as go
import plotly.express as px
from astropy.time import Time
import seaborn as sns

import numpy as np
import pandas as pd

import kpfcc.io as io

palette = sns.color_palette("colorblind")
named_colors = ['blue', 'red', 'green', 'gold', 'maroon', 'gray', 'orange', 'magenta', 'purple']
color_scales = {'blue':[[0, 'rgba(0,0,0,0)'],[1, 'rgb(0,0,255)']],
                'red':[[0, 'rgba(0,0,0,0)'],[1, 'rgb(255,0,0)']],
                'green':[[0, 'rgba(0,0,0,0)'],[1, 'rgb(0,255,0)']],
                'gold':[[0, 'rgba(0,0,0,0)'],[1, 'rgb(255, 215, 0)']],
                'maroon':[[0, 'rgba(0,0,0,0)'],[1, 'rgb(128,0,0)']],
                'gray':[[0, 'rgba(0,0,0,0)'],[1, 'rgb(128,128,128)']],
                'orange':[[0, 'rgba(0,0,0,0)'],[1, 'rgb(255,165,0)']],
                'magenta':[[0, 'rgba(0,0,0,0)'],[1, 'rgb(255,0,255)']],
                'purple':[[0, 'rgba(0,0,0,0)'],[1, 'rgb(127,0,255)']],}

def run_plot_suite(star_tracker, manager):
    """
    A one line call to kick off all the plotting/reporting routines

    Args:
        star_tracker (object): a object from the StarTracker class
        outputdir (str): directory to save outputs
        build_starmaps (boolean): if true, spend the time to build the star maps

    Returns:
        None
    """
    all_program_cofs = []
    all_program_info = {}
    all_program_first_forecast = {}

    all_star_maps_all_programs = []
    all_star_colors_all_programs = []
    all_star_nObs_all_programs = []
    all_stars_in_program = []

    all_stars = []
    all_completions = []

    prog_colors = {}
    for p, item in enumerate(star_tracker.programs):
        prog_colors[p] = palette[p%len(palette)]

        program_COF, program_total_obs_request, program_first_forecast, stars_in_program, \
            stars_completion = generate_single_program_plot_suite(star_tracker, star_tracker.programs[p],
                                                                build_cadence_plots=manager.build_starmaps)
        all_program_cofs.append(program_COF)
        all_stars.append(stars_in_program)
        all_completions.append(stars_completion)

        all_program_info[star_tracker.programs[p]] = program_total_obs_request
        all_program_first_forecast[star_tracker.programs[p]] = program_first_forecast

        all_star_maps, all_star_cols, all_star_nObs, stars_in_program = \
            create_single_program_birds_eye(star_tracker, star_tracker.programs[p], prog_colors[p])
        all_star_maps_all_programs.append(all_star_maps)

        # program_color = named_colors[p%len(named_colors)]
        # all_star_colors_all_programs += [program_color]*len(all_star_maps)
        all_star_colors_all_programs += [prog_colors[p]]*len(all_star_maps)
        all_star_nObs_all_programs.append(all_star_nObs)

        all_stars_in_program.append(stars_in_program)

        generate_birds_eye(star_tracker, all_star_maps, all_star_cols, all_star_nObs,
                             star_tracker.programs[p], stars_in_program)

    flat_all_star_maps_all_programs = np.concatenate(all_star_maps_all_programs)
    flat_all_star_nObs_all_programs = np.concatenate(all_star_nObs_all_programs)
    flat_all_stars_in_program = np.concatenate(all_stars_in_program)

    generate_admin_view_plot_suite(star_tracker, all_program_cofs, all_program_info,
                                       all_program_first_forecast, flat_all_star_maps_all_programs,
                                       all_star_colors_all_programs,
                                       flat_all_star_nObs_all_programs, flat_all_stars_in_program,
                                       prog_colors)
    all_stars = np.array(all_stars).flatten()
    all_completions = np.array(all_completions).flatten()
    complete_frame = pd.DataFrame({"Starname":all_stars, 'CompletetionPercent':all_completions})
    complete_frame.to_csv(manager.semester_directory + 'reports/admin/' + str(star_tracker.manager.current_day) + \
                                                    '/cofs/final_completion_rates.csv', index=False)
    print("All plots generated. Saved to " + manager.semester_directory + "reports/")

def write_cadence_plot_files(manager):
    print("Writing cadence plot files.")
    with open(os.path.join(manager.output_directory, "raw_combined_semester_schedule_Round2.txt"), 'r', encoding='utf-8') as file:
        lines = file.readlines()
        combined_semester_schedule_stars = [line.strip().split(',') for line in lines]

    turn_on_off_frame = pd.read_csv(manager.turn_on_off_file)
    all_starmaps = {}
    for i in range(len(manager.requests_frame)):
        if manager.database_info_dict != {}:
            starmap = io.build_observed_map_past( \
                manager.database_info_dict[manager.requests_frame['Starname'][i]], manager.starmap_template_filename)
        else:
            starmap = io.build_observed_map_past([[],[],[],[]], manager.starmap_template_filename)

        starmap_updated = io.build_observed_map_future(manager, combined_semester_schedule_stars,
                            manager.requests_frame['Starname'][i], starmap)

        all_starmaps[manager.requests_frame['Starname'][i]] = starmap_updated
        future_unique_days_forecasted = 0
        for k, item in enumerate(manager.combined_semester_schedule_stars):
            if manager.requests_frame['Starname'][i] in manager.combined_semester_schedule_stars[k]:
                future_unique_days_forecasted += 1

        try:
            past_unique_dates_for_star = manager.database_info_dict[manager.requests_frame['Starname'][i]][1]
        except:
            past_unique_dates_for_star = []
        write_one_cadence_plot_file(manager.requests_frame['Starname'][i], starmap_updated,
                                    turn_on_off_frame, manager.requests_frame,
                                    future_unique_days_forecasted,
                                    past_unique_dates_for_star,
                                    manager.current_day, manager.output_directory)


def write_one_cadence_plot_file(starname, starmap, turn_frame, requests_frame, future_obs,
                            unique_hst_dates_observed, current_day, outputdir):
    """
    Write the file from which later we can produce the cadence plot

    Args:
        starname (int): the string of the star's name
        starmap (dataframe): contains information on observation history and forecast
        turn_frame (dataframe): contains information on the first and last day a target
                                 is observable in each quarter of the night. Pre-computed.
        requests_frame (dataframe): the information on the request strategies
        unique_hst_dates_observed (array): list of dates previously observed
        current_day (str): today's date, format YYYY-MM-DD
        outputdir (str): the path to the save directory
    Returns:
        None
    """
    request_id = requests_frame.index[requests_frame['Starname']==str(starname)][0]
    program_code = requests_frame.loc[request_id,'Program_Code']
    save_path = outputdir + "/cadences/" + str(program_code) + "/"
    os.makedirs(save_path, exist_ok = True)

    n_obs_desired = requests_frame.loc[request_id,'# of Nights Per Semester']
    n_obs_taken = len(unique_hst_dates_observed)
    n_obs_scheduled = future_obs #np.sum(starmap['N_obs']) #n_obs_desired - n_obs_taken #
    cadence = requests_frame.loc[request_id,'Minimum Inter-Night Cadence']
    turn_index = turn_frame.index[turn_frame['Starname']==str(starname)][0]
    turns = [[turn_frame['Q1_on_date'][turn_index], turn_frame['Q1_off_date'][turn_index]],
             [turn_frame['Q2_on_date'][turn_index], turn_frame['Q2_off_date'][turn_index]],
             [turn_frame['Q3_on_date'][turn_index], turn_frame['Q3_off_date'][turn_index]],
             [turn_frame['Q4_on_date'][turn_index], turn_frame['Q4_off_date'][turn_index]]]

    commentsfile = open(save_path + starname + "_Cadence_Interactive.csv", 'w')
    commentsfile.write('#starname:' + str(starname) + '\n')
    commentsfile.write('#programcode:' + str(program_code) + '\n')
    commentsfile.write('#Nobs_scheduled:' + str(n_obs_scheduled) + '\n')
    commentsfile.write('#Nobs_desired:' + str(n_obs_desired) + '\n')
    commentsfile.write('#Nobs_taken:' + str(n_obs_taken) + '\n')
    commentsfile.write('#cadence:' + str(cadence) + '\n')
    commentsfile.write('#q1_start:' + str(turns[0][1]) + '\n')
    commentsfile.write('#q1_end:' + str(turns[0][0]) + '\n')
    commentsfile.write('#q2_start:' + str(turns[1][1]) + '\n')
    commentsfile.write('#q2_end:' + str(turns[1][0]) + '\n')
    commentsfile.write('#q3_start:' + str(turns[2][1]) + '\n')
    commentsfile.write('#q3_end:' + str(turns[2][0]) + '\n')
    commentsfile.write('#q4_start:' + str(turns[3][1]) + '\n')
    commentsfile.write('#q4_end:' + str(turns[3][0]) + '\n')
    commentsfile.write('#current_day:' + str(current_day) + '\n')
    starmap = starmap[starmap.Allocated == True]
    starmap.to_csv(commentsfile, index=False)
    commentsfile.close()

def single_request_cof(star_tracker, star):
    """
    Compute cumulative observation function for a single request

    Args:
        star_tracker (object): a object from the StarTracker class
        star (str): the name of the star in question

    Returns:
        running_total_obs_counter (array): of length semester_length, the running number of
                                           observations over the course of the semester
        running_total_obs_percent (array): of length semester_length, the running percent of
                                           request completion over the course of the semester
        total_observations_requested (int): the number of exposure to complete the request
        first_forecast (dataframe): contains information on the first forecast of the target
    """
    requested_nobs_per_night, total_observations_requested, exposure_time, slots_per_night, \
                 program = star_tracker.get_star_stats(star)

    # to add an absolute time accounting COF, need:
    # parse past open shutter time + overhead per night
    # define expected time (incl. overhead) per night
    # define total time ask for request
    first_forecast = star_tracker.get_star_first_forecast(program, star)

    observations_past = star_tracker.get_star_past(star)
    observations_future = star_tracker.get_star_forecast(star)
    running_total_obs_counter = []
    running_total_obs_percent = []
    counter = 0
    for i in range(len(star_tracker.manager.all_dates_array)):
        if star_tracker.manager.all_dates_array[i] in list(observations_past.keys()):
            counter += observations_past[star_tracker.manager.all_dates_array[i]]
        if star_tracker.manager.all_dates_array[i] in list(observations_future.keys()):
            counter += observations_future[star_tracker.manager.all_dates_array[i]]
        running_total_obs_counter.append(counter)
        running_total_obs_percent.append(np.round((counter/total_observations_requested)*100))

    return running_total_obs_counter, running_total_obs_percent, total_observations_requested, \
            first_forecast

def programmatic_cof(star_tracker, list_of_request_cofs, total_program_observations):
    """
    Build the cumulative obesrvation function for a group of requests

    Args:
        star_tracker (object): a object from the StarTracker class
        list_of_request_cofs (array): an array of COF arrays
        total_program_observations (int): sum of all exposures for the set of requests

    Returns:
        program_cof_counter (array): of length semester_length, the running number of
                                      observations over the course of the semester
        program_cof_percent (array): of length semester_length, the running percent of
                                     request completion over the course of the semester
    """
    program_cof_counter = []
    program_cof_percent = []
    all_counter = 0
    for i in range(len(star_tracker.manager.all_dates_array)):
        day_counter = 0
        for j in range(len(list_of_request_cofs)):
            day_counter += list_of_request_cofs[j][i]
        all_counter += day_counter
        program_cof_counter.append(day_counter)
        program_cof_percent.append(np.round((day_counter/total_program_observations)*100,2))
    return program_cof_counter, program_cof_percent

def build_multi_request_cof(star_tracker, all_cofs_array, star_info, all_first_forecasts,
                            even_burn_rate=True, flag=False, color_palette={}):
    """
    Build the COF plot

    Args:
        star_tracker (object): a object from the StarTracker class
        all_cofs_array (array): an array of COF arrays
        star_info (dict): keys are the star names, values are
        all_first_forecasts (array): array of the first forecast arrays for each star in set
        even_burn_rate (boolean): if True, show an even burn rate line
        flag (boolean): if True, the plot is only for the admin view

    Returns:
        None
    """
    if flag and color_palette == {}:
        print("Must input a color_palette value when flag==True.")

    starnames = list(star_info.keys())
    lines = []
    fig = go.Figure()
    if even_burn_rate:
        burn_line = np.linspace(0, 100, len(star_tracker.manager.all_dates_array))
        for b in range(len(burn_line)):
            burn_line[b] = np.round(burn_line[b],2)
        fig.add_trace(go.Scatter(
            x=star_tracker.manager.all_dates_array,
            y=burn_line,
            mode='lines',
            line=dict(color='black', width=2, dash='dash'),
            name="Even Burn Rate",
            hovertemplate= 'Date: %{x}' + '<br>% Complete: %{y}'
        ))
    for i, y_values in enumerate(all_cofs_array):
        if flag:
            color_line = f"rgb({int(color_palette[i][0]*255)}, {int(color_palette[i][1]*255)}, {int(color_palette[i][2]*255)})"
        else:
            color_line = np.random.choice(named_colors)
        star = starnames[i]
        lines.append(str(star) + "," + str(np.round(y_values[-1],2)))
        label = star_info[star]
        fig.add_trace(go.Scatter(
            x=star_tracker.manager.all_dates_array,
            y=y_values,
            mode='lines',
            line=dict(color=color_line, width=2),
            name=star,
            hovertemplate= 'Date: %{x}' + '<br>% Complete: %{y}' + '<br># Obs Requested: ' + \
                str(label) + '<br>'
        ))
        try:
            fig.add_trace(go.Scatter(
                x=star_tracker.manager.all_dates_array,
                y=all_first_forecasts[star]['pct_complete'],
                mode='lines',
                line=dict(color=color_line, width=2, dash='dash'),
                name=star + "_First_Forecast",
                hovertemplate= 'Date: %{x}' + '<br>% Complete: %{y}' + '<br># Obs Requested: ' + \
                    str(label) + '<br>'
            ))
        except:
            continue
    fig.add_vrect(
            x0=star_tracker.manager.current_day,
            x1=star_tracker.manager.current_day,
            annotation_text="Today",
            line_dash="dash",
            fillcolor=None,
            line_width=2,
            line_color='black',
            annotation_position="bottom left"
        )
    fig.update_layout(
        title="Cumulative Observation Function (COF)",
        xaxis_title="Calendar Date",
        yaxis_title="Request % Complete",
        showlegend=True
    )
    admin_path = star_tracker.manager.reports_directory + "admin/" + star_tracker.manager.current_day + "/cofs/"
    if flag:
        fileout_path = admin_path
        fileout_name = 'all_programs_pct_by_program.html'
        fileout_name_report = 'all_programs_pct_by_program_report.txt'
    else:
        fileout_path = star_tracker.manager.reports_directory + list(star_info.keys())[0] + "/" + \
            star_tracker.manager.current_day + "/cofs/"
        fileout_name = list(star_info.keys())[0] + '_COF_pct_by_star.html'
        fileout_name_report =  list(star_info.keys())[0] + "_completion_Report.txt"
        # write copy to admin folder
        fig.write_html(admin_path + fileout_name)
        file = open(admin_path + fileout_name_report, "w")
        for l in lines:
            file.write(l + "\n")
        file.close()

    fig.write_html(fileout_path + fileout_name)
    file = open(fileout_path + "completion_Report.txt", "w")
    for l in lines:
        file.write(l + "\n")
    file.close()


def generate_cadence_view_plot(star_tracker, program, star):
    """
    Build cadence overview plot

    Args:
        program (str): the program code
        starname (str): the name of the star in question

    Returns:
        None
    """
    # Read the file and process data to generate the plot
    # Example: Read data from filename and generate a scatter plot
    # Replace this with your actual data processing code
    # For this example, we'll generate a dummy scatter plot
    # files should always have the following format (with {} indicating values):
    # line1 --- "starname:{string}"
    # line2 --- "programcode:{string}"
    # line3 --- "Nobs_scheduled:{int}"
    # line4 --- "Nobs_desired:{int}"
    # line5 --- "Nobs_taken:{int}"
    # line6 --- "cadence:{int}"
    # line7 --- "q1_start:{string}"
    # line8 --- "q1_end:{string}"
    # line9 --- "q2_start:{string}"
    # line10 --- "q1_end:{string}"
    # line11 --- "q3_start:{string}"
    # line12 --- "q3_end:{string}"
    # line13 --- "q4_start:{string}"
    # line14 --- "q4_end:{string}"
    # line15 --- "current_day:{string}"
    # then a csv with columns: ['Date', 'Observed', 'Quarter', 'Weathered', '']
    file_in = star_tracker.cadence_files_dir + program + "/" + star + "_Cadence_Interactive.csv"
    data_fig = pd.read_csv(file_in, comment='#')
    parameter_dict = {}
    with  open(file_in,'r') as cmt_file:
        for line in cmt_file:
            if line[0] == '#':
                line = line[1:]
                para = line.split(':')
                if len(para) == 2:
                    parameter_dict[ para[0].strip()] = para[1].strip()

    # use data file to determine which unique nights got observed, for "6th" quarter row up top
    observed_this_night = []
    holders = []
    dates_on_sky = data_fig['Date'].unique().tolist()
    for d in range(len(dates_on_sky)):
        today = data_fig[data_fig['Date'] == dates_on_sky[d]]
        summation = np.sum(today['Observed'])
        if summation > 0:
            observed_this_night.append(True)
        else:
            observed_this_night.append(False)
        holders.append(5.5)
    sumholder = [np.sum(observed_this_night)]*len(dates_on_sky)
    outs = pd.DataFrame({'Dates':dates_on_sky, "Gotten":observed_this_night, "Holders":holders, \
                            "SumHolder":sumholder})

    # plot the dots of the data
    # for some reason if the first element is "True", then the "True" values are plotted in the
    # color of the first color of the list
    # So I just check if this is satisfied and if so, swap the order of the colors
    if data_fig['Observed'][0]:
        fig = px.scatter(data_fig, x="Date", y="Quarter", color='Observed',
            color_discrete_sequence=["red", "black"], hover_data=["N_obs"],
            title='Cadence Plot for ' + str(parameter_dict['starname']) + ' in Program ' + \
            str(parameter_dict['programcode']) + ", scheduled for " + \
            str(int(parameter_dict['Nobs_scheduled']) + int(parameter_dict['Nobs_taken'])) + \
                    " obs out of " + \
            str(parameter_dict['Nobs_desired']) + " requested (" + \
            str(parameter_dict['Nobs_taken']) + " already taken), at " + \
            str(parameter_dict['cadence']) + " day cadence.")
    else:
        fig = px.scatter(data_fig, x="Date", y="Quarter", color='Observed',
            color_discrete_sequence=["black", "red"], hover_data=["N_obs"],
            title='Cadence Plot for ' + str(parameter_dict['starname']) + ' in Program ' + \
            str(parameter_dict['programcode']) + ", scheduled for " + \
            str(int(parameter_dict['Nobs_scheduled']) + int(parameter_dict['Nobs_taken'])) + \
                    " obs out of " + \
            str(parameter_dict['Nobs_desired']) + " requested (" + \
            str(parameter_dict['Nobs_taken']) + " already taken), at " + \
            str(parameter_dict['cadence']) + " day cadence.")
    fig.update_traces(marker={'size': 10})
    fig.update_layout(yaxis = dict(tickmode = 'linear',tick0 = 0, dtick = 1))

    # plot the lines atop the dots
    for j in range(len(data_fig)):
        if data_fig['Observed'][j]:
            fig.add_shape(type="line", x0=data_fig['Date'][j], x1=data_fig['Date'][j],
                            y0=data_fig['Quarter'][j]-0.5, y1=data_fig['Quarter'][j]+0.5,
                            line=dict(color="red",width=3))
        else:
            fig.add_shape(type="line", x0=data_fig['Date'][j], x1=data_fig['Date'][j],
                            y0=data_fig['Quarter'][j]-0.5, y1=data_fig['Quarter'][j]+0.5,
                            line=dict(color="black",width=3))

    shown = False
    for k in range(len(data_fig)):
        if data_fig['Weathered'][k]:
            if shown == False:
                fig.add_shape(showlegend=True, type="line", x0=data_fig['Date'][k],
                    x1=data_fig['Date'][k], y0=data_fig['Quarter'][k]-0.5,
                    y1=data_fig['Quarter'][k]+0.5, line=dict(color="yellow",width=3),
                    label=dict(text="Weathered", font=dict(size=1, color="yellow")))
                shown = True
            else:
                fig.add_shape(type="line", x0=data_fig['Date'][k],
                    x1=data_fig['Date'][k], y0=data_fig['Quarter'][k]-0.5,
                    y1=data_fig['Quarter'][k]+0.5, line=dict(color="yellow",width=3))

    # plot the "6th" quarter row
    for j in range(len(outs)):
        if outs['Gotten'][j]:
            fig.add_shape(type="line", x0=outs['Dates'][j], x1=outs['Dates'][j],
                y0=outs['Holders'][j]-0.5, y1=outs['Holders'][j]+0.5,
                line=dict(color="red",width=3))
        else:
            fig.add_shape(type="line", x0=outs['Dates'][j], x1=outs['Dates'][j],
                y0=outs['Holders'][j]-0.5, y1=outs['Holders'][j]+0.5,
                line=dict(color="black",width=3))

    # Add the green boxes of accessibility
    fig.add_shape(type="rect", x0=parameter_dict['q1_start'], x1=parameter_dict['q1_end'],
        y0=0, y1=1, fillcolor='lime', opacity=0.3, showlegend=True, name='Yes',
        legendgrouptitle=dict(text='Observable'))
    fig.add_shape(type="rect", x0=parameter_dict['q2_start'], x1=parameter_dict['q2_end'],
        y0=1, y1=2, fillcolor='lime', opacity=0.3)
    fig.add_shape(type="rect", x0=parameter_dict['q3_start'], x1=parameter_dict['q3_end'],
        y0=2, y1=3, fillcolor='lime', opacity=0.3)
    fig.add_shape(type="rect", x0=parameter_dict['q4_start'], x1=parameter_dict['q4_end'], y0=3,
        y1=4, fillcolor='lime', opacity=0.3)

    # Add the dashed line for current day
    fig.add_vline(x=parameter_dict['current_day'], line_dash='dash',
        line=dict(color="black",width=1,))
    fig.add_annotation(x=parameter_dict['current_day'], y=0.1, text="Today",
        showarrow=False,yshift=10)

    admin_path = star_tracker.manager.reports_directory + "admin/" + star_tracker.manager.current_day + "/cadences/"
    program_path = star_tracker.manager.reports_directory + program + "/" + star_tracker.manager.current_day + "/cadences/"
    fileout_name = star + '_cadence_overview.html'
    fig.write_html(admin_path + fileout_name)
    fig.write_html(program_path + fileout_name)

def process_single_target_for_birds_eye(star_tracker, starname):
    """
    Collect the data for a birds eye view plot

    Args:
        star_tracker (object): a object from the StarTracker class
        starname (str): the name of the star in question

    Returns:
        starmap (array): n_nights_in_semester by n_slots_in_night, 1 when target is scheduled
        nobs (int): the sum of all future observations
    """
    schedule_array = np.array(star_tracker.forecast)
    consolidate = []
    for i in range(len(schedule_array)):
        if schedule_array[i][0] != 'Past':
            consolidate.append(np.array(schedule_array[i]))
    star_tracker.forecast_no_Past = np.array(consolidate)
    starmap = np.zeros(np.shape(star_tracker.forecast_no_Past))
    indices = np.where(star_tracker.forecast_no_Past == starname)
    if len(indices[0]) != 0:
        coords = list(zip(indices[0], indices[1]))
        for i in range(len(coords)):
            starmap[coords[i][0], coords[i][1]] = 1
    nobs = np.sum(np.sum(starmap,axis=1))
    return starmap, nobs

def create_single_program_birds_eye(star_tracker, program, color):
    """
    Prepare data for a birds eye view plot of one program's requests

    Args:
        star_tracker (object): a object from the StarTracker class
        program (str): the program code
        color (object): seaborn color palette value for this program

    Returns:
        all_star_maps (array): all of the starmaps of requests in the program
        all_star_cols (array): the colors to plot each request on the birds eye
        all_star_nobs (array): the number of observations in each scheduled window
        stars_in_program (array): the list of the stars within the program
    """
    all_star_maps = []
    all_star_cols = []
    all_star_nobs = []
    stars_in_program = []
    color_counter = 0
    for i in range(len(star_tracker.manager.requests_frame)):
        if star_tracker.manager.requests_frame['Program_Code'][i] == program:
            starname = star_tracker.manager.requests_frame['Starname'][i]
            stars_in_program.append(starname)
            requested_nobs_per_night, total_observations_requested, exposure_time, slots_per_night,\
                program = star_tracker.get_star_stats(starname)
            starmap, nobs = process_single_target_for_birds_eye(star_tracker, starname)
            all_star_maps.append(starmap)
            all_star_cols.append(color)
            # all_star_cols.append(named_colors[color_counter%len(named_colors)])
            all_star_nobs.append(nobs/slots_per_night)
            color_counter += 1

    return all_star_maps, all_star_cols, all_star_nobs, stars_in_program

def generate_birds_eye(star_tracker, all_star_maps, all_star_cols, all_star_nobs, program,
                    stars_in_program, nonqueue=False):
    """
    Plot a birds eye view of one program's requests

    Args:
        star_tracker (object): a object from the StarTracker class
        all_star_maps (array): all of the starmaps of requests in the program
        all_star_cols (array): the colors to plot each request on the birds eye
        all_star_nobs (array): the number of observations in each scheduled window
        program (str): the program code
        stars_in_program (array): the list of the stars within the program
        nonqueue (boolean): if True, plot the non-queue times on the birds eye view
    Returns:
        None
    """
    fig = go.Figure()
    for i in range(len(all_star_maps)):
        all_star_maps_transpose = np.array(all_star_maps[i]).T
        # col = all_star_cols[i]
        col = f"rgb({int(all_star_cols[i][0]*255)}, {int(all_star_cols[i][1]*255)}, {int(all_star_cols[i][2]*255)})"
        nobscounter = [str(all_star_nobs[i])]
        fig.add_trace(go.Heatmap(
            z = all_star_maps_transpose,
            colorscale=[[0, 'rgba(0,0,0,0)'],[1, col]],
            # colorscale=color_scales[col],
            zmin=0, zmax=1,
            opacity=1.0,
            showscale=False,
            name=stars_in_program[i],
            text=nobscounter,
            hovertemplate='<b>' + str(stars_in_program[i]) +
                '</b><br><b>Date: %{y}</b><br><b>Slot: %{x}</b><br>Forecasted N_Obs: ' + \
                nobscounter[0] + '<extra></extra>',
            showlegend=True,
        ))
    if nonqueue:
        for i in range(len(star_tracker.nonqueue)):
            starname = 'RM__' + star_tracker.nonqueue['Starname'][i]
            starmap, nobs = process_single_target_for_birds_eye(star_tracker, starname)
            starttime = Time(star_tracker.nonqueue['Start'][i])
            stoptime = Time(star_tracker.nonqueue['Stop'][i])

            starmap_transpose = np.array(starmap).T
            fig.add_trace(go.Heatmap(
                z = starmap_transpose,
                colorscale=[[0, 'rgba(0,0,0,0)'], [1, 'rgb(0, 0, 0)']],
                zmin=0, zmax=1,
                opacity=1.0,
                showscale=False,
                name=star_tracker.nonqueue['Starname'][i],
                text=[str(nobs)],
                hovertemplate='<b>' + str(star_tracker.nonqueue['Starname'][i]) +
                    '</b><br><b>Date: %{y}</b><br><b>Slot: %{x}</b><br>Window: ' + \
                    str(starttime) + " to " + str(stoptime) + '<extra></extra>',
                showlegend=True,
            ))
    fig.update_layout(
        title="Birds Eye View Semester Schedule",
        yaxis_title="Slot in Night",
        xaxis_title="Night in Semester",
        yaxis=dict(tickvals=np.arange(0, np.shape(all_star_maps)[2], 10)),
        xaxis=dict(tickvals=np.arange(0, np.shape(all_star_maps)[1], 7)),
        template="plotly",
        showlegend=True,
    )
    admin_path = star_tracker.manager.reports_directory + "admin/" + star_tracker.manager.current_day + "/birds_eyes/"
    fileout_path = star_tracker.manager.reports_directory + program + "/" + star_tracker.manager.current_day+ "/birds_eyes/"
    if program == 'admin':
        fileout_name = 'all_programs_birds_eyes.html'
    else:
        fileout_name = program + '_birds_eyes.html'
        # write copy to admin folder
        fig.write_html(admin_path + fileout_name)
    fig.write_html(fileout_path + fileout_name)

def generate_admin_view_plot_suite(star_tracker, set_of_program_cofs, program_info,
                                   all_program_first_forecasts, all_star_maps_all_programs,
                                   all_star_colors_all_programs, all_star_nobs_all_programs,
                                   all_stars_in_program, prog_colors):
    """
    Construct all plots for admin view

    Args:
        star_tracker (object): a object from the StarTracker class
        set_of_program_cofs (array): each element is an array containing the cof of a program
        program_info (dict): keys are program codes and values are total obs for program
        all_program_first_forecasts (array): each element is an array containing the first forecast
                                            for the program, cof.
        all_star_maps_all_programs (array): all of the programmatic running total percentages
        all_star_colors_all_programs (array): the colors to plot each program on the birds eye
        all_star_nobs_all_programs (array): the number of observations in each scheduled window
        all_stars_in_program (array): the list of the stars within the program
    Returns:
        None
    """
    print("Building admin view of COFs.")
    build_multi_request_cof(star_tracker, set_of_program_cofs, program_info,
            all_first_forecasts=all_program_first_forecasts, flag=True, color_palette=prog_colors)
    print("Building admin view of Birds Eye")
    generate_birds_eye(star_tracker, all_star_maps_all_programs, all_star_colors_all_programs,
            all_star_nobs_all_programs, 'admin', all_stars_in_program, nonqueue=True)

def generate_single_program_plot_suite(star_tracker, program_code, build_cadence_plots=False):
    """
    Plot the birds eye view for the single program PI view

    Args:
        star_tracker (object): a object from the StarTracker class
        program_code (str): the program code
    Returns:
        program_cof_percent (array): the calculated cof array for each star in program
        total_program_requested_nobs (int): the total requested observations within the program
        program_first_forecast (array): the first forecast for each star in program
    """
    stars_completion = []
    stars_in_program = []
    stars_cofs_counters = []
    stars_cofs_percents = []
    stars_cofs_first_forecasts = {}
    stars_info = {}
    print("Building individual program report for ", program_code)
    for i in range(len(star_tracker.manager.requests_frame['Program_Code'])):
        if star_tracker.manager.requests_frame['Program_Code'][i] == program_code:
            star = star_tracker.manager.requests_frame['Starname'][i]
            if star_tracker.manager.requests_frame['# of Nights Per Semester'][i] != 0:
                stars_in_program.append(star)
                running_total_obs_counter, running_total_obs_percent, total_observations_requested,\
                    first_forecast = single_request_cof(star_tracker, star)
                stars_completion.append(running_total_obs_percent[-1])
                stars_cofs_counters.append(running_total_obs_counter)
                stars_cofs_percents.append(running_total_obs_percent)
                stars_cofs_first_forecasts[star] = first_forecast
                stars_info[star] = total_observations_requested
                if build_cadence_plots:
                    generate_cadence_view_plot(star_tracker, program_code, star)
    total_program_requested_nobs = np.sum(list(stars_info.values()))
    program_cof_counter, program_cof_percent = programmatic_cof(star_tracker, stars_cofs_counters,
        total_program_requested_nobs)
    stars_cofs_counters.append(program_cof_counter)
    stars_cofs_percents.append(program_cof_percent)
    stars_info[program_code] = np.sum(list(stars_info.values()))
    program_first_forecast = star_tracker.get_star_first_forecast(program_code, program_code)
    stars_cofs_first_forecasts[program_code] = program_first_forecast

    # reverse ordering so the program code always appears at the top of the plotly figure
    stars_cofs_counters.reverse()
    stars_cofs_percents.reverse()
    stars_info = dict(reversed(list(stars_info.items())))
    build_multi_request_cof(star_tracker, stars_cofs_percents, stars_info,
        all_first_forecasts=stars_cofs_first_forecasts)

    return program_cof_percent, total_program_requested_nobs, program_first_forecast, \
            stars_in_program, stars_completion
