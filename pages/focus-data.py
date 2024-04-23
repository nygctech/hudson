import dash
from dash import dcc, html, callback
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
import plotly.express as px
import pandas as pd
import os
from pathlib import Path
import glob
import numpy as np
import re
from datetime import datetime
from dash import dash_table
from dash_canvas import DashCanvas

dash.register_page(__name__,name='Auto Focus Data')
       
# Define the layout of the app
layout = html.Div([
    html.H2("Auto Focus Data"),
    html.H3("Autofocus Info"),
    html.Div(id = "autofocus-time"),
    dash_table.DataTable(
        id='autofocus-pts',
        fill_width=False,
        style_table={'overflowX': 'auto'},
        style_cell={
            'height': 'auto',
            # all three widths are needed
            'minWidth': '200px', 'width': '200px', 'maxWidth': '200px',
            'whiteSpace': 'normal'},
    ),
    "Section",
    dcc.Dropdown(id = 'focus-section',placeholder="Choose a section"),
    html.P("Double click a series number to isolate a curve."), 
    dcc.Graph(id='focus-data-plot'),
    html.Div([
        DashCanvas(id='autofocus-image',
            tool='line',
            lineWidth=5,
            lineColor='red',
            ),
    

])])

@callback(Output("focus-section","options"),
          Input("exp_dir","data"))
def focus_sections(exp_dir):
    sections = [f.name[0:12] for f in os.scandir(Path(exp_dir)/'logs'/'focus_data') if f.name[-4:] == '.txt']
    return list(dict.fromkeys(sections))

# Function to read data from multiple focus files
def read_data_from_files(file_pattern):
    all_data = []
    focus_nps = {}
    
    i = 1
    
    for file_path in glob.glob(file_pattern):
        df = pd.read_csv(file_path, delim_whitespace=True, names=['Objective_Position', 'Focus_Metric'])
        #section_id = file_path.split('_')[0]
        df['Series'] = str(i)
        #df['Series'] = pd.Categorical(df['Series'])
        
        # make np array from each focus file 
        with open(file_path, 'r') as file:
            # Read lines from the file
            lines = file.readlines()

        # Initialize an empty list to store the data
        data = []

        # Iterate over each line
        for line in lines:
            # Split the line into columns based on the delimiter (e.g., space, comma, tab)
            columns = line.split(" ")  # Change the split delimiter if needed
            # Convert the columns to floats and append to the data list
            data.append([float(columns[0]), float(columns[1])])

        # Convert the data list to a numpy array
        data_array = np.array(data)
        
        focus_nps[i] = data_array
        i += 1
        all_data.append(df)
    return pd.concat(all_data, ignore_index=True),focus_nps



### Kunal's code to plot FOVs on round 1 image ###

def get_fovs(exp_dir, section_name):
    files = os.listdir(Path(exp_dir) / 'exp_conf' / section_name)
    fovs_ = []
    for f in files: 
        id_, r, c = _.split('_')
        r = int(r[1:])/scale
        c = (int(c.split('.')[0][1:]))/scale
        fovs_.append((r,c))
    return fovs_

def focus_position_map(exp_dir, section_name):
    #fig2 = px.imshow(image_small**gamma, width=len(image_small.col), height = len(image_small.row), aspect = 'equal', color_continuous_scale = 'magma', zmin = 125**gamma, zmax = 3000**gamma)
    # TODO
    # get preview image from 1st cycle and plot
    # get how much the preview image was downscaled

    scale = 4
    width = 2048 / scale
    height = 16 / scale
    r_offset = 64 / scale
    c_offset = 0 / scale
    # Add rectangles to the figure
    for r, c in get_fovs(exp_dir, section_name):

        fig2.add_shape(
            type="rect",
            x0=c+c_offset, y0=r+r_offset, x1=(c+c_offset+width), y1=(r+r_offset+height),
            #line=rect['line']
        )
    
        # Configure layout
        fig2.update_layout(
            title=section_name,
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            showlegend=False
    )


### End of Kunal's code #####

# define mixed gaussian func
from scipy.optimize import least_squares

def fit_mixed_gaussian(data):
    """Fit focus data & return optimal objective focus step.Focus objective step vs frame JPEG file size is fit to a mixed
       gaussian model. The optimal objective focus step is returned at step
       of the max fit JPEG file. If the objective focus step is not
       returned, False is returned.

       **Parameters:**
        - data (array Nx2): Focus data where the 1st column is the objective
          step and the 2nd column is the corresponding file size.

       **Returns:**
        - int: Optimal focus objective step if found (if not, False).

    """
    focus_start = 2000
    focus_stop = 62000
    nyquist_obj = 235


    name_ = 'FitMixedGaussian::'
    #hs = self.hs

    # initialize values
    max_peaks = 4
    # Initialize varibles
    amp = []; amp_lb = []; amp_ub = []
    cen = []; cen_lb = []; cen_ub = []
    sigma = []; sigma_lb = []; sigma_ub = []
    y = data[:,1]
    SST = np.sum((y-np.mean(y))**2)
    R2 = 0
    tolerance = 0.9                                                             # Tolerance for R2
    xfun = data[:,0]; yfun = data[:,1]

    # Add peaks until fit reaches threshold
    while len(amp) <= max_peaks and R2 < tolerance:
        # set initial guesses
        max_y = np.max(y)
        amp.append(max_y*10000)
        index = np.argmax(y)
        y = np.delete(y, index)
        index = np.where(data[:,1] == max_y)[0][0]
        cen.append(data[index,0])
        sigma.append(np.sum(data[:,1]**2)**0.5*10000)
        p0 = np.array([amp, cen, sigma])
        p0 = p0.flatten()

        # set bounds
        amp_lb.append(0); amp_ub.append(np.inf)
        cen_lb.append(np.min(data[:,0])); cen_ub.append(np.max(data[:,0]))
        sigma_lb.append(0); sigma_ub.append(np.inf)
        lo_bounds = np.array([amp_lb, cen_lb, sigma_lb])
        up_bounds = np.array([amp_ub, cen_ub, sigma_ub])
        lo_bounds = lo_bounds.flatten()
        up_bounds = up_bounds.flatten()

        # Optimize parambeters
        results = least_squares(res_gaussian, p0, bounds=(lo_bounds,up_bounds),
                                args=(data[:,0],data[:,1]))

        if results.success:
            R2 = 1 - np.sum(results.fun**2)/SST
            #self.message(False,name_,'R2=',R2,'with',len(amp),'peaks')


        if results.success and R2 > tolerance:
            _objsteps = range(focus_start, focus_stop,
                              int(nyquist_obj/2))
            _focus = gaussian(_objsteps, results.x)
            optobjstep = int(_objsteps[np.argmax(_focus)])
            if optobjstep in (focus_start, focus_stop):
                #self.message(False, name_, 'Peak at endpoint: ', optobjstep)
                optobjstep = False
        else:
            optobjstep = False
            #if len(amp) == max_peaks:
                #self.message(False, name_, 'Bad fit')
                #break

    return optobjstep

# define gaussian
def gaussian(x, *args):
    """Gaussian function for curve fitting."""

    name_ = 'Gaussian::'
    if len(args) == 1:
        args = args[0]

    n_peaks = int(len(args)/3)


    if len(args) - n_peaks*3 != 0:
        print('Unequal number of parameters')
    else:
        for i in range(n_peaks):
            amp = args[0:n_peaks]
            cen = args[n_peaks:n_peaks*2]
            sigma = args[n_peaks*2:n_peaks*3]

        g_sum = 0
        for i in range(len(amp)):
            g_sum += amp[i]*(1/(sigma[i]*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen[i])/sigma[i])**2)))

        return g_sum
    
def res_gaussian(args, xfun, yfun):
    """Gaussian residual function for curve fitting."""

    g_sum = gaussian(xfun, args)

    return yfun-g_sum

@callback(Output("focus-data-plot","figure"),
          Input("exp_dir","data"),
          Input("focus-section","value"))
def focus_plot(exp_dir,section):
    if section == None:
        raise PreventUpdate
        
    # Specify the directory path and file pattern
    directory_path = f'{exp_dir}/logs/focus_data/'
    file_pattern = f'{directory_path}{section}_r*_c*.txt'

    # Read data from multiple files
    df_all,focus_files = read_data_from_files(file_pattern)
    
    # Find the mixed Gaussian peaks

    for key in focus_files:
        peak = fit_mixed_gaussian(focus_files[key]) 

        # New row to be added
        new_row = pd.DataFrame({'Series': [str(key)], 'Objective_Position': [peak], 'Focus_Metric': [0.085]})
        # Concatenate the existing DataFrame with the new row
        df_all = pd.concat([df_all, new_row], ignore_index=True)
        

    fig = px.scatter(df_all, x='Objective_Position', y='Focus_Metric', color = 'Series',title='Focus Metric vs Objective Position')
    
    # Add a solid line at y=0.085
    fig.add_hline(y=0.085, line_dash="solid", line_color="grey")

    
    return fig

    
@callback(Output('autofocus-time','children'),
          Output('autofocus-pts','data'),
          Output('autofocus-pts','columns'),
          Input('exp_dir','data'))
def autofocus_time(exp_dir):
    log_file_path = f'{exp_dir}/logs/20211207_AlternatingSABER.log'
    with open(log_file_path, 'r', encoding='iso-8859-1') as file:
        lines = file.readlines()

    data = []
    total_minutes = 0
    autofocus_completed_count = 0
    for line in lines:
        match_time = re.search(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}),.*?Completed in (\d+) minutes', line)
        match_pt = re.search(r'Autofocus::GetFocusData:: (Good|Bad)\s+point:: \[\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\]', line)
    
        if match_time:
            minutes = int(match_time.group(2))
            if minutes > 0:
                total_minutes += minutes
                autofocus_completed_count += 1
        elif match_pt:
            point_type = match_pt.group(1)
            x = int(match_pt.group(2))
            y = int(match_pt.group(3))
            median_focus_step = int(match_pt.group(4))
            index = int(match_pt.group(5))
            data.append({'Index': index, 'Point Type': point_type, 'X': x, 'Y': y, 'Median Focus Step': median_focus_step})
    # Convert dictionary to DataFrame
    df = pd.DataFrame.from_dict(data)
    
    # Drop duplicates to keep only the latest entry for each index
    df = df.drop_duplicates(subset='Index', keep='last')
    
    return f'Total time taken to autofocus: {total_minutes} minutes. \n\nNumber of autofocus attempts: {autofocus_completed_count}', df.to_dict('records'), [{"name": i, "id": i} for i in df.columns]
