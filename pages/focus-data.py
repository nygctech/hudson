import dash
from dash import dcc, html, callback
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
import plotly.express as px
import pandas as pd
import os
from pathlib import Path
import glob

dash.register_page(__name__,name='Focus Data')
       

# Define the layout of the app
layout = html.Div([
    html.H2("Focus Data"),
    "Section",
    dcc.Dropdown(id = 'focus-section',placeholder="Choose a section"),
    dcc.Graph(id='focus-data-plot'),

])

@callback(Output("focus-section","options"),
          Input("exp_dir","data"))
def focus_sections(exp_dir):
    sections = [f.name[0:12] for f in os.scandir(Path(exp_dir)/'logs'/'focus_data') if f.name[-4:] == '.txt']
    return list(dict.fromkeys(sections))

# Function to read data from multiple files
def read_data_from_files(file_pattern):
    all_data = []
    i = 1
    for file_path in glob.glob(file_pattern):
        df = pd.read_csv(file_path, delim_whitespace=True, names=['Objective_Position', 'Focus_Metric'])
        #section_id = file_path.split('_')[0]
        df['Series'] = str(i)
        #df['Series'] = pd.Categorical(df['Series'])
        i += 1
        all_data.append(df)
    return pd.concat(all_data, ignore_index=True)

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
    df_all = read_data_from_files(file_pattern)
    
    fig = px.scatter(df_all, x='Objective_Position', y='Focus_Metric', color = 'Series',title='Focus Metric vs Objective Position')
    return fig
