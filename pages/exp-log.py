import dash
from dash import dcc, html, callback
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from dash import dash_table
import pandas as pd
import re
from datetime import datetime
import plotly.graph_objects as go
import plotly.express as px

dash.register_page(__name__,name='Experiment Log')
       

# Define the layout of the app
layout = html.Div([
    html.H2("Log Data Summary"),
    html.H3("Log File Path"),
    html.Div(id = "log-file"),
    #html.H3("Autofocus Info"),
    #html.Div(id = "autofocus-time"),
    #dash_table.DataTable(
    #    id='autofocus-pts',
     #   fill_width=False,
     #   style_table={'overflowX': 'auto'},
      #  style_cell={
      #      'height': 'auto',
            # all three widths are needed
      #      'minWidth': '200px', 'width': '200px', 'maxWidth': '200px',
      #      'whiteSpace': 'normal'},
    #),
    html.H3("Experiment Priming"),
    dash_table.DataTable(
        id='priming',
        fill_width=False,
        style_table={'overflowX': 'auto'},
        style_cell={
            'height': 'auto',
            # all three widths are needed
            'minWidth': '300px', 'width': '300px', 'maxWidth': '300px',
            'whiteSpace': 'normal'},
    ),
    html.H3("Laser Powers"),
    dcc.Graph(id='lasers'),
    html.H3("Experiment Temperature"),
    dcc.Graph(id='temperature'),
    
])

## Store log file path
@callback(Output("log-file","children"),
          Input("exp_dir","data"))
def log_file_path(exp_dir):
    return f'{exp_dir}/logs/20211207_AlternatingSABER.log'


## Display priming data
@callback(Output('priming', 'data'),
          Output('priming', 'columns'),
          Input('log-file','children'),
)
def update_table(log_file_path):
    # Read the log file line by line
    with open(log_file_path, 'r', encoding='iso-8859-1') as file:
        lines = file.readlines()

    # Define a regular expression pattern to extract timestamp and reagent information
    pattern = re.compile(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}) - Priming (\w+)')

    # Create a list to store extracted data
    log_data = []

    # Iterate over each log entry and extract information
    for line in lines:
        match = pattern.match(line)

        if match:
            timestamp_str, reagent = match.groups()
            timestamp = datetime.strptime(timestamp_str, "%Y-%m-%d %H:%M:%S,%f")
            log_data.append({"Timestamp": timestamp, "Reagent": reagent})

    df = pd.DataFrame(log_data)                 
    return df.to_dict('records'), [{"name": i, "id": i} for i in df.columns]

## Line graphs of laser power vs time
@callback(Output('lasers','figure'),
          Input('log-file','children'))
def laser_power_graph(log_file_path):
    with open(log_file_path, 'r', encoding='iso-8859-1') as file:
        lines = file.readlines()
    # Define a regular expression pattern
    pattern = r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}) - (\w+)::rcvd:: (\d+)mW'

    powers = {}

    for line in lines:

        # Match the pattern in the log entry
        match = re.match(pattern, line)

        if match:
            # Extract groups from the match
            timestamp_str, color, power_str = match.groups()

            # Convert timestamp string to datetime object
            timestamp = datetime.strptime(timestamp_str, "%Y-%m-%d %H:%M:%S,%f")
            # Convert power string to integer
            power = int(power_str)

            if color not in powers:
                powers[color]=[]
            powers[color].append({"Timestamp": timestamp, "Power (mW)": power})

    fig = go.Figure()
        
    for color in powers:
        df = pd.DataFrame(powers[color])
        fig.add_trace(go.Scatter(x=df['Timestamp'], y=df['Power (mW)'], mode='lines', name=color,line=dict(color=color[:-5])))
        
    # Customize layout
    fig.update_layout(title='Laser Power Over Time',
                      xaxis_title='Time',
                      yaxis_title='Power (mW)')

    return fig

## Temperature log
@callback(Output('temperature','figure'),
          Input('log-file','children'))
def laser_power_graph(log_file_path):
    with open(log_file_path, 'r', encoding='iso-8859-1') as file:
        lines = file.readlines()
    # Initialize lists to store data for each cycle
    cycle_data = {}

    # Define a regular expression pattern to match the cycle number
    cycle_pattern = r"cycle(\d+)"

    # Regular expression pattern to match timestamps and temperatures
    pattern = r"(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}) - PySeq::A::Temperature:: (\d+\.\d+) °C"

    # Iterate through lines to extract data
    for line in lines:
        # Check cycle
        match = re.search(cycle_pattern, line)

        # Check if a match is found
        if match:
            cycle = match.group(1)
            if cycle not in cycle_data and cycle is not None:
                cycle_data[cycle]= {'Timestamp': [], 'Temperature': []}      

        # Add temperature data for each cycle
        match = re.search(pattern, line)
        if match:
            timestamp_str = match.group(1)
            temperature = float(match.group(2))
            
            timestamp = datetime.strptime(timestamp_str, "%Y-%m-%d %H:%M:%S,%f")
            cycle_data[cycle]['Timestamp'].append(timestamp)
            cycle_data[cycle]['Temperature'].append(temperature)
            

    fig = go.Figure()
    
    for cycle in cycle_data:
        df = pd.DataFrame(cycle_data[cycle])
        fig.add_trace(go.Scatter(x=df['Timestamp'], y=df['Temperature'], mode='lines', name=f'Cycle {cycle}'))
    # Customize layout
    fig.update_layout(title='Temperature Log',
                      xaxis_title='Time',
                      yaxis_title='Temperature (°C)')

    return fig


'''    
@callback(Output('autofocus-time','children'),
          Output('autofocus-pts','data'),
          Output('autofocus-pts','columns'),
          Input('log-file','children'))
def autofocus_time(log_file_path):
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
    '''
    