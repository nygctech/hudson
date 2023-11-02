import dash
from dash import dcc, ctx, html, callback, Input, Output, State 
import plotly.express as px
import numpy as np
import plotly.express as px
import plotly.graph_objs as go
from pathlib import Path
import os
from dash.exceptions import PreventUpdate

dash.register_page(__name__,name = 'Channel Mixing Analysis')

import yaml

    
def wavelength_to_rgb(wavelength):
    """
    Convert a wavelength in the visible spectrum to an RGB color.

    :param wavelength: Wavelength in nanometers (approximately 380 to 750 nm).
    :return: RGB color as a tuple of three values between 0 and 255.
    """
    gamma = 0.8
    intensity_max = 255
    factor = 0

    if 380 <= wavelength < 440:
        red = -(wavelength - 440) / (440 - 380)
        green = 0.0
        blue = 1.0
    elif 440 <= wavelength < 490:
        red = 0.0
        green = (wavelength - 440) / (490 - 440)
        blue = 1.0
    elif 490 <= wavelength < 510:
        red = 0.0
        green = 1.0
        blue = -(wavelength - 510) / (510 - 490)
    elif 510 <= wavelength < 580:
        red = (wavelength - 510) / (580 - 510)
        green = 1.0
        blue = 0.0
    elif 580 <= wavelength < 645:
        red = 1.0
        green = -(wavelength - 645) / (645 - 580)
        blue = 0.0
    elif 645 <= wavelength <= 750:
        red = 1.0
        green = 0.0
        blue = 0.0
    else:
        red = 0.0
        green = 0.0
        blue = 0.0

    # Adjust intensity
    if 380 <= wavelength < 420:
        factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380)
    elif 420 <= wavelength < 645:
        factor = 1.0
    elif 645 <= wavelength <= 750:
        factor = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)

    r = int(intensity_max * (red * factor) ** gamma)
    g = int(intensity_max * (green * factor) ** gamma)
    b = int(intensity_max * (blue * factor) ** gamma)

    return r, g, b

layout = html.Div([
            html.H3('Channel Mixing Analysis'),
            dcc.Dropdown(
                id="section-yaml",
                placeholder='Section',
            ),
            dcc.Dropdown(
                id="cycle-dropdown",
                placeholder='Cycle',
            ),
            dcc.Graph(id='square-plot'),
])

#dropdown to select section of interest
@callback(Output('section-yaml','options'),
          Input('exp_dir','data'))
def update_section_yaml(exp_dir):
    yaml_files = [f.name for f in os.scandir(Path(exp_dir)/'final_zarr') if f.name[-5:] == '.yaml']
    
    return yaml_files

@callback(Output('cycle-dropdown','options'),
          Input('section-yaml','value'),
          Input('exp_dir','data'))
def update_cycle_dropdown(section,exp_dir):
    if section is None:
        raise PreventUpdate
        
    yaml_file_path = f'{exp_dir}/final_zarr/{section}'
    with open(yaml_file_path, 'r') as file:
        yaml_data = yaml.safe_load(file)
    return list(yaml_data.keys())

@callback(
    Output('square-plot', 'figure'),
    Input('cycle-dropdown', 'value'),
    Input('section-yaml','value'),
    Input('exp_dir','data'),
)
def update_square_plot(cycle,section,exp_dir):
    if section is None:
        raise PreventUpdate
    elif cycle is None:
        raise PreventUpdate
    
    # Specify the path to your YAML file
    yaml_file_path = Path(exp_dir)/'final_zarr'/section
    
    # Open and read the YAML file
    with open(yaml_file_path, 'r') as file:
        yaml_data = yaml.safe_load(file)
    
    # Read how many channels are in the selected cycle 
    channels = len(yaml_data[cycle]['channels'])
    
    # Create different cases for 2 or 4 channels 
  
    if channels == 2:

        # Create traces for the two squares
        overlap = yaml_data[cycle]['alpha'][0][0]
        offset = 1-abs(overlap)
        
        bg = yaml_data[cycle]['background'][0][0] * overlap / 100
        # proxy for % area to subtract from square 
        

        square1_color = 'rgb'+str(wavelength_to_rgb(yaml_data[cycle]['channels'][0]))
        square2_color = 'rgb'+str(wavelength_to_rgb(yaml_data[cycle]['channels'][1]))

        square1 = go.Scatter(
            x=[-0.5, 0.5, 0.5, -0.5, -0.5],
            y=[-0.5, -0.5, 0.5, 0.5, -0.5],
            mode='lines',
            text='Channel 1',
            name='Channel 1',
            line=dict(color=square1_color),
            fill='toself', 
            fillcolor=square1_color,
            opacity=0.8
        )
        square2 = go.Scatter(
            x=[-0.5 + offset, 0.5 + offset, 0.5 + offset, -0.5 + offset, -0.5 + offset],
            y=[-0.5, -0.5, 0.5, 0.5, -0.5],
            mode='lines',
            text='Channel 2',
            name='Channel 2',
            line=dict(color=square2_color),
            fill='toself', 
            fillcolor=square2_color,
            opacity=0.8
        )

        # Create the figure
        figure = {
            'data': [square1, square2],
            'layout': go.Layout(
                title=f'Channel overlap: {overlap:.2f}',
                xaxis=dict(range=[-8, 8]),
                yaxis=dict(range=[-2, 2]),
                showlegend=True
            )
        }
        return figure
    elif channels == 4:
        
        # Create traces for the first two channels
        overlap = yaml_data[cycle]['alpha'][0][0]
        offset = 1-abs(overlap)

        square1_color = 'rgb'+str(wavelength_to_rgb(yaml_data[cycle]['channels'][0]))
        square2_color = 'rgb'+str(wavelength_to_rgb(yaml_data[cycle]['channels'][1]))

        square1 = go.Scatter(
            x=[-0.5, 0.5, 0.5, -0.5, -0.5],
            y=[-0.5, -0.5, 0.5, 0.5, -0.5],
            mode='lines',
            text='Channel 1',
            name='Channel 1',
            line=dict(color=square1_color),
            fill='toself', 
            fillcolor=square1_color,
            opacity=0.8
        )
        square2 = go.Scatter(
            x=[-0.5 + offset, 0.5 + offset, 0.5 + offset, -0.5 + offset, -0.5 + offset],
            y=[-0.5, -0.5, 0.5, 0.5, -0.5],
            mode='lines',
            text='Channel 2',
            name='Channel 2',
            line=dict(color=square2_color),
            fill='toself', 
            fillcolor=square2_color,
            opacity=0.8
        )
        
        # Create traces for the third and fourth channels
        overlap2 = yaml_data[cycle]['alpha'][2][1]
        offset2 = 1-abs(overlap2)

        square3_color = 'rgb'+str(wavelength_to_rgb(yaml_data[cycle]['channels'][2]))
        square4_color = 'rgb'+str(wavelength_to_rgb(yaml_data[cycle]['channels'][3]))

        square3 = go.Scatter(
            x=[2, 3, 3, 2, 2],
            y=[-0.5, -0.5, 0.5, 0.5, -0.5],
            mode='lines',
            text='Channel 3',
            name='Channel 3',
            line=dict(color=square3_color),
            fill='toself', 
            fillcolor=square3_color,
            opacity=0.8
        )
        square4 = go.Scatter(
            x=[2+offset2, 3+offset2, 3+offset2, 2+offset2, 2+offset2],
            y=[-0.5, -0.5, 0.5, 0.5, -0.5],
            mode='lines',
            text='Channel 4',
            name='Channel 4',
            line=dict(color=square4_color),
            fill='toself', 
            fillcolor=square4_color,
            opacity=0.8
        )


        # Create the figure
        figure = {
            'data': [square1, square2, square3, square4],
            'layout': go.Layout(
                title=f'Channel 1 & 2 overlap: {overlap:.2f}; Channel 3 & 4 overlap: {overlap2:.2f}',
                xaxis=dict(range=[-8, 8]),
                yaxis=dict(range=[-2, 2]),
                showlegend=True
            )
        }

    return figure

      