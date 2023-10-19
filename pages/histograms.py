import dash
from dash import dcc, ctx, html, callback, Input, Output, State 
import plotly.express as px
import numpy as np
import plotly.express as px
import cv2
import os
from pathlib import Path
from dash.exceptions import PreventUpdate

dash.register_page(__name__)

def calculate_histogram(image_path):
    image = cv2.imread(image_path)
    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    hist = cv2.calcHist([gray_image], [0], None, [256], [0, 256])
    return hist.ravel()


layout = html.Div([
            html.H3('Preview images with histograms'),
            "Select the image you want to display",        
            dcc.Dropdown(
                id="section-histogram-dropdown",
                placeholder="Tissue Section",
            ),
            dcc.Dropdown(
                id="marker-histogram-dropdown",
                placeholder="Marker",
            ),
            html.Button('Previous Image', id='prev-button',n_clicks=0),
            html.Button('Next Image', id='next-button',n_clicks=0),
            dcc.Graph(id="histogram"),
        ])


@callback(Output("section-histogram-dropdown","options"), 
          Input('exp_dir', 'data'))
def set_sections(exp_dir):
    return [f.name for f in os.scandir(Path(exp_dir)/'dashboard'/'preview') if f.is_dir()]


@callback(Output("marker-histogram-dropdown","options"),
          Output("marker-histogram-dropdown","value"),
          Input('exp_dir', 'data'),
          Input('section-histogram-dropdown', 'value'), 
          Input('marker-histogram-dropdown', 'value'),
         )
def set_markers(exp_dir, section, m):

    if section is not None:
        markers = [f.name[0:-4] for f in os.scandir(Path(exp_dir)/'dashboard'/'preview'/section) if f.name[-4:] == '.jpg']

        if m in markers:
            ind = markers.index(m)
        else:
            ind = 0

        if len(markers) == 0:
            m = ['No Markers']
        else:
            m = markers[ind]
      
        return markers, m
    else:
        return [], None

# @callback(Output("dropdown","options"), 
#           Input('exp_dir', 'data'))
# def set_sections(exp_dir):
#     return [f for f in os.listdir(Path(exp_dir)/'dashboard') if '.jpg' in f]

# @callback(
#     Output('histogram', 'figure'),
#     Input('dropdown', 'value'),
#     Input('exp_dir','data'),
#     prevent_initial_call = True
    
# )
# def update_histogram(image,exp_dir):
#     image_path = exp_dir+'/dashboard/'+image
   
#     hist = calculate_histogram(image_path)
    
#     return {
#         'data': [
#             {'x': list(range(256)), 'y': hist, 'type': 'bar', 'name': 'Intensity'},
#         ],
#         'layout': {
#             'title': 'Histogram of Intensity Values',
#             'xaxis': {'title': 'Intensity Value'},
#             'yaxis': {'title': 'Frequency'}
#         }
#     }

# Update histogram based on image selected
@callback(
    Output("histogram", "figure"), 
    Output("marker-histogram-dropdown","value", allow_duplicate = True),
    Input("section-histogram-dropdown", "value"),
    Input("marker-histogram-dropdown", "value"),
    Input("marker-histogram-dropdown","options"),
    Input("next-button","n_clicks"),
    Input("prev-button","n_clicks"),
    Input('exp_dir', 'data'),
    prevent_initial_call = True
)
def display_image(section,marker,marker_,next,prev, exp_dir):
    
    if marker is None:
        raise PreventUpdate

    if "marker-histogram-dropdown" == ctx.triggered_id and marker is not None:
        m = marker
    elif "section-histogram-dropdown" == ctx.triggered_id and marker is not None:
        m = marker
    
    # !! add cases for pressing next/prev when no image is selected    
    elif "next-button" == ctx.triggered_id:
        index = marker_.index(marker)
        if index < len(marker_)-1:
            m = marker_[index + 1]
        else:
            m = marker_[0]
    elif "prev-button" == ctx.triggered_id:
        index = marker_.index(marker)
        if index > 0:
            m = marker_[index - 1]
        else:
            m = marker_[-1]
    else:
        m = marker
        
    image_path = f'{exp_dir}/dashboard/{section}_{m}.jpg'
   
    hist = calculate_histogram(image_path)
    
    return {
        'data': [
            {'x': list(range(256)), 'y': hist, 'type': 'bar', 'name': 'Intensity'},
        ],
        'layout': {
            'title': 'Histogram of Intensity Values',
            'xaxis': {'title': 'Intensity Value'},
            'yaxis': {'title': 'Frequency'}
        }
    }, m




      