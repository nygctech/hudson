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

dash.register_page(__name__,name='Summary Report')
       

# Define the layout of the app
layout = html.Div([
    html.H2("Experiment Summary Report"),
    html.H3("Autofocus Evaluation"),
    html.H3("Imaging Time"),
    html.H3("Number of Cycles"),
    
])
