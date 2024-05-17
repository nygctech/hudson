import dash
from dash import dcc, html, callback
from dash.dependencies import Input, Output
import plotly.express as px
import csv
import plotly.graph_objs as go
import os
from pathlib import Path
from dash.exceptions import PreventUpdate
from dash_canvas import DashCanvas
from dash import dash_table
import pandas as pd

dash.register_page(__name__,name='Cell Type Analysis and Segmentation')

# Read the sections dicts created by snakemake and displays cell type analysis for the sections 
    
app = dash.Dash(__name__)

layout = html.Div([
    html.H3("Aggregate Tissue Section Cell Types Analysis"),
    dcc.Graph(id='stacked-bar-plot'),
    html.H3("Individual Section Cell Type Analysis and Segmentation"),
    dcc.Dropdown(id='sections',placeholder='Tissue Section'),
    dcc.Graph(id='section-pie-chart'),
    html.Div([
        DashCanvas(id='canvas-section',
            tool='line',
            lineWidth=5,
            lineColor='red',
            ),
    ]),
    'Key for cell types:',
    dash_table.DataTable(
        id='table',
        fill_width=False,
        style_table={'overflowX': 'auto'},
        style_cell={
            'height': 'auto',
            # all three widths are needed
            'minWidth': '200px', 'width': '200px', 'maxWidth': '200px',
            'whiteSpace': 'normal'},
    ),
])

# Make a stacked bar plot of all the sections
@callback(
    Output('stacked-bar-plot', 'figure'),
    Input('exp_dir', 'data'),
    [Input('stacked-bar-plot', 'relayoutData')]
)
def update_stacked_bar_plot(exp_dir,relayoutData):
    sections = [f.name[0:-18] for f in os.scandir(Path(exp_dir)/'tables') if f.name[-18:] == '_celltypecount.csv']
    section_dfs = {}

    # Combine all the section csv
    for section in sections:
        data = pd.read_csv(Path(exp_dir)/'tables'/f'{section}_celltypecount.csv')
        section_dfs[section] = data
    
    # Combine DataFrames
    combined_df = pd.concat(section_dfs.values(), keys=section_dfs.keys())
    
    # Reset index to remove hierarchical index created by concat
    combined_df.reset_index(level=0, inplace=True)
    combined_df.rename(columns={'level_0': 'Section'}, inplace=True)

    fig = px.bar(combined_df, x="Section", y="Fractional Count", color="Cell Type", title="Cell Types in Tissue Sections")

    return fig

#Section dropdown - select a section to generate a pie chart
@callback(Output('sections','options'),
          Input('exp_dir','data'))
def update_section_dropdown(exp_dir):
    return [f.name[0:-18] for f in os.scandir(Path(exp_dir)/'tables') if f.name[-18:] == '_celltypecount.csv']

#Generate a cell type pie chart for specified section 
@callback(Output('section-pie-chart','figure'),
          Input('exp_dir','data'),
          Input('sections','value'))
def update_pie_chart(exp_dir,section):
    if section is None:
        raise PreventUpdate
        
    data = pd.read_csv(Path(exp_dir)/'tables'/f'{section}_celltypecount.csv') 
    fig = px.pie(data, values='Count', names='Cell Type', title=f'{section} Cell Types')
    
    return fig

@callback(
    Output("canvas-section","image_content"),
    Input("sections","value"),
)
def section_image(section):
    if section is None:
        raise PreventUpdate
        
    return f'/static/{section}_seg_tab20c.jpeg'

# Callback to update the DataTable based on the dropdown selection
@callback(Output('table', 'data'),
          Output('table', 'columns'),
          Input('sections', 'value'),
          Input('exp_dir','data'),
          prevent_initial_call = True,
)
def update_table(section,exp_dir):
    if section is None:
        raise PreventUpdate
        
    df = pd.read_csv(f'{exp_dir}/dashboard/{section}_cell_type_key.csv')                  
    return df.to_dict('records'), [{"name": i, "id": i} for i in df.columns]
    
