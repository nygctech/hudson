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

# Create dictionary where keys are sections and values are a dictionary of cell types & their fractional counts
def cell_types(exp_dir):
    table_path = f'{exp_dir}/table'
    
    sections = [f.name[0:-4] for f in os.scandir(Path(exp_dir)/'tables') if f.name[-4:] == '.csv']
    
    cell_type_counts = {}
    cell_type_counts_frac = {}
    
    for section in sections:
        # Create a dictionary to count category appearances
        category_count_dict = {}

        # Open the CSV file and read its content
        with open(f'{exp_dir}/tables/{section}.csv', 'r') as csv_file:
            csv_reader = csv.reader(csv_file)

            # Skip the header row if it exists
            next(csv_reader, None)

            # Iterate through the rows in the CSV file
            for row in csv_reader:
                category = row[1]  # Assuming the category is in the second column (index 1)

                if category in category_count_dict:
                    category_count_dict[category] += 1
                else:
                    category_count_dict[category] = 1                
        
        cell_type_counts[section] = category_count_dict
        
        #create dict with fractional values
        sum_of_cells = sum(category_count_dict.values())
        for category in category_count_dict.keys():
            category_count_dict[category] = category_count_dict[category]/sum_of_cells
        cell_type_counts_frac[section] = category_count_dict    
        
    
    return cell_type_counts, cell_type_counts_frac

def unique_types(cell_type_counts):
    return list({cell_type for counts in cell_type_counts.values() for cell_type in counts.keys()})
    
    
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

#stacked bar plot
@callback(
    Output('stacked-bar-plot', 'figure'),
    Input('exp_dir', 'data'),
    [Input('stacked-bar-plot', 'relayoutData')]
)
def update_stacked_bar_plot(exp_dir,relayoutData):
    data = []
    _,cell_type_counts = cell_types(exp_dir)
    unique_cell_types = unique_types(cell_type_counts)

    for cell_type in unique_cell_types:
        counts = [cell_type_counts[tissue_section].get(cell_type, 0) for tissue_section in cell_type_counts.keys()]
        data.append(go.Bar(x=list(cell_type_counts.keys()), y=counts, name=cell_type, width=0.2))

    layout = go.Layout(
        barmode='stack',
        xaxis=dict(title='Tissue Section'),
        yaxis=dict(title='Cell Type Count'),
        title='Cell Types in Tissue Sections',
        width=0.3 * 1920,
    )

    return {'data': data, 'layout': layout}

#section dropdown
@callback(Output('sections','options'),
          Input('exp_dir','data'))
def update_section_dropdown(exp_dir):
    return [f.name[0:-4] for f in os.scandir(Path(exp_dir)/'tables') if f.name[-4:] == '.csv']

#pie chart
# figure out a way to display absolute value of each cell type?
@callback(Output('section-pie-chart','figure'),
          Input('exp_dir','data'),
          Input('sections','value'))
def update_pie_chart(exp_dir,section):
    if section is None:
        raise PreventUpdate
        
    cell_type_counts,_ = cell_types(exp_dir)
    data = cell_type_counts[section]
    fig = px.pie(values=list(data.values()), names=list(data.keys()), title=f'{section} Cell Types')
    fig.update_layout(width=0.3 * 1920)  
    
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

    