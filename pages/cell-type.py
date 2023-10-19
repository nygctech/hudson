import dash
from dash import dcc, html, callback
from dash.dependencies import Input, Output
import plotly.express as px
import csv
import plotly.graph_objs as go
import os
from pathlib import Path

dash.register_page(__name__)

def cell_types(exp_dir):
    table_path = f'{exp_dir}/table'
    
    sections = [f.name[0:-4] for f in os.scandir(Path(exp_dir)/'tables') if f.name[-4:] == '.csv']
    
    cell_type_counts = {}
    
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
    
    return cell_type_counts

def unique_types(cell_type_counts):
    return list({cell_type for counts in cell_type_counts.values() for cell_type in counts.keys()})
    
    
app = dash.Dash(__name__)


layout = html.Div([
    html.H1("Stacked Bar Plot of Cell Types in Tissue Sections"),
    dcc.Graph(id='stacked-bar-plot'),
])

@callback(
    Output('stacked-bar-plot', 'figure'),
    Input('exp_dir', 'data'),
    [Input('stacked-bar-plot', 'relayoutData')]
)
def update_stacked_bar_plot(exp_dir,relayoutData):
    data = []
    cell_type_counts = cell_types(exp_dir)
    unique_cell_types = unique_types(cell_type_counts)

    for cell_type in unique_cell_types:
        counts = [cell_type_counts[tissue_section].get(cell_type, 0) for tissue_section in cell_type_counts.keys()]
        data.append(go.Bar(x=list(cell_type_counts.keys()), y=counts, name=cell_type))

    layout = go.Layout(
        barmode='stack',
        xaxis=dict(title='Tissue Section'),
        yaxis=dict(title='Cell Type Count'),
        title='Cell Types in Tissue Sections',
    )

    return {'data': data, 'layout': layout}

if __name__ == '__main__':
    app.run_server(debug=True)

    