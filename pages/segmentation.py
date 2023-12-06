import dash
from dash import dcc, html, callback
from dash.dependencies import Input, Output
import csv
import matplotlib.pyplot as plt
import os
from pathlib import Path
from dash.exceptions import PreventUpdate
import imageio.v2 as imageio
import numpy as np
import plotly.graph_objects as go
from dash_canvas import DashCanvas
from dash import dash_table
import pandas as pd

dash.register_page(__name__,name='Cell Type Segmentation')
       
layout = html.Div([
    html.H3("Tissue Section Cell Types Segmentation"),
    dcc.Dropdown(id='section-seg',placeholder='Tissue Section'),
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

#section dropdown
@callback(Output('section-seg','options'),
          Input('exp_dir','data'))
def update_section_dropdown(exp_dir):
    return [f.name[0:-4] for f in os.scandir(Path(exp_dir)/'tables') if f.name[-4:] == '.csv']


@callback(
    Output("canvas-section","image_content"),
    Input("section-seg","value"),
)
def section_image(section):
    if section is None:
        raise PreventUpdate
        
    return f'/static/{section}_seg_tab20c.jpeg'

# Callback to update the DataTable based on the dropdown selection
@callback(Output('table', 'data'),
          Output('table', 'columns'),
          Input('section-seg', 'value'),
          Input('exp_dir','data'),
          prevent_initial_call = True,
)
def update_table(section,exp_dir):
    if section is None:
        raise PreventUpdate
        
    df = pd.read_csv(f'{exp_dir}/dashboard/{section}_cell_type_key.csv')                  
    return df.to_dict('records'), [{"name": i, "id": i} for i in df.columns]
'''    

@callback(
    Output("section-img", "figure"),
    Input("section-seg","value"),
    Input('exp_dir', 'data'),
    prevent_initial_call = True
)
def display_image(section,exp_dir):
    if section == None:
        raise PreventUpdate
    
    # find section mask tiff file
    tiff_path = f'{exp_dir}/masks/{section}.tiff'
    matrix = read_tiff_as_matrix(tiff_path)
    
    # read csv file with cell types
    csv_path = f'{exp_dir}/tables/{section}.csv'
    #cell_type_id, cell_types = read_csv_for_cell_types(csv_path)
    cell_type_id = {}
    cell_types={}
    color = 1
    with open(csv_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)

        # Skip the header row if it exists
        next(csv_reader, None)

        # Iterate through the rows in the CSV file
        for row in csv_reader:
            category = row[1]  # Assuming the category is in the second column (index 1)
            index = int(row[0])

            if category in cell_type_id:
                cell_type_id[category].append(index)
            else:
                cell_type_id[category]=[index]
            if category not in cell_types:
                cell_types[category] = color 
                color +=1
        
    for category in cell_types:
        condition = cell_type_id[category]
        new_value = cell_types[category]
        matrix = replace_values_condition_np_where(matrix, condition, new_value)
    
    threshold = max(cell_types.values(), key=lambda x: x)
    mask = matrix > threshold

    # Replace values based on the mask
    matrix[mask] = 0
    
    # Define the file path for the TIFF file
    file_path = f'{exp_dir}/dashboard/output_image.tiff'

    # Save the matrix as a TIFF file
    imageio.imwrite(file_path, matrix)
    tiff_image = imageio.imread(file_path)

    # Create a figure using plotly.express.imshow
    fig = px.imshow(tiff_image, color_continuous_scale='gray')

    # Update layout
    fig.update_layout(coloraxis_showscale=False)

    #normed_matrix = new_matrix / np.max(new_matrix)
    
    #fig = px.imshow(normed_matrix, color_continuous_scale='Viridis')

    # Update layout
    #fig.update_layout(coloraxis_showscale=False)
    
    # resize matrix 
    rows, cols = matrix.shape
    quarter_matrix = matrix[:rows//5, :cols//5]
    print(quarter_matrix)
    print(np.max(quarter_matrix))
    
    # Create a colormap (you can choose a different one based on your preference)
    #cmap = plt.get_cmap('viridis')

    # Normalize the matrix values to the range [0, 1] for proper color mapping
    normed_matrix = quarter_matrix / np.max(quarter_matrix)
    
    fig = px.imshow(normed_matrix, color_continuous_scale='Viridis')

    # Update layout
    #fig.update_layout(coloraxis_showscale=False)
   
    # Apply the colormap to the normalized matrix
    #colored_matrix = cmap(normed_matrix)
    
    #fig = px.imshow(colored_matrix)

    
    return fig
'''


'''
def read_tiff_as_matrix(file_path):
    # Read the TIFF file using imageio
    image = imageio.imread(file_path)

    # Convert the image to a NumPy array (matrix)
    image_matrix = np.array(image)

    return image_matrix

def read_csv_for_cell_types(csv):
    cell_type_id = {}
    cell_types={}
    color = 1
    with open(csv, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)

        # Skip the header row if it exists
        next(csv_reader, None)

        # Iterate through the rows in the CSV file
        for row in csv_reader:
            category = row[1]  # Assuming the category is in the second column (index 1)
            index = int(row[0])

            if category in cell_type_id:
                cell_type_id[category].append(index)
            else:
                cell_type_id[category]=[index]
            if category not in cell_types:
                cell_types[category] = color 
                color +=1
    return cell_type_id, cell_types

def replace_values_condition_np_where(matrix, condition_values, new_value):
    # Use np.where to create a new matrix with replaced values
    new_matrix = np.where(np.isin(matrix, condition_values), new_value, matrix)

    return new_matrix

def resize_matrix(matrix, new_shape):
    if not matrix.size or new_shape[0] <= 0 or new_shape[1] <= 0:
        raise ValueError("Invalid input")

    # Calculate the block size for each dimension
    block_size = (matrix.shape[0] // new_shape[0], matrix.shape[1] // new_shape[1])

    # Reshape the matrix into non-overlapping blocks
    blocks = matrix[:new_shape[0] * block_size[0], :new_shape[1] * block_size[1]].reshape(
        new_shape[0], block_size[0], new_shape[1], block_size[1]
    )

    # Calculate the average value for each block
    resized_matrix = blocks.mean(axis=(1, 3))

    return resized_matrix

'''
          

    