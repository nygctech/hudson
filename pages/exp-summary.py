import dash
from dash import dcc, html, callback
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from dash import dash_table
import pandas as pd
import cv2
import os
import plotly.graph_objects as go
import plotly.express as px
import yaml
import numpy as np

dash.register_page(__name__,name='Summary Report')

def compute_contrast(image):
    # Convert the image to grayscale
    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    # Compute the gradient magnitude
    sobelx = cv2.Sobel(gray_image, cv2.CV_64F, 1, 0, ksize=5)
    sobely = cv2.Sobel(gray_image, cv2.CV_64F, 0, 1, ksize=5)
    gradient_magnitude = np.sqrt(sobelx**2 + sobely**2)

    # Return the sum of gradient magnitudes as contrast
    return np.sum(gradient_magnitude)

def max_contrast_section(image):
# Divide the image into 9 equal sections
    sections = []
    height, width = image.shape[:2]
    section_height = height // 3
    section_width = width // 3

    for y in range(0, height, section_height):
        for x in range(0, width, section_width):
            section = image[y:y+section_height, x:x+section_width]
            sections.append(section)

    # Compute contrast for each section
    contrasts = [compute_contrast(section) for section in sections]

    # Find the index of the section with the highest contrast
    max_contrast_index = np.argmax(contrasts)
    
    # Display the section with the highest contrast
    fig = px.imshow(cv2.cvtColor(sections[max_contrast_index], cv2.COLOR_BGR2RGB))

    return fig
       

# Define the layout of the app
layout = html.Div([
    html.H2("Experiment Summary Report"),
    html.H3("Autofocus Evaluation"),
    html.H3("Imaging Time"),
    html.H3("Cycle Summary"),
    html.H4("Markers imaged in each cycle"),
    html.Div(id='cycle-markers'),
    html.Div(id = 'marker-thumbnails'),
    
])

@callback(Output("cycle-markers","children"),
          Input("exp_dir","data"))
def get_cycle_markers(exp_dir):
    yaml_file = exp_dir+"/mouse_genotype2.yaml"
          

    with open(yaml_file, 'r') as file:
        data = yaml.safe_load(file)
        markers_section = data.get("markers", {})
        output = []
        for key, value in markers_section.items():
            if isinstance(value, dict):
                cycle_number = key
                markers_list = ", ".join(list(value.values()))
                output.append(f"In cycle {cycle_number}, the markers imaged were: "+ markers_list)
                
    items_html = [html.Li(item) for item in output]
    return html.Ul(items_html)

@callback(Output("marker-thumbnails","children"),
          Input("exp_dir","data"))
def get_thumbnails(exp_dir):
    # List all files in the directory
    files = os.listdir(exp_dir+"/dashboard")

    # Filter JPEG files
    jpeg_files = [file for file in files if file.endswith('.jpg')]

    # Open each JPEG file with cv2
    graph_components = []
    i = 1
    for jpeg_file in jpeg_files:
        # Construct the full path to the image
        image_path = os.path.join(exp_dir+'/dashboard',jpeg_file)

        # Read the image using cv2
        image = cv2.imread(image_path)

        # Append the graph component to the list
        graph_components.append(
            dcc.Graph(
                id=f'graph{i}',
                figure=max_contrast_section(image),
                style={'display': 'inline-block','width': '70%', 'height': '300px'},
                config={'displayModeBar': False}
            )
        )
        i += 1

    
    return graph_components



