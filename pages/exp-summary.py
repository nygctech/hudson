import dash
from dash import dcc, html, callback
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from dash import dash_table
import pandas as pd
import re
import cv2
import os
from datetime import datetime
import plotly.graph_objects as go
import plotly.express as px
import yaml
import numpy as np

dash.register_page(__name__,order =2, name='Summary Report')

'''
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
 '''      

# Define the layout of the app
layout = html.Div([
    html.H2("Experiment Summary Report"),
    html.H3("Experiment Duration"),
    html.Div(id = "exp-duration",style={'whiteSpace': 'pre-line'}),
    html.H3("Imaging Time"),
    dash_table.DataTable(
        id='img-time',
        fill_width=False,
        style_table={'overflowX': 'auto'},
        style_cell={
            'height': 'auto',
            # all three widths are needed
            'minWidth': '300px', 'width': '300px', 'maxWidth': '300px',
            'whiteSpace': 'normal'},
    ),
    html.H3("Cycle Summary"),
    html.H4("Markers imaged in each cycle"),
    html.Div(id='cycle-markers'),
    html.H4("Thumbnails of markers"),
    html.Div(id = 'marker-thumbnails'),
    
])


## Get start time, end time, and total duration of exp
def extract_timestamp(line):
    # Define a regular expression pattern to match the timestamp
    timestamp_pattern = r"(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3})"
    
    # Use re.search to find the timestamp in the log entry
    match = re.search(timestamp_pattern, line)
    
    # Return the timestamp if found, otherwise return None
    return match.group(1) if match else None

# Get duration
def calculate_duration(start_timestamp, end_timestamp, timestamp_format='%Y-%m-%d %H:%M:%S,%f'):
    # Convert timestamps to datetime objects
    start_time = datetime.strptime(start_timestamp, timestamp_format)
    end_time = datetime.strptime(end_timestamp, timestamp_format)

    # Calculate the duration
    duration = end_time - start_time

    # Extract days, hours, minutes, and seconds
    days, seconds = duration.days, duration.seconds
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    # Return a formatted string
    duration_string = f"{days} days, {hours} hours, {minutes} minutes, {seconds} seconds"
    return duration_string

@callback(Output("exp-duration","children"),
          Input("exp_dir","data"))
def get_log_duration(exp_dir):
    log_file_path = f'{exp_dir}/logs/20211207_AlternatingSABER.log'

    # Get number of lines
    try:
        with open(log_file_path, 'r',encoding='iso-8859-1') as file:
            number_of_lines = len(file.readlines())
            #print(f"Number of lines in the log file: {number_of_lines}")
    except FileNotFoundError:
        print(f"The file '{log_file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    try:
        with open(log_file_path, 'r', encoding='iso-8859-1') as file:
            line_id = 1

            for line in file:
                line.strip()

                ## Get start time from line 1
                if line_id == 1:
                    start_time = extract_timestamp(line)
                    #print(f'Experiment start time: {start_time}')
                    break

            # Iterate through the lines in reverse order
            lines = file.readlines()
            for line in reversed(lines):
                end_time = extract_timestamp(line)

                if end_time:
                    #print(f"Experiment end time: {end_time}")
                    break  # Stop iterating once a line with a timestamp is found

    except FileNotFoundError:
        print(f"The file '{log_file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # Call the function to calculate the duration
    duration = calculate_duration(start_time, end_time)
    
    return f'Experiment start time: {start_time}\nExperiment end time: {end_time}\nExperiment duration: {duration}\n'


def parse_imgtime_log(log_file_path):
    with open(log_file_path, 'r', encoding='iso-8859-1') as file:
        lines = file.readlines()

    data = []
    for line in lines:
        match = re.search(r'(\d{4}-\d{2}-\d{2}) (\d{2}:\d{2}:\d{2}),.*?cycle(\d+)::m(\d+)::Imaging completed in (\d+) minutes', line)
        if match:
            cycle_number = match.group(3)
            m_number = match.group(4)
            imaging_time = int(match.group(5))
            data.append({'Cycle': cycle_number,'Imaging Time (minutes)': imaging_time})

    df = pd.DataFrame(data)
    return df

@callback(Output('img-time','data'),
          Output('img-time','columns'),
          Input('exp_dir','data'))
def img_time(exp_dir):
    # Function to parse log file and extract information
    log_file_path = f'{exp_dir}/logs/20211207_AlternatingSABER.log'
    if log_file_path is not None:
        df = parse_imgtime_log(log_file_path)
        avg_imaging_time = df.groupby('Cycle')['Imaging Time (minutes)'].mean().reset_index().round(2)
        
       # avg_imaging_time.drop(columns=['M Number'])
        overall_avg = df['Imaging Time (minutes)'].mean()
        
        # Row to append
        new_row = pd.DataFrame([{'Cycle': 'Average', 'Imaging Time (minutes)': overall_avg}])

        # Append the row
        avg_imaging_time = pd.concat([avg_imaging_time, new_row], ignore_index=True).round(2)

        return avg_imaging_time.to_dict('records'), [{"name": i, "id": i} for i in avg_imaging_time.columns]

    
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
    # Get list of image files in the folder
    image_files = [f.name for f in os.scandir(f'{exp_dir}/dashboard') if f.name[-13:] == 'thumbnail.jpg']
    #images = []
    image_files.sort()

    # Open each JPEG file with cv2
    graph_components = []
    i = 1
    for jpeg_file in image_files:
        # Construct the full path to the image
        image_path = os.path.join('/static',jpeg_file)

        # Append the graph component to the list
        graph_components.append(
            html.Div([
                html.Img(
                    id=f'img{i}',
                    src = image_path,
                    style={'width':'30%','height':'30%'},
                ),
                html.P(jpeg_file),
            ],style={'display': 'inline-block','width': '30%','height':'30%','whiteSpace': 'normal'}),
            
        )
        i += 1

    return graph_components

    

'''
    
'''
'''
    for img_file in image_files:
        # Create img element for each image
        img_path = f'static/{img_file}'
        images.append(html.Img(src=img_path, style={'width': '300px', 'height': '300px'}))

        # Create paragraph element for image caption
        #caption_element = html.P(img_file, style={'text-align': 'center'})
        
        # Append both img and caption elements to images list
        images.append(caption_element)
        
    return images

'''