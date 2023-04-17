# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import pandas as pd
import os
import dash
from dash import dcc, ctx
from dash import html 
from dash.dependencies import Input, Output, State

app = dash.Dash(__name__)
server = app.server

# locate folder of interest
### !!! folder naming? can it be named something other than "assets"
folders = ["assets"]

# get paths of files in folder
### !!! only read .jpg or .jpeg
file_names = os.listdir("assets")
image_paths = []
for image in file_names:
    image_paths.append("assets/"+ image)

### !!! add dropdown menu to select image to view
### !!! Prev Image button 
### !!! use the canvas component to display 
### !!! look into css/html, find templates? 

### !!! store every step in a tab
### !!! for the future: test data dir - write automated tests to test components
### !!! ^ pytest

controls = [
    dcc.Dropdown(
        id="dropdown",
        options=[{"label": x, "value": x} for x in folders],
        value=folders[0],
    )
]

app.layout = html.Div(
    [html.H1("Image Thumbnails - Test"), html.Div(controls), html.Div(id="folder-files"),
    html.Button('Next Image', id='next',n_clicks=0),
    html.Div(html.Img(src=image_paths[0]),id='container-button-click')              
    ]
)


@app.callback(Output("folder-files", "children"), Input("dropdown", "value"))
def list_all_files(folder_name):
    # This is relative, but you should be
    # able to provide the absolute path too
    file_names = os.listdir(folder_name)

    file_list = html.Ul([html.Li(file) for file in file_names])

    return file_list


@app.callback(Output("container-button-click","children"),
              Input("next", "n_clicks"))
def updateIndex(n_clicks):
    index = n_clicks%len(file_names) 

    """
    When the next or previous buttons are clicked, update the image index
    ensuring it stays within 0 <= index < len(contents)
    """

    # update data if necessary and return it
    '''
    if "next" == ctx.triggered_id and index < len(file_names)-1:
        index = index + 1
    elif "next" == ctx.triggered_id and index >= len(file_names) - 1:

    '''
    
    return html.Img(src=image_paths[index])

'''

def display_all_images(folder_name):
    file_names = os.listdir(folder_name)

    image_list = html.Ul([html.Li(html.Img(src = folder_name + '/'+ file for file in file_names))])

    return image_list

'''

if __name__ == "__main__":
    app.run_server(debug=True)

