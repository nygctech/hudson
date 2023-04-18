# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import pandas as pd
import os
import dash
from dash import dcc, ctx
from dash import html 
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

app = dash.Dash(__name__)
server = app.server

# locate folder of interest
### !!! folder naming? can it be named something other than "assets"
folders = [f for f in os.listdir('.') if os.path.isdir(f)]

#image_paths=["assets/" + x for x in os.listdir("assets") if ".jpg" in x or ".jpeg" in x]

### !!! use the canvas component to display 
### !!! look into css/html, find templates? 
### !!! store every step in a tab
### !!! for the future: test data dir - write automated tests to test components
### !!! ^ pytest

app.layout = html.Div([
    html.H1("Image Thumbnails - 17 Apr"), 
    html.Div([
        "Select the folder with your images:",
        dcc.Dropdown(
            id="folder-dropdown",
            options=[{"label": x, "value": x} for x in folders],
            #value=folders[0],
            placeholder="Select the folder with your images"
        ),
    ]), 
    #list files in selected folder in a bulleted list
    #html.Div(id="folder-files"),
    html.Div([
        "Select the image you want to display",
        dcc.Dropdown(
            id="image-dropdown",
            #options=[{"label": x, "value": x} for x in image_paths],
            #value=image_paths[0],
            placeholder="Select the image you want to display",
        ),
    ]),
    html.Button('Previous Image', id='prev',n_clicks=0),
    html.Button('Next Image', id='next',n_clicks=0),
    html.Div(
        #html.Img(src=image_paths[0]),
        id='container-image'
    ),              
])

# Update image-dropdown depending on folder selected 
@app.callback(Output("image-dropdown", "options"), Input("folder-dropdown", "value"))
def update_image_dropdown(folder_name):
    images = [x for x in os.listdir(folder_name) if ".jpg" in x or ".jpeg" in x]
    #list all files in bulleted list
    #file_list = html.Ul([html.Li(file) for file in file_names])
    return images

# Display image selected in image-dropdown, utilize buttons to click through
@app.callback(
    Output("container-image", "children"), 
    Output("image-dropdown","value"),
    Input("folder-dropdown","value"),
    Input("image-dropdown", "value"),
    Input("image-dropdown","options"),
    Input("next","n_clicks"),
    Input("prev","n_clicks"),
)
def display_image(folder,image,image_names,next,prev):
    if image is None:
        raise PreventUpdate

    if "image-dropdown" == ctx.triggered_id:
        return html.Img(src=folder+"/"+image),image
    #add cases for pressing next/prev when no image is selected    
    elif "next" == ctx.triggered_id:
        index = image_names.index(image)
        if index < len(image_names)-1:
            index += 1
            return html.Img(src=folder+"/"+image_names[index]),image_names[index]
        else:
            return html.Img(src=folder+"/"+image_names[0]),image_names[0]
    elif "prev" == ctx.triggered_id:
        index = image_names.index(image)
        if index > 0:
            index -= 1
            return html.Img(src=folder+"/"+image_names[index]),image_names[index]
        else:
            return html.Img(src=folder+"/"+image_names[len(image_names)-1]),image_names[len(image_names)-1]

if __name__ == "__main__":
    app.run_server(debug=True)

'''
OLD CODE FOR BACKUP 
# Click through images with next/prev buttons
@app.callback(
    Output("container-image","children"),

    Input("folder-dropdown","value"),

)
def updateIndex(next,prev,folder):
    
    #index = n_clicks%len(image_paths) 
    """
    When the next or previous buttons are clicked, update the image index
    ensuring it stays within 0 <= index < len(contents)
    """

        if trig == "prevImage.n_clicks" and index > 0:
        data['index'] = index - 1
    elif trig == "nextImage.n_clicks" and index < len(contents) - 1:
        data['index'] = index + 1
    else:
        # index should not change, it would cause undefined index
        raise PreventUpdate
    return data

    # update data if necessary and return it

    if "next" == ctx.triggered_id and index < len(file_names)-1:
        index = index + 1
    elif "next" == ctx.triggered_id and index >= len(file_names) - 1:

    return html.Img(src=image_paths[index])
'''

'''

def display_all_images(folder_name):
    file_names = os.listdir(folder_name)

    image_list = html.Ul([html.Li(html.Img(src = folder_name + '/'+ file for file in file_names))])

    return image_list

'''

