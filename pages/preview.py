import dash
from dash import dcc, ctx, html, callback, Input, Output, State 
from dash_canvas import DashCanvas
from dash.exceptions import PreventUpdate
# from dash.dependencies import callback, Input, Output, State
import os
from pathlib import Path

dash.register_page(__name__)

# point Dash to image directory of interest here
# image_directory = '/nethome/kpandit/hudson/images'
# list_of_images = [f for f in os.listdir(image_directory) if '.png' in f]

# point Dash to image directory of interest here
# image_directory = f'{args.directory}/dashboard/preview/m387ntga1'
# list_of_images = [f for f in os.listdir(image_directory) if '.jpg' in f]
# static_image_route = '/static/'

#image_directory = '/commons/groups/nygcfaculty/PySeq/20210505_mouse_genotype_3/tiffs'
#list_of_images = [f for f in os.listdir(image_directory) if '.tiff' in f]
static_image_route = '/static/'


layout = html.Div([
    html.H1("Image Thumbnails - 2 May"), 
    #html.Div(["Folder being served: ", image_directory]),
    html.Div([
        "Select the image you want to display",
        dcc.Dropdown(
            id="section-dropdown",
            placeholder="Tissue Section",
        ),
        dcc.Dropdown(
            id="marker-dropdown",
            placeholder="Marker",
        ),
    ]),
    html.Button('Previous Image', id='prev',n_clicks=0),
    html.Button('Next Image', id='next',n_clicks=0),
    html.Div(
        #html.Img(src=image_paths[0]),
        id='container-image'
    ),
    html.Div([
        DashCanvas(id='canvas-image',
            tool='line',
            lineWidth=5,
            lineColor='red',
  #          image_content=static_image_route+list_of_images[0]
            ),
    ]),           
])
'''
    # folder dropdown 
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
'''

'''
# Update image-dropdown depending on folder selected 
@app.callback(Output("image-dropdown", "options"), Input("folder-dropdown", "value"))
def update_image_dropdown(folder_name):
    images = [x for x in os.listdir(folder_name) if ".jpg" in x or ".jpeg" in x]
    #list all files in bulleted list
    #file_list = html.Ul([html.Li(file) for file in file_names])
    return images
'''
# # Serve static images
# @layout.server.route('{}<image_path>.png'.format(static_image_route))
# def serve_image(image_path):
#     image_name = '{}.png'.format(image_path)
#     if image_name not in list_of_images:
#         raise Exception('"{}" is excluded from the allowed static files'.format(image_path))
#     return flask.send_from_directory(image_directory, image_name)


@callback(Output("section-dropdown","options"), 
          Input('exp_dir', 'data'))
def set_sections(exp_dir):
    return [f.name for f in os.scandir(Path(exp_dir)/'dashboard'/'preview') if f.is_dir()]


@callback(Output("marker-dropdown","options"),
          Output("marker-dropdown","value"),
          Input('exp_dir', 'data'),
          Input('section-dropdown', 'value'), 
          Input('marker-dropdown', 'value'),
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

# Update canvas image based on image selected in image-dropdown
@callback(
    Output("canvas-image", "image_content"), 
    Output("marker-dropdown","value", allow_duplicate = True),
    #Input("folder-dropdown","value"),
    Input("section-dropdown", "value"),
    Input("marker-dropdown", "value"),
    Input("marker-dropdown","options"),
    Input("next","n_clicks"),
    Input("prev","n_clicks"),
    Input('exp_dir', 'data'),
    prevent_initial_call = True
)
def display_image(section,marker,marker_,next,prev, exp_dir):
    
    if marker is None:
        raise PreventUpdate

    if "marker-dropdown" == ctx.triggered_id and marker is not None:
        m = marker
    elif "section-dropdown" == ctx.triggered_id and marker is not None:
        m = marker
    
    # !! add cases for pressing next/prev when no image is selected    
    elif "next" == ctx.triggered_id:
        index = marker_.index(marker)
        if index < len(marker_)-1:
            m = marker_[index + 1]
        else:
            m = marker_[0]
    elif "prev" == ctx.triggered_id:
        index = marker_.index(marker)
        if index > 0:
            m = marker_[index - 1]
        else:
            m = marker_[-1]
    else:
        m = marker

    return f'/static/{section}_{m}.jpg', m