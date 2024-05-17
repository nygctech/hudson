import dash
from dash import html

dash.register_page(__name__, name="Home",order=1)

layout = html.Div([
    html.H2("This is a Visual Dashboard for the Hudson."),
    html.Div([
        "The overarching purpose of the dashboard is to quickly summarize and assess the quality of a Hudson experiment.", 
        html.Br(),
        "The dashboard has pages that visualize different aspects of the experiment in a user friendly manner.",
        html.Br(),
        "This dashboard is run on Plotly Dash -- a low-code framework for rapidly building data apps in Python.",
    ]),
    html.Img(src="https://i.redd.it/hldaw5pz8xcb1.png"),
])
