import dash
from dash import html

dash.register_page(__name__, name="Home",order=1)

layout = html.Div([
    html.H1("This is a Visual Dashboard for the Hudson."),
    html.Img(src="https://i.redd.it/hldaw5pz8xcb1.png"),
])
