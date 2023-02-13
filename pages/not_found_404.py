from dash import html
import dash
import dash_bootstrap_components as dbc

dash.register_page(__name__)

# !!! the filename has to be 'not_found_404.py'

layout = html.Div(
    children=[
        html.H1("Oops, that page doesn't exist"),
        html.Hr(),
        dbc.Button("Go back to Home", size="lg", id="home_btn_404", href="/"),
    ],
    className="notfound404",
)
