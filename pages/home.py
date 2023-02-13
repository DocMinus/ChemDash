import dash
from dash import html
from . import navigation

dash.register_page(__name__, path="/")

layout = html.Div(
    children=[
        navigation.navbar,
        html.Br(),
        html.H3(children="Welcome to the Smart Chemistry Dashboard!"),
        html.Label("Please use above menus to navigate."),
        html.Br(),
        html.Label("For details, check out the README or have a look into the code."),
        html.Br(),
        html.Label("Written early January 2023 by DocMinus"),
    ]
)
