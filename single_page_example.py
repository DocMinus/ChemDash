"""
Chemical Dashboard for reactions
Version 0.1.1 (Jan 03, 14:15:00 2023)
Update: Jan 20, 2023.
Oxy balance included

@author: Alexander Minidis (DocMinus)

license: MIT
Copyright (c) 2023 DocMinus
"""


import dash_bootstrap_components as dbc
from dash import Dash, Input, Output, State, ctx, dcc, html
from rdkit import RDLogger

from src.calculations.oxybalance import n_count, oxy_balance
from src.ids import *
from src.mol_img import mol_image

RDLogger.logger().setLevel(RDLogger.CRITICAL)

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div(
    [
        html.H1("Safety Calculations"),
        html.Hr(),
        html.Div(
            [
                html.H3("Enter compound as smiles:"),
                html.Label(
                    "here TNT as test example: C1(C=C([N+](=O)[O-])C(C)=C([N+]([O-])=O)C=1)[N+]([O-])=O"
                ),
                html.Br(),
                dcc.Input(
                    id="input",
                    type="text",
                    value="C1(C=C([N+](=O)[O-])C(C)=C([N+]([O-])=O)C=1)[N+]([O-])=O",
                    debounce=True,
                ),
                html.Button("Submit", id="button0"),
            ]
        ),
        html.Br(),
        html.Div(id="image"),
        html.Span(id="o_balance"),
        html.Div(id="num_of_n"),
        html.Br(),
        html.Div(
            [
                html.Button("Clear Form", id="reset", n_clicks=0),
            ]
        ),
    ]
)


@app.callback(
    [
        Output("o_balance", "children"),
        Output("num_of_n", "children"),
        Output("image", "children"),
    ],
    [
        Input("button0", "n_clicks"),
        State("input", "value"),
        Input("reset", "n_clicks"),
    ],
)
def update_image_0(button0, input, reset):
    button_triggered = ctx.triggered_id

    num_of_n = 0
    o_balance = 0
    image = mol_image(input)
    alt_text = "No or invalid smiles. "

    if len(image):
        o_balance = oxy_balance(input)
        num_of_n = n_count(input)

    _i0 = html.Img(
        src="data:image/svg+xml;base64,{}".format(image),
        alt=alt_text,
    )

    if button_triggered == "reset":
        o_balance = 0
        image = ""
        _i0 = html.Img(src="", alt=alt_text)

    return (
        f"Oxygene Balance = {o_balance} !Safe values > - 200!",
        f"Nitrogen Count = {num_of_n}",
        _i0,
    )


@app.callback(
    Output("input", "value"),
    [Input("reset", "n_clicks")],
)
def clear(reset):
    return ""


if __name__ == "__main__":
    debug = True
    app.run_server(debug=debug, port=8069)
