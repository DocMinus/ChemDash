import dash
from dash import html, Input, Output, dcc
import dash_bio as dashbio
import plotly.graph_objects as go
import plotly.subplots as sp

import pandas as pd

from . import navigation

from src.calculations.oxybalance import oxy_balance
from src.calculations.mol_properties import molecule
from src.mol_img import mol_image
from src.ids import *

dash.register_page(__name__, path=MENU1SUB1_HREF)

layout = html.Div(
    [
        navigation.navbar,
        html.H2("Some Property Calculations"),
        html.Hr(),
        html.Div(
            [
                dashbio.Jsme(
                    width="25%",
                    height="36vh",
                    id="jsme",
                    smiles="",
                ),
            ]
        ),
        html.Br(),
        html.Span(id="o_balance"),
        html.Div(id="sum_formula"),
        html.Div(id="smiles-output"),
        html.Div(id="calc_image"),
        html.Div(
            dcc.Graph(id="properties-table"),
            style={"width": "50%", "margin": "1", "display": "inline-block"},
        ),
    ]
)


@dash.callback(
    [
        Output("o_balance", "children"),
        Output("sum_formula", "children"),
        Output("smiles-output", "children"),
        Output("calc_image", "children"),
        Output("properties-table", "figure"),
    ],
    [
        Input("jsme", "eventSmiles"),
    ],
)
def calc_and_display(jsme_smiles):
    if jsme_smiles == None:
        jsme_smiles = ""  # else error in image module

    o_balance = 0
    sum_formula = "C0H0"
    calc_image = mol_image(jsme_smiles)
    alt_text = "No/invalid smiles. "
    descriptors = pd.DataFrame()
    fig = sp.make_subplots(rows=1, cols=1, specs=[[{"type": "table"}]])
    show_first_x_props = 5
    emtpy_table = pd.DataFrame(
        dict.fromkeys(range(show_first_x_props), None), index=[0]
    )

    if len(calc_image):
        o_balance = oxy_balance(jsme_smiles)
        _mol = molecule(jsme_smiles)
        sum_formula = _mol.sumformula()
        descriptors = _mol.moldescriptors()
        # only selected few properties shown; again, only proof of concept.
        # slicing from 1:x instead of 0:x removes the ID
        descriptors = descriptors.iloc[:1, 1:show_first_x_props]
        # returning a table within a callback, leads to some additional lines of code
        # instead of desired simple fig.show().
        fig.add_trace(
            go.Table(
                header=dict(values=descriptors.columns),
                cells=dict(values=[descriptors[col] for col in descriptors.columns]),
            )
        )
    else:
        fig.add_trace(
            go.Table(
                header=dict(values=emtpy_table.columns),
                cells=dict(values=[emtpy_table[col] for col in emtpy_table.columns]),
            )
        )

    # can't get displayModeBar False to work, probably doing it wrong.
    #fig.update_layout(displayModeBar=False, width=800, height=400)
    fig.update_layout(width=800, height=400)

    _image = html.Img(
        src="data:image/svg+xml;base64,{}".format(calc_image),
        alt=alt_text,
    )

    return (
        f"Oxygene balance = {o_balance} !Safe values > - 200!",
        f"SumFormula = {sum_formula}",
        f"Smiles = {jsme_smiles}",
        _image,
        fig,
    )
