# Update small refactor
import dash
from dash import Input, Output, State, ctx, dcc, html

from src.calculations.molecule import Molecule
from src.calculations.oxybalance import n_count, oxy_balance
from src.ids import *
from src.mol_img import mol_image

from . import navigation

dash.register_page(__name__, path=MENU1SUB2_HREF)

layout = html.Div(
    [
        navigation.navbar,
        html.H1("Example Property Calculations via smiles-input"),
        html.Hr(),
        html.Div(
            [
                html.H3("Enter compound"),
                html.Label("Smiles: "),
                dcc.Input(id="calc_smiles_", type="text", value="", debounce=True),
                html.Button("Submit", id="calc_button_"),
            ]
        ),
        html.Br(),
        html.Span(id="o_balance_"),
        html.Div(id="num_of_n_"),
        html.Div(id="sum_formula_"),
        html.Div(id="mol_weight_"),
        html.Div(id="smiles-output_"),
        html.Div(id="calc_image_"),
        html.Br(),
        html.Div(
            [
                html.Button("Clear Form", id="reset", n_clicks=0),
            ]
        ),
    ]
)


@dash.callback(
    [
        Output("o_balance_", "children"),
        Output("num_of_n_", "children"),
        Output("sum_formula_", "children"),
        Output("mol_weight_", "children"),
        Output("smiles-output_", "children"),
        Output("calc_image_", "children"),
    ],
    [
        Input("calc_button_", "n_clicks"),
        State("calc_smiles_", "value"),
        Input("reset", "n_clicks"),
    ],
)
def calc_and_display(calc_button_, calc_smiles_, reset):
    button_triggered = ctx.triggered_id

    num_of_n = 0
    o_balance = 0
    sum_formula = "C0H0"
    mol_weight = 0
    calc_image = mol_image(calc_smiles_)
    alt_text = "No/invalid smiles. "

    if len(calc_image):
        o_balance = oxy_balance(calc_smiles_)
        num_of_n = n_count(calc_smiles_)
        _mol = Molecule(calc_smiles_)
        sum_formula = _mol.sumformula()
        mol_weight = _mol.molwt()

    _image = html.Img(
        src="data:image/svg+xml;base64,{}".format(calc_image),
        alt=alt_text,
    )

    if button_triggered == "reset":
        o_balance = 0
        sum_formula = "C0H0"
        mol_weight = 0
        calc_image = ""
        _image = html.Img(src="", alt=alt_text)

    return (
        f"OxyBalance = {o_balance} !Safe values > - 200!",
        f"Nitrogen Count = {num_of_n}",
        f"SumFormula = {sum_formula}",
        f"Mol.Weight = {round(mol_weight,2)}",
        f"Smiles = {calc_smiles_}",
        _image,
    )


@dash.callback(
    Output("calc_smiles_", "value"),
    [Input("reset", "n_clicks")],
)
def clear(reset):
    return ""
