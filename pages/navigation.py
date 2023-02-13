import dash
from dash import html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

from src.ids import *

navbar = dbc.Navbar(
    dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.Img(
                                src=dash.get_asset_url(MENU_LOGO),
                                height="50px",
                            ),
                            dbc.NavbarBrand("Chemistry Dashboard", className="ms-2"),
                        ],
                        width={"size": "auto"},
                    )
                ],
                align="center",
                className="ms-2",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Nav(
                                [
                                    dbc.NavItem(dbc.NavLink("Home", href="/")),
                                    dbc.NavItem(
                                        dbc.DropdownMenu(
                                            children=[
                                                dbc.DropdownMenuItem(
                                                    MENU1SUB1_NAME,
                                                    href=MENU1SUB1_HREF,
                                                ),
                                                dbc.DropdownMenuItem(
                                                    MENU1SUB2_NAME,
                                                    href=MENU1SUB2_HREF,
                                                ),
                                            ],
                                            nav=True,
                                            in_navbar=True,
                                            label=MENU1_NAME,
                                        )
                                    ),
                                    dbc.NavItem(
                                        dbc.DropdownMenu(
                                            children=[
                                                dbc.DropdownMenuItem(
                                                    MENU2SUB1_NAME,
                                                    href=MENU2SUB1_HREF,
                                                ),
                                            ],
                                            nav=True,
                                            in_navbar=True,
                                            label=MENU2_NAME,
                                        )
                                    ),
                                ],
                                navbar=True,
                            )
                        ],
                        width={"size": "auto"},
                    )
                ],
                align="center",
            ),
            dbc.Col(dbc.NavbarToggler(id="navbar-toggler", n_clicks=0)),
        ],
        fluid=True,
    ),
    color="primary",
    dark=True,
)


@dash.callback(
    Output("navbar-collapse", "is_open"),
    [Input("navbar-toggler", "n_clicks")],
    [State("navbar-collapse", "is_open")],
)
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open
