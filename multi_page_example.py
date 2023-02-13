#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Chemical Dashboard for reactions
Version 0.3.7 (Jan 03, 14:15:00 2023)
Update: Feb 10, 2023.

@author: Alexander Minidis (DocMinus)

license: MIT
Copyright (c) 2023 DocMinus
"""

import dash
from dash import html
import dash_bootstrap_components as dbc

app = dash.Dash(
    __name__,
    use_pages=True,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    suppress_callback_exceptions=True,
)

debug = True
server = app.server
app.layout = html.Div(children=[dash.page_container])

if __name__ == "__main__":
    app.run(debug=debug, port=8069)
