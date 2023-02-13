import dash
from dash import html

from . import navigation
from src.ids import *


dash.register_page(__name__, path=MENU2SUB1_HREF)


hyperlinks = [
    ("ASKCOS MIT", "https://askcos.mit.edu/"),
    (
        "Organic Chemistry Data - Collections, Principles, Tables (e.g. pka), etc.",
        "https://organicchemistrydata.org",
    ),
    (
        "ACS Tools for Innovative Chemistry",
        "https://www.acsgcipr.org/tools-for-innovation-in-chemistry/",
    ),
    ("Chemix - editor for labdiagrams, etc", "https://chemix.org"),
    ("Docminus's Gibthub:", "https://github.com/DocMinus"),
]

links = [
    html.Li([html.A(link[0], href=link[1], target="_blank")]) for link in hyperlinks
]

layout = html.Div(
    [
        navigation.navbar,
        html.H1("Collection of some more (or less) interesting resources on the web"),
        html.Hr(),
        html.Ul(links),
    ]
)
