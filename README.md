[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License](https://img.shields.io/github/license/tbd)](https://github.com/tbd/blob/master/LICENSE)
[![version](https://img.shields.io/github/v/release/tbd)](https://github.com/tbd/releases)

# ChemDash
## Proof of concept small molecule (web) GUI
Here is a proof of concept of a "small molecule chemistry centric gui", running in a web-browser<br>
It's build in Pyton and uses the Plotly Dash components.<br>
For installation, see the environment folder.<br>

Files:<br>
```
python single_page_example.py 
```
A simple input form (one one single page) to enter a smiles which calculates some values that you can easily extend with additional numbers.
```
python multi_page_example.py
```
A concept for a dashboard where you could have multiple pages for different purposes, accessible via a menu. This example also contains a structural editor version of the above example.
<br>
### Comments and acknowledgements
Dash is a somewhat straightforward way to make a GUI, especially in context of data-(analysis). The callbacks though aren't always that straightforward imho to use. And although html is used, dash uses it's own wrapper which can make things a bit confusing. Then again, it could just be me.<br>
A number of tutorials helped me, among them I want to thank and highlight:<br> 
* [ArjanCodes dash tutorial](https://github.com/ArjanCodes/2022-dash), including the Youtube videos
* [Algorithms](https://github.com/siddharthajuprod07/algorithms), including the Youtube videos
