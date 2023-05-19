# define diverse permanent variables accessible/used in multiple places.
import os
import sys

CURRENT_SRC_DIR = os.path.abspath(os.path.abspath(sys.argv[0]))
TMP_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_SRC_DIR), "assets"))
MENU_LOGO = "navbar_logo.png"


# menu names and hrefs
MENU1_NAME = "Calculations"
MENU1SUB1_NAME = "Single Molecule (editor)"
MENU1SUB1_HREF = "/calculations_jsme"
MENU1SUB2_NAME = "Single Molecule (smiles)"
MENU1SUB2_HREF = "/calculations_smi"
MENU2_NAME = "Other"
MENU2SUB1_NAME = "Links"
MENU2SUB1_HREF = "/links_of_interest"
