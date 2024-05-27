from ase.io import read
from pathlib import Path
from ase.db import connect

""" This script is used to extract the transition states
of the requested trajectory files and make new directories for them.
"""

# Define the names of the trajectory files to be read
def main(**kwargs):
    names = ["Ba7Nb4MoO20"] #kwargs["names"]
    dopant_sites = ["p1", "p2", "p3"]
    systems = [name + "_" + site for name in names for site in dopant_sites]

    for system in systems:
        newdir = Path.cwd() / "systems" / system
        newdir.mkdir(exist_ok=True)
    
    return True, None



    