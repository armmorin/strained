from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from ase.dft.bandgap import bandgap
import pandas as pd
from ase.io import read,write
from ase.db 

plt.rcParams.update({'font.size': 11})

# Define the starting position
here = Path(__file__).resolve().parent
