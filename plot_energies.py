import matplotlib.pyplot as plt
import numpy as np
from  pathlib import Path
from ase.io import read

# Plot the total energies of the systems in the optimized geometries
optimized_dir = Path('optimized')

optimized_dirs = [d for d in optimized_dir.iterdir() if d.is_dir()]

# Descend on each subdirectory in the optimized directory
for sub_dir in optimized_dirs:
    # Get the name of the subdirectory
    sub_dir_name = sub_dir.name
    
    # Descend on each subdirectory
    