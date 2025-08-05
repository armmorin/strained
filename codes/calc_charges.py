from pathlib import Path
from typing import Optional
from ase import Atoms
from ase.io import read, write
from ase.db.sqlite import SQLite3Database
from ase.calculators.calculator import ReadError
from ase.spacegroup import get_spacegroup
#from rich import Console
from os import environ
from ase.db import connect
import shutil
from herculestools.dft import (
    RunConfiguration as RC, 
    create_Vasp_calc,
    set_magnetic_moments
    )

#c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

#here = Path.cwd() # Where the code is located
here = Path(__file__).resolve().parent.parent
home = RC.home # Where the calculations are stored
struc_dir = RC.structures_dir # Where the structures are stored

# Connect to the database
db_path = f"{struc_dir}/hexag_perovs_strained.db"

#e_config_dict = {'Hf': 10, 'Zr': 12, 'O': 6}

# Define a function to generate and calculate the pristine structures. 
# The output will be used to generate the alloy structures.
def main(db_id:int=1, vasp: dict = {}, **kwargs):
    """Calculate the partial charges of the given systems"""
    
    # Read from the database the given structure
    db = connect(db_path)
    row = db.get(db_id)
    atoms = row.toatoms()
    #electrons = sum([e_config_dict[atom.symbol] for atom in atoms])
    #bands = int(electrons / 2)
    # Create a list of the last 8 bands
    #ibands = sorted([bands - i for i in range(16)])
    row_dir = row.dir
    sys_dir = Path(home / row_dir)
    try:
        element = kwargs['element']
        chg_dir = sys_dir / element
        atoms = read(chg_dir / 'POSCAR')
    except KeyError as e:
        print(f"Error: {e}")
        chg_dir = sys_dir / 'charges'
        chg_dir.mkdir(parents=True, exist_ok=True)
    
    # Remove the vasprun on the charges directory
    vasprun = chg_dir / 'vasprun.xml'
    if vasprun.exists():
        vasprun.unlink()

    # Copy the WAVECAR file to the charges directory
    wavecar = sys_dir / 'WAVECAR'
    if wavecar.exists():
        shutil.copy(wavecar, chg_dir)
        
    # Setup the vasp parameters
    vasp_settings = {   
        'ngxf': 240,
        'ngyf': 240,
        'ngzf': 240,
        'lcharg' : True,
        'lwave' : False,
        'laechg' : True,
        'ncore' : ncore,
    }

    vasp_settings.update(vasp)
    
    # Calculate the energy of the structure
    set_magnetic_moments(atoms)    
    calc = create_Vasp_calc(atoms, 'PBEsol', chg_dir)
    calc.set(
        **vasp_settings
        )
    
    # Calculate the bands
    atoms.get_potential_energy()
    
    
    # Calculate the energy of the structure
    # try:
    #     atoms.get_potential_energy()
    # except ReadError:
    #     # The vasprun is faulty, so we need to rerun the calculation
    #     vasprun = chg_dir / 'vasprun.xml'
    #     vasprun.unlink()
        
    #c.log(f"Calculated energy for {name} in directory {subdir}")
    
    # Move the important files from scratch to home
    for file in chg_dir.glob("*"):
        if file.is_file() and file.name != 'WAVECAR':
            # Take the current path and move it to the home directory
            file_parts = list(file.parts)
            new_path = file_parts[:2] + ['energy'] + file_parts[3:]
            new_file = Path(*new_path)
            new_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(file, new_file)

    return True, {'db_id': db_id}

if __name__ == '__main__':
    from sys import argv
    index = int(argv[1])
    args = [float(arg) for arg in argv[2:]]

    print(main(index, *args))
    