from os import environ
from rich.console import Console
from typing import Optional, Tuple
from ase import Atoms
from ase.constraints import FixAtoms
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from pathlib import Path
from ase.io import read

from herculestools.dft import (
    RunConfiguration,
    create_VASP_calc,
    set_magnetic_moments,
)

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

def main(rattle_std: float=0.03,
         distance: float=0.5, vasp: dict={},
         **kwargs) -> Tuple[bool, dict]:
    # Ask the PQ for database id
    
    db_path = RunConfiguration.structures_dir / "hexag_perovs_wdiscards.db"
    db_new_path = "structures/hexag_perovs_strained.db"
    
    db: SQLite3Database
    
    # Load the structure and extract fields
    with connect(db_path) as db:
        entry = db.get(kwargs['sys_id'])
        atoms: Atoms = entry.toatoms()
        name = kwargs['name']
    
    # Take the dopant position from the last character of the name
    dops = int(name[-1])
    
    # Set up calculator and auxilliary things
    lets_run = True
    direc: Path = RunConfiguration.home / f"relaxations/{name}"
    direc.mkdir(parents=True, exist_ok=True)
    
    if ((outcar := direc/"vasp.out")).exists():
        outcar_txt = outcar.read_text()
        if "reached required accuracy - stopping structural energy minimisation" in outcar_txt:
            lets_run = False
    
    if not ((fl := direc/'vasprun.xml')).exists():
        # Create the supercell, and rattle all atoms but the specific oxygen
        atoms = atoms.repeat([2, 2, 1])        
        cn = FixAtoms([31])
        atoms.set_constraint(cn)
        atoms.rattle(rattle_std)  # <--- Tolerance #1
        atoms.set_constraint()
        # \/ Parameter #3
        atoms[31].position += distance * atoms.cell[0] / atoms.cell.lengths()[0]
        atoms.wrap()

    else:
        try:
            atoms = read(fl)
        except:
            pass

    calc = create_VASP_calc(atoms, 'PBEsol', direc.name, direc)
    set_magnetic_moments(atoms)
    calc.set(
        kpar = nnodes,
        ncore = ncore,
        ibrion = 2,
        ediffg = -0.05,
        nsw = 250,
        **vasp)
    
    # Does the relaxation as well
    if lets_run:
        # Apply calculator to atoms
        atoms.calc = calc
        atoms.get_potential_energy()

    # Save the result to the database
    with connect(db_new_path) as db:
        DB_ID = update_or_write(db, atoms, name=name+"_r", dopant=dops, dir=direc.as_posix())

    return True, {'db_id': DB_ID}

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)

if __name__ == '__main__':
    from sys import argv
    index = int(argv[1])
    args = [float(arg) for arg in argv[2:]]

    print(main(index, *args))
