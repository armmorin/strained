from os import environ
from rich.console import Console
from ase.optimize import FIRE
from typing import Optional, Tuple
from ase import Atoms
from ase.constraints import FixAtoms
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from pathlib import Path
from ase.io import read
import shutil
from xml.etree.ElementTree import ParseError
from ase.calculators.calculator import ReadError
from herculestools.dft import (
    RunConfiguration,
    create_VASP_calc,
    set_magnetic_moments,
)

here = Path(__file__).parent.parent
home = RunConfiguration.home
struc_dir = RunConfiguration.structures_dir
c = Console()
#nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

def main(rattle_std: float=0.03,
         distance: float=0.5, vasp: dict={},
         **kwargs) -> Tuple[bool, dict]:
    # Ask the PQ for database id
    
    ref_db_path = "structures/hexag_perovs_wdiscards.db"
    db_path = struc_dir/ "hexag_perovs_strained.db"
    
    db: SQLite3Database
    
    # Load the structure and extract fields
    with connect(ref_db_path) as db:
        entry = db.get(kwargs['sys_id'])
        atoms: Atoms = entry.toatoms()
        nameparts = entry.name.split("_")
        name = "_".join(nameparts[:2])
    
    # Take the dopant position from the last character of the name
    dops = int(name[-1])
    c.log(f"Relaxing db entry {entry.id} with name {nameparts[0]} with the dopant in position {dops}.")
    
    # Set up calculator and auxilliary things
    direc: Path = home / f"relaxations/{name}"
    direc.mkdir(parents=True, exist_ok=True)
    c.log(f"Working in {direc}")
    direc_string = direc.as_posix()

    # Check if the atoms object is a supercell.
    #lets_run = True
    #if len(atoms) == 32:
        # Generate the supercell
        #c.log(f"Generating supercell for {name}")
        #atoms = atoms.repeat([2, 2, 1])

    # Check if the relaxation has been done for the same structure
    # if ((outcar := direc/"vasp.out")).exists():
    #     try:
    #         # Check if the number of atoms in the vasprun is the same as the number of atoms in the atoms object
    #         if len(atoms) == len(read(vasprun)):
    #             # Check if the relaxation has been done
    #                 outcar_txt = outcar.read_text()
    #                 if "reached required accuracy - stopping structural energy minimisation" in outcar_txt:
    #                     c.log(f"Relaxation for {name} has already been done")
    #                     lets_run = False
    #     except (ReadError, ParseError):
    #         c.log("Error reading the vasprun.xml file. Running the relaxation again.")
            
    # Set up the calculator
    
    trajs_list = list(direc.glob("*relax.traj"))
    trajs_list = [traj for traj in trajs_list if traj.is_file() and traj.stat().st_size > 0]
    
    vasprun = direc/'vasprun.xml'
    if len(trajs_list) > 0:
        last_traj = sorted(trajs_list, key=lambda x: int(x.stem.split("_")[-2]))[-1]
        c.log(f"Attempting restart from trajectory: {last_traj.name}")
        atoms = read(last_traj)
        counter = int(last_traj.stem.split("_")[-2]) + 1
        c.log(f"Restarting from step {counter}")
    else:
        counter = 0
        cn = FixAtoms([31])
        # If vasprun.xml is exists, try to read the last structure from it
        if vasprun.exists():
           try:
               atoms = read(vasprun)
               c.log("Read the last structure from vasprun.xml")
           except (ReadError, ParseError):
               c.log("Error reading the vasprun.xml file. Starting from the initial structure.")
        else:
            c.log("No vasprun.xml file found. Starting from the initial structure.")
            atoms = read(direc/"POSCAR")
            atoms.set_constraint(cn)
            atoms.rattle(rattle_std)
            atoms.set_constraint()
            atoms[31].position += distance * atoms.cell[0] / atoms.cell.lengths()[0]
            atoms.wrap()
            c.log(f"Rattled the structure with std {rattle_std} and moved the dopant by {distance} along the x-axis")


    c.log(f"Setting up the calculator for {name}")
    calc = create_VASP_calc(atoms, 'PBEsol', direc.name, direc)
    set_magnetic_moments(atoms)
    vasp_settings = {
            'ibrion' : 2,
            'nsw' : 1,
            #'ediffg' : -0.05,
            'nsw' : 0,
            'ncore' : ncore,
    }
    vasp_settings.update(vasp)
    calc.set(**vasp_settings)
    
    # Does the relaxation as well
    #if lets_run:
    # Apply calculator to atoms
    c.log(f"Relaxing {name}")
    atoms.calc = calc
    atoms.get_potential_energy()
    opt = FIRE(atoms, trajectory=f"{direc_string}/{name}_{counter}_relax.traj", logfile=f"{direc_string}/{name}_relax.log")
    opt.run(fmax=0.05)
    
    # Move the important files from scratch to home
    for file in direc.glob("*"):
        if file.is_file() and file.name != 'WAVECAR':
            # Take the current path and move it to the home directory
            file_parts = list(file.parts)
            new_path = file_parts[:2] + ['energy'] + file_parts[3:]
            new_file = Path(*new_path)
            new_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(file, new_file)
    
    c.log(f"Relaxation for {name} has been done. Files have been moved to {new_file.parent}")
    
    # Save the result to the database
    with connect(db_path) as db:
        DB_ID = update_or_write(db, atoms, name=name+"_r", dopant=dops, dir=direc.relative_to(home).as_posix())
    
    c.log(f"The database has been updated. The new id is {DB_ID}")
    # Copy the database back to the home directory
    home_db = here / "structures/hexag_perovs_strained.db"
    shutil.copy(db_path, home_db)
    
    return True, {'db_id': DB_ID}

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs) -> int:
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
