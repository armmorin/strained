from ase.io import write
#from pathlib import Path
from os import environ
from rich.console import Console
from ase.atoms import Atoms
from ase.db import connect
from ase.db.sqlite import SQLite3Database
import numpy as np
from perqueue.constants import INDEX_KW
from herculestools.dft import (
    RunConfiguration,
    create_Vasp_calc,
    set_magnetic_moments
)

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

def main(vasp:dict = {}, **kwargs):

    """
    This function applies a strain to a structure and saves the new structure in a new database.
    The strain is applied to the in-plane and out-of-plane directions of the structure.
    """

    # Take the in_plane and out_of_plane strain values
    in_plane = kwargs["in_plane"]
    out_of_plane_idx, *_ = kwargs[INDEX_KW]
    out_of_plane = np.linspace(-2, 2, 5)[out_of_plane_idx]
    
    # Load the configuration
    RunConfiguration.load()
    db_path = RunConfiguration.structures_dir / "hexag_perovs_re-nebs.db"
    db: SQLite3Database

    # Connect to the database
    db = connect(db_path)
    
    # We use the id of the entry in the database to get the structure
    entry = db.get(kwargs["system_id"])
    en_name = entry.name
    # Split the name of the entry to get the components of the name
    name_components = en_name.split("_")
    # Take the dopant position from the last character of the second component of the name
    dops = int(name_components[1][-1])
    
    # Get the new name of the job with the first two components of the name
    db_name = "_".join(name_components[0:2])
    atoms = entry.toatoms()

    # Create new variables depending on the values of the strain.
    ip_distortion = (1 + in_plane/100)
    if ip_distortion < 1:
        name_ip = f"c{abs(in_plane)}"
    elif ip_distortion > 1:
        name_ip = f"s{in_plane}"
    else:
        name_ip = f"e{in_plane}"
    
    op_distortion = (1 + out_of_plane/100)

    if op_distortion < 1:
        name_op = f"c{abs(out_of_plane)}"
    elif op_distortion > 1:
        name_op = f"s{out_of_plane}"
    else:
        name_op = f"e{out_of_plane}"
  
    # Create a new directory to save the new structures
    job_dir = RunConfiguration.home / 'supercell' / db_name / name_ip / name_op
    job_dir.mkdir(parents=True, exist_ok=True)
    
    # Make a supercell of the structure
    supercell = atoms.repeat((2, 2, 1))
    supercell.set_cell(supercell.get_cell() * [ip_distortion, ip_distortion, op_distortion], scale_atoms=True)
    print(supercell)
    supercell.set_pbc([True, True, True])
    set_magnetic_moments(supercell)
    
    # Create the VASP calculator
    calc = create_Vasp_calc(supercell, 'PBEsol', job_dir, job_dir)
    calc.set(**vasp,
            #nsw = 250,
            #ibrion=2,
            #eddifg = -0.05,
            kpar = nnodes,
            ncore = ncore)
    
    # Get the potential energy of the new structure
    traj_name = f"{job_dir}/{db_name}_{name_ip}_{name_op}.traj"
    print(traj_name)
    write(traj_name, supercell)
    supercell.get_potential_energy()
    
    # Save the new structure in the new database
    db_new_path = RunConfiguration.structures_dir / "hexag_perovs_strained.db"
    db_new = connect(db_new_path)
    db_id = update_or_write(db_new, supercell, name=db_name, sys_id=kwargs["system_id"], dopant=dops, in_plane=ip_distortion, out_of_plane=op_distortion, dir=str(job_dir))
    
    return True, {"db_id": db_id}

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)