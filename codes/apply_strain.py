from ase.io import write
#from pathlib import Path
from os import environ
from rich.console import Console
from ase.atoms import Atoms
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
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
    For each strain value applied to the structure, a mask is created to allow the structure to relax in the out-of-plane direction.
    """

    # Load the configuration
    RunConfiguration.load()
    db_path = RunConfiguration.structures_dir / "hexag_perovs_re-nebs.db"
    db: SQLite3Database

    # Connect to the database
    db = connect(db_path)
    
    # We use the id of the entry in the database to get the structure
    sys_id = kwargs["system_id"]
    entry = db.get(sys_id)
    en_name = entry.name
    # Split the name of the entry to get the components of the name
    name_components = en_name.split("_")
    # Take the dopant position from the last character of the second component of the name
    dops = int(name_components[1][-1])
    
    # Get the new name of the job with the first two components of the name
    db_name = "_".join(name_components[0:2])
    atoms = entry.toatoms()

    # Take the in_plane
    in_plane = kwargs["in_plane"]
    
    # Create new variables depending on the values of the strain. c = compressive, s = tensile, e = no strain
    ip_distortion = (1 + in_plane/100)
    in_plane = int(in_plane)
    if ip_distortion < 1:
        name_ip = f"c{abs(in_plane)}"
    elif ip_distortion > 1:
        name_ip = f"s{in_plane}"
    else:
        name_ip = f"e{in_plane}"
    
    # Create a new directory to save the new structures
    job_dir = RunConfiguration.home / 'best_fit' / db_name / name_ip
    job_dir.mkdir(parents=True, exist_ok=True)
    direc = job_dir.relative_to(RunConfiguration.home)
    
    # Apply the strain to the structure
    atoms.set_cell(atoms.get_cell() * [ip_distortion, ip_distortion, 1], scale_atoms=True)
    print(atoms)
    atoms.set_pbc([True, True, True])
    set_magnetic_moments(atoms)
    
    # Apply a mask to the structure to allow the relaxation in the out-of-plane direction
    mask = (0,0,1,0,0,0)
    UnitCellFilter(atoms, mask=mask)
  
    # Create the VASP calculator
    calc = create_Vasp_calc(atoms, 'PBEsol', direc, direc)
    calc.set(**vasp,
            kpar = nnodes,
            ncore = ncore)
    
    # Get the potential energy of the new structure
    traj_name = f"{direc}/{db_name}_{name_ip}.traj"
    print(traj_name)
    write(traj_name, atoms)
    atoms.get_potential_energy()
    opt = BFGS(atoms, logfile=f"{direc}/opt.log")
    opt.run(fmax=0.05)
    
    # Save the new structure in the new database
    db_new_path = "structures/hexag_perovs_strained.db"
    db_new = connect(db_new_path)
    db_id = update_or_write(db_new, atoms, name=f"{db_name}_{name_ip}", dopant=dops, in_plane=ip_distortion, dir=direc.as_posix())
    
    return True, {"db_id": db_id}

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)