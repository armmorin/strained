from ase.io import write
from os import environ
from rich.console import Console
from ase.atoms import Atoms
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from ase.optimize import FIRE
from ase.constraints import UnitCellFilter
from pathlib import Path
import shutil
from ase.io.trajectory import Trajectory
from herculestools.dft import (
    RunConfiguration,
    create_Vasp_calc,
    set_magnetic_moments
)

here = Path(__file__).parent.parent

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

def main(vasp:dict = {}, **kwargs):

    """
    This function takes an entry from the database, creates a supercell with a given in-plane strain, 
    and relaxes the structure in the out-of-plane direction by applying a mask to the structure.
    """

    # Load the configuration
    db_path = RunConfiguration.structures_dir / "hexag_perovs_strained.db"
    db: SQLite3Database

    # Connect to the database
    db = connect(db_path)
    
    # We use the id of the entry in the database to get the structure
    db_id = kwargs["db_id"]
    entry = db.get(db_id)
    en_name = entry.name
    dops = entry.dopant
    
    # Split the name of the entry to get the components of the name
    name_components = en_name.split("_")
    
    # Get the new name of the job with the first two components of the name
    db_name = "_".join(name_components[0:2])
    atoms = entry.toatoms()
    
    # Take the in_plane from the INDEX_KW so that it matches the index of the strain_range
    in_plane = kwargs['in_plane']#[INDEX_KW]
        
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
    job_dir = RunConfiguration.home / 'distorted' / db_name / kwargs["mask"] / name_ip
    job_dir.mkdir(parents=True, exist_ok=True)
    direc = job_dir.relative_to(RunConfiguration.home)
    
    # Check if there is not an already relaxed structure.
    counter = len(trajectories := [traj for traj in (job_dir.glob("*.traj")) if traj.stat().st_size > 0])
    
    # Apply a mask to the structure to relax. The mask can have different shapes.
    mask_dict = {'biaxial': (0,0,1,0,0,0),
            'uniaxial': (0,1,1,0,0,0)}
    
    # Depending on the type of mask we apply the strain to the structure
    distortion_dict = {'biaxial': (ip_distortion, ip_distortion, 1),
                        'uniaxial': (ip_distortion, 1, 1)}
    
    # From the reading of the mask keyword, we get the mask to apply to the structure.
    mask_name = kwargs['mask']
    distortion = distortion_dict[mask_name]
    
    # If there are no previous trajectory files generated, we apply the strain from the beginning
    if  counter == 0:
        atoms.set_cell(atoms.get_cell() * distortion, scale_atoms=True)
        atoms.set_pbc([True, True, True])
        set_magnetic_moments(atoms)

    else:
        # Sort the trajectories by the most recent
        trajectories.sort(key=lambda x: x.stat().st_mtime)
        # Get the last atoms object from the most recent trajectory
        atoms = Trajectory(trajectories[-1], 'r')[-1]
    
    # Lastly, we apply the mask to the structure
    ucf = UnitCellFilter(atoms, mask=mask_dict[mask_name])

    # Create the VASP calculator
    calc = create_Vasp_calc(atoms, 'PBEsol', direc)
    calc.set(**vasp,
            kpar = nnodes,
            ncore = ncore)
    
    # Get the potential energy of the new structure
    traj_name = f"{direc}/{db_name}_{name_ip}_{counter}.traj"
    print(traj_name)
    traj = Trajectory(traj_name, 'w', atoms)
    atoms.get_potential_energy()
    opt = FIRE(ucf, logfile=f"{direc}/{db_name}_{name_ip}.log",
               #dt=0.01, maxstep=0.05, dtmax=0.2, Nmin=15, finc=1.03, fdec=0.6
               )
    opt.attach(traj)
    opt.run(fmax=0.03)
    
    # Move the important files from scratch to home
    for file in job_dir.glob("*"):
        if file.is_file() and file.name != 'WAVECAR':
            # Take the current path and move it to the home directory
            file_parts = list(file.parts)
            new_path = file_parts[:2] + ['energy'] + file_parts[3:]
            new_file = Path(*new_path)
            new_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(file, new_file)
    
    # Save the new structure in the new database
    new_id = update_or_write(db, atoms, name=f"{db_name}_{kwargs['mask']}_{name_ip}", dopant=dops, in_plane=ip_distortion, mask=kwargs['mask'], dir=direc.as_posix())
    
    # Copy the database back to the home directory
    home_db = here / "structures/hexag_perovs_strained.db"
    shutil.copy(db_path, home_db)
    
    return True, {"db_id": new_id, "in_plane": in_plane}

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)
