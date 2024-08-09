from os import environ
from pathlib import Path
from typing import Optional, Tuple
from xml.etree.ElementTree import ParseError
from rich.console import Console
from ase.calculators.calculator import ReadError
from ase import Atoms
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from ase.io import read, write
from typing import List

from herculestools.dft import (
    RunConfiguration,
    create_Vasp_calc,
    set_magnetic_moments,
)

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

def read_poscar(poscar: Path) -> Atoms | None:
    try:
        atoms: Atoms = read(poscar)
    except Exception as e:
        raise POSCARError(f"Error in {poscar}") from e

    return atoms if atoms is not None else None

class POSCARError(Exception):
    pass

## DEFINING A FUNCTION TO REGULATE THE PROCESS.
def start_run(atoms: Atoms | List[Atoms], direc: Path,
            vacancy: int, vasp: dict={}) -> Atoms:

    # Making a dictionary of the paths to the files.
    files = {
        'contcar' : direc / "CONTCAR",
        'poscar'  : direc / "POSCAR",
        'vasprun' : direc / "vasprun.xml"
    }

    try:
        atoms = read_poscar(files['poscar'])
        atoms = read(files['contcar'])
        atoms = read(files['vasprun'])
        print(f"A structure with {len(atoms)} atoms will start from {direc.name}")
        setup_run(atoms, direc, vasp)
        

    except (POSCARError, IndexError):
        # Add a check to see if the system already has vacancies from the atoms object.
        if len(atoms) == 32:
            raise ValueError("The system is not a supercell")
           
        else:
            print(f"Restarting from supercell structure in {direc.name}")

        setup_run(atoms, direc, vasp)

        return atoms

    except (ReadError, FileNotFoundError, ParseError):
        print(f"No CONTCAR or vasprun.xml found in {direc.name}")
        setup_run(atoms, direc, vasp)

    return atoms   

def setup_run(atoms: Atoms | List[Atoms], direc: Path, vasp: dict = {}) -> Atoms:
    set_magnetic_moments(atoms)
    calc = create_Vasp_calc(atoms, 'PBEsol', direc)
    calc.set(**vasp,
            ediffg=-0.05,
            ibrion=2,
            nsw=250,
            kpar=nnodes,
            ncore=ncore)

    return atoms

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID

    return db.write(atoms, name=name, **kwargs)

# UPDATED for newer versions of HERCULEStools and PerQueue
def main(**kwargs) -> Tuple[bool, Optional[dict]]:
    """
    Perform the preNEB calculation for a given structure.

    """
    
    # Connect to the reference database
    db_path= RunConfiguration.structures_dir / "hexag_perovs_strained.db"
    db = connect(db_path)
    db: SQLite3Database

    # Get the name of the structure from the database.
    print(f"DB id is: {(db_id:=kwargs['db_id'])}")
    entry = db.get(db_id)
    atoms = entry.toatoms()
    en_name = entry.name
    name_components = en_name.split("_")
    dops = int(name_components[1][-1])
    db_name = "_".join(name_components[0:2])
    
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
        
    # Set up directories:
    name = db_name + "_" + name_ip
    i_direc = Path(RunConfiguration.home / f"preNEB/{name}_init")
    i_direc.mkdir(parents=True, exist_ok=True)
    f_direc = Path(i_direc.parent / f"{name}_final")
    f_direc.mkdir(parents=True, exist_ok=True)
    
    initial_vac = 30
    init = start_run(atoms=atoms, direc=i_direc, vacancy=initial_vac)
    if not (traj := f"{i_direc / name_ip}.traj").exists():
        write(traj, init)
    i_energy = init.get_potential_energy()
    print(f"The energy of the initial structure is: {i_energy:.3f} eV")

    final_vac = 31
    final = start_run(atoms=atoms, direc=f_direc, vacancy=final_vac)
    if not (traj := f"{f_direc / name_ip}.traj").exists():
        write(traj, final)
    f_energy = final.get_potential_energy()
    print(f"The energy of the final structure is: {f_energy:.3f} eV")

    # Get the energy difference from the two positions.
    dE = abs(i_energy - f_energy)

    # Save the result to the database
    iID = update_or_write(db, init,  name+"_vi", dopant=dops, dir=i_direc.as_posix(), in_plane=in_plane, delta_e = dE)
    fID = update_or_write(db, final, name+"_vf", dopant=dops, dir=f_direc.as_posix(),  in_plane=in_plane, delta_e = dE)

    print(f"The energy difference between images is :{dE:.3f} eV")

    return True, {
        'initial_id': iID,
        'final_id': fID,
        'initial_vac': initial_vac,
        'final_vac': final_vac
        }

if __name__ == '__main__':
    from sys import argv
    index = int(argv[1])
    args = [arg for arg in argv[2:]]

    print(main(index, *args))
