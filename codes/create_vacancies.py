from os import environ
from pathlib import Path
from typing import Optional, Tuple
from xml.etree.ElementTree import ParseError
from rich.console import Console
from ase.calculators.calculator import ReadError
from ase import Atoms
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from ase.io import read
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
            print(f"Supercell is generated in {direc.name}")
            sc = atoms.repeat((2, 2, 1))
            atoms = sc.copy()
            atoms.pop(vacancy) # type: ignore
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

    Args:
        db_id (Optional[int]): The database id of the structure.
        op_min (Optional[float]): The minimum value of the out-of-plane strain.
        ip_dist (Optional[float]): The in-plane strain value.
        **kwargs: Additional keyword arguments.

    Returns:
        Tuple[bool, Optional[dict]]: A tuple containing a boolean indicating whether the calculation was successful,
        and an optional dictionary with the initial and final ids, initial and final vacancies.
    """
    
    RunConfiguration.load()
    db_path = RunConfiguration.structures_dir / "hexag_perovs_strained.db"
    db: SQLite3Database
    db = connect(db_path)

    # Get the name of the structure from the database.
    entry = db.get(kwargs['db_id'])
    en_name = entry.name
    name_components = en_name.split("_")
    name = "_".join(name_components[0:2])
    dops = int(name_components[1][-1])
    atoms = entry.toatoms()

    # Set up directories:
    i_direc = Path(RunConfiguration.home / f"preNEB/{name}_init")
    f_direc = Path(i_direc.parent / f"{name}_final")
    f_direc.mkdir(parents=True, exist_ok=True)

    initial_vac = 30
    init = start_run(atoms=atoms, direc=i_direc, vacancy=initial_vac)
    i_energy = init.get_potential_energy()
    print(f"The energy of the initial structure is: {i_energy:.3f} eV")

    final_vac = 31
    final = start_run(atoms=atoms, direc=f_direc, vacancy=final_vac)
    f_energy = final.get_potential_energy()
    print(f"The energy of the final structure is: {f_energy:.3f} eV")

    # Get the energy difference from the two positions.
    dE = abs(i_energy - f_energy)

    # Save the result to the database
    iID = update_or_write(db, init,  name+"_vi", dopant=dops, dir=i_direc.as_posix(), delta_e = dE)
    fID = update_or_write(db, final, name+"_vf", dopant=dops, dir=f_direc.as_posix(), delta_e = dE)

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
