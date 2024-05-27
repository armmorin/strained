from os import environ
from rich.console import Console
import numpy as np
#from pathlib import Path
from ase.io import write
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from ase.atoms import Atoms
from ase.dft.bandgap import bandgap
from herculestools.dft import RunConfiguration, create_Vasp_calc, set_magnetic_moments

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

# This function will take the energy values of the out-of-plane strain and apply a polynomial fit to them. These values can be read from the database.
def main(vasp= {}, **kwargs):
    
    """
    This function will apply a polynomial fit to the energy values of the out-of-plane strain.
    The minimum of the polynomial will be taken as the out-of-plane strain value and
    a new structure will be created which will be relaxed and saved in the new database.
    """
    
    # Load the configuration
    RunConfiguration.load()
    old_db_path = RunConfiguration.structures_dir / "hexag_perovs_re-nebs.db"
    old_db: SQLite3Database
    old_db = connect(old_db_path)
    db_path = RunConfiguration.structures_dir / "hexag_perovs_strained.db"
    db: SQLite3Database
    db = connect(db_path)

    sys_id = kwargs["system_id"]
    entries = db.select(sys_id=sys_id)
    values = []
    for entry in entries:
        ip_dist = entry.ip_distortion
        direc = entry.direc.parent
        values.append([entry.op_distortion, entry.energy])

    values = np.array(values)
    # Sort the values by the out-of-plane strain from smallest to largest.
    values = values[values[:,0].argsort()]

    # Apply a polynomial fit to the energy values and get the minimum value of the polynomial.
    poly = np.polyfit(values[:,0], values[:,1], 3)
    # Get the zero of the polynomial
    op_min = -poly[1]/(2*poly[0])

    row_db = old_db.get(id=sys_id)
    atoms_row = row_db.toatoms()
    min_atoms = atoms_row.copy()
    min_atoms.repeat(2,2,1)
    get_cell = min_atoms.cell
    new_cell = get_cell * [ip_dist, ip_dist, op_min]
    min_atoms.set_cell(new_cell, scale_atoms=True)

    # Split the name of the entry to get the components of the name
    en_name = row_db.name
    name_components = [en_name.split("_")]
    dops = name_components[1][-1]
    db_name = "_".join(name_components[0:2])
    min_name = f"{db_name}_{ip_dist}"
    # Write the new structure as a file
    atoms_dir = f"{direc}/{min_name}.traj"
    write(atoms_dir, min_atoms)
    
    calc = create_Vasp_calc(min_atoms, 'PBEsol', f"{db_name}_{ip_dist}_min.traj", direc)
    set_magnetic_moments(min_atoms)
    min_atoms.set_pbc([True, True, True])
    calc.set(**vasp,
            #nsw = 250,
            #ibrion=2,
            #eddifg = -0.05,
            kpar = nnodes,
            ncore = ncore)
    min_atoms.get_potential_energy()
    gap, *_ = bandgap(min_atoms.calc, direct=True)
    
    # Save the new structure in the new database
    min_id = update_or_write(db_path, min_atoms, min_name, dopant=dops, ip_distortion=ip_dist, op_distortion=op_min, gap=gap, dir=direc)

    return True, {"db_id": min_id, "op_min": op_min, "ip_dist": ip_dist}
    
def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)
    