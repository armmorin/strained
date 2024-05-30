from os import environ
from rich.console import Console
import numpy as np
from pathlib import Path
from ase.io import write
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from ase.atoms import Atoms
from ase.dft.bandgap import bandgap
from herculestools.dft import RunConfiguration, create_Vasp_calc, set_magnetic_moments
import matplotlib.pyplot as plt
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

    db_ids = list(kwargs['db_id'].values())
    values = []
    
    for db_id in db_ids: # It fails here. AttributeError: 'int' object has no attribute 'count'
        entry = db.get(id=db_id)
        ip_dist = entry.in_plane
        values.append([entry.out_of_plane, entry.energy])

    print(f"The in-plane strain is {ip_dist}.")
    values = np.array(values)
    # Sort the values by the out-of-plane strain from smallest to largest.
    values = values[values[:,0].argsort()]
    print(f"Values: \n {values}")
    # Apply a polynomial fit to the energy values and get the minimum value of the polynomial.
    poly = np.polyfit(values[:,0], values[:,1], 2)
    # Get the zero of the polynomial
    op_min = -poly[1]/(2*poly[0])

    row_db = old_db.get(id=kwargs['system_id'])
    print(f"System name: {row_db.name}")
    atoms = row_db.toatoms()
    print(f"Size of the supercell: \n {atoms.get_cell()}")
    get_cell = atoms.cell
    print(f"A distortion of {ip_dist} in the in-plane direction and {op_min} in the out-of-plane direction will be applied.")
    new_cell = get_cell * [ip_dist, ip_dist, op_min]
    print(f"New size of the supercell: \n {new_cell}")
    op_min = round(op_min, 2)
    atoms.set_cell(new_cell, scale_atoms=True)

    # Create new variables depending on the values of the strain.
    ip_distortion = f"{(ip_dist - 1)*100:.2f}"
    print("ip_distortion=(ip_dist - 1)*100")
    ip_distortion = round(float(ip_distortion))
    print(f"{ip_distortion}=({ip_dist})-1*100")      
    if ip_dist < 1:
        name_ip = f"c{abs(ip_distortion)}"
    elif ip_dist > 1:
        name_ip = f"s{ip_distortion}"
    else:
        name_ip = f"e{ip_distortion}"

    print(name_ip)
    
    # Split the name of the entry to get the components of the name
    en_name = row_db.name
    name_components = en_name.split("_")
    dops = int(name_components[1][-1])
    db_name = "_".join(name_components[0:2])
    print(f"System name: {db_name}")
    direc = Path(f"supercell/{db_name}/{name_ip}/fit")
    print(direc)
    direc.mkdir(parents=True, exist_ok=True)
    min_name = f"{db_name}_{name_ip}_fit"
    print(min_name)
    # Write the new structure as a file
    atoms_dir = f"{direc}/{min_name}.traj"
    write(atoms_dir, atoms)

    # Generate a plot of the polynomial fit. UNTESTED
    x = np.linspace(values[0,0], values[-1,0], 100)
    y = np.polyval(poly, x)
    plt.plot(x, y, label="Polynomial fit")
    plt.scatter(values[:,0], values[:,1], label="Data points")
    plt.xlabel("Out-of-plane strain")
    plt.ylabel("Energy")
    plt.title(f"Polynomial fit for {min_name}")
    plt.savefig(f"{direc}/{min_name}.png")
    
    calc = create_Vasp_calc(atoms, 'PBEsol', direc, direc)
    set_magnetic_moments(atoms)
    atoms.set_pbc([True, True, True])
    calc.set(**vasp,
            nsw = 250,
            ibrion=2,
            ediffg = -0.05,
            kpar = nnodes,
            ncore = ncore)
    atoms.get_potential_energy()
    print(f"Potential energy of the new structure: {atoms.get_potential_energy()}")
    gap, *_ = bandgap(atoms.calc, direct=True)
    
    # Save the new structure in the new database
    min_id = update_or_write(db, atoms, name=min_name, dopant=dops, in_plane=ip_dist, out_of_plane=op_min, gap=gap, dir=direc.as_posix())

    return True, {"db_id": min_id}
    
def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)
