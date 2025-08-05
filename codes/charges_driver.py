from perqueue.queue import PersistentQueue
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup
from perqueue.task_classes import Task
from pathlib import Path
from ase.db import connect
from itertools import product
from ase.io import write
from ase.visualize import view
from ase import Atoms
from herculestools.dft import RunConfiguration as RC

# Extract the rows from the database
db_path = 'structures/hexag_perovs_strained.db'
db = connect(db_path)

# Define the run configuration
home = RC.home

# Define the different parameters for the calculations
masks = ['x', 'y']
distortions = [1.0, 0.95]


# Create a function to extract the different elements from each of the crystal structures and calculate the charges in separate directories
def extract_elements(db_id: int = 1, db = db):
    db_row = db.get(db_id)
    atoms = db_row.toatoms()
    cell = atoms.get_cell()
    elements = list(set([atom.symbol for atom in atoms]))
    dir = Path(db_row.dir)
    # Generate separate sublists for each of the atomic species
    for element in elements:
        element_dir = home / dir / element
        element_dir.mkdir(parents=True, exist_ok=True)
        #print(element_dir)
        # Remove all other species from the structure
        element_atoms = [atom for atom in atoms if atom.symbol == element]
        element_atoms = Atoms(element_atoms)
        element_atoms.set_cell(cell)
        element_atoms.set_pbc(atoms.get_pbc())
        #view(element_atoms)
        # Store the new structure in the new directory
        write(element_dir / 'POSCAR', element_atoms)
    return elements

# for systems in symmetries and strain and directions and types:
all_ids = []
for mask,dist in product(masks,distortions):
    if dist == 1.0:
        #Get the rows from the database
        row = db.get(
            f"in_plane={dist},climb=True"
            )
    else:
        #print(f"symmetry: {sym}, distortion: {dist}, type: {typ}")
        row = db.get(
                f"mask={mask}_axis,in_plane={dist},climb=True"
            )
    all_ids.append(row.id)
    
          
# From the contents of the list, create a task for each id
for id in all_ids:
    
    totchg = Task(name=f'charges_{id}', code='codes/calc_charges.py', args={'db_id':id}, resources="24:1:xeon24el8:8h")
    
    elements = extract_elements(id)
    
    element_charges = []
    for element in elements:
    
        elmchg = Task(name=f'{element}_charges_{id}', code='codes/calc_charges.py', args={'db_id':id, 'element':element}, resources="24:1:xeon24el8:8h")
        element_charges.append(elmchg)

    swg = StaticWidthGroup(element_charges, width=1)
    wf = Workflow([totchg, swg])

    with PersistentQueue() as pq:
        pq.submit(wf)
    
        #break
    #break
