from ase.db import connect
from pathlib import Path
from sys import argv
import numpy as np
from herculestools.dft import RunConfiguration
from perqueue import PersistentQueue
from perqueue.task_classes.task import Task
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup
import itertools

CORES = 96 
NODE = "epyc96"
NEB_RSC = f"{CORES*2}:1:{NODE}:50h"
REL_RSC = f"{CORES}:1:{NODE}:50h"
TEST_RSC = f"48:1:xeon24el8_test:10m"
WIDTH = 5

def find_id_from_name(name:str) -> list[int]:
    """ Find the id of the entry in the database from the name of the entry. """
    db = connect(RunConfiguration.structures_dir / "hexag_perovs_wdiscards.db")
    # If the name is found, return the id.
    if name:
        row = db.get(name=f"{name}_r")
        return [row.id]
    # If the name is not found, return None
    else:    
        print(f"Entry with name {name} not found in the database.")
        return []

current_dir = Path(__file__).resolve().parent
sys_name = argv[1].split(",")

strain_range = np.linspace(-3, 3, WIDTH)
masks = ['biaxial', 'uniaxial']

# Calculate the NEB for the initial and final structures. Do one internal image, use climbing image. Save the barriers.
neb = Task(current_dir / "neb.py", args=None, resources=NEB_RSC)
preneb = Task(current_dir / "preneb.py", args=None, resources=REL_RSC)

swf = Workflow([preneb, neb])
        
for name in sys_name:
    db_id = find_id_from_name(name)[0]
    args = {"sys_id": db_id}
    args['name'] = name
    relax = Task(current_dir / "relax.py", args=args.copy(), resources=REL_RSC)  # Assign the system id to the task.


for mask, i in itertools.product(masks, strain_range):
    args = {}
    args["in_plane"]= i
    args["mask"] = mask    
    strain = Task(current_dir / "apply_strain.py", args=args.copy(), resources=REL_RSC)    
    
    wf = Workflow({relax: [], strain: [relax], swf: [strain]})

    with PersistentQueue() as pq:
        pq.submit(wf)