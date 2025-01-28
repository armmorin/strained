from ase.db import connect
from pathlib import Path
from sys import argv
import numpy as np
from herculestools.dft import RunConfiguration as RC
from perqueue import PersistentQueue
from perqueue.task_classes.task import Task
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup
import math

NCORES = 56
NNODES = 1 
if NCORES == 96:
    NODE = "epyc96"
elif NCORES == 56:
    NODE = "xeon56"
else:
    NODE = f"xeon{NODE}el8"
NODE = "xeon56"
REL_RSC = f"{NCORES*NNODES}:1:{NODE}:50h"
TEST_RSC = "48:1:xeon24el8_test:30m"

def find_id_from_name(name:str) -> list[int]:
    """ Find the id of the entry in the database from the name of the entry. """
    db = connect("structures/hexag_perovs_wdiscards.db")
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

# Calculate the NEB for the initial and final structures. After the conventional NEB has converged, use climbing image. Save the barriers.    
args = {}
for name in sys_name:
    db_id = find_id_from_name(name)[0]
    args = {"sys_id":db_id, "name":name}
    args['vasp'] = {'symprec': 1e-09, 'isym': 0, 'algo': 'VeryFast'}
    relax = Task(name=f'relax_{name}', code=current_dir / "relax.py", args=args, resources=TEST_RSC)  # Assign the system id to the task.
    
    args['strain_list'] = [-.75, 0.75, 2.5]
    #args['strain_list'] = [-5, -3, -1.5, 1.5, 3]#np.linspace(-5, 5, 2)
    args['mask_list'] = ['biaxial']
    #args['mask_list'] = ['biaxial', 'x_axis', 'y_axis']
    args['shape'] = (len(args['strain_list']), len(args['mask_list']))
    strain = Task(name=f'strain_{name}',code=current_dir / "apply_strain.py", args=args, resources=REL_RSC)    
    
    swg = StaticWidthGroup([strain], width=math.prod(args['shape']))
    
    wf = Workflow({relax: [], swg: [relax]})
    
    with PersistentQueue() as pq:
        pq.submit(wf)