from ase.db import connect
from pathlib import Path
from sys import argv
import numpy as np
from herculestools.dft import RunConfiguration as RC
from perqueue import PersistentQueue
from perqueue.task_classes.task import Task
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup
import math

CORES = 96 
NODE = "epyc96"
NEB_RSC = f"{CORES*2}:1:{NODE}:50h"
REL_RSC = f"{CORES}:1:{NODE}:50h"
TEST_RSC = f"48:1:xeon24el8_test:10m"

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
    args['strain_list'] = [-3, -1.5, 1.5, 3, 5]#np.linspace(-5, 5, 2)
    args['mask_list'] = ['biaxial', 'x_axis', 'y_axis']
    args['shape'] = (len(args['strain_list']), len(args['mask_list']))
    db_id = find_id_from_name(name)[0]
    #args["sys_id"] = db_id
    
    relax = Task(name=f'relax_{name}', code=current_dir / "relax.py", args={"sys_id":db_id}, resources=REL_RSC)  # Assign the system id to the task.
    strain = Task(name=f'strain_{name}',code=current_dir / "apply_strain.py", args=args, resources=REL_RSC)    
    
    preneb = Task(name= f'preNEB_{name}', code = current_dir / "preneb.py", args=None, resources=REL_RSC)
    neb = Task(name= f'NEB_{name}', code = current_dir / "neb.py", args=None, resources=NEB_RSC)
    cineb = Task(name = f'CINEB_{name}', code = current_dir / "neb.py", args={'climb':True}, resources=NEB_RSC)
    swf = Workflow([preneb, neb, cineb])

    swg = StaticWidthGroup([strain], width=math.prod(args['shape']))
    
    wf = Workflow({relax: [], swg: [relax]})
    
    with PersistentQueue() as pq:
        pq.submit(wf)