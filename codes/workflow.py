from ase.db import connect
from pathlib import Path
from sys import argv
import numpy as np
from herculestools.dft import RunConfiguration as RC
from perqueue import PersistentQueue
from perqueue.task_classes.task import Task
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup
import math

# Define the resources for the different tasks.
# The resources are defined based on the number of cores and nodes.
NCORES = 56
NNODES = 1 
if NCORES == 96:
    NODE = "epyc96"
elif NCORES == 56:
    NODE = "xeon56"
else:
    NODE = f"xeon{NODE}el8"
NODE = "xeon56"
NEB_RSC = f"{NCORES*NNODES}:1:{NODE}:50h"
REL_RSC = f"{NCORES*NNODES}:1:{NODE}:50h"
TEST_RSC = "48:1:xeon24el8_test:30m"

# Function to find the id of an entry in the database from its name.
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

# Calculate the NEB for the initial and final structures. 
# After the conventional NEB has converged, use climbing image. 
# Save the barriers.    
args = {}
for name in sys_name:
    # Set up the run configuration.
    args['strain_list'] = [-5, -3, -1.5, 1.5, 3] # Strain values in percent.
    args['mask_list'] = ['biaxial', 'x_axis', 'y_axis'] # Mask types for the strain.
    args['shape'] = (len(args['strain_list']), len(args['mask_list'])) # Shape of the strain list.
    db_id = find_id_from_name(name)[0]
    
    # Define the tasks for the workflow.
    relax = Task(name=f'relax_{name}', code=current_dir / "relax.py", args={"sys_id":db_id, "name":name}, resources="40:1:xeon40el8:50h")  # Assign the system id to the task.
    strain = Task(name=f'strain_{name}',code=current_dir / "apply_strain.py", args=args, resources=REL_RSC)  # Apply strain to the relaxed structure.
    preneb = Task(name= f'preNEB_{name}', code = current_dir / "preneb.py", args=None, resources=REL_RSC) # Prepare the NEB calculation.
    neb = Task(name= f'NEB_{name}', code = current_dir / "neb.py", args=None, resources=NEB_RSC) # Perform the NEB calculation.
    cineb = Task(name = f'CINEB_{name}', code = current_dir / "neb.py", args={'climb':True}, resources=NEB_RSC) # Perform the climbing image NEB calculation.
    swf = Workflow([preneb, neb, cineb]) # Define a subworkflow for the NEB calculations.

    # Define a StaticWidthGroup for the different arguments of the system.
    # This will ensure that the tasks are executed in parallel for each strain and mask combination.
    swg = StaticWidthGroup([strain], width=math.prod(args['shape']))
    
    # Define the workflow for the system.
    wf = Workflow({relax: [], swg: [relax]})
    
    with PersistentQueue() as pq:
        pq.submit(wf)