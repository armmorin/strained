from ase.db import connect
from pathlib import Path
from sys import argv
import numpy as np
from herculestools.dft import RunConfiguration
from perqueue import PersistentQueue
from perqueue.task_classes.task import Task
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup

def find_id_from_name(name:str) -> list[int]:
    """ Find the id of the entry in the database from the name of the entry. """
    db = connect(RunConfiguration.structures_dir / "hexag_perovs_re-nebs.db")
    # If the entry comes with a separator, then it is a list of names.
    if "," in name:
        names = name.split(",")
        db_ids = []
        for n in names:
            for row in db.select(name=f"{n}_r"):
                db_ids.append(row.id)
        return db_ids

    # If the entry is a single name, then find the id of the entry.    
    elif name:
        row = db.get(name=f"{name}_r")
        return [row.id]

    # If the name is not found, return None
    else:    
        print(f"Entry with name {name} not found in the database.")
        return []

current_dir = Path(__file__).resolve().parent
sys_name = argv[1]

WIDTH=5
strain_range = np.linspace(-2, 2, WIDTH)
#RESOURCES = "48:1:xeon24el8_test:30m"

# Generate the vacancies on the new structure and calculate the energy.
preneb = Task(current_dir / "create_vacancies.py", args={}, resources="56:1:xeon56:50h")
#preneb = Task(current_dir / "create_vacancies.py", args={}, resources=RESOURCES)

# Calculate the NEB for the initial and final structures. Do one internal image, use climbing image. Save the barriers.
neb = Task(current_dir / "neb.py", args={}, resources="112:1:xeon56:50h")

for db_id in find_id_from_name(sys_name):
    # For each in-plane strain, apply a range of out-of-plane strains contained in a StaticWidthGroup.
    args = {"system_id": db_id}
    for i in strain_range:
        args["in_plane"]=i
        strain = Task(current_dir / "apply_strain.py", args=args.copy(), resources="40:1:xeon40:50h")
        #strain = Task(current_dir / "apply_strain.py", args=args.copy(), resources=RESOURCES)
        #strain_swg = StaticWidthGroup(task, width=WIDTH)
        
        wf = Workflow([strain, preneb, neb])

        with PersistentQueue() as pq:
            pq.submit(wf)
        #break