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
    db = connect(RunConfiguration.home / "../discarded/structures/hexag_perovs_wdiscards.db")
    # if "," in name:
    #     names = name.split(",")
    #     db_ids = []
    #     for n in names:
    #         for row in db.select(name=f"{n}_vi"):
    #             db_ids.append(row.id)
    #    return db_ids
    # elif name:
    #     row = db.get(name=f"{name}_vi")
    #     return [row.id]

    if name:
        row = db.get(name=f"{name}_vi")
        return [row.id]
    # If the name is not found, return None
    else:    
        print(f"Entry with name {name} not found in the database.")
        return []

current_dir = Path(__file__).resolve().parent
sys_name = argv[1].split(",")

WIDTH=5
strain_range = np.linspace(-2, 2, WIDTH)
#RESOURCES = "48:1:xeon24el8_test:30m"
CORES = 96 
NODE = "epyc96"

#RESOURCES="192:1:epyc96:50h"

# Calculate the NEB for the initial and final structures. Do one internal image, use climbing image. Save the barriers.
for name in sys_name:
    print(name)
    neb = Task(current_dir / "neb.py", args={'name':name}, resources=f"{CORES*2}:1:{NODE}:50h")

    for db_id in find_id_from_name(name):
        # For each in-plane strain, apply a range of out-of-plane strains contained in a StaticWidthGroup.
        args = {"system_id": db_id}
        for i in strain_range:
            args["in_plane"]=i
            strain = Task(current_dir / "strained_preneb.py", args=args.copy(), resources=f"{CORES}:1:{NODE}:50h")

            wf = Workflow([strain, neb])

            with PersistentQueue() as pq:
                pq.submit(wf)
            