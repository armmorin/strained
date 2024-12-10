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
NEB_RSC = "400:1:xeon40el8:50h" #f"{CORES*2}:1:{NODE}:50h"
REL_RSC = f"{CORES}:1:{NODE}:50h"
TEST_RSC = f"48:1:xeon24el8_test:10m"

current_dir = Path(__file__).resolve().parent
#sys_name = argv[1].split(",")

# def find_id_from_name(id:int) -> list[int]:
#     """ Find the id of the entry in the database from the name of the entry. """
#     db = connect("structures/hexag_perovs_strained.db")
#     # If the name is found, return the id.
#     if id:
#         row = db.get(id=id)
#         return [row.name]
#     # If the name is not found, return None
#     else:    
#         print(f"Entry with d id {id} not found in the database.")
#         return []

# Calculate the NEB for the initial and final structures. After the conventional NEB has converged, use climbing image. Save the barriers.    
#args = {}
db_id = 43
neb = Task(name= f'NEB_{db_id}', code = current_dir / "neb.py", args={'db_id':db_id}, resources=NEB_RSC)
cineb = Task(name = f'CINEB_{db_id}', code = current_dir / "neb.py", args={'climb':True}, resources=NEB_RSC)

wf = Workflow({neb: [], cineb: [neb]})

with PersistentQueue() as pq:
    pq.submit(wf)