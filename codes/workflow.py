from ase.db import connect
from pathlib import Path
from sys import argv
import numpy as np
from herculestools.dft import RunConfiguration
from perqueue import PersistentQueue
from perqueue.task_classes.task import Task
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup

def find_id_from_name(name:str):
    """ Find the id of the entry in the database from the name of the entry. """
    db = connect(RunConfiguration.structures_dir / "hexag_perovs_re-nebs.db")
    # If the entry comes with a separator, then it is a list of names.
    if "," in name:
        names = name.split(",")
        db_ids = []
        for n in names:
            for row in db.select(name=n):
                db_id = row.id
                db_ids.append(db_id)
        return db_ids
    
    # If the entry is a single name, then find the id of the entry.    
    elif name:
        row = db.get(name=name)
        db_id = row.id
        return [db_id]
    
    # If the name is not found, return None
    else:    
        print(f"Entry with name {name} not found in the database.")
        return None

args = {}
current_dir = Path(__file__).resolve().parent
sys_name = argv[1]

# The workflow is based on the following steps:
# 1. For each name in the database, get the id of the entry.
#    The id can be a single id or a list of ids.
db_ids = [db_id for db_id in find_id_from_name(sys_name)]
#RESOURCES = "24:1:xeon24el8_test:10m"
RESOURCES = "56:1:xeon56:50h"
# 2. For each id, apply a range of in-plane strains.
for db_id in db_ids:
    args["system_id"] = db_id
    WIDTH=5
    strain_range = np.linspace(-2, 2, WIDTH)
    # 2.1 For each in-plane strain, apply a range of out-of-plane strains contained in a StaticWidthGroup.
    for i in strain_range:
        args["in_plane"]=i
        task = Task(current_dir / "apply_strain.py", args=args.copy(), resources=RESOURCES)
        strain_swg = StaticWidthGroup(task, width=WIDTH)
        #print(strain_swg)
    
    # 3 Read into the energies for each in-plane strain and apply a polynomial fit to every out-of-plane strains.
    # Relax the structure and calculate the energy and the band gap.
    op_min = Task(current_dir / "apply_fit.py", args=args.copy(), resources="1:xeon24el8:10m")
    #print(op_min)
    # 4. Generate the vacancies on the new structure and calculate the energy.
    # This is done for the initial and the final structures of the NEB.
    preneb = Task(current_dir / "create_vacancies.py", args=args.copy(), resources=RESOURCES)
    #print(preneb)
    # 5. Calculate the NEB for the initial and final structures. Do one internal image, use climbing image. Save the barriers.
    neb = Task(current_dir / "neb.py", args=args.copy(), resources="112:1:xeon56:50h")
    #print(neb)

    wf = Workflow([strain_swg, op_min, preneb, neb])
    
    with PersistentQueue() as pq:
        pq.submit(wf)
    
    #break    
# 6.* Plot the barriers against the strain for each structure.




