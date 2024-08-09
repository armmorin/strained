from perqueue import PersistentQueue
from perqueue.selection import Selection
from ase.db import connect
from pathlib import Path

with PersistentQueue() as pq:
    entries = pq.get_entries()

s = Selection(states="s")
targets = s.filter(entries)
db = connect('structures/hexag_perovs_strained.db')
        
def rm_wavecar(id: int):
    entry = db.get(selection=f'id={id}')
    name = entry.name
    jobdir = Path(entry.dir)
    direc_list = [direc.parent for direc in jobdir.rglob("vasprun.xml")]
    for direc in direc_list:
        print(direc)
        wavecar = direc / "WAVECAR"
        chg = direc / "CHG"
        procar = direc / "PROCAR"
        chgcar = direc / "CHGCAR"
        files = [wavecar, chg, procar, chgcar]
        for file in files:
            if file.exists():
                file.unlink()
                print(f"Removing {file.name} for {name}")
 
for item in [en.data for en in targets]:
    if 'initial_id' in item:
        initial_id = item['initial_id']
        rm_wavecar(initial_id)

    if 'final_id' in item:
        final_id = item['final_id']
        rm_wavecar(final_id)
    
    if 'neb_id' in item:
        neb_id = item['neb_id']
        rm_wavecar(neb_id)
    
    if 'db_id' in item:
        db_id  = item['db_id']
        rm_wavecar(db_id)

    else:
        pass
