from ase.db import connect
from perqueue.queue import PersistentQueue
from sys import argv
from herculestools.dft import RunConfiguration as RC

home = RC.home

keys = argv[1].split(',')
db_path = home / 'structures' / 'hexag_perovs_strained.db'

codes_dict = {
    'relax' : 'relax.py',
    'strain' : 'apply_strain.py',
    'preneb' : 'preneb.py',
    'neb' : 'neb.py',
}

direcs = {
'neb' : home / 'NEB',
'preneb' : home / 'preNEB',
'strain' : home / 'distorted',
'relax' : home / 'relaxations'
}

# Given a key, or keys, return the corresponding directory of that entry.
for key in keys:
    with PersistentQueue() as pq:
        entry = pq.get_entry(int(key))
        print(f"Getting the directory for entry {entry.key}.")
        code = entry.task.code.name
        print(f"Code name is {code}.")
        if code != codes_dict['neb']:
            try:
                #db_id = entry._task.args['db_id']
                db_id = entry.data['db_id']
                print(f"db_id is {db_id}.")
            except KeyError:
                print(f"db_id not found in the task arguments for entry {entry.key}.")
                continue
        else:        
            db_id = entry.data['neb_id']
            print(f"db_id is {db_id}.")
        
        with connect(db_path) as db:
            row = db.get(id=db_id)
            row_dir = row.dir
            print((home / row_dir).resolve())
        #     split_name = row.name.split('_')
        #     joint_name = '_'.join(split_name[:-1])
        #     # Get the name of the entry from the database if possible
        #     if row.dir:
        #         job_dir = row.dir
        #     else:
        #         job_dir = direcs[code + '_path'] / joint_name
                
           