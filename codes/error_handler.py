from perqueue.queue import PersistentQueue
from perqueue.selection import Selection
from pathlib import Path
from ase.db import connect
from typing import Optional

calc_failed = "CalculationFailed"
time_out = "DUE TO TIME LIMIT"
error_string = "I REFUSE TO CONTINUE WITH THIS SICK JOB"
matching_str = "EEEEEEE  R     R" # String from the bottom row of the ERROR message.
common_errors = {
    "eddav": "Call to ZHEGV failed",
    "zbrent" : "ZBRENT: fatal error in bracketing",
    "bravais": "Inconsistent Bravais lattice" ,
    "fexcf": "ERROR FEXCF: supplied exchange-correlation",
    "ibrion0" : "Fatal error! IBRION=0",
    "wavecar" : "ERROR: while reading WAVECAR",
    "lapack" : "LAPACK: Routine ZPOTRF failed",
    "subrot": "ERROR in subspace rotation",
}

# Common errors and their solutions by modifying the INCAR file
common_solutions = {
    "eddav":    {'algo': 'VeryFast'},
    "zbrent":   {'potim': 0.25,
                'addgrid' : True},
    "bravais":  {'symprec': 1e-03,
                'algo': 'VeryFast'},
    "fexcf":    {'potim': 0.25},
    "ibrion0":  {'ibrion' : 2},
    "wavecar":  {'algo': 'VeryFast',
                'istart': 1},
    "lapack":   {'algo': 'VeryFast',
                'potim': 0.25},
    "subrot":   {'algo': 'Normal'}
}

cwd = Path.cwd()
db = connect('structures/hexag_perovs_strained.db')

with PersistentQueue() as pq:
    entries = pq.get_entries()

s = Selection(states='f')
targets = s.filter(entries)
codes = [pq.get_code(en.key) for en in targets]
codes = list(set(codes))
print(f"Codes in the queue are {codes}")

# PQ assigns the newest MQ id to its entry while still keeping the last state before the run was completed.
errors_dict = dict()
timeout_jobs = []

# Function to read the entry data and extract the relevant information from the ASE database
def extract_data(entry) -> Optional[str]:
    data = entry.data
    try:
        if 'initial_id' in data:
            initial_id = data['initial_id']
            initial = db.get(initial_id)
            initial_name = initial.name[:-3]
            return initial_name
        elif 'final_id' in data:
            final_id = data['final_id']
            final = db.get(final_id)
            final_name = final.name[:-3]
            return final_name
        elif 'neb_id' in data:
            neb_id = data['neb_id']
            neb = db.get(neb_id)
            neb_name = neb.name[:-4]
            return neb_name
    except TypeError:
        return None

# Get the timed out and calculation failed jobs
for en in targets:
    mq_id = en.mq_id
    pq_key = en.key
    code = pq.get_code(pq_key)
    status = en.state.serialize()
    data = en.data
    # Getting the name of the last successful entry
    base_name = extract_data(en)

    # Use the mq_id to get the directory of the calculation.
    error_file = Path(f"perqueue.runner.{mq_id}.err")
    backtrace = "Backtrace for this error:"
    
    if error_file.is_file() and error_file.stat().st_size > 0:
        readlines = error_file.read_text().split('\n')
        print(f"Checking {error_file} for entries {en.key} with name {base_name} and status {status}")
        if (time_out_line := [line for line in readlines if time_out in line]):
            timeout_jobs.append(pq_key)
        elif (backtrace in readlines):
            timeout_jobs.append(pq_key)
        elif (calc_failed_dir := [line.split(' ') for line in readlines if calc_failed in line][-1][-5]):
            calc_failed_dir = Path(calc_failed_dir)
            vaspout = calc_failed_dir / "vasp.out"
            if vaspout.is_file() and vaspout.stat().st_size > 0:
                vaspout_lines = vaspout.read_text().split('\n')
                next_lines = [line.strip().strip('|').strip(' ') for idx, line in enumerate(vaspout_lines) if matching_str in line for line in vaspout_lines[idx+2:-5]]
                error_msg = " ".join(next_lines)
                errors_dict[pq_key] = error_msg
        else:
            pass

# Resubmit the jobs that timed out
with PersistentQueue() as pq:
    print(f"Resubmitting {len(timeout_jobs)} jobs because of time out")
    for job in timeout_jobs:
        en = pq.get_entry(job)
        pq.resubmit(en)

# Now match the error messages with the common errors and resubmit the jobs
    print(f"Resubmitting {len(errors_dict)} jobs because of calculation failed")
    for pq_key, error_msg in errors_dict.items():
        en = pq.get_entry(pq_key)
        new_args = pq.get_args(pq_key) 
        for key, value in common_errors.items():
            if value in error_msg:
                new_args['vasp'] = common_solutions[key]
                print(f"Resubmitting {pq_key} with new args {new_args}")
                # Save the new arguments
                pq.save_args(pq_key, new_args, False)
                # Resubmit the job
                pq.resubmit(en)