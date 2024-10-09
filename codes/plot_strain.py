from ase.db import connect
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from perqueue import PersistentQueue
from perqueue.selection import Selection

def get_lowest_energy(initial_id: int, final_id: int):
    """ Get the lowest energy id from two entries in the database."""
    db = connect('structures/hexag_perovs_strained.db')
    with db:
        initial = db.get(initial_id)
        final = db.get(final_id)
        if initial.energy < final.energy:
            return initial_id
        else:
            return final_id

working_dir = Path.cwd()

db = connect('structures/hexag_perovs_strained.db')

df = pd.DataFrame(columns=['name', 'mask', 'x_strain', 'y-strain', 'z_strain', 'energy', 'barrier'])

# Get the reference system information from the database
ref_dict = {}
with db:
    ref = db.get(id=1) # The original reference should already have the NEB step done. It is running ATM.
    ref_atoms = ref.toatoms()
    ref_dict['ref_volume'] = ref.volume
    ref_dict['ref_a_lattice'] = ref_atoms.get_cell()[0,0]
    ref_dict['ref_b_lattice'] = ref_atoms.get_cell()[1,1]
    ref_dict['ref_c_lattice'] = ref_atoms.get_cell()[2,2]
    
db_ids = []

s = Selection(states='srq', names='neb*')
with PersistentQueue() as pq:
    entries = pq.get_entries()

targets = s.filter(entries)

for entry in targets:
    key = entry.key
    #print(f"Processing entry {key}")
    try:
        #entry.data['neb_id']
        neb_id = entry.data['neb_id']
        db_ids.append(neb_id)
    except:
        entry = pq.get_entry(key-1)
        initial_id = entry.data['initial_id']
        final_id = entry.data['final_id']
        neb_id = get_lowest_energy(initial_id, final_id)
        db_ids.append(neb_id)
        #print(f"NEB step not done for entry {entry.name}")

# Get the data for each of the completed entries and compare it against the reference values.        
for id in db_ids:
    try:
        db_id = db.get(id)
    except KeyError:
        continue
    mask = db_id.mask
    energy = db_id.energy
    volume = db_id.volume
    d_vol = (volume - ref_dict['ref_volume']) / ref_dict['ref_volume'] * 100
    a_lattice = db_id.toatoms().get_cell()[0,0]
    d_a = (a_lattice - ref_dict['ref_a_lattice']) / ref_dict['ref_a_lattice'] * 100
    b_lattice = db_id.toatoms().get_cell()[1,1]
    d_b = (b_lattice - ref_dict['ref_b_lattice']) / ref_dict['ref_b_lattice'] * 100
    c_lattice = db_id.toatoms().get_cell()[2,2]
    d_c = (c_lattice - ref_dict['ref_c_lattice']) / ref_dict['ref_c_lattice'] * 100
    if hasattr(db_id, 'barrier'):
        barrier = db_id.barrier
    else:
        barrier = np.nan
    df = pd.concat([df, pd.DataFrame({'name': db_id.name, 'mask': mask, 'x_strain': d_a, 'y_strain': d_b, 'z_strain': d_c, 'energy': energy, 'barrier': barrier}, index=[0])], ignore_index=True)
 
df.to_csv(working_dir / 'strained_values.csv', index=False)          

axis_dict = {
    0: 'x',
    1: 'y',
    2: 'z'
}
# Now we can loop over the different masks and plot different values.
masks = df['mask'].unique()
for mask in masks:
    mask_df = df[df['mask'] == mask]
    # Plot all values that have the same mask
    x_strain = mask_df['x_strain']
    y_strain = mask_df['y_strain']
    z_strain = mask_df['z_strain']
    energy = mask_df['energy']
    barrier = mask_df['barrier']
    
    #y = energy
    y  = barrier
        
    # Plot all three axes of strain against energy in the same plot
    fig, ax = plt.subplots(1, 3, figsize=(15,5))
    for i, strain in enumerate([x_strain, y_strain, z_strain]):
        ax[i].plot(strain, y, 'o', label=f'{axis_dict[i]} strain')
        ax[i].set_xlabel(f'{mask} {axis_dict[i]} strain (%)')
        x_fit = np.linspace(strain.min(), strain.max(), 100)
        y_fit = np.polyval(np.polyfit(strain, y, 2), x_fit)
        ax[i].plot(x_fit, y_fit, label=f'{axis_dict[i]} fit')
        ax[i].grid()
        ax[i].legend()
        
    if y.name == 'energy':
        ax[0].set_ylabel('Total Energy (eV)')
    else:
        ax[0].set_ylabel('Migration Barrier (eV)')
    
    plt.tight_layout()
    
    if y.name == 'energy':
        plt.savefig(working_dir / f'{mask}_energy.png')
    else:
        plt.savefig(working_dir / f'{mask}_barrier.png')
    plt.close()