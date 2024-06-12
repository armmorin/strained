from ase.db import connect
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from perqueue import PersistentQueue

working_dir = Path.cwd()

db = connect('structures/hexag_perovs_strained.db')

db_ids = []
with PersistentQueue() as pq:
    entries = pq.get_entries()
        
for entry in db.select(selection='db_id'):
    db_id = entry.id

# Get the entries with the db_id
    sys_id = pq.get_args()
    db_ids.append((db_id))

# Create a dictionary with the system_id as key and the db_ids as values
ids_dict = {}
for db_id in db_ids:
    
# Now we can loop over the sys_ids and get a fit curve for each
df = pd.DataFrame(columns=['name','ip_strain','op_strain', 'energy'])
values_dict = {}
for sys_id in ids_dict.items():
    db_ids = sys_id[1]
    out_of_plane = {}
    energies = {}
    for db_id in db_ids:
        entry = db.get(id=db_id)
        name_components = entry.name.split('_')
        name = '_'.join(name_components[0:2])
        ip_strain = entry.in_plane
        energies[ip_strain] = entry.energy
        out_of_plane[ip_strain] = entry.toatoms().get_cell()[2,2]
    
        #out_of_plane.append((ip_strain,entry.toatoms().get_cell()[2,2]))
    print(out_of_plane)
    print(energies)
    
    # Use the out_of_plane value when ip_strain is 0 as a reference
    ref_op_strain = [op_strain for ip_strain, op_strain in out_of_plane.items() if ip_strain == 0]
    #ref_op_strain = [op_strain for ip_strain, op_strain in out_of_plane if ip_strain == 0][0]
    # Get the out-of-plane strain by comparing the strained and unstrained lattice parameters
    for ip_strain, op_strain in out_of_plane:
        op_strain = (op_strain - ref_op_strain) / ref_op_strain * 100
        print(f'{name} {ip_strain} {op_strain} {energies[ip_strain]}')
        # Store the values in a dictionary with the name as key
        #if name in values_dict:
        #    values_dict[name].append((ip_strain, op_strain, energy))
        #else:
        #    values_dict[name] = [(ip_strain, op_strain, energy)]
        
        # Store the values in a DataFrame
        df = df._append({'name': name, 'ip_strain': ip_strain, 'op_strain': op_strain, 'energy': f"{energy:.4f}"}, ignore_index=True)
        df.sort_values(by=['name', 'ip_strain', 'op_strain'], inplace=True)
        df.reset_index(drop=True, inplace=True)

    #break

plots_dir = working_dir / 'plots'
plots_dir.mkdir(parents=True, exist_ok=True)

# Fit a curve to the data
for name in df['name'].unique():
    for ip_strain in df['ip_strain'].unique():
        values = []
        for index, row in df.iterrows():
            if row['name'] == name and row['ip_strain'] == ip_strain:
                values.append([row['op_strain'], row['energy']])
    
        values = np.array(values)
        x = values[:,0].astype(float)
        y = values[:,1].astype(float)
        poly = np.polyfit(x, y, 2)
        en_min = np.min(np.polyval(poly, x))
        op_min = x[np.argmin(np.polyval(poly, x))]
        df = df._append({'name': name, 'ip_strain': ip_strain, 'op_strain': f"{op_min:.2f}", 'energy': f"{en_min:.4f}"}, ignore_index=True)          
        df.sort_values(by=['name', 'ip_strain', 'energy'], inplace=True)
        df.reset_index(drop=True, inplace=True)
        
        # Plot the fit
        plt.plot(x, y, 'o', label = f'ip={ip_strain}')
        x_fit = np.linspace(x[0], x[-1], 100)
        y_fit = np.polyval(poly, x_fit) 
        plt.plot(x_fit, y_fit, label=f'ip={ip_strain} fit')
    
    plt.title(f'{name} fit')
    plt.xlabel('Out-of-plane strain')
    plt.ylabel('Energy') 
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{plots_dir}/{name}_fit.png')
    
    plt.close()  

df.to_csv(working_dir / 'fit_values.csv', index=False)