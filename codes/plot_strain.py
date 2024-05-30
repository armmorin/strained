from ase.db import connect
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

working_dir = Path.cwd()

db = connect('structures/hexag_perovs_strained.db')

# Get the entries with the sys_id value
sys_ids = []
for entry in db.select(selection='sys_id'):
    sys_id = entry.sys_id
    db_id = entry.id
    sys_ids.append((sys_id, db_id))

# Create a dictionary with the sys_id as key and the db_ids as values
ids_dict = {}
for sys_id, db_id in sys_ids:
    if sys_id in ids_dict:
        ids_dict[sys_id].append(db_id)
    else:
        ids_dict[sys_id] = [db_id]

# Now we can loop over the sys_ids and get a fit curve for each
df = pd.DataFrame(columns=['name','ip_strain', 'op_strain', 'energy'])
values = {}
for sys_id in ids_dict.items():
    db_ids = sys_id[1]
    for db_id in db_ids:
        entry = db.get(id=db_id)
        name_components = entry.name.split('_')
        name = '_'.join(name_components[0:2])
        ip_strain = entry.in_plane
        op_strain = entry.out_of_plane
        energy = entry.energy
        # Store the values in a dictionary with the name as key
        if name in values:
            values[name].append((ip_strain, op_strain, energy))
        else:
            values[name] = [(ip_strain, op_strain, energy)]
        
        # Store the values in a DataFrame
        df = df._append({'name': name, 'ip_strain': ip_strain, 'op_strain': op_strain, 'energy': energy}, ignore_index=True)
        df.sort_values(by=['name', 'ip_strain', 'op_strain'], inplace=True)
        df.reset_index(drop=True, inplace=True)

    #break

plots_dir = working_dir / 'plots'
plots_dir.mkdir(parents=True, exist_ok=True)

# Fit a curve to the data
for name in df['name'].unique():
    min_set = []
    for ip_strain in df['ip_strain'].unique():
        values = []
        for index, row in df.iterrows():
            if row['name'] == name and row['ip_strain'] == ip_strain:
                values.append([row['op_strain'], row['energy']])
    
        values = np.array(values)
        x = values[:,0]
        y = values[:,1]
        poly = np.polyfit(x, y, 2)
        op_min = -poly[1]/(2*poly[0])
        en_min = np.polyval(poly, op_min)
        
        # Plot the data
        plt.plot(x, y, 'o')
        x_fit = np.linspace(x[0], x[-1], 100)
        y_fit = np.polyval(poly, x_fit)
        plt.plot(x_fit, y_fit, label='Polynomial fit')
        plt.xlabel('Out-of-plane strain')
        plt.ylabel('Energy')
        plt.legend()
        plt.savefig(f'{plots_dir}/{name}_{ip_strain}_fit.png')
        plt.close()  
        # Store the minimum value of the polynomial together with the in-plane strain
        min_set.append((ip_strain, op_min))
        df = df._append({'name': name, 'ip_strain': ip_strain, 'op_strain': f"{op_min:.2f}", 'energy': f"{en_min:.6f}"}, ignore_index=True)
    
    df.sort_values(by=['name', 'ip_strain', 'energy'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    # Plot the minimum values
    x = np.array(min_set)[:,0]
    y = np.array(min_set)[:,1]
    plt.title(f'{name} minimum energy')
    plt.plot(x, y, 'o-')
    plt.xlabel('In-plane strain')
    plt.ylabel('Out-of-plane strain')
    plt.savefig(f'{plots_dir}/{name}_strain.png')
    plt.close()

print(df)
df.to_csv(working_dir / 'fit_values.csv', index=False)