from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import  MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments  import LightStructureEnvironments
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import json
from perqueue.queue import PersistentQueue
from perqueue.selection import Selection
from ase.db import connect
import pandas as pd
from ase.neighborlist import NeighborList
from ase import data
from ase.neighborlist import PrimitiveNeighborList
# from typing import Optional, List
# from ase.atoms import Atoms
from ase.visualize import view

# Function to identify the closest atom to a given set of coordinates
def get_most_similar_atom(atoms, coords, symbol):
    diff_positions = atoms.positions - coords
    distances = np.linalg.norm(diff_positions, axis=1)
    # print the index, distances and the symbol of the atom as a tuple
    symbol_distances = [(i, d, atoms.symbols[i]) for i, d in enumerate(distances)]
    # Remove the rows that contain an 'O' or the same symbol as the first atom in the structure
    symbol_distances = [sd for sd in symbol_distances if sd[2] == symbol]
    # Sort the distances and get the index of the closest atom
    symbol_distances.sort(key=lambda x: x[1])
    # Get the index of the closest atom
    return symbol_distances[0][0], symbol_distances[0][1]

# Defining the function to plot the coordination environments
def plot_coord_envs(coord_environments: dict, name: str) -> None:
    steps = range(5)
    cations = ['donor', 'acceptor', 'moving O']
    colors = ['Blues_r', 'Reds_r', 'Greens_r']
    markers = ['s', 'D', 'o']
    color_map = {'donor': 'blue', 'acceptor': 'red', 'moving O': 'green'}

    fig, ax = plt.subplots(figsize=(10, 6))

    handles = [plt.Line2D([], [], marker='o', color='w', label=cation,
                    markerfacecolor=color_map[cation], markersize=16)
        for cation in cations]
    
    for cation, color in zip(cations, colors):
        ce_numbers = [coord_environments[f'img{i}'][cation]['ce number'] for i in steps]
        csm = [coord_environments[f'img{i}'][cation]['csm'] for i in steps]
        
        scatter = ax.scatter(steps, ce_numbers, c=csm, cmap=plt.get_cmap(color), 
                            marker=markers[cations.index(cation)], edgecolor='k',
                            label=cation, s=250, alpha=0.75)
        # Annotate the ce symbol for each cation
        for i, txt in enumerate(ce_numbers):
            if cation == 'donor':
                ax.annotate(coord_environments[f'img{i}'][cation]['ce symbol'], (steps[i], ce_numbers[i]), fontsize=12,
                            xytext=(-15, 10), textcoords='offset points')
            elif cation == 'acceptor':
                ax.annotate(coord_environments[f'img{i}'][cation]['ce symbol'], (steps[i], ce_numbers[i]), fontsize=12,
                            xytext=(-20, -15), textcoords='offset points')
            else:
                ax.annotate(coord_environments[f'img{i}'][cation]['ce symbol'], (steps[i], ce_numbers[i]), fontsize=12,
                            xytext=(5, 10), textcoords='offset points')
            
        if cation == 'moving O':
            ax.plot(steps, ce_numbers, color=color_map[cation], 
                linestyle='--', alpha=0.25)
        else:
            ax.plot(steps, ce_numbers, color=color_map[cation], linewidth=2,
                linestyle='-', alpha=0.5)
        
        # Add a colorbar for each cation with its name
        cbar = plt.colorbar(scatter, ax=ax, pad=0.01, shrink=0.4, label='')
        # Add a label to the colorbar
        cbar.ax.set_ylabel(f'{cation} CSM', rotation=270, labelpad=-50)
        cbar.ax.yaxis.set_label_position('left')
        cbar.ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    
    ax.grid(True)
    # Add a spacing for the legend 
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1, 1), labelspacing=1.2)
    ax.set_xlabel('NEB Step')
    ax.set_ylabel('Coordination Number')
    ax.set_title(f'Coordination Environments for {name}')

    plt.tight_layout()
    plt.savefig(f'neb_geometries/{name}/{name}_coord_envs.pdf')
    #plt.show()
    plt.close()
    return None

# Function to flatten a dictionary
def flatten_dict(d, parent_key='', sep=' '):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


# Define the path to the database

db_path = 'structures/hexag_perovs_strained.db'

# Create a directory to store the coordination environment plots

geoms_dir = Path('neb_geometries')
geoms_dir.mkdir(exist_ok=True, parents=True)

# Get the entries from the queue

with PersistentQueue() as pq:
    entries = pq.get_entries()

selection = Selection(names='*CINEB*', states='s')
targets = selection.filter(entries)

# Generate a list of the database entries from the selection
db_entries = []

for entry in targets:
    db_id = entry.data['neb_id']
    traj = entry.data['trajectory']
    db_entries.append((db_id, traj))

# The following variables can be extracted in an automated way.
a_site = 'Ba' # The A site atom
oop_o_idx = 126 # The index of the out-of-the-plane migrating oxygen atom
ionic_r_oxy = 1.4

# Create the dataframe to store the coordination environment information
coord_env_df = pd.DataFrame()

strain_dict = {'s3': 1.030, 's2': 1.025, 's1': 1.015, 's0': 1.0075, 'e0': 1.000, 'c0':0.9925, 'c1': 0.985, 'c3': 0.970}

# Loop over the database entries and extract the geometries
db = connect(db_path)
for db_id, traj in db_entries:
    row = db.get(id=db_id)
    name = row.name
    print(f'Processing {db_id}: {name}')
    try:
        MASK = row.mask
    except AttributeError:
        rows = db.select(name=name)
        for row in rows:
            if hasattr(row, 'mask'):
                MASK = row.mask
                break
    neb_geoms = geoms_dir / name
    neb_geoms.mkdir(exist_ok=True, parents=True)
    neb_dir = Path(row.dir)
    traj_path = f'{neb_dir}/{traj}'
    neb_traj = read(traj_path, index='-5:')
    
    # Translate all images to center the migrating oxygen atom in the z direction
    for neb_step in neb_traj:
        cell = neb_step.cell
        neb_step.translate([0, 0, cell[2,2]/2])
        neb_step.wrap()

    name_parts = name.split('_')
    basename = name_parts[0]
    position = name_parts[1]
    strain = row.in_plane
    #strain = strain_dict[applied_strain]


    # Read the TS atoms. I would extract the TS atoms from the stored atoms element in the DB, but they should not be too different from the saddle point image.
    ts_atoms = neb_traj[2]

    # Get the TS cation indices from the closest atoms to the moving oxygen atom
    #ts_cation_ids = []
    
    # Get the nearest atoms to the moving oxygen atom
    nl = NeighborList([ionic_r_oxy]*len(ts_atoms), self_interaction=False, bothways=True, primitive=PrimitiveNeighborList)
    nl.update(ts_atoms)
    nearest_indices, _ = nl.get_neighbors(oop_o_idx)
    # Remove the atoms that belong to the A site or are oxygen atoms
    ts_cation_ids = [i for i in nearest_indices if ts_atoms.symbols[i] != a_site and ts_atoms.symbols[i] != 'O']


    moving_o_dict = {}
    b1_site_dict = {}
    b2_site_dict = {}

    # Set up the local geometry finder
    lgf = LocalGeometryFinder()
    lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True, structure_refinement="none")

    coord_environments = {}
    coord_environments['name'] = name
    coord_environments['position'] = position
    coord_environments['mask'] = MASK
    coord_environments['strain'] = strain
    coord_environments['delta_e'] = row.delta_e
    coord_environments['forward_e'] = row.forward_e
    coord_environments['reverse_e'] = row.reverse_e
    coord_environments['barrier'] = row.barrier
    
    ### START OF THE CSM SUBROUTINE:
    # Read the structures.
    for i, neb_step in enumerate(neb_traj):
        cell = neb_step.cell
  
        # Create a new neb_step object to eliminate the a site atoms
        neb_step = neb_step[[atom.index for atom in neb_step if atom.symbol != a_site]] #type: ignore

        # Atoms of interest will be reindexed once the a site atoms are removed. 
        # Use the find most similar atom function to get the indices of the moving oxygen and the b sites.
        moving_o_idx, _ = get_most_similar_atom(neb_step, ts_atoms.positions[oop_o_idx], 'O')
        b1_site, _ = get_most_similar_atom(neb_step, ts_atoms.positions[ts_cation_ids[0]], ts_atoms.symbols[ts_cation_ids[0]])
        b2_site, _ = get_most_similar_atom(neb_step, ts_atoms.positions[ts_cation_ids[1]], ts_atoms.symbols[ts_cation_ids[1]])

        # Add a check to see if the b1 and b2 sites are identical
        if b1_site == b2_site:
            print(f"b1 and b2 sites are identical in {name} for image {i}")
            break

        structure = AseAtomsAdaptor.get_structure(neb_step)
        lgf.setup_structure(structure)

        # Get the structure environments only for the specified indices.
        # This selection is done to avoid the calculation of the coordination environments for all atoms in the structure.
        se = lgf.compute_structure_environments(
            only_indices= [
                        b1_site, 
                        b2_site,
                        moving_o_idx
                    ],
            maximum_distance_factor= ionic_r_oxy
        )
        
        #strategy= SimplestChemenvStrategy(distance_cutoff=1.4,angle_cutoff=0.3)
        strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()

        lse = LightStructureEnvironments.from_structure_environments(strategy=strategy,structure_environments=se)
        
        b1_site_dict[i] = lse.coordination_environments[b1_site][0]
        b2_site_dict[i] = lse.coordination_environments[b2_site][0]
        moving_o_dict[i] = lse.coordination_environments[moving_o_idx][0]
        
        # generate dictionaries for each of the coordination environment properties.
        coord_environments[f'img{i}'] = {
            # The expected output for each of the dictionaries is:
            # b1 site coordination: {'ce_symbol': 'T:4', 'ce_fraction': 0.9784470107138163, 'csm': 0.046576132532119084, 'permutation': [0, 1, 2, 3]}
            'acceptor': { 
                'ce symbol': b1_site_dict[i]['ce_symbol'].split(':')[0],
                'ce number': int(b1_site_dict[i]['ce_symbol'].split(':')[1]),
                'ce fraction': b1_site_dict[i]['ce_fraction'],
                'csm': b1_site_dict[i]['csm']
                },
            'donor': {
                'ce symbol': b2_site_dict[i]['ce_symbol'].split(':')[0],
                'ce number': int(b2_site_dict[i]['ce_symbol'].split(':')[1]),
                'ce fraction': b2_site_dict[i]['ce_fraction'],
                'csm': b2_site_dict[i]['csm']
                },
            'moving O': {
                'ce symbol': moving_o_dict[i]['ce_symbol'].split(':')[0],
                'ce number': int(moving_o_dict[i]['ce_symbol'].split(':')[1]),
                'ce fraction': moving_o_dict[i]['ce_fraction'],
                'csm': moving_o_dict[i]['csm']
                }
        }

    # Plot the coordination environments
    plot_coord_envs(coord_environments, name)

    ## END OF THE CSM SUBROUTINE
        
    # Flatten the dictionaries to store the coordination environment information in a dataframe
    coord_environments = flatten_dict(coord_environments)
    coord_env_df = pd.concat([coord_env_df, pd.DataFrame(coord_environments, index=[name])], ignore_index=True)

    #break
pd_to_csv = coord_env_df.to_csv('coord_envs.csv')

