from os import environ
from pathlib import Path
from typing import Optional, Tuple
from rich.console import Console
from ase import Atoms
from ase.io import read,write
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from ase.optimize import FIRE, BFGS
from ase.neb import NEB, NEBTools

from herculestools.dft import (
    RunConfiguration,
    create_Vasp_calc,
    set_magnetic_moments,
    create_neb_path
)

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

def main(db_id: Optional[int] = None,
         N_images: int=3, climb: bool=True, fmax: float=0.03, parallel: bool= False,
         fire: bool=True, vasp: dict={}, **kwargs) -> Tuple[bool, Optional[dict]]:

    """
    Perform the main NEB (Nudged Elastic Band) calculation.

    Args:
        db_id (Optional[int]): The ID of the NEB in the database. Default is None.
        N_images (int): The number of images in the NEB calculation. Default is 3.
        climb (bool): Whether to use the climbing image NEB method. Default is False.
        fmax (float): The maximum force tolerance for the optimization. Default is 0.05.
        parallel (bool): Whether to run the calculation in parallel. Default is False.
        fire (bool): Whether to use the FIRE optimizer. Default is False.
        vasp (dict): Additional parameters for the VASP calculator. Default is an empty dictionary.

    Returns:
        Tuple[bool, Optional[dict]]: A tuple containing a boolean indicating whether the calculation was successful,
        and an optional dictionary with the trajectory file path and NEB ID.
    """

    # Connect to the database
    db_path = RunConfiguration.get().structures_dir / "hexag_perovs_strained.db"
    db: SQLite3Database
    db = connect(db_path)

    iID = kwargs['initial_id']
    fID = kwargs['final_id']
    vi = kwargs['initial_vac']
    vf = kwargs['final_vac']

    # Get the name of the structure from the database. 
    # If no initial or final ID is given, then the db_id is used.
    if iID == 0 or fID == 0:
        #print(db_id)
        db_name = db.get(db_id).name
        # Keep the name and the dopant position, and remove the type of job.
        names_list = db_name.split('_')
        name = '_'.join(names_list[:-1])
        #print(name)
        iID = db.get(selection=f"name={name}_vi").id
        fID = db.get(selection=f"name={name}_vf").id

    else:
        db_name = db.get(iID).name
        names_list = db_name.split('_')
        name = '_'.join(names_list[:-1])
        #print(name)

    neb_dir = Path(RunConfiguration.get().home / 'NEB' / f"{name}")
    neb_dir.mkdir(parents=True, exist_ok=True)
    
    # Get the initial and final structures from the converged results. XXX: This is a temporary solution.
    initial_row = db.get(iID)
    initial_dir = initial_row.dir
    initial= read(f"{initial_dir}/vasprun.xml", index=-1)
    #write(f"{neb_dir}/{name}_initial.traj", initial)
    #print(initial.get_atomic_numbers())
    final_dir: Atoms = db.get(fID).dir
    final = read(f"{final_dir}/vasprun.xml", index=-1)
    #write(f"{neb_dir}/{name}_final.traj", final)
    #print(final.get_atomic_numbers())

    # Order the endpoints correctly
    initial.append(initial.pop(vf-1))
    #print(initial.get_atomic_numbers())
    final.append(final.pop(vi))
    #print(final.get_atomic_numbers())
    
    # Create images and master directory
    neb = create_neb_path(initial, final, N_images, climb=climb, parallel = parallel)
    neb_dir = Path(RunConfiguration.get().home / 'NEB' / f"{name}")
    neb_dir.mkdir(parents=True, exist_ok=True)

    # The trajectory file must be written only once.
    start_traj = neb_dir / f"{name}_start.traj"
    if not start_traj.exists():
        write(f"{neb_dir}/{name}_start.traj", neb.images)
    else:
        pass

    # If parallel is True, then allow_shared_calculator = False
    if parallel:
        shared_calc = False
    else:
        shared_calc = True

    # Restart a calculation. For the case when new files are created, 
    # a timestamp will enable to read the most recent of these files.
    # Optimizers don't allow append mode, so this is a work-around
    path = ".traj"  
    for i in neb_dir.glob(f"*{path}"):
        if i.stat().st_size == 0:
           i.unlink()

    trajs = [i for i in neb_dir.glob(f"*{path}") if i.is_file() and i.stat().st_size > 0]
    counter = len(trajs)
    if counter > 0:
        traj = max(trajs, key= lambda a: a.stat().st_mtime)
        traj_old = neb_dir / traj.name
        neb = NEB(read(str(traj_old) + f"@-{N_images+2}:"), climb=climb, parallel=parallel,
                  method="improvedtangent", allow_shared_calculator=shared_calc)
        initial = neb.images[0]
        final = neb.images[-1]

    counter += 1
    traj_file = neb_dir / f"{neb_dir.name}_{counter}{path}"

    if counter > 10:
        fmax = fmax + 0.005
    
    print(f"Starting NEB calculation for {name}")
    
    # Single point calculations for endpoints
    i_dir = neb_dir / 'init'
    initial = setup_run(initial, i_dir, vasp)
    initial.get_potential_energy()

    f_dir = neb_dir / 'final'
    final = setup_run(final, f_dir, vasp)
    final.get_potential_energy()

    # Create calculator for internal images
    for i, img in enumerate(neb.images[1:-1]):
        direc = neb_dir / f"iimg{i+1:02}"
        img = setup_run(img, direc, vasp)
        img.get_potential_energy()
        
    # Run path relaxation
    # Add a boolean for FIRE, and call it from the function's variables.
    if fire:
        opt = FIRE(neb, trajectory=str(traj_file), logfile=f"{neb_dir}/{neb_dir.stem}_{counter}.log", 
                   # Suggested by chatGPT:
                   dt=0.05, 
                   maxstep=0.1, 
                   dtmax=0.2,
                   Nmin=10,
                   finc=1.05,
                   fdec=0.4,
                   astart=0.1, 
                   fa=0.98, 
                   a=0.1, 
                   downhill_check=True, force_consistent=False)
    else:
        opt = BFGS(neb, trajectory=str(traj_file), logfile=f"{neb_dir}/{neb_dir.stem}_{counter}.log")
    opt.run(fmax=fmax)

    # Once the calculation has completed:
    converged_traj = read(f"{Path(traj_file)}@-{N_images+2}:")
    tool = NEBTools(converged_traj)
    Ef, dE = tool.get_barrier()
    Er = Ef-dE
    Emax = max(Ef, Er)
    
    # Here, I will add a check to see if the energy barrier given is actually a maximum. If not, the next few lines will correct it.
    energies = [trajectory.get_potential_energy() for trajectory in converged_traj]
    neg_barrier = False
    if energies[0] != min(energies) or energies[-1] != min(energies):
        neg_barrier = True
        print(f"The energy barrier at {direc.name} is not a maximum. Correcting it.")
        min_idx = energies.index(min(energies))
        fwd_barrier = abs(energies[0] - energies[min_idx])
        bwd_barrier = abs(energies[-1] - energies[min_idx])    
        max_barrier = max(energies) - energies[min_idx]
    
        if fwd_barrier > Ef:
            Ef = fwd_barrier
        if bwd_barrier > Er:
            Er = bwd_barrier
        if max_barrier > Emax:
            Emax = max_barrier

    db_id = update_or_write(db, initial, name=f"{name}_neb", dir = str(neb_dir), forward_e = Ef, delta_e = dE, reverse_e = Er, barrier = Emax, neg_barrier = neg_barrier)

    print(f"The energy barrier is :{Emax:.3f}")
    return True, {'trajectory': str(traj_file.relative_to(neb_dir)), 'db_id': int(db_id)}

def setup_run(atoms: Atoms, direc: Path, vasp: dict = {}) -> Atoms:
    atoms.calc = create_Vasp_calc(atoms, 'PBEsol', direc)
    set_magnetic_moments(atoms)
    atoms.calc.set(**vasp,
                   ibrion=-1,
                   nsw=1,
                   kpar=nnodes,
                   ncore=ncore)
    return atoms

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)