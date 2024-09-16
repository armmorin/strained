from os import environ
from pathlib import Path
from typing import Optional, Tuple
from rich.console import Console
from ase import Atoms
from ase.io import read, write
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from ase.optimize import FIRE, BFGS
from ase.neb import NEB, NEBTools
from ase.io.trajectory import Trajectory
import shutil
from herculestools.dft import (
    RunConfiguration,
    create_Vasp_calc,
    set_magnetic_moments,
    create_neb_path
)

here = Path(__file__).resolve().parent.parent

c = Console()
nnodes = int(environ['SLURM_NNODES'])
ncore = int(environ['NCORE'])

def main(db_id: Optional[int] = None,
         N_images: int = 3, climb: bool = False, fmax: float = 0.05, parallel: bool = False,
         fire: bool = True, vasp: dict = {}, **kwargs) -> Tuple[bool, Optional[dict]]:
    """
    Perform the main NEB (Nudged Elastic Band) calculation.

    Args:
        db_id (Optional[int]): The ID of the NEB in the database. Default is None.
        N_images (int): The number of images in the NEB calculation. Default is 3.
        climb (bool): Whether to use the climbing image NEB method. Default is True.
        fmax (float): The maximum force tolerance for the optimization. Default is 0.05.
        parallel (bool): Whether to run the calculation in parallel. Default is False.
        fire (bool): Whether to use the FIRE optimizer. Default is True.
        vasp (dict): Additional parameters for the VASP calculator. Default is an empty dictionary.

    Returns:
        Tuple[bool, Optional[dict]]: A tuple containing a boolean indicating whether the calculation was successful,
        and an optional dictionary with the trajectory file path and NEB ID.
    """

    # Set parallel execution
    shared_calc = not parallel

    # Database connection
    db_path = RunConfiguration.structures_dir / "hexag_perovs_strained.db"
    db: SQLite3Database = connect(db_path)

    iID = kwargs.get('initial_id', 0)
    fID = kwargs.get('final_id', 0)
    vi = kwargs.get('initial_vac',127)
    vf = kwargs.get('final_vac',126)

    # Get the name of the structure from the database.
    if iID == 0 or fID == 0:
        c.log(f"No initial or final ID provided. Using db_id: {db_id}")
        db_name = db.get(db_id).name
        names_list = db_name.split('_')
        name = '_'.join(names_list[:-1])
        mask = names_list[2]
        iID = db.get(selection=f"name={name}_vi").id
        fID = db.get(selection=f"name={name}_vf").id
    else:
        db_name = db.get(iID).name
        names_list = db_name.split('_')
        name = '_'.join(names_list[:-1])
        mask = names_list[2]

    neb_dir = Path(RunConfiguration.home / 'NEB' / f"{name}")
    neb_dir.mkdir(parents=True, exist_ok=True)
        
    # The trajectory file must be written only once.
    trajs = []
    path = "traj"
    if not (start_traj := neb_dir / f"{name}_start.{path}").exists():
    
        initial_row = db.get(iID)
        print(iID, initial_row.formula, initial_row.natoms)
        initial: Atoms = initial_row.toatoms()
        final_row = db.get(fID)
        print(fID, initial_row.formula, initial_row.natoms)
        final: Atoms = final_row.toatoms()

        # Order the endpoints correctly
        print(initial.numbers[vf], initial.positions[vf])
        initial.append(initial.pop(vf))
        print(final.numbers[vi-1], final.positions[vi-1])
        final.append(final.pop(vi-1))

        # Create images and master directory
        neb = create_neb_path(initial, final, N_images, climb=climb, parallel=parallel)

        write(str(start_traj), neb.images)

    # Restart a calculation if possible
    else:
        trajs = [i for i in neb_dir.glob(f"*{path}") if i.is_file() and i.stat().st_size > 0]
        traj = max(trajs, key=lambda a: a.stat().st_mtime)
        # Use this trajectory to check for any newer WAVECAR files created after the trajectory's last modification
        traj_mtime = traj.stat().st_mtime
        for wavecar in neb_dir.rglob("WAVECAR"):
            if wavecar.stat().st_mtime > traj_mtime:
                c.log(f"Found a newer WAVECAR file in {wavecar.parent}.")
                wavecar.unlink()
            else:
                continue
        
        initial = read(traj.as_posix(), index=0)
        final = read(traj.as_posix(), index=-1)
        neb = NEB(read(f"{traj}@-{N_images+2}:"), climb=climb, parallel=parallel,
                method="improvedtangent", allow_shared_calculator=shared_calc)

    counter = len(trajs) + 1
    traj_file = neb_dir / f"{neb_dir.name}_{counter}.{path}"

    if counter > 10:
        fmax += 0.005

    c.log(f"Starting NEB calculation for {name}")

    # Single point calculations for endpoints
    initial = setup_run(initial, neb_dir / 'init', vasp)
    initial.get_potential_energy()

    final = setup_run(final, neb_dir / 'final', vasp)
    final.get_potential_energy()

    # Create calculator for internal images
    for i, img in enumerate(neb.images[1:-1]):
        img = setup_run(img, neb_dir / f"iimg{i+1:02}", vasp)
        img.get_potential_energy()

    # Run path relaxation with checker
    traj_log = neb_dir / f"{neb_dir.stem}.log"
    run_neb_with_checker(neb, traj_file, traj_log, fire, fmax)

    # Post-calculation processing
    
    #converged_traj =  write(f"{name}_converged.traj", neb.images)
    converged_traj = read(f"{Path(traj_file)}@-{N_images+2}:")
    #converged_traj = Trajectory(f"{name}_converged.traj",'w',[initial]+neb.images[1:-1]+[final])
    # Get the energies and find the transition state
    try:
        energies = [image.get_potential_energy() for image in converged_traj]
    except Exception as e:
        converged_traj = read(f"{Path(traj_file)}@-{N_images+2}:")
        print(f"Could not write converged trajectory. Error: {e}")
        energies = [image.get_potential_energy() for image in converged_traj]
    
    dE = energies[-1] - energies[0]
    # barrier = max(energies) - min(energies)

    # Move the important files from scratch to home
    for file in neb_dir.glob("*"):
        if file.is_file() and file.name != 'WAVECAR':
            # Take the current path and move it to the home directory
            file_parts = list(file.parts)
            new_path = file_parts[:2] + ['energy'] + file_parts[3:]
            new_file = Path(*new_path)
            new_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(file, new_file)
    
    if energies[0] == min(energies) or energies[-1] == min(energies):
        neg_barrier = 'no'
        Ef = max(energies) - energies[0]
        Er = max(energies) - energies[-1]
        barrier = max(Ef, Er)
        ts_atoms = converged_traj[energies.index(max(energies))]
    else:
        c.log(f"The energy barrier at {neb_dir.name} is not a maximum. Correcting it.")
        neg_barrier = 'yes'
        Ef = energies[0] - min(energies)
        Er = energies[-1] - min(energies)
        barrier = max(energies) - min(energies)
        ts_atoms = converged_traj[energies.index(min(energies))]

    # Emax = max(Ef, Er)
    nebtool = NEBTools(converged_traj)
    nebtool.plot_bands(label=str(neb_dir / neb_dir.stem))

    db_id = update_or_write(db, ts_atoms, name=f"{name}_neb", dir=str(neb_dir), mask=mask,
                            forward_e=Ef, delta_e=dE, reverse_e=Er, barrier=barrier, neg_barrier=neg_barrier)

    # Copy the database back to the home directory
    home_db = here / "structures/hexag_perovs_strained.db"
    shutil.copy(db_path, home_db)
    
    c.log(f"The energy barrier is : {barrier:.3f}")
    return True, {'trajectory': str(traj_file.relative_to(neb_dir)), 'neb_id': int(db_id)}

def monitor_log(logfile: Path, n_last=10, threshold=1e-3) -> bool:
    """
    Monitor the last `n_last` iterations in the log file to check for convergence.
    If the change in energy or forces is below `threshold`, return True for stagnation.
    """
    energies, forces = [], []

    with open(logfile, 'r') as f:
        lines = f.readlines()
        for line in lines[-n_last:]:
            if 'Step' in line:
                continue
            parts = line.split()
            # Some forces add a consistency check at the end
            energy = parts[-2][:-1] if not parts[-2][-1].isnumeric() else parts[-2] 
            energies.append(float(energy))
            forces.append(float(parts[-1]))

    if len(energies) < n_last:
        return False  # Not enough data yet

    energy_change = max(abs(energies[i] - energies[i+1]) for i in range(len(energies) - 1)) 
    force_change = max(abs(forces[i] - forces[i+1]) for i in range(len(forces) - 1))

    return energy_change < threshold and force_change < threshold

def run_neb_with_checker(neb: NEB, traj_file: Path, logfile: Path, fire: bool, fmax: float,
                         n_iter_check: int = 10, n_checks: int = 0, max_restarts: int = 3) -> None:
    """
    Run the NEB calculation with a convergence checker.
    The checker monitors the last `n_iter_check` iterations in the log file.
    It adjusts the parameters or switches the optimizer if stagnation is detected.
    """
    def setup_optimizer(fire, neb):
        dt = 0.05 #/ 2
        maxstep = 0.1 #/ 2
        dtmax = 0.3 #/ 2
        Nmin = 10
        finc = 1.05 #* 0.9
        fdec = 0.4 #* 1.1
        return FIRE(neb, trajectory=str(traj_file), logfile=logfile.as_posix(), dt=dt, maxstep=maxstep, dtmax=dtmax, 
                Nmin=Nmin, finc=finc, fdec=fdec, force_consistent=False) if fire else BFGS(neb, trajectory=str(traj_file), logfile=logfile.as_posix())

    optimizer = setup_optimizer(fire, neb)
    optimizer.run(fmax=fmax)

    check_counter = 0
    while check_counter < max_restarts:
        if monitor_log(logfile, n_last=n_iter_check):#, threshold=fmax / 10):
            c.log("Stagnation detected, adjusting optimizer parameters.")
            fire = not fire
            optimizer = setup_optimizer(fire, neb)
            optimizer.run(fmax=fmax)
            check_counter += 1
        else:
            break

    if check_counter >= max_restarts:
        c.log(f"Calculation stagnated after {max_restarts} restarts.")

def setup_run(atoms: Atoms, directory: Path, vasp_params: dict) -> Atoms:
    """
    Setup and run the DFT calculation for the given atoms object.

    Args:
        atoms (Atoms): The atomic structure to be calculated.
        directory (Path): The directory where the calculation will be performed.
        vasp_params (dict): Additional parameters for the VASP calculator.

    Returns:
        Atoms: The atoms object with the calculator attached.
    """
    directory.mkdir(parents=True, exist_ok=True)
    calc = create_Vasp_calc(atoms, 'PBEsol', directory)
    calc.set(**vasp_params,
            ibrion=-1,
            nsw=1,
            kpar=nnodes,
            ncore=ncore)
    set_magnetic_moments(atoms)
    return atoms

def update_or_write(db: SQLite3Database, atoms: Atoms, name: str, **kwargs):
    if db.count(name=name) > 0:
        ID = next(db.select(name=name)).id
        db.update(ID, atoms, name=name, **kwargs)
        return ID
    return db.write(atoms, name=name, **kwargs)
