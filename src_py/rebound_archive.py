"""
Utilities for working with Rebound archives.

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

# Core
import numpy as np

# Astronomy
import rebound

# Utilities
from datetime import datetime, timedelta
from collections import namedtuple
from itertools import chain
from pathlib import Path
import os
from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from utils import hash_id_crc32
from astro_utils import mjd_to_datetime, OrbitalElement

# Typing
from typing import List, Tuple, Dict, Set, Optional

# ********************************************************************************************************************* 
# Directory for archives; make if missing
dir_archive: str = '../data/rebound/archive'
Path(dir_archive).mkdir(parents=True, exist_ok=True)

# ********************************************************************************************************************* 
def calc_fname_archive(sim: rebound.simulation, mjd0: int, mjd1: int) -> str:
    """Generate filename for the simulation archive from a simulation"""
    # Tokens for frequency and hash of bodies
    freq_str = f'sf{sim.steps_per_day}'
    body_hash_str: str = f'bh{hash_id_crc32(tuple(sim.body_names))}'
    fname: str =  f'{sim.body_collection}_{mjd0}_{mjd1}_{freq_str}_{body_hash_str}.bin'
    return fname

# ********************************************************************************************************************* 
def make_archive_impl(sim_epoch: rebound.Simulation, 
                      mjd0: int, 
                      mjd1: int, 
                      time_step: float, 
                      save_step: int,
                      save_elements: bool,
                      progbar: bool) -> None:
    """
    Create a rebound simulation archive and save it to disk.
    INPUTS:
        fname_archive:  the file name to save the archive to
        sim_epoch:      rebound simulation object as of the epoch time; to be integrated in both directions
        mjd0:           the earliest MJD to simulate back to
        mjd1:           the latest MJD to simulate forward to
        time_step:      the time step in days for the simulation
        save_step:      the interval for saving snapshots to the simulation archive
        save_elements:  flag indiciting whether to save orbital elements
        progbar:        flag - whether to display a progress bar
    OUTPUTS:
        None.  saves the assembled simulation archive to disk.
    """
    # Get archive filename from simulation
    fname_archive: str = calc_fname_archive(sim_epoch, mjd0, mjd1)
    # Path of archive file
    path_archive: Path = Path(dir_archive, fname_archive)
    pathstr_archive: str = path_archive.as_posix()

    # Look up the epoch from the base simulation
    epoch: int = sim_epoch.epoch

    # Look up body IDs and names
    body_ids: np.ndarray = sim_epoch.body_ids
    body_names: np.ndarray = sim_epoch.body_names

    # Convert epoch, start and end times relative to a base date of the simulation start
    # This way, time is indexed from t0=0 to t1 = (dt1-dt0)
    t0: float = 0.0
    t1: float = mjd1 - mjd0
    t_epoch: float = float(epoch) - mjd0
    
    # Create copies of the simulation to integrate forward and backward
    sim_fwd: rebound.Simulation = sim_epoch.copy()
    sim_back: rebound.Simulation = sim_epoch.copy()

    # Set the time counter on both simulation copies to the epoch time
    sim_fwd.t = t_epoch
    sim_back.t = t_epoch
    # Flip sign of dt on sim_back
    sim_back.dt *= -1.0

    # Set the times for snapshots in both directions;
    ts: np.array = np.arange(t0, t1+time_step, time_step)
    idx: int = np.searchsorted(ts, t_epoch, side='left')
    ts_fwd: np.array = ts[idx:]
    ts_back: np.array = ts[:idx][::-1]

    # The epochs corresponding to the times in ts; as MJD and as datetime
    epochs: np.array = ts + (epoch - t_epoch)
    dt0: datetime = mjd_to_datetime(mjd0)
    epochs_dt: np.array = np.array([dt0 + timedelta(t) for t in ts])

    # File names for forward and backward integrations
    path_fwd: Path = Path(dir_archive, fname_archive.replace('.bin', '_fwd.bin'))
    path_back: Path = Path(dir_archive, fname_archive.replace('.bin', '_back.bin'))
    pathstr_fwd: str = path_fwd.as_posix()
    pathstr_back: str = path_back.as_posix()

    # Number of snapshots
    M_back: int = len(ts_back)
    M_fwd: int = len(ts_fwd)
    M: int = M_back + M_fwd
    # Number of particles
    N: int = sim_epoch.N

    # Initialize arrays for the position and velocity
    shape_qv: Tuple[int] = (M, N, 3)
    q: np.array = np.zeros(shape_qv, dtype=np.float64)
    v: np.array = np.zeros(shape_qv, dtype=np.float64)
    
    # Initialize arrays for orbital elements if applicable  
    if save_elements:
        # Arrays for a, e, inc, Omega, omega, f
        shape_elt: Tuple[int] = (M, N)
        orb_a: np.array = np.zeros(shape_elt, dtype=np.float64)
        orb_e: np.array = np.zeros(shape_elt, dtype=np.float64)
        orb_inc: np.array = np.zeros(shape_elt, dtype=np.float64)
        orb_Omega: np.array = np.zeros(shape_elt, dtype=np.float64)
        orb_omega: np.array = np.zeros(shape_elt, dtype=np.float64)
        orb_f: np.array = np.zeros(shape_elt, dtype=np.float64)

        # Wrap these into a named tuple
        elts: OrbitalElement = \
            OrbitalElement(a=orb_a, e=orb_e, inc=orb_inc,
                           Omega=orb_Omega, omega=orb_omega, f=orb_f)
    else:
        elts = OrbitalElement(a=0.0, e=0.0, inc=0.0, Omega=0.0, omega=0.0, f=0.0)

    # Subfunction: process one row of the loop    
    def process_row(sim: rebound.Simulation, fname: str, t: float, row: int):
        # Integrate to the current time step with an exact finish time
        sim.integrate(t, exact_finish_time=1)
        
        # Serialize the position and velocity
        sim.serialize_particle_data(xyz=q[row])
        sim.serialize_particle_data(vxvyvz=v[row])
        
        # Save a snapshot on multiples of save_step or the first / last row
        if (row % save_step == 0) or (row in (0, M-1)):
            sim.simulationarchive_snapshot(filename=fname)

        # Save the orbital elements if applicable
        if save_elements:
            # Compute the orbital elements vs. the sun as primary
            primary = sim.particles['Sun']
            orbits = sim.calculate_orbits(primary=primary, jacobi_masses=False)
            # Save the elements on this date as an Nx6 array
            # This approach only traverses the (slow) Python list orbits one time
            elt_array = np.array([[orb.a, orb.e, orb.inc, orb.Omega, orb.omega, orb.f] \
                                  for orb in orbits])
            # Save the elements into the current row of the named orbital elements arrays
            # The LHS row mask starts at 1 b/c orbital elements are not computed for the first object (Sun)
            orb_a[row, 1:] = elt_array[:, 0]
            orb_e[row, 1:] = elt_array[:, 1]
            orb_inc[row, 1:] = elt_array[:, 2]
            orb_Omega[row, 1:] = elt_array[:, 3]
            orb_omega[row, 1:] = elt_array[:, 4]
            orb_f[row, 1:] = elt_array[:, 5]

    # Integrate the simulation forward in time
    idx_fwd: List[Tuple[int, float]] = list(enumerate(ts_fwd))
    if progbar:
        idx_fwd = tqdm_auto(idx_fwd)
    for i, t in idx_fwd:
        # Row index for position data
        row: int = M_back + i
        # Process this row
        process_row(sim=sim_fwd, fname=pathstr_fwd, t=t, row=row)
        
    # Integrate the simulation backward in time
    idx_back: List[Tuple[int, float]] = list(enumerate(ts_back))
    if progbar:
        idx_back = tqdm_auto(idx_back)
    for i, t in idx_back:
        # Row index for position data
        row: int = M_back - (i+1)
        # Process this row
        process_row(sim=sim_back, fname=pathstr_back, t=t, row=row)

    # Load the archives with the forward and backward snapshots
    sa_fwd: rebound.SimulationArchive = rebound.SimulationArchive(pathstr_fwd) if idx_fwd else []
    sa_back: rebound.SimulationArchive = rebound.SimulationArchive(pathstr_back) if idx_back else []

    # Filename for numpy arrays of position and velocity
    path_np = Path(dir_archive, fname_archive.replace('.bin', '.npz'))
    pathstr_np = path_np.as_posix()

    # Save the object name hashes
    hashes: np.array = np.zeros(N, dtype=np.uint32)
    sim_epoch.serialize_particle_data(hash=hashes)

    # Combine the forward and backward archives in forward order from t0 to t1
    sims = chain(reversed(sa_back), sa_fwd)
    # Process each simulation snapshot in turn
    for i, sim in enumerate(sims):
        # Save a snapshot on multiples of save_step
        sim.simulationarchive_snapshot(pathstr_archive)   

    # Save the numpy arrays with the object hashes, position and velocity
    np.savez(pathstr_np, 
             q=q, v=v, elts=elts,
             ts=ts, epochs=epochs, epochs_dt=epochs_dt,
             hashes=hashes, body_ids=body_ids, body_names=body_names)
    
    # Delete the forward and backward archives
    path_fwd.unlink(missing_ok=True)
    path_back.unlink(missing_ok=True)
    
# ********************************************************************************************************************* 
def make_archive(sim_epoch: rebound.Simulation, 
                 mjd0: int, mjd1: int, 
                 time_step: int, save_step: int,
                 save_elements: bool,
                 progbar: bool) -> rebound.SimulationArchive:
    """
    Load a rebound archive if available; otherwise generate it and save it to disk.
    INPUTS:
        fname_archive:  the file name to save the archive to
        sim_epoch:      rebound simulation object as of the epoch time; to be integrated in both directions
        mjd0:           the earliest MJD to simulate back to
        mjd1:           the latest MJD to simulate forward to
        time_step:      the time step in days for the simulation
        save_step:      the interval for saving snapshots to the simulation archive
        save_elements:  flag indiciting whether to save orbital elements
        progbar:        flag - whether to display a progress bar
    OUTPUTS:
        SimulationArchive: Rebound SimulationArchive object loaded from disk
    """
    # Get archive filename from simulation
    fname_archive: str = calc_fname_archive(sim_epoch, mjd0, mjd1)
    # Path of archive
    path_archive = os.path.join(dir_archive, fname_archive)
    try:
        # First try to load the named archive
        sa = rebound.SimulationArchive(filename=path_archive)
        print(f'Loaded archive {fname_archive}.')
    except:
        # If the archive is not on disk, save it to disk
        print(f'Generating archive {fname_archive}\n'
              f'from mjd {mjd0} to {mjd1}, time_step={time_step}, save_step={save_step}...')
        make_archive_impl(sim_epoch=sim_epoch, mjd0=mjd0, mjd1=mjd1, time_step=time_step, 
                          save_step=save_step, save_elements=save_elements, progbar=progbar)
        # Load the new archive into memory
        sa = rebound.SimulationArchive(filename=path_archive)

    # Bind the archive filename to the archive
    sa.fname_archive: str = fname_archive
    sa.fname_np: str = fname_archive.replace('.bin', '.npz')
    sa.body_collection = sim_epoch.body_collection
    # Bind various dates and time steps
    sa.epoch = sim_epoch.epoch
    sa.mjd0: int = mjd0
    sa.mjd1: int = mjd1
    sa.time_step: int = time_step
    # Bind the simulation at the epoch
    sa.sim_epoch = sim_epoch

    return sa
