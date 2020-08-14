"""
Harvard IACS Masters Thesis
Utilities for working with Rebound simulations and archives

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

# Core
import numpy as np
import pandas as pd

# Astronomy
import rebound

# Utilities
from datetime import datetime, timedelta
from collections import namedtuple
from pathlib import Path
import os
from itertools import chain

# UI
import matplotlib.pyplot as plt
# from tqdm import tqdm as tqdm_console
from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from utils import hash_id_crc32
from astro_utils import datetime_to_mjd, mjd_to_datetime, datetime_to_year, cart_to_sph
from db_config import db_engine
from horizons import make_sim_horizons, extend_sim_horizons, extend_sim_horizons_ast
from utils import rms

# Typing
from typing import List, Tuple, Dict, Set, Optional

# ********************************************************************************************************************* 
# Named tuple data type for orbital elements
OrbitalElement = namedtuple('OrbitalElement', 'a e inc Omega omega f')

# Directories for single simulations and archives
dir_sim: str = '../data/rebound/sim'
dir_archive: str = '../data/rebound/archive'

# ********************************************************************************************************************* 
def make_sim(body_collection: str, body_names_add: Optional[List[str]], epoch: int, add_as_test: bool,
             integrator: str = 'ias15', steps_per_day: int = 4, epsilon: float = 2.0**-30,
             save_file: bool = True) -> rebound.Simulation:
    """
    Create or load simulation with the specified objects at the specified time
    INPUTS:
        body_collection:    Collection of massive objects used to initialize simulation
        body_names_add:     Additional bodies added to the starting body collection
        add_as_test:        Flag indicating whether additional bodies are massless test particlesnot
        epoch:              MJD as of which state vectors are taken; INTEGER in range 40400 - 77600
        integrator:         String identifier for the Rebound integrator, e.g. 'ias15'
        steps_per_day:      Number of integration steps per day used to set dt and min_dt
        epsilon:            Tolerance for ias15 adaptive integrator; rebound default is 1.0E-9
        save_file:          Flag - whether to save the simulation to disk
    """
    # Filename for archive
    file_date: str = f'{epoch}'
    fname_sim: str = f'{body_collection}_{file_date}.bin'
    path_sim: Path = Path(dir_sim, f'{body_collection}_{file_date}.bin').as_posix()
    path_npz: Path = Path(dir_sim, f'{body_collection}_bodies.npz').as_posix()

    # If this file already exists, load it and check for both extra and missing bodies
    sim: rebound.Simulation
    try:
        # Attempt to load the archive file
        sim = rebound.Simulation(path_sim)
        # print(f'Loaded {fname_sim}')

        # Add body_ids and body_names to sim
        with np.load(path_npz, allow_pickle=True) as npz:
            sim.body_ids = npz['body_ids']
            sim.body_names = npz['body_names']
    except:           
        # Initialize simulation
        print(f'Unable to load {fname_sim}, building from Horizons data...')
        sim = make_sim_horizons(body_collection=body_collection, epoch=epoch)
    
    # Move to center of momentum
    # sim.move_to_com()

    # Set integrator and time step
    sim.integrator = integrator
    dt: float = 1.0 / steps_per_day if steps_per_day > 0 else 0
    sim.dt = dt
    # Special handling for the ias15 integrator
    if integrator == 'ias15':
        # the ias15 integrator object
        ias15 = sim.ri_ias15
        # Apply the time step dt as the minimum adaptive time step
        ias15.min_dt = dt
        # Set tolerance epsilon for adaptive time steps
        sim.ri_ias15.epsilon = epsilon
        # Compute error by max(acc_err / acc) when flag=0 (more demanding)
        # Compute error by max(acc_err) / max(acc) when flag=1 (less demanding)
        sim.ri_ias15.epsilon_global = 0

    # Save a snapshot to the archive file if requested
    if save_file:
        sim.simulationarchive_snapshot(filename=path_sim, deletefile=True)
        np.savez(path_npz, body_ids=sim.body_ids, body_names=sim.body_names)

    # Save additional data attributes to the simulation for use downstream
    sim.body_collection = body_collection
    sim.body_names_add = body_names_add
    sim.epoch = epoch
    sim.steps_per_day: int = steps_per_day   

    # Extend the simulation for additional bodies if they were requested
    if body_names_add:
        extend_sim_horizons(sim=sim, body_names=body_names_add, add_as_test=add_as_test)

    # Return the simulation
    return sim

# ********************************************************************************************************************* 
def extend_sim(sim: rebound.Simulation, body_names: List[str], add_as_test: bool):
    """
    Extend an existing simulation to include a list of named bodies
    INPUTS:
        sim:         A rebound simulation object
        body_names:  List of string body names; references DB table KS.Body, field name BodyName
                     The state vectors of these bodies will be added to the simulation, sim.
        add_as_test: Flag indicating whether bodies added as massless test particles or not
    RETURNS:
        None:        Modifies sim in place
    """
    # Generate list of missing object names
    body_names_present: Set[str] = set(sim.body_names)
    body_names_missing: List[str] = [nm for nm in body_names if nm not in body_names_present]
    objects_missing: bool = len(body_names_missing) > 0

    # Extend the simulation and save it with the augmented bodies
    if objects_missing:
        extend_sim_horizons(sim, body_names=body_names_missing, add_as_test=add_as_test)

    return sim

# ********************************************************************************************************************* 
def extend_sim_ast(sim: rebound.Simulation, n_ast: int, add_as_test: bool):
    """
    Extend an existing simulation to include a list of named bodies
    INPUTS:
        sim:         A rebound simulation object
        n_ast:       The number of asteroids to add to the simulation; adds the first n_ast asteroids
        add_as_test: Flag indicating whether bodies added as massless test particles or not
    RETURNS:
        None:        Modifies sim in place
    """
    # Extend the simulation and save it with the augmented bodies
    extend_sim_horizons_ast(sim, n_ast=n_ast, add_as_test=add_as_test)

    return sim

# ********************************************************************************************************************* 
def make_sim_planets(epoch: int, integrator: str ='ias15', epsilon: float = 2.0**-32, steps_per_day: int = 16):
    """Create a simulation with the sun and 8 planets at the specified time"""
    # Arguments for make_sim
    body_collection: str = 'Planets'
    body_names_add: Optional[List[str]] = None
    add_as_test: bool = True
    save_file = True

    # Build a simulation with the selected objects
    sim = make_sim(body_collection=body_collection, body_names_add=body_names_add, epoch=epoch,
                   add_as_test=add_as_test, integrator=integrator, 
                   epsilon=epsilon, steps_per_day=steps_per_day, save_file=save_file)
    return sim

# ********************************************************************************************************************* 
def make_sim_de435(epoch: int, integrator: str ='ias15', epsilon: float = 2.0**-32, steps_per_day: int = 16):
    """Create a simulation with all the massive objects used in the DE435 integration"""
    # Arguments for make_sim
    body_collection: str = 'DE435'
    body_names_add: Optional[List[str]] = None
    add_as_test: bool = True
    save_file = True

    # Build a simulation with the selected objects
    sim = make_sim(body_collection=body_collection, body_names_add=body_names_add, epoch=epoch,
                   add_as_test=add_as_test, integrator=integrator, 
                   epsilon=epsilon, steps_per_day=steps_per_day, save_file=save_file)
    return sim

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
    
# ********************************************************************************************************************* 
def load_sim_np(fname_np: str) -> Tuple[np.array, np.array, Dict[str, np.array]]:
    """Load numpy arrays for position, velocity, and catalog data from the named file"""
    # Path of numpy data file
    path_np= os.path.join(dir_archive, fname_np)
    # Load the numpy data file
    with np.load(path_np, allow_pickle=True) as npz:
        # Extract position, velocity and hashes
        q = npz['q']
        v = npz['v']
        elts_np = npz['elts']
        ts = npz['ts']
        epochs = npz['epochs']
        epochs_dt = npz['epochs_dt']
        hashes = npz['hashes']
        body_ids = npz['body_ids']
        body_names = npz['body_names']
        # body_names_list: List[str] = [nm for nm in body_names]

    # Wrap the catalog into a dictionary
    catalog = {
        'ts': ts,
        'epochs': epochs,
        'epochs_dt': epochs_dt,
        'hashes': hashes,
        'body_ids': body_ids,
        'body_names': body_names,
        }

    # For some reason, np.save() squishes a namedtuple into an ND array.  Restore it to a named tuple
    elts = OrbitalElement(a=elts_np[0], e=elts_np[1], inc=elts_np[2], 
                          Omega=elts_np[3], omega=elts_np[4], f=elts_np[5])

    # Return the position, velocity, and catalog        
    return q, v, elts, catalog

# ********************************************************************************************************************* 
def get_asteroids() -> pd.DataFrame:
    """Return list of known asteroid names and IDs"""
    with db_engine.connect() as conn:
        sql = 'CALL KS.GetAsteroids()'
        ast = pd.read_sql(sql, con=conn)
    return ast
