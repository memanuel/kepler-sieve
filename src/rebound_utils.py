"""
Utilities for working with Rebound simulations.

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

# Core
import numpy as np
import pandas as pd
import sqlalchemy

# Astronomy
import rebound

# Utilities
import time
from datetime import datetime, timedelta
from pathlib import Path
import os

# UI
import matplotlib.pyplot as plt
from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from utils import hash_id_crc32, rms
from astro_utils import mjd_to_datetime, datetime_to_mjd, OrbitalElement
from db_config import db_engine
from horizons import make_sim_horizons, extend_sim_horizons, extend_sim_horizons_ast

# Typing
from typing import List, Tuple, Dict, Set, Optional

# ********************************************************************************************************************* 
# Directory for simulations; make if missing
dir_sim: str = '../data/rebound/sim'
Path(dir_sim).mkdir(parents=True, exist_ok=True)

# ********************************************************************************************************************* 
def get_asteroids() -> pd.DataFrame:
    """Return list of known asteroid names and IDs"""
    with db_engine.connect() as conn:
        sql = 'CALL KS.GetAsteroids()'
        ast = pd.read_sql(sql, con=conn)
    return ast

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
        # print(f'Saved simulation snapshot to {path_sim} and body names and IDs to {path_npz}.')

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
def integrate_numpy(sim_epoch: rebound.Simulation, 
                    mjd0: int, 
                    mjd1: int, 
                    time_step: float, 
                    save_elements: bool,
                    progbar: bool) -> None:
    """
    Perform an integration in rebound and save to Numpy arrays.
    INPUTS:
        sim_epoch:      rebound simulation object as of the epoch time; to be integrated in both directions
        mjd0:           the earliest MJD to simulate back to
        mjd1:           the latest MJD to simulate forward to
        time_step:      the time step in days for the simulation
        save_elements:  flag indiciting whether to save orbital elements
        progbar:        flag - whether to display a progress bar
    OUTPUTS:
        Numpy arrays with output of integration
        body_ids:       N integer body IDs of the bodies that were integrated
        body_names:     N body names
        epochs:         M times as of which the integration was saved; MJDs
        q:              MxNx3 array of positions in AU
        v:              MxNx3 array of velocities in AU / day
    """
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
    def process_row(sim: rebound.Simulation, t: float, row: int):
        # Integrate to the current time step with an exact finish time
        sim.integrate(t, exact_finish_time=1)
        
        # Serialize the position and velocity
        sim.serialize_particle_data(xyz=q[row])
        sim.serialize_particle_data(vxvyvz=v[row])
        
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
        process_row(sim=sim_fwd, t=t, row=row)
        
    # Integrate the simulation backward in time
    idx_back: List[Tuple[int, float]] = list(enumerate(ts_back))
    if progbar:
        idx_back = tqdm_auto(idx_back)
    for i, t in idx_back:
        # Row index for position data
        row: int = M_back - (i+1)
        # Process this row
        process_row(sim=sim_back, t=t, row=row)

    return body_ids, body_names, epochs, q, v, elts

# ********************************************************************************************************************* 
def integration_np2df(body_ids: np.array, body_names: np.array, epochs: np.array, 
                      q: np.array, v: np.array):
    """
    Arrange Numpy arrays with integration output into a Pandas DataFrame with one row per observation.\
    INPUTS:
        body_ids:       N integer body IDs of the bodies that were integrated
        body_names:     N body names
        epochs:         M times as of which the integration was saved; MJDs
        q:              MxNx3 array of positions in AU
        v:              MxNx3 array of velocities in AU / day
    OUTPUTS:
        df: DataFrame with columns
        BodyID:       N integer body IDs of the bodies that were integrated
        BodyName:     N body names
        MJD:          M times as of which the integration was saved; MJDs
        qx, qy, qx:   Positions in AU in the BME
        vx, vy, vz:   Velocities in AU / day
    """
    # Array sizes
    M: int = epochs.shape[0]
    N: int = body_ids.shape[0]

    # The time stamps
    MJD = epochs.repeat(N)
    TimeID = np.array(MJD*24*60).astype(np.int32)

    # The ID and name of each body
    BodyID = np.tile(body_ids, M)
    BodyName = np.tile(body_names, M)
    
    # The positions
    qx = q[:,:,0].flatten()
    qy = q[:,:,1].flatten()
    qz = q[:,:,2].flatten()

    # The velocities
    vx = v[:,:,0].flatten()
    vy = v[:,:,1].flatten()
    vz = v[:,:,2].flatten()    

    # Wrap into a dictionary
    data_dict = {
        'TimeID': TimeID,
        'BodyID': BodyID,
        'BodyName': BodyName,
        'MJD': MJD,
        'qx': qx,
        'qy': qy,
        'qz': qz,
        'vx': vx,
        'vy': vy,
        'vz': vz,
    }

    # Convert to a DataFrame
    df = pd.DataFrame(data_dict)
    return df

# ********************************************************************************************************************* 
def integrate_df(sim_epoch: rebound.Simulation, 
                 mjd0: int, 
                 mjd1: int, 
                 time_step: float, 
                 save_elements: bool,
                 progbar: bool) -> None:
    """
    Perform an integration in rebound and save to Pandas DataFrame.
    INPUTS:
        sim_epoch:      rebound simulation object as of the epoch time; to be integrated in both directions
        mjd0:           the earliest MJD to simulate back to
        mjd1:           the latest MJD to simulate forward to
        time_step:      the time step in days for the simulation
        save_elements:  flag indiciting whether to save orbital elements
        progbar:        flag - whether to display a progress bar
    OUTPUTS:
        df: DataFrame with columns
        BodyID:       N integer body IDs of the bodies that were integrated
        BodyName:     N body names
        MJD:          M times as of which the integration was saved; MJDs
        qx, qy, qx:   Positions in AU in the BME
        vx, vy, vz:   Velocities in AU / day
    """
    # Delegate to integrate_numpy
    body_ids, body_names, epochs, q, v, elts = \
            integrate_numpy(sim_epoch=sim_epoch, mjd0=mjd0, mjd1=mjd1, time_step=time_step, 
            save_elements=save_elements, progbar=progbar)    
    # Delegate to integration_np2df
    df = integration_np2df(body_ids=body_ids, body_names=body_names, epochs=epochs, q=q, v=v)

    return df

