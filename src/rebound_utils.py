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
from utils import arange_inc
from astro_utils import datetime_to_mjd, mjd_to_datetime, datetime_to_year, cart_to_sph
from horizons import make_sim_horizons, extend_sim_horizons
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
def make_sim(body_collection: str, body_names_add: List[str], epoch: int, add_as_test: bool,
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
    # fname_sim: str = os.path.join(dir_sim, f'{body_collection}_{file_date}.bin')
    # fname_npz: str = os.path.join(dir_sim, f'{body_collection}_bodies.npz')
    fname_sim: str = f'{body_collection}_{file_date}.bin'
    path_sim: Path = Path(dir_sim, f'{body_collection}_{file_date}.bin')
    path_npz: Path = Path(dir_sim, f'{body_collection}_bodies.npz')    

    # If this file already exists, load it and check for both extra and missing bodies
    sim: rebound.Simulation
    try:
        # Attempt to load the archive file
        sim = rebound.Simulation(path_sim.as_posix())
        # print(f'Loaded {fname_sim}')

        # Add body_ids and body_names to sim
        with np.load(fname_npz, allow_pickle=True) as npz:
            sim.body_ids = npz['body_ids']
            sim.body_names = npz['body_names']
    except:           
        # Initialize simulation
        print(f'Unable to load {fname_sim}, building from Horizons data...')
        sim = make_sim_horizons(body_collection=body_collection, epoch=epoch)
    
    # Extend the simulation for additional bodies if they were requested
    if body_names_add:
        extend_sim_horizons(sim=sim, body_names=body_names_add, add_as_test=add_as_test)

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

    # Save a snapshot to the archive file if requested
    if save_file:
        sim.simulationarchive_snapshot(filename=path_sim.as_posix(), deletefile=True)
        np.savez(path_npz.as_posix(), body_ids=sim.body_ids, body_names=sim.body_names)

    # Save additional data attributes to the simulation for use downstream
    sim.body_collection = body_collection
    sim.body_names_add = body_names_add
    sim.epoch = epoch
    sim.steps_per_day: int = steps_per_day
    sim.archive_filename = f'{body_collection}_<mjd0>_<mjd1>_sf{steps_per_day}.bin'

    # Return the simulation
    return sim

# ********************************************************************************************************************* 
def extend_sim(sim: rebound.Simulation, 
               body_names: List[str], 
               add_as_test: bool):
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
    # hashes_present: Set[int] = set(p.hash.value for p in sim.particles)
    # body_names_missing: List[str] = [nm for nm in body_names if rebound.hash(nm).value not in hashes_present]
    body_names_present: Set[str] = set(sim.body_names)
    body_names_missing: List[str] = [nm for nm in body_names if nm not in body_names_present]

    # Extend the simulation and save it with the augmented bodies
    if objects_missing:
        extend_sim_horizons(sim, body_names=body_names_missing, add_as_test=add_as_test)

    return sim

# ********************************************************************************************************************* 
def calc_fname_archive(sim: rebound.simulation, mjd0: int, mjd1: int) -> str:
    """Generate filename for the simulation archive from a simulation"""
    fname_archive: str = sim.archive_filename.replace('<mjd0>', f'{mjd0}').replace('<mjd1>', f'{mjd1}')
    return fname_archive

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
    # fname_fwd: str = path_archive.replace('.bin', '_fwd.bin')
    # fname_back: str = path_archive.replace('.bin', '_back.bin')
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
        elts = OrbitalElement()

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
def sim_cfg_array(sim: rebound.Simulation, body_names: Optional[List[str]]=None) -> np.array:
    """Extract the Cartesian configuration of each body in the simulation"""
    # Allocate array of configurations
    num_objects: int = sim.N if body_names is None else len(body_names)
    cfgs: np.array = np.zeros(shape=(num_objects, 6))

    # Iterate through particles
    ps = sim.particles

    # Create list of (i, p) pairs
    if body_names is None:
        ips = enumerate(ps)
    else:
        ips = [(i, ps[body_name]) for i, body_name in enumerate(body_names)]

    # Extract the configuration of each particle
    for i, p in ips:
        cfgs[i] = np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz])

    return cfgs

# ********************************************************************************************************************* 
def sim_elt_array(sim: rebound.Simulation, body_names=None) -> np.array:
    """Extract the orbital elements of each body in the simulation"""
    # Allocate array of elements
    num_objects: int = sim.N-1 if body_names is None else len(body_names)
    elts: np.array = np.zeros(shape=(num_objects, 9))
    
    # Iterate through particles AFTER the primary
    ps: List[rebound.Particles] = sim.particles
    primary: rebound.Particle = ps[0]
    if body_names is None:
        for i, p in enumerate(ps[1:]):
            orb = p.calculate_orbit(primary=primary)
            elts[i] = np.array([orb.a, orb.e, orb.inc, orb.Omega, orb.omega, 
                                orb.f, orb.M, orb.pomega, orb.l])
    else:
        for i, body_name in enumerate(body_names):
            # look up the particle with this name
            p = ps[body_name]
            orb = p.calculate_orbit(primary=primary)
            elts[i] = np.array([orb.a, orb.e, orb.inc, orb.Omega, orb.omega, 
                                orb.f, orb.M, orb.pomega, orb.l])
    return elts

# ********************************************************************************************************************* 
def report_sim_difference(sim0: rebound.Simulation, sim1: rebound.Simulation, verbose: bool=False) -> \
                          Tuple[np.array, np.array]:
    """Report the difference between two simulations on a summary basis"""
    # Get the body names
    body_names = sim0.body_names

    # Extract configuration arrays for the two simulations
    cfg0: np.array = sim_cfg_array(sim0, body_names)
    cfg1: np.array = sim_cfg_array(sim1, body_names)
    
    # Displacement of each body to earth
    earth_idx: int = np.argmax(body_names == 'Earth')
    q0: np.array = cfg0[:, 0:3] - cfg0[earth_idx, 0:3]
    q1: np.array = cfg1[:, 0:3] - cfg1[earth_idx, 0:3]
    
    # Right Ascension and Declination
    r0, asc0, dec0 = cart_to_sph(q0)
    r1, asc1, dec1 = cart_to_sph(q1)

    # Error in asc and dec; convert from radians to arcseconds
    asc_err: np.array = np.degrees(np.abs(asc1-asc0)) * 3600
    dec_err: np.array = np.degrees(np.abs(dec1-dec0)) * 3600

    # Take differences
    cfg_diff: np.array = (cfg1 - cfg0)
    pos_diff: np.array = cfg_diff[:, 0:3]
    vel_diff: np.array = cfg_diff[:, 3:6]

    # Error in position and velocity in heliocentric coordinates; skip the sun
    pos_err: np.array = np.linalg.norm(pos_diff, axis=1)
    pos_err_den: np.array = np.linalg.norm(cfg0[:, 0:3], axis=1)
    pos_err_den[0] = 1.0
    pos_err_rel: np.array = pos_err / pos_err_den
    vel_err: np.array = np.linalg.norm(vel_diff, axis=1)
    vel_err_den: np.array = np.linalg.norm(cfg0[:, 3:6], axis=1)
    vel_err_den[0] = 1.0
    vel_err_rel: np.array = vel_err / vel_err_den

    if verbose:
        print(f'\nPosition difference - absolute & relative')
        print(f'(Angle errors in arcseconds, position in AU)')
        print(f'Body       : Phi     : Theta   : Pos AU  : Pos Rel : Vel Rel')
        body_names_short: List[str] = [nm.replace(' Barycenter', '') for nm in body_names]
        for i, nm in enumerate(body_names_short):
            print(f'{nm:10} : {asc_err[i]:5.2e}: {dec_err[i]:5.2e}: {pos_err[i]:5.2e}: '
                  f'{pos_err_rel[i]:5.2e}: {vel_err_rel[i]:5.2e}')
        print(f'Overall    : {rms(asc_err):5.2e}: {rms(dec_err):5.2e}: {rms(pos_err):5.2e}: '
              f'{rms(pos_err_rel):5.2e}: {rms(vel_err_rel):5.2e}')

    # Extract orbital element arrays from the two simulations
    elt0: np.array = sim_elt_array(sim0, body_names[1:])
    elt1: np.array = sim_elt_array(sim1, body_names[1:])

    # Take differences
    elt_diff: np.array = (elt1 - elt0)
    # Angle differences are mod two pi
    two_pi: float = 2.0 * np.pi
    elt_diff[:,2:] = (elt_diff[:,2:] +np.pi ) % two_pi - np.pi
    
    # Compute RMS difference by orbital element
    elt_rms: np.array = rms(elt_diff, axis=0)
    elt_err: np.array = np.abs(elt_diff)

    # Names of selected elements
    elt_names: List[str] = ['a', 'e', 'inc', 'Omega', 'omega', 'f', 'M', 'pomega', 'long']

    # Report RMS orbital element differences
    if verbose:
        print(f'\nOrbital element errors:')
        print(f'elt    : RMS      : worst      : max_err  : HRZN        : REB')
        for j, elt in enumerate(elt_names):
            idx = np.argmax(elt_err[:, j])
            worse = body_names_short[idx+1]
            print(f'{elt:6} : {elt_rms[j]:5.2e} : {worse:10} : {elt_err[idx, j]:5.2e} : '
                  f'{elt0[idx, j]:11.8f} : {elt1[idx, j]:11.8f}')
        print(f'RMS (a, e, inc) =          {rms(elt_diff[:,0:3]):5.2e}')
        print(f'RMS (f, M, pomega, long) = {rms(elt_diff[:,5:9]):5.2e}')      

    # One summary error statistic
    ang_err: np.array = rms(np.array([asc_err, dec_err]))
    
    # Return the RMS position error and angle errors
    return pos_err, ang_err
    
# ********************************************************************************************************************* 
def test_integration(sa: rebound.SimulationArchive,
                     test_bodies: List[str], 
                     sim_name: str, 
                     test_name: str, 
                     verbose: bool = False,
                     make_plot: bool = False) -> \
                     Tuple[np.array, np.array]:
    """Test the integration of the planets against Horizons data"""
    # Extract body_collection from the simulation archive
    body_collection: str = sa.body_collection

    # Test bodies is input manually in case it's a strict subset of the bodies in the collection
    if test_bodies is None:
        test_bodies = sa.sim_epoch.body_names

    # Extract mjd0, mjd1 from simulation archive
    mjd0: int = sa.mjd0
    mjd1: int = sa.mjd1
    # Set the time step for comparing snapshots
    step: int = 20

    # Times to be tested - every step days in simulation range
    t_test = np.arange(np.int32(sa.tmin), np.int32(sa.tmax), np.int32(step)).astype(np.int32)

    # Epochs to be tested - generate from test times
    test_epochs = t_test + mjd0
    sz: int = test_epochs.shape[0]
    # Generate verbose output on the last test epoch if requested
    verbose_epochs = test_epochs[sz-1:sz]

    # Convert epochs to dates and years for plotting
    test_dates = [mjd_to_datetime(mjd) for mjd in test_epochs]
    test_years = [datetime_to_year(test_date) for test_date in test_dates]
    
    # Errors on these dates
    ang_errs: List[np.array] = []
    pos_errs: List[np.array] = []
    
    # Test the specified times
    for t in t_test:
        # The epoch of this test time
        epoch: int = mjd0 + t
        # The reference simulation from Horizons
        sim0: rebound.Simulation = make_sim_horizons(body_collection=body_collection, epoch=epoch)
        # The test simulation from the simulation archive
        sim1: rebound.Simulation = sa.getSimulation(t=t, mode='exact')
        # Verbosity flag and screen print if applicable
        report_this_date: bool = (epoch in verbose_epochs) and verbose
        if report_this_date:
            dt_t: datetime = mjd_to_datetime(epoch)
            print(f'\nDifference on {dt_t}:')
        # Run the test
        pos_err, ang_err = report_sim_difference(sim0=sim0, sim1=sim1, verbose=report_this_date)
        # Save position and angle errors
        pos_errs.append(pos_err)
        ang_errs.append(ang_err)
    
    # Plot error summary
    pos_err_rms: np.array = np.array([rms(x) for x in pos_errs])
    ang_err_rms: np.array = np.array([rms(x) for x in ang_errs])
        
    if make_plot:
        # Chart titles
        sim_name_chart = sim_name
        test_name_chart = test_name

        # Error in the position
        fig, ax = plt.subplots(figsize=[16,10])
        ax.set_title(f'Position Error of {test_name_chart} in {sim_name_chart} Integration')
        ax.set_ylabel(f'RMS Position Error in AU')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0,))
        # ax.plot(test_years, pos_err_rms, marker='o', color='red')
        ax.plot(test_epochs, pos_err_rms, marker='o', color='red')
        ax.grid()
        fname: str = f'../figs/integration_test/sim_error_{sim_name}_{test_name}_pos.png'
        fig.savefig(fname=fname, bbox_inches='tight')

        # Error in the angle
        fig, ax = plt.subplots(figsize=[16,10])
        ax.set_title(f'Angle Error of {test_name_chart} in {sim_name_chart} Integration')
        ax.set_ylabel(f'RMS Angle Error vs. Earth in Arcseconds')
        # ax.plot(test_years, ang_err_rms, marker='o', color='blue')
        ax.plot(test_epochs, ang_err_rms, marker='o', color='blue')
        ax.grid()
        fname: str = f'../figs/integration_test/sim_error_{sim_name}_{test_name}_angle.png'
        fig.savefig(fname=fname, bbox_inches='tight')
    
    if verbose:
        print(f'\nError by Date:')
        print('DATE       : ANG   : AU  ')
        for i, dt_t in enumerate(test_dates):
            print(f'{dt_t.date()} : {ang_err_rms[i]:5.3f} : {pos_err_rms[i]:5.3e}')
    
    # Compute average error
    mean_ang_err = np.mean(ang_err_rms)
    mean_pos_err = np.mean(pos_err_rms)
    print(f'\nMean RMS error in {sim_name} integration of {test_name} test objects:')
    print(f'AU   : {mean_pos_err:5.3e}')
    print(f'angle: {mean_ang_err:5.3f}')
    
    # Return summary of errors in position and angles
    return pos_err_rms, ang_err_rms
