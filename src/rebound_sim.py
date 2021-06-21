"""
Utilities for creating and extending Rebound simulations.

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

# Core
import numpy as np

# Astronomy
import rebound

# Utilities
from pathlib import Path
import os

# MSE imports
from horizons import make_sim_horizons, extend_sim_horizons

# Typing
from typing import List, Optional

# ********************************************************************************************************************* 
# Directory for simulations; make if missing
dir_sim: str = '../data/rebound/sim'
Path(dir_sim).mkdir(parents=True, exist_ok=True)

# ********************************************************************************************************************* 
def make_sim(body_collection: str, body_names_add: Optional[List[str]], epoch: int, add_as_test: bool,
             integrator: str = 'ias15', steps_per_day: int = 4, epsilon: float = 2.0**-30,
             save_file: bool=True, load_file: bool=False, verbose: bool=False) -> rebound.Simulation:
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
        load_file:          Flag - whether to use a loaded simulation
    """
    # Filename for archive
    file_date: str = f'{epoch}'
    fname_sim: str = f'{body_collection}_{file_date}.bin'
    path_sim: Path = Path(dir_sim, f'{body_collection}_{file_date}.bin').as_posix()
    path_npz: Path = Path(dir_sim, f'{body_collection}_bodies.npz').as_posix()

    # If this file already exists, load it and check for both extra and missing bodies
    sim: rebound.Simulation
    sim_is_loaded: bool = False
    if load_file:
        try:
            # Attempt to load the archive file
            sim = rebound.Simulation(path_sim)
            # print(f'Loaded {fname_sim}')

            # Add body_ids and body_names to sim
            with np.load(path_npz, allow_pickle=True) as npz:
                sim.body_ids = npz['body_ids']
                sim.body_names = npz['body_names']
            # Set sim_is_loaded to indicate success
            sim_is_loaded=True
        except:
            # Initialize simulation
            if verbose:
                print(f'Unable to load {fname_sim}, building from Horizons data...')
    
    # If simulation was not loaded from disk (either by choice or a failure) generate it from Horizons DB
    if not sim_is_loaded:
        sim = make_sim_horizons(body_collection=body_collection, epoch=epoch, verbose=verbose)

    # DO NOT Move to center of momentum with sim.move_to_com() method !!!
    # The Horizons data is already in the barycentric frame

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
    sim.steps_per_day = steps_per_day   

    # Extend the simulation for additional bodies if they were requested
    if body_names_add:
        extend_sim_horizons(sim=sim, body_names=body_names_add, add_as_test=add_as_test)

    # Return the simulation
    return sim

# ********************************************************************************************************************* 
def make_sim_planets(epoch: int, integrator: str ='ias15', epsilon: float = 2.0**-32, 
                     steps_per_day: int = 48, load_file: bool=False):
    """Create a simulation with the sun and 8 planets at the specified time"""
    # Arguments for make_sim
    body_collection: str = 'Planets'
    body_names_add: Optional[List[str]] = None
    add_as_test: bool = True
    save_file = True

    # Build a simulation with the selected objects
    sim = make_sim(body_collection=body_collection, body_names_add=body_names_add, epoch=epoch,
                   add_as_test=add_as_test, integrator=integrator, 
                   epsilon=epsilon, steps_per_day=steps_per_day, save_file=save_file, load_file=load_file)
    return sim

# ********************************************************************************************************************* 
def make_sim_de435(epoch: int, integrator: str ='ias15', epsilon: float = 2.0**-32, 
                   steps_per_day: int=48, load_file: bool=False):
    """Create a simulation with all the massive objects used in the DE435 integration"""
    # Arguments for make_sim
    body_collection: str = 'DE435'
    body_names_add: Optional[List[str]] = None
    add_as_test: bool = True
    save_file = True

    # Build a simulation with the selected objects
    sim = make_sim(body_collection=body_collection, body_names_add=body_names_add, epoch=epoch,
                   add_as_test=add_as_test, integrator=integrator, 
                   epsilon=epsilon, steps_per_day=steps_per_day, save_file=save_file, load_file=load_file)
    return sim

# ********************************************************************************************************************* 
def reverse_velocity(sim):
    """Reverse the velocities in a simulation for backwards time integration; modifies sim in place"""
    for p in sim.particles:
        vx, vy, vz = p.vx, p.vy, p.vz
        p.vx = -vx
        p.vy = -vy
        p.vz = -vz

# # ********************************************************************************************************************* 
# def load_sim_np(fname_np: str) -> Tuple[np.array, np.array, Dict[str, np.array]]:
#     """Load numpy arrays for position, velocity, and catalog data from the named file"""
#     # Path of numpy data file
#     path_np= os.path.join(dir_archive, fname_np)
#     # Load the numpy data file
#     with np.load(path_np, allow_pickle=True) as npz:
#         # Extract position, velocity and hashes
#         q = npz['q']
#         v = npz['v']
#         elts_np = npz['elts']
#         ts = npz['ts']
#         epochs = npz['epochs']
#         epochs_dt = npz['epochs_dt']
#         hashes = npz['hashes']
#         body_ids = npz['body_ids']
#         body_names = npz['body_names']
#         # body_names_list: List[str] = [nm for nm in body_names]

#     # Wrap the catalog into a dictionary
#     catalog = {
#         'ts': ts,
#         'epochs': epochs,
#         'epochs_dt': epochs_dt,
#         'hashes': hashes,
#         'body_ids': body_ids,
#         'body_names': body_names,
#         }

#     # For some reason, np.save() squishes a namedtuple into an ND array.  Restore it to a named tuple
#     elts = OrbitalElement(a=elts_np[0], e=elts_np[1], inc=elts_np[2], 
#                           Omega=elts_np[3], omega=elts_np[4], f=elts_np[5], M=elts_np[6])

#     # Return the position, velocity, and catalog        
#     return q, v, elts, catalog
