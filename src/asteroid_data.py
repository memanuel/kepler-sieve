"""
Harvard IACS Masters Thesis
Generate TensorFlow datasets for asteroid trajectories.

Michael S. Emanuel
Sat Sep 21 10:38:38 2019
"""

# Core
import tensorflow as tf
import scipy, scipy.interpolate
import numpy as np
import pandas as pd

# Astronomy
import astropy
from astropy.units import au, day

# Utility
from datetime import datetime
from tqdm.auto import tqdm

# Local imports
from astro_utils import datetime_to_mjd
from asteroid_integrate import load_ast_elt
from asteroid_dataframe import spline_ast_vec_dir
from rebound_utils import load_sim_np

# Type names
from typing import Tuple, Dict

# ********************************************************************************************************************* 
# DataFrame of asteroid snapshots
ast_elt = load_ast_elt()

# Alias number of dimensions in space to avoid magic number 3
space_dims = 3

# Speed of light; express this in AU / day
light_speed_au_day = astropy.constants.c.to(au / day).value

# Data type
dtype = tf.float32

# ********************************************************************************************************************* 
def get_earth_sun_pos_file():
    """Get the position and velocity of earth and sun consistent with the asteroid data; look up from data file"""
    # selected data type for TF tensors
    dtype = np.float32
    
    # Name of the numpy archive
    n0, n1 = 0, 1000
    fname_np: str = f'../data/asteroids/sim_asteroids_n_{n0:06}_{n1:06}.npz'

    # The full array of positions and velocities
    q, v, elts, catalog = load_sim_np(fname_np=fname_np)

    # The object names
    object_names = catalog['object_names']

    # The snapshot times; offset to start time t0=0
    ts = catalog['ts'].astype(dtype)
    # Convert ts from relative time vs. t0 to MJD
    dt0 = datetime(2000, 1, 1)
    t_offset = datetime_to_mjd(dt0)
    ts += t_offset
   
    # Extract position of earth and sun in Barycentric coordinates
    earth_idx = object_names.index('Earth')
    sun_idx = object_names.index('Sun')
    q_earth = q[:, earth_idx, :]
    q_sun = q[:, sun_idx, :]

    # Velocity of earth and sun in Barycentric coordinates
    v_earth = v[:, earth_idx, :]
    v_sun = v[:, sun_idx, :]

    return ts, q_earth.astype(dtype), q_sun.astype(dtype), v_earth.astype(dtype), v_sun.astype(dtype)

# ********************************************************************************************************************* 
# Create 1D linear interpolator for earth positions; just need one instance for the whole module

# Get position and velocity of earth and sun in barycentric frame at reference dates from file
ts, q_earth, q_sun, v_earth, v_sun = get_earth_sun_pos_file()

# Build interpolator for barycentric position of earth and sun
earth_interp = scipy.interpolate.interp1d(x=ts, y=q_earth, kind='cubic', axis=0)
sun_interp = scipy.interpolate.interp1d(x=ts, y=q_sun, kind='cubic', axis=0)

# Build interpolator for barycentric velocity of earth and sun
earth_interp_v = scipy.interpolate.interp1d(x=ts, y=v_earth, kind='cubic', axis=0)
sun_interp_v = scipy.interpolate.interp1d(x=ts, y=v_sun, kind='cubic', axis=0)

# Build interpolator for heliocentric earth position (old; don't use this anymore!)
q_earth_helio = q_earth - q_sun
earth_interp_helio = scipy.interpolate.interp1d(x=ts, y=q_earth_helio, kind='cubic', axis=0)

# Alias earth_interp to earth_interp_bary for consistent naming scheme
earth_interp_bary = earth_interp

# ********************************************************************************************************************* 
def get_earth_pos(ts: np.ndarray, dtype = np.float64) -> np.array:
    """
    Get position of earth consistent with asteroid data at the specified times (MJDs) in barycentric frame
    INPUTS:
        ts: Array of times expressed as MJDs
        dtype: Desired datatype, e.g. np.float64 or np.float32
    """
    # Compute interpolated position at desired times
    q_earth = earth_interp(ts).astype(dtype)

    return q_earth

# ********************************************************************************************************************* 
def get_earth_pos_helio(ts: np.ndarray) -> np.array:
    """
    Get position of earth consistent with asteroid data at the specified times (MJDs) in heliocentric frame
    INPUTS:
        ts: Array of times expressed as MJDs
    """
    # Compute interpolated position at desired times
    q_earth = earch_interp_helio(ts)

    return q_earth

# ********************************************************************************************************************* 
def get_earth_pos_flex(ts: np.ndarray, heliocentric: bool) -> np.array:
    """
    Get position of earth consistent with asteroid data at the specified times (MJDs) in barycentric frame
    INPUTS:
        ts: Array of times expressed as MJDs
        heliocentric: flag indicating whether to use heliocentric coordinates;
                      default false uses barycentric coordinates
    """
    # Use selected interpolator
    interpolator = earth_interp_helio if heliocentric else earth_interp_bary

    # Compute interpolated position at desired times
    q_earth = interpolator(ts)

    return q_earth

# ********************************************************************************************************************* 
def get_sun_pos_vel(ts: np.ndarray, dtype = np.float64) -> np.array:
    """
    Get barycentric position and velocity of sun consistent with asteroid data at the specified times (MJDs)
    INPUTS:
        ts:    Array of times expressed as MJDs
        dtype: Desired datatype, e.g. np.float64 or np.float32
    """
    # Compute interpolated position at desired times
    q_sun = sun_interp(ts).astype(dtype)
    v_sun = sun_interp_v(ts).astype(dtype)
    return q_sun, v_sun

# ********************************************************************************************************************* 
# def get_earth_sun_pos(ts) -> np.array:
#     """
#     Get barycentric position of earth and sun consistent with asteroid data at the specified times (MJDs)
#     INPUTS:
#         ts: Array of times expressed as MJDs
#     """
#     # Choose barycentric or heliocentric interpolator
#     interpolator_earth = earth_interp

#     # Compute interpolated position at desired times
#     q_earth = earth_interp(ts)
#     q_sun = sun_interp(ts)
#     return q_earth, q_sun

# ********************************************************************************************************************* 
# Dataset of predicted asteroid directions
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_data_one_file(n0: int, n1: int) -> Tuple[Dict[str, np.array], Dict[str, np.array]]:
    """
    Wrap the data in one file of asteroid trajectory data into dictionary of numpy arrays
    INPUTS:
        n0: the first asteroid in the file, e.g. 0
        n1: the last asteroid in the file (exclusive), e.g. 1000
    OUTPUTS:
        inputs: a dict of numpy arrays
        outputs: a dict of numpy arrays
    """
    # selected data type for TF tensors
    dtype = np.float32
    
    # Name of the numpy archive
    fname_np: str = f'../data/asteroids/sim_asteroids_n_{n0:06}_{n1:06}.npz'
    
    # The full array of positions and velocities
    q, v, elts, catalog = load_sim_np(fname_np=fname_np)

    # The object names
    object_names = catalog['object_names']

    # The snapshot times; offset to start time t0=0
    ts = catalog['ts'].astype(dtype)
    # Convert ts from relative time vs. t0 to MJD
    dt0 = datetime(2000, 1, 1)
    t_offset = datetime_to_mjd(dt0)
    ts += t_offset
    
    # mask for selected asteroids
    mask = (n0 <= ast_elt.Num) & (ast_elt.Num < n1)
    
    # count of selected asteroids
    # asteroid_name = list(ast_elt.Name[mask].to_numpy())
    N_ast: int = np.sum(mask)
    # offset for indexing into asteroids; the first [10] objects are sun and planets
    ast_offset: int = len(object_names) - N_ast

    # Extract position of the sun in barycentric coordinates
    sun_idx = 0
    q_sun = q[:, sun_idx, :]
    v_sun = q[:, sun_idx, :]
    # Position of the earth in barycentric or heliocentric coordinates
    earth_idx = 3
    q_earth = q[:, earth_idx, :]
    # v_earth = v[:, earth_idx :]

    # shrink down q and v to slice with asteroid data only; 
    q = q[:, ast_offset:, :]
    v = v[:, ast_offset:, :]
    
    # swap axes for time step and body number; TF needs inputs and outputs to have same number of samples
    # this means that inputs and outputs must first be indexed by asteroid number, then time time step
    q = np.swapaxes(q, 0, 1)
    v = np.swapaxes(v, 0, 1)
    
    # Compute relative displacement to earth; instantaneous, before light time adjustment
    q_rel_inst = q - q_earth
    # Distance to earth
    r_earth_inst = np.linalg.norm(q_rel_inst, axis=2, keepdims=True)
    # Light time in days from asteroid to earth in days (time units is days)
    light_time = (r_earth_inst / light_speed_au_day)
    # Adjusted relative position, accounting for light time; simulation velocity units are AU /day
    dq_lt = light_time * v
    q_rel = q_rel_inst - dq_lt
    # Adjusted distance to earth, accounting for light time
    r = np.linalg.norm(q_rel, axis=2, keepdims=True)
    # Direction from earth to asteroid as unit vectors u = (ux, uy, uz)    
    u = q_rel / r

    # dict with inputs   
    inputs = {
        'a': ast_elt.a[mask].to_numpy().astype(dtype),
        'e': ast_elt.e[mask].to_numpy().astype(dtype),
        'inc': ast_elt.inc[mask].to_numpy().astype(dtype),
        'Omega': ast_elt.Omega[mask].to_numpy().astype(dtype),
        'omega': ast_elt.omega[mask].to_numpy().astype(dtype),
        'f': ast_elt.f[mask].to_numpy().astype(dtype),
        'epoch': ast_elt.epoch_mjd[mask].to_numpy().astype(dtype),
        'asteroid_num': ast_elt.Num[mask].to_numpy(),
        'asteroid_name': ast_elt.Name[mask].to_numpy(),
        'ts': np.tile(ts, reps=(N_ast,1,)),
    }
    
    # dict with outputs
    outputs = {
        'q': q.astype(dtype),
        'v': v.astype(dtype),
        'u': u.astype(dtype),
        'r': r.astype(dtype)
    }
    
    return inputs, outputs
