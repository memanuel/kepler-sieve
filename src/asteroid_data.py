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
from asteroid_integrate import load_data as load_data_asteroids
from rebound_utils import load_sim_np

# Type names
from typing import Tuple, Dict

# ********************************************************************************************************************* 
# DataFrame of asteroid snapshots
ast_elt = load_data_asteroids()

# Alias number of dimensions in space to avoid magic number 3
space_dims = 3

# Speed of light; express this in AU / day
light_speed_au_day = astropy.constants.c.to(au / day).value

# ********************************************************************************************************************* 
def make_ragged_tensors(t_np: np.array, u_np: np.array, ast_num_np: np.array, batch_size: int):
    """
    Convert t, u, ast_num into ragged tensors
    INPUTS:
        t_np: A numpy array with observation times as MJDs; size (N,)
        u_np: A numpy array with directions as MJDs; size (N,3)
        ast_num_np: A numpy array with the asteroid numbers; size (N,)
        batch_size: Used to pad end of data set with times having zero observations
    """
    # Unique times and their indices
    t_unq, inv_idx, = np.unique(t_np, return_inverse=True)
    
    # Pad t_unq with extra times so it has an even multiple of batch_size
    data_size: int = t_unq.shape[0]
    if batch_size is not None:
        num_pad: int = -data_size % batch_size
        t_unq = np.pad(t_unq, pad_width=(0,num_pad), mode='constant')
        t_unq[data_size:] = t_unq[data_size-1] + np.arange(num_pad)

    # The row IDs for the ragged tensorflow are what numpy calls the inverse indices    
    value_rowids = inv_idx 

    # Tensor with distinct times
    t = tf.convert_to_tensor(value=t_unq)
    # Ragged tensors for direction u and asteroid number ast_num
    u = tf.RaggedTensor.from_value_rowids(values=u_np, value_rowids=inv_idx)
    ast_num = tf.RaggedTensor.from_value_rowids(values=ast_num_np, value_rowids=value_rowids)

    # Return the tensors for t, u, ast_num
    return t, u, ast_num
    
# ********************************************************************************************************************* 
def make_ztf_dataset(ztf: pd.DataFrame, batch_size: int= None):
    """
    Create a tf.Dataset from ZTF DataFrame
    INPUTS:
        ztf: A DataFrame of ZTF data.  Columns must include mjd, ux, uy, uz.
        batch_size: Batch size used for TensorFlow DataSet.  None wraps all data into one big batch.
    Outputs:
        ds: A TensorFlow DataSet object.
    """
    # Extract inputs for make_ragged_tensors from ztf DataFrame
    t_np = ztf.mjd.values.astype(np.float32)
    cols_u = ['ux', 'uy', 'uz']
    u_np = ztf[cols_u].values.astype(np.float32)
    ast_num_np = ztf.nearest_ast_num.values

    # Sort arrays in increasing order of observation time; otherwise make_ragged_tensors won't work
    sort_idx = np.argsort(t_np)
    t_np = t_np[sort_idx]
    u_np = u_np[sort_idx]
    ast_num_np = ast_num_np[sort_idx]

    # Convert to tensors (regular for t, ragged for u and ast_num)
    t, u_r, ast_num_r = make_ragged_tensors(t_np=t_np, u_np=u_np, ast_num_np=ast_num_np, batch_size=batch_size)

    # Set batch_size to all the times if it was not specified
    if batch_size is None:
        batch_size = u_r.shape[0]

    # Number of entries to pad at the end so batches all divisible by batch_size    
    num_pad = t.shape[0] - u_r.shape[0]

    # Get row lengths
    row_len = tf.cast(u_r.row_lengths(), tf.int32)
    # Pad row len with zeros for the last batch
    row_len = tf.pad(row_len, paddings=[[0,num_pad]], mode='CONSTANT', constant_values=0)
    
    # Pad u_r into a regular tensor
    pad_default = np.array([0.0, 0.0, 65536.0])
    u = u_r.to_tensor(default_value=pad_default)
    
    # Pad u with entries for the last batch
    pad_shape_u = [[0,num_pad], [0,0], [0,0]]
    u = tf.pad(u, paddings=pad_shape_u, mode='CONSTANT', constant_values=65536.0)
        
    # Index associated with time points
    idx = tf.range(t.shape[0], dtype=tf.int32)
    
    # Pad ast_num_r into a regular tensor; use default -1 to indicate padded observation
    ast_num = ast_num_r.to_tensor(default_value=-1)
    
    # Pad ast_num with entries for the last batch
    pad_shape_ast_num = [[0,num_pad],[0,0]]
    ast_num = tf.pad(ast_num, paddings=pad_shape_ast_num, mode='CONSTANT', constant_values=-1)
    
    # Wrap into tensorflow Dataset
    inputs = {
        't': t,
        'idx': idx,
        'row_len': row_len,
        'u_obs': u, 
    }
    outputs = {
        'ast_num': ast_num,
    }
    ds = tf.data.Dataset.from_tensor_slices((inputs, outputs))
    
    # Batch the dataset
    drop_remainder: bool = False
    buffer_size = 16384
    ds = ds.batch(batch_size, drop_remainder=drop_remainder).shuffle(buffer_size).repeat()
    # Return the dataset as well as tensors with the time snapshots and row lengths
    return ds, t, row_len
    
# ********************************************************************************************************************* 
def orbital_element_batch(ast_nums: np.ndarray):
    """
    Return a batch of orbital elements
    INPUTS:
        ast_nums: Numpy array of asteroid numbers to include in batch
    OUTPUTS:
        elts: Dictionary with seven keys for a, e, inc, Omega, omega, f, epoch
    """
    # The orbital elements and epoch
    dtype = np.float32
    a = ast_elt.a[ast_nums].to_numpy().astype(dtype)
    e = ast_elt.e[ast_nums].to_numpy().astype(dtype)
    inc = ast_elt.inc[ast_nums].to_numpy().astype(dtype)
    Omega = ast_elt.Omega[ast_nums].to_numpy().astype(dtype)
    omega = ast_elt.omega[ast_nums].to_numpy().astype(dtype)
    f = ast_elt.f[ast_nums].to_numpy().astype(dtype)
    epoch = ast_elt.epoch_mjd[ast_nums].to_numpy().astype(dtype)
    
    # Wrap into dictionary
    elts = {
        'a': a,
        'e': e,
        'inc': inc,
        'Omega': Omega,
        'omega': omega,
        'f': f,
        'epoch': epoch
    }
    
    return elts

# ********************************************************************************************************************* 
def orbital_element_batch_v1(n0: int, batch_size: int=64):
    """Return a batch of orbital elements"""
    # Get start and end index location of this asteroid number
    i0: int = ast_elt.index.get_loc(n0)
    i1: int = i0 + batch_size
    ast_nums = np.arange(i0, i1+1, dytpe=np.int32)

    # Delegate to orbital_element_batch
    return orbital_element_batch(ast_nums)

# ********************************************************************************************************************* 
# OLD STUFF - PROBABLY GET RID OF THIS LATER
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def get_earth_pos_file(heliocentric: bool):
    """Get the position of earth consistent with the asteroid data; look up from data file"""
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
   
    # Extract position of sun and earth in Barycentric coordinates
    sun_idx = object_names.index('Sun')
    earth_idx = object_names.index('Earth')
    q_sun_bary = q[:, sun_idx, :]
    q_earth_bary = q[:, earth_idx, :]
    # Compute earth position in heliocentric coordinates
    q_earth_helio = q_earth_bary - q_sun_bary
    # Selected flavor (barycentric vs. heliocentric)
    q_earth = q_earth_helio if heliocentric else q_earth_bary
    # Convert to selected data type
    return q_earth.astype(dtype), ts

# ********************************************************************************************************************* 
# Create 1D linear interpolator for earth positions; just need one instance for the whole module

# Get position of earth at reference dates from file
q_earth_bary, t_bary = get_earth_pos_file(heliocentric=False)
q_earth_helio, t_helio = get_earth_pos_file(heliocentric=True)
# Build the interpolators for barycentric and heliocentric
earth_interp_bary = scipy.interpolate.interp1d(x=t_bary, y=q_earth_bary, kind='cubic', axis=0)
earth_interp_helio = scipy.interpolate.interp1d(x=t_helio, y=q_earth_helio, kind='cubic', axis=0)

# ********************************************************************************************************************* 
def get_earth_pos(ts, heliocentric: bool = False) -> np.array:
    """
    Get position of earth consistent with asteroid data at the specified times (MJDs)
    INPUTS:
        ts: Array of times expressed as MJDs
        heliocentric: flag indicating whether to use heliocentric coordinates;
                      default false uses barycentric coordinates
    """
    # Choose barycentric or heliocentric interpolator
    interpolator = earth_interp_helio if heliocentric else earth_interp_bary

    # Compute interpolated position at desired times
    q_earth = interpolator(ts)
    return q_earth

# ********************************************************************************************************************* 
def make_data_one_file(n0: int, n1: int, heliocentric: bool = False) \
    -> Tuple[Dict[str, np.array], Dict[str, np.array]]:
    """
    Wrap the data in one file of asteroid trajectory data into a TF Dataset
    INPUTS:
        n0: the first asteroid in the file, e.g. 0
        n1: the last asteroid in the file (exclusive), e.g. 1000
        heliocentric: flag indicating whether to use heliocentric coordinates
                      default false uses barycentric coordinates
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
    if heliocentric:
        q_earth -= q_sun
        # v_earth -= v_sun

    # shrink down q and v to slice with asteroid data only; 
    q = q[:, ast_offset:, :]
    v = v[:, ast_offset:, :]
    
    # swap axes for time step and body number; TF needs inputs and outputs to have same number of samples
    # this means that inputs and outputs must first be indexed by asteroid number, then time time step
    q = np.swapaxes(q, 0, 1)
    v = np.swapaxes(v, 0, 1)
    if heliocentric:
        q -= q_sun
        v -= v_sun
    
    # Compute relative displacement to earth
    q_rel = q - q_earth
    # Distance to earth
    r_earth = np.linalg.norm(q_rel, axis=2, keepdims=True)
    # Light time in days from asteroid to earth in days (time units is days)
    light_time = (r_earth / light_speed_au_day)
    # adjusted relative position, accounting for light time; simulation velocity units are AU /day
    dq_lt = light_time * v
    q_rel_lt = q_rel - dq_lt
    # adjusted distance to earth, accounting for light time
    r_lt = np.linalg.norm(q_rel_lt, axis=2, keepdims=True)
    # Direction from earth to asteroid as unit vectors u = (ux, uy, uz)    
    u = q_rel_lt / r_lt

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
        # 'q_earth' : q_earth.astype(dtype),
    }
    
    return inputs, outputs

# ********************************************************************************************************************* 
def make_dataset_pos_file(n0: int, n1: int, include_vel: bool) -> tf.data.Dataset:
    """
    Wrap the data in one file of asteroid trajectory data into a TF Dataset
    INPUTS:
        n0: the first asteroid in the file, e.g. 0
        n1: the last asteroid in the file (exclusive), e.g. 1000
    OUTPUTS:
        ds: a tf.data.Dataset object for this 
    """
    # Load data
    inputs_all, outputs_all = make_data_one_file(n0, n1)
    
    # Use all inputs
    inputs = inputs_all
    
    # Wrap up selected outputs
    outputs ={'q': outputs_all['q']}
    if include_vel:
        outputs['v'] = outputs_all['v']
    
    # Wrap into a dataset
    ds = tf.data.Dataset.from_tensor_slices((inputs, outputs))    
    return ds

# ********************************************************************************************************************* 
def make_dataset_pos_n(n: int, include_vel: bool):
    """Convert the nth file, with asteroids (n-1)*1000 to n*1000, to a tf Dataset"""
    batch_size: int = 1000
    n0 = batch_size*n
    n1 = n0 + batch_size
    return make_dataset_pos_file(n0, n1, include_vel)

# ********************************************************************************************************************* 
def combine_datasets_pos(n0: int, n1: int, include_vel: bool = False, 
                     batch_size: int = 64, progbar: bool = False):
    """Combine datasets in [n0, n1) into one large dataset"""
    # Iteration range for adding to the initial dataset
    ns = range(n0+1, n1)
    ns = tqdm(ns) if progbar else ns
    # Initial dataset
    ds = make_dataset_pos_n(n0, include_vel)
    # Extend initial dataset
    for n in ns:
        try:
            ds_new = make_dataset_pos_n(n, include_vel)
            ds = ds.concatenate(ds_new)
        except:
            pass
    # Batch the dataset
    drop_remainder: bool = True    
    return ds.batch(batch_size, drop_remainder=drop_remainder)

# ********************************************************************************************************************* 
def make_dataset_ast_pos(n0: int, num_files: int, include_vel: bool = False):
    """Create a dataset spanning files [n0, n1); user friendly API for combine_datasets"""
    n1: int = n0 + num_files
    progbar: bool = False
    return combine_datasets_pos(n0=n0, n1=n1, include_vel=include_vel, progbar=progbar)

# ********************************************************************************************************************* 
# Datasets where the output is the asteroid direction from earth (ux, uy, uz)
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_dataset_dir_file(n0: int, n1: int) -> tf.data.Dataset:
    """
    Wrap the data in one file of asteroid trajectory data into a TF Dataset
    INPUTS:
        n0: the first asteroid in the file, e.g. 0
        n1: the last asteroid in the file (exclusive), e.g. 1000
    OUTPUTS:
        ds: a tf.data.Dataset object for this 
    """
    # Load data
    inputs_all, outputs_all = make_data_one_file(n0, n1)
    
    # Use all inputs
    inputs = inputs_all
    
    # Wrap up selected outputs
    outputs ={'u': outputs_all['u']}

    # Wrap into a dataset
    ds = tf.data.Dataset.from_tensor_slices((inputs, outputs))    
    return ds

# ********************************************************************************************************************* 
def make_dataset_dir_n(n: int):
    """Convert the nth file, with asteroids (n-1)*1000 to n*1000, to a tf Dataset with direction from earth"""
    batch_size: int = 1000
    n0 = batch_size*n
    n1 = n0 + batch_size
    return make_dataset_dir_file(n0, n1)

# ********************************************************************************************************************* 
def combine_datasets_dir(n0: int, n1: int, batch_size: int = 64, progbar: bool = False):
    """Combine datasets in [n0, n1) into one large dataset"""
    # Iteration range for adding to the initial dataset
    ns = range(n0+1, n1)
    ns = tqdm(ns) if progbar else ns
    # Initial dataset
    ds = make_dataset_dir_n(n0)
    # Extend initial dataset
    for n in ns:
        try:
            ds_new = make_dataset_dir_n(n)
            ds = ds.concatenate(ds_new)
        except:
            pass
    # Batch the dataset
    drop_remainder: bool = True
    return ds.batch(batch_size, drop_remainder=drop_remainder)

# ********************************************************************************************************************* 
def make_dataset_ast_dir(n0: int, num_files: int):
    """Create a dataset spanning files [n0, n1); user friendly API for combine_datasets"""
    n1: int = n0 + num_files
    progbar: bool = False
    ds = combine_datasets_dir(n0=n0, n1=n1, progbar=progbar)
    return ds

