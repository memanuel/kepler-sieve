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

# ********************************************************************************************************************* 
def make_ragged_tensors(t_np: np.array, u_np: np.array, element_id_np: np.array, batch_size: int):
    """
    Convert t, u, element_id Numpy arrays into ragged tensors
    INPUTS:
        t_np: A numpy array with observation times as MJDs; size (N,)
        u_np: A numpy array with directions as MJDs; size (N,3)
        element_id_np: A numpy array with the element_id's; size (N,)
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
    # Ragged tensors for direction u and element_id
    u = tf.RaggedTensor.from_value_rowids(values=u_np, value_rowids=inv_idx)
    element_id = tf.RaggedTensor.from_value_rowids(values=element_id_np, value_rowids=value_rowids)

    # Return the tensors for t, u, element_id
    return t, u, element_id
    
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
    # ast_num_np = ztf.nearest_ast_num.values
    element_id_np = ztf.element_id.values

    # Sort arrays in increasing order of observation time; otherwise make_ragged_tensors won't work
    sort_idx = np.argsort(t_np)
    t_np = t_np[sort_idx]
    u_np = u_np[sort_idx]
    # ast_num_np = ast_num_np[sort_idx]
    element_id_np = element_id_np[sort_idx]

    # Convert to tensors (regular for t, ragged for u and element_id)
    # t, u_r, ast_num_r = make_ragged_tensors(t_np=t_np, u_np=u_np, element_id_np=element_id_np, batch_size=batch_size)
    t, u_r, element_id_r = make_ragged_tensors(t_np=t_np, u_np=u_np, element_id_np=element_id_np, batch_size=batch_size)

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
    
    # Pad element_id_r into a regular tensor; use default -1 to indicate padded observation
    # ast_num = ast_num_r.to_tensor(default_value=-1)
    element_id = element_id_r.to_tensor(default_value=-1)
    
    # Pad ast_num with entries for the last batch
    # pad_shape_ast_num = [[0,num_pad],[0,0]]
    # ast_num = tf.pad(ast_num, paddings=pad_shape_ast_num, mode='CONSTANT', constant_values=-1)
    pad_shape_elt_id = [[0,num_pad],[0,0]]
    element_id = tf.pad(element_id, paddings=pad_shape_elt_id, mode='CONSTANT', constant_values=-1)
    
    # Wrap into tensorflow Dataset
    inputs = {
        't': t,
        'idx': idx,
        'row_len': row_len,
        'u_obs': u, 
    }
    outputs = {
        # 'ast_num': ast_num,
        'element_id': element_id,
    }
    ds = tf.data.Dataset.from_tensor_slices((inputs, outputs))
    
    # Batch the dataset
    drop_remainder: bool = False
    buffer_size = 16384
    ds = ds.batch(batch_size, drop_remainder=drop_remainder).shuffle(buffer_size).repeat()
    # Return the dataset as well as tensors with the time snapshots and row lengths
    return ds, t, row_len
    
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
# Build TensorFlow data sets with directions of asteroids - using spline_ast_vec_dir
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_dataset_ast_dir_spline(n0: int, n1: int, site_name: str = 'geocenter', 
                                N_t: int = 1000, batch_size: int = 64) -> tf.data.Dataset:
    """
    Create a tf.Dataset for testing asteroid direction model.
    Data comes from output of spline_ast_vec_dir in asteroid_dataframe.
    INPUTS:
        n0: First asteroid number; inclusive, e.g. 1
        n1: Last asteroid number; exclusive, e.g. 65
        site_name: Name of site, e.g. 'palomar'
        N_t: Number of times to sample, e.g. 1000
        batch_size: batch_size for this dataset
    """
    # DataFrame of asteroid snapshots
    ast_elt_all = load_ast_elt()    

    # Elements of selected asteroids
    # ast_elt = ast_elt_all.loc[ast_nums]
    mask = (n0 <= ast_elt_all.Num) & (ast_elt_all.Num < n1)
    ast_elt = ast_elt_all[mask]
    N_ast = ast_elt.shape[0]

    # Range of times for sampling
    dt0 = datetime(2000, 1, 1)
    dt1 = datetime(2040, 1, 1)
    mjd0 = datetime_to_mjd(dt0)
    mjd1 = datetime_to_mjd(dt1)

    # Data type for this dataset
    dtype = np.float32

    # Select a random subset of times; these must be sorted
    np.random.seed(42)
    ts = np.sort(np.random.uniform(low=mjd0, high=mjd1, size=N_t).astype(dtype=dtype))

    # dict with inputs
    inputs = {
        'a': ast_elt.a.values.astype(dtype),
        'e': ast_elt.e.values.astype(dtype),
        'inc': ast_elt.inc.values.astype(dtype),
        'Omega': ast_elt.Omega.values.astype(dtype),
        'omega': ast_elt.omega.values.astype(dtype),
        'f': ast_elt.f.values.astype(dtype),
        'epoch': ast_elt.epoch_mjd.values.astype(dtype),
        'asteroid_num': ast_elt.Num.values.astype(np.int32),
        'ts': np.tile(ts, reps=(N_ast,1,)),
    }

    # Build splined direction at selected times
    df_ast, df_earth, df_dir = spline_ast_vec_dir(n0=n0, n1=n1, mjd=ts, site_name=site_name)

    # Extract outputs: u, r
    cols_u = ['ux', 'uy', 'uz']
    u_flat = df_dir[cols_u].values
    r_flat = df_dir.delta.values

    # Reshape outputs into rectangular arrays; shape is (ast_num, time_idx, space)
    u = u_flat.reshape((N_ast, N_t, space_dims))
    r = r_flat.reshape((N_ast, N_t))

    # Outputs as Python dict
    outputs = {
        'u': u, 
        'r': r,
    }

    # Wrap into a dataset
    ds = tf.data.Dataset.from_tensor_slices((inputs, outputs))

    # Batch the dataset and return it
    drop_remainder = True
    ds = ds.batch(batch_size=batch_size, drop_remainder=drop_remainder)

    return ds

# ********************************************************************************************************************* 
# Build TensorFlow data sets with positions and velocities of asteroids - manually from integration outputs
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_data_one_file(n0: int, n1: int) -> Tuple[Dict[str, np.array], Dict[str, np.array]]:
    """
    Wrap the data in one file of asteroid trajectory data into a TF Dataset
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
    """
    # Load data
    inputs_all, outputs_all = make_data_one_file(n0, n1)
    
    # Use all inputs
    inputs = inputs_all
    
    # Wrap up selected outputs
    outputs ={
        'u': outputs_all['u'],
        'r': outputs_all['r']
    }

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
