"""
Interpolated position of earth and sun.

Michael S. Emanuel
2021-05-26
"""

# Core
import scipy, scipy.interpolate
import numpy as np
import pandas as pd

# ********************************************************************************************************************* 
# Data type
dtype_default = np.float64

# ********************************************************************************************************************* 
def get_earth_sun_file(dtype: np.ScalarType = dtype_default):
    """Get the position and velocity of earth and sun consistent with the asteroid data; look up from data file"""
    # Name of the HDF5 archive
    fname_h5: str = '../data/planets/StateVectors_HiRes.h5'
    
    # Extract the two data frames
    df_earth = pd.read_hdf(fname_h5, key='df_earth')
    df_sun  = pd.read_hdf(fname_h5, key='df_earth')

    # The snapshot times as MJDs
    ts = df_earth.mjd.astype(dtype).values
   
    # Columns with the position and velocity
    cols_pos = ['qx', 'qy', 'qz']
    cols_vel = ['vx', 'vy', 'vz']

    # Extract position of earth and sun in Barycentric coordinates
    q_earth = df_earth[cols_pos].astype(dtype).values
    v_earth = df_earth[cols_vel].astype(dtype).values
    q_sun = df_sun[cols_pos].astype(dtype).values
    v_sun = df_sun[cols_vel].astype(dtype).values
    return ts, q_earth, q_sun, v_earth, v_sun

# ********************************************************************************************************************* 
# Create 1D linear interpolator for earth positions; just need one instance for the whole module

# Get position and velocity of earth and sun in barycentric frame at reference dates from file
ts, q_earth, q_sun, v_earth, v_sun = get_earth_sun_file()

# Build interpolator for barycentric position of earth and sun
earth_interp_q = scipy.interpolate.interp1d(x=ts, y=q_earth, kind='cubic', axis=0)
sun_interp_q = scipy.interpolate.interp1d(x=ts, y=q_sun, kind='cubic', axis=0)

# Build interpolator for barycentric velocity of earth and sun
earth_interp_v = scipy.interpolate.interp1d(x=ts, y=v_earth, kind='cubic', axis=0)
sun_interp_v = scipy.interpolate.interp1d(x=ts, y=v_sun, kind='cubic', axis=0)

# ********************************************************************************************************************* 
def get_earth_pos(ts: np.ndarray, dtype = dtype_default) -> np.array:
    """
    Get position of earth consistent with asteroid data at the specified times (MJDs) in barycentric frame
    INPUTS:
        ts: Array of times expressed as MJDs
        dtype: Desired datatype, e.g. np.float64 or np.float32
    """
    # Compute interpolated position at desired times
    q_earth = earth_interp_q(ts).astype(dtype)

    return q_earth

# ********************************************************************************************************************* 
def get_sun_pos(ts: np.ndarray, dtype = dtype_default) -> np.array:
    """
    Get position of sun consistent with asteroid data at the specified times (MJDs) in barycentric frame
    INPUTS:
        ts: Array of times expressed as MJDs
        dtype: Desired datatype, e.g. np.float64 or np.float32
    """
    # Compute interpolated position at desired times
    q_sun = sun_interp_q(ts).astype(dtype)

    return q_sun

# ********************************************************************************************************************* 
def get_earth_vectors(ts: np.ndarray, dtype = dtype_default) -> np.array:
    """
    Get barycentric position and velocity of earth consistent with asteroid data at the specified times (MJDs)
    INPUTS:
        ts:    Array of times expressed as MJDs
        dtype: Desired datatype, e.g. np.float64 or np.float32
    """
    # Compute interpolated position at desired times
    q_earth = earth_interp_q(ts).astype(dtype)
    v_earth = earth_interp_v(ts).astype(dtype)
    return q_earth, v_earth

# ********************************************************************************************************************* 
def get_sun_vectors(ts: np.ndarray, dtype = dtype_default) -> np.array:
    """
    Get barycentric position and velocity of sun consistent with asteroid data at the specified times (MJDs)
    INPUTS:
        ts:    Array of times expressed as MJDs
        dtype: Desired datatype, e.g. np.float64 or np.float32
    """
    # Compute interpolated position at desired times
    q_sun = sun_interp_q(ts).astype(dtype)
    v_sun = sun_interp_v(ts).astype(dtype)
    return q_sun, v_sun
