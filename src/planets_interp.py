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
def load_earth_sun_vectors(interval: str, dtype: np.ScalarType = dtype_default):
    """
    Load the position and velocity of earth and sun consistent with the asteroid data from an HDF5 data file
    INPUTS:
        interval: One of 'M' (1 minute), 'H' (1 hour), or 'D' (1 day); interval vectors are saved with
        dtype:    The data type used; should be one of np.float64 or np.float32
    """

    # Check interval was valid
    if interval not in ('M', 'H', 'D'):
        raise ValueError("Bad interval! Must be one of 'M' (1 minute), 'H' (1 hour) or 'D' (1 day)")

    # Name of the HDF5 archive
    suffix_tbl = {
        'M': 'Minute',
        'H': 'Hour', 
        'D': 'Day',
    }
    suffix = suffix_tbl[interval]
    fname_h5: str = f'../data/planets/StateVectors_{suffix}.h5'

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
def make_earth_sun_interpolators(interval: str, dtype = dtype_default) -> np.array:
    """
    Construct interpolators for the earth and sun at the desired interval: minute (M), hour (H), or day (D)
    INPUTS:
        interval: One of 'M' (1 minute), 'H' (1 hour), or 'D' (1 day); interval vectors are saved with
        dtype: Desired datatype, e.g. np.float64 or np.float32
    OUTPUTS:
        tuple of four scipy.interpolate instances:
        earth_interp_q, sun_interp_q, earth_interp_v, sun_interp_v
    """
    # Load position and velocity of earth and sun in barycentric frame at hourly intervals from file
    ts, q_earth, q_sun, v_earth, v_sun = load_earth_sun_vectors(interval=interval)

    # Build interpolator for barycentric position of earth and sun
    earth_interp_q = scipy.interpolate.interp1d(x=ts, y=q_earth, kind='cubic', axis=0)
    sun_interp_q = scipy.interpolate.interp1d(x=ts, y=q_sun, kind='cubic', axis=0)

    # Build interpolator for barycentric velocity of earth and sun
    earth_interp_v = scipy.interpolate.interp1d(x=ts, y=v_earth, kind='cubic', axis=0)
    sun_interp_v = scipy.interpolate.interp1d(x=ts, y=v_sun, kind='cubic', axis=0)

    return earth_interp_q, sun_interp_q, earth_interp_v, sun_interp_v

# ********************************************************************************************************************* 
# Create 1D linear interpolator for earth positions; just need one instance for the whole module

# Load position and velocity of earth and sun in barycentric frame at hourly intervals from file
ts, q_earth, q_sun, v_earth, v_sun = load_earth_sun_vectors(interval='H')

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
        ts: Array of times at which to interpolate expressed as MJDs
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
