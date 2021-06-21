"""
Load integrated asteroid trajectories as Pandas DataFrames.
Add earth or sun position / state vectors to asteroid DataFrame.

Functions in this module:
get_df_shape(df, id_col)
make_spline_df(df, cols_spline, id_col, time_col)
spline_data(df, ts, cols_spline, id_col, time_col)
spline_ast_data(df_ast, ts, cols_spline)
spline_ast_vec(vec, ts)
spline_ast_elt(elt, ts)
test_ast_spline_elt()

Michael S. Emanuel
Sat Sep 21 10:38:38 2019
"""

# Core
import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline

# Local imports
from planets_interp import get_sun_vectors, get_sun_pos
from orbital_element import unpack_elt_np
from asteroid_data import load_ast_dir
from asteroid_element import get_ast_elts_ts, get_ast_data_ts
from orbital_element import elt2vec, elt2pos, anomaly_M2f

# Types
from typing import List, Tuple, Optional, Callable
spline_type_ast = Callable[[np.ndarray, np.ndarray], np.ndarray]

# ********************************************************************************************************************* 
# Back end for building splines from DataFrames - general
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def get_df_shape(df: pd.DataFrame, id_col: str) -> Tuple[int, int]:
    """
    Get the shape of a DataFrame used to build a spline
    INPUTS:
        df:             DataFrame with key column, time column, and data columns to be splined
        id_col:         Name of the column identifying entities (e.g. AsteroidID)
    OUTPUTS:
        N_obj:          Number of distinct objects in input DataFrame (e.g. asteroids)
        N_t_in:         Number of time snaps for each object
    """
    # The ID values and the unique set of values
    id_val_in = df[id_col].values
    id_val_unq = np.unique(id_val_in)

    # Number of objects (e.g. asteroids)
    N_obj = id_val_unq.size

    # Number of rows in input and splined output
    N_row_in = df.shape[0]

    # Calculated number of times in input
    N_t_in = N_row_in // N_obj

    return N_obj, N_t_in

# ********************************************************************************************************************* 
def make_spline_df(df: pd.DataFrame, cols_spline: List[str],
                   time_col: Optional[str]=None, id_col: Optional[str]=None) -> spline_type_ast:
    """
    Build a splining interpolator from a DataFrame
    INPUTS:
        df:             DataFrame with key column, time column, and data columns to be splined
        cols_spline:    List of column names to be splined
        id_col:         Name of the column identifying entities (e.g. AsteroidID)
        time_col:       Name of the column with the spline times
    OUTPUTS:
        spline_func:    An interpolation function that accepts array inputs x and y, e.g.
                        z = spline_func(x, y)
                        x and y should be 1D arrays of the same size, sz.
                        z is an array of shape (sz, k) where k is the number of columns
    """
    # If the id column was not provided, default to the first column
    id_col = id_col or df.columns[0]
    # If the time column was not provided, default to the second column
    time_col = time_col or df.columns[1]

    # Get the size of the input DataFrame
    N_obj, N_t_in = get_df_shape(df=df, id_col=id_col)
    N_row = N_obj * N_t_in

    # Dimension of the output is the number of columns to spline
    k = len(cols_spline)

    # Data to be splined: x axis is time; column name in the asteroids DataFrame is mjd
    # Note that this is assumed to be shared by all the different IDs
    x = df[time_col].values[0:N_t_in]

    # Data to be splined: y axis is the ID, e.g. AsteroidID; 
    # This column isn't *really* splined, it just uniquely identifies a row
    y = df[id_col].values[0:N_row:N_t_in]

    # Set splining order: cubic on x axis (time), linear on dummy AsteroidID axis
    kx=3
    ky=1

    # Construct a splining function for each splined column
    spline_funcs = []
    for col in cols_spline:
        zi = df[col].values.reshape((N_obj, N_t_in)).T
        spline_func_i = RectBivariateSpline(x=x, y=y, z=zi, kx=kx, ky=ky)
        spline_funcs.append(spline_func_i)

    # Wrap the column splines into a single spline function that will output one array of shape (sz, k)
    def spline_func(x: np.ndarray, y: np.ndarray):
        z = np.zeros(shape=(x.size, k))
        for i in range(k):
            z[:, i] = spline_funcs[i].ev(x, y)
        return z

    # Return the assembled spline function
    return spline_func

# ********************************************************************************************************************* 
# End to end functions to build splines for 
# (1) asteroid orbital elements (2) asteroid positions (3) asteroid state vectors
# where the positions / vectors are computed via an intermediate spline on the orbital elements
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def pad_dates(mjd0: int, mjd1: int):
    """Extend dates on either end of a spline to ensure it can handle all dates between mjd0 and mjd1."""
    pad: int = 32
    mjd0 = mjd0 - pad
    mjd1 = mjd1 + pad
    return mjd0, mjd1

# ********************************************************************************************************************* 
def make_spline_ast_elt(n0: int, n1: int, mjd0: int, mjd1: int) -> spline_type_ast:
    """
    Build a splining interpolator that returns splined orbital elements for asteroids
    Assumes that M is monotonic (e.g. winding number is incorporated)
    INPUTS:
        n0:         First asteroid number to return (inclusive)
        n1:         Last asteroid number to return (exclusive)
        mjd0:       First date on which to return data (inclusive)
        mjd1:       Last date on which to return data (exclusive)
    OUTPUTS:
        spline_elt: An interpolation function that accepts array inputs ts and asteroid_id, e.g.
                    elt = spline_elt(ts, element_id)
    """
    # Pad dates
    mjd0, mjd1 = pad_dates(mjd0=mjd0, mjd1=mjd1)
    # Get the orbital elements for these asteroids; these are the inputs to build the spline
    elt_in = get_ast_elts_ts(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Delegate to make_spline_df
    cols_spline = ['a', 'e', 'inc', 'Omega', 'omega', 'f', 'M']
    time_col = 'mjd'
    id_col = 'AsteroidID'
    spline_elt_raw = make_spline_df(df=elt_in, cols_spline=cols_spline, time_col=time_col, id_col=id_col)

    # Wrap up a function that returns the position vector
    def spline_elt(ts: np.ndarray, asteroid_id: np.ndarray):
        """Splined elements as a function of time and asteroid_id"""
        # Caclulate the splined orbital elements for the input times and asteroid_ids
        elt_out = spline_elt_raw(ts, asteroid_id)
        # Unpack the elements
        a, e, inc, Omega, omega, f, M = unpack_elt_np(elt_out)
        # Recover f from M and e
        f = anomaly_M2f(M=M, e=e)
        # Overwrite the f column in the splined element array; avoid allocating and copying a new array
        elt_out[:, 5] = f
        # elt_out[:, 6] = M
        # Return the splined orbital elements with corrected f
        return elt_out

    # Return the function that calculates orbital elements
    return spline_elt

# ********************************************************************************************************************* 
def make_spline_ast_elt_no_winding(n0: int, n1: int, mjd0: int, mjd1: int) -> spline_type_ast:
    """
    Build a splining interpolator that returns splined orbital elements for asteroids.
    Version suitable to use when M is not monotonic (e.g. winding number is not incorporated)
    INPUTS:
        n0:         First asteroid number to return (inclusive)
        n1:         Last asteroid number to return (exclusive)
        mjd0:       First date on which to return data (inclusive)
        mjd1:       Last date on which to return data (exclusive)
    OUTPUTS:
        spline_elt: An interpolation function that accepts array inputs ts and asteroid_id, e.g.
                    elt = spline_elt(ts, element_id)
    """
    # Pad dates
    mjd0, mjd1 = pad_dates(mjd0=mjd0, mjd1=mjd1)
    # Get the orbital elements for these asteroids; these are the inputs to build the spline
    elt_in = get_ast_elts_ts(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)
    # Add cosine and sine of M
    elt_in['Mx'] = np.cos(elt_in.M)
    elt_in['My'] = np.sin(elt_in.M)

    # Delegate to make_spline_df
    cols_spline = ['a', 'e', 'inc', 'Omega', 'omega', 'Mx', 'My']
    time_col = 'mjd'
    id_col = 'AsteroidID'
    spline_elt_raw = make_spline_df(df=elt_in, cols_spline=cols_spline, time_col=time_col, id_col=id_col)

    # Wrap up a function that returns the position vector
    def spline_elt(ts: np.ndarray, asteroid_id: np.ndarray):
        """Splineed elements as a function of time and asteroid_id"""
        # Caclulate the splined orbital elements for the input times and asteroid_ids
        elt_out = spline_elt_raw(ts, asteroid_id)
        # Unpack the elements; abuse API unpack_elt; really it just unpacks seven arrays in order
        a, e, inc, Omega, omega, Mx, My = unpack_elt_np(elt_out)
        # Recover M from Mx and My
        M = np.arctan2(My, Mx)
        # Recover f from M and e
        f = anomaly_M2f(M=M, e=e)
        # Overwrite the f and M columns in the splined element array; avoid allocating and copying a new array
        elt_out[:, 5] = f
        elt_out[:, 6] = M
        # Return the splined orbital elements with corrected f and M
        return elt_out

    # Return the function that calculates orbital elements
    return spline_elt

# ********************************************************************************************************************* 
def make_spline_ast_pos(n0: int, n1: int, mjd0: int, mjd1: int) -> spline_type_ast:
    """
    Build a splining interpolator that returns asteroid positions by way of splined orbital elements
    INPUTS:
        n0:         First asteroid number to return (inclusive)
        n1:         Last asteroid number to return (exclusive)
        mjd0:       First date on which to return data (inclusive)
        mjd1:       Last date on which to return data (exclusive)
    OUTPUTS:
        spline_q:   An interpolation function that accepts array inputs ts and asteroid_id, e.g.
                    q_ast = spline_elt(ts, element_id)
    """
    # Delegate to make_spline_ast_elt to build an interpolator for the orbital elements
    spline_elt = make_spline_ast_elt(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Wrap up a function that returns the position vector
    def spline_q(ts: np.ndarray, asteroid_id: np.ndarray):
        """Splineed position as a function of time and asteroid_id"""
        # Caclulate the splined orbital elements for the input times and asteroid_ids
        elt = spline_elt(ts, asteroid_id)
        # Unpack the elements
        a, e, inc, Omega, omega, f, M = unpack_elt_np(elt)        
        # Compute q in the heliocentric frame
        q_hel: np.ndarray = elt2pos(a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)
        # Position  of the sun
        q_sun: np.ndarray = get_sun_pos(ts=ts)
        # The position of the asteroid in the barycentric frame
        q: np.ndarray = q_hel + q_sun
        return q

    # Return the function that calculates positions
    return spline_q

# ********************************************************************************************************************* 
def make_spline_ast_vec(n0: int, n1: int, mjd0: int, mjd1: int) -> spline_type_ast:
    """
    Build a splining interpolator that returns asteroid state vectors by way of splined orbital elements
    INPUTS:
        n0:         First asteroid number to return (inclusive)
        n1:         Last asteroid number to return (exclusive)
        mjd0:       First date on which to return data (inclusive)
        mjd1:       Last date on which to return data (exclusive)
    OUTPUTS:
        spline_vec: An interpolation function that accepts array inputs ts and asteroid_id, e.g.
                    q_ast, v_ast = spline_elt(ts, element_id)
    """
    # Delegate to make_spline_ast_elt to build an interpolator for the orbital elements
    spline_elt = make_spline_ast_elt(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)
    # Wrap up a function that returns the position vector
    def spline_vec(ts: np.ndarray, asteroid_id: np.ndarray):
        """Splineed position as a function of time and asteroid_id"""
        # Caclulate the splined orbital elements for the input times and asteroid_ids
        elt = spline_elt(ts, asteroid_id)
        # Unpack the elements
        a, e, inc, Omega, omega, f, M = unpack_elt_np(elt)        
        # Compute q, v in the heliocentric frame
        q_hel, v_hel = elt2vec(a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)
        # Vectors of the sun
        q_sun, v_sun = get_sun_vectors(ts=ts)
        # The position of the asteroid in the barycentric frame
        q: np.ndarray = q_hel + q_sun
        v: np.ndarray = v_hel + v_sun
        return q, v

    # Return the function that calculates state vectors
    return spline_vec

# ********************************************************************************************************************* 
def make_spline_ast_pos_direct(n0: int, n1: int, mjd0: int, mjd1: int) -> spline_type_ast:
    """
    Build a splining interpolator that returns splined position directly.
    INPUTS:
        n0:         First asteroid number to return (inclusive)
        n1:         Last asteroid number to return (exclusive)
        mjd0:       First date on which to return data (inclusive)
        mjd1:       Last date on which to return data (exclusive)
    OUTPUTS:
        spline_q:   An interpolation function that accepts array inputs ts and asteroid_id, e.g.
                    q_ast = spline_elt(ts, element_id)
    """
    # Pad dates
    mjd0, mjd1 = pad_dates(mjd0=mjd0, mjd1=mjd1)
    # Get the orbital elements for these asteroids; these are the inputs to build the spline
    df_in = get_ast_data_ts(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Delegate to make_spline_df
    cols_spline = ['qx', 'qy', 'qz']
    time_col = 'mjd'
    id_col = 'AsteroidID'
    spline_pos_np = make_spline_df(df=df_in, cols_spline=cols_spline, time_col=time_col, id_col=id_col)

    # Wrap up a function that returns the position vector with named arguments
    def spline_pos(ts: np.ndarray, asteroid_id: np.ndarray):
        """Splined elements as a function of time and asteroid_id"""
        return spline_pos_np(ts, asteroid_id)

    # Return the function that calculates position directly
    return spline_pos

# ********************************************************************************************************************* 
def make_spline_ast_dir(n0: int, n1: int, mjd0: int, mjd1: int) -> spline_type_ast:
    """
    Build a splining interpolator that returns splined direction.
    INPUTS:
        n0:         First asteroid number to return (inclusive)
        n1:         Last asteroid number to return (exclusive)
        mjd0:       First date on which to return data (inclusive)
        mjd1:       Last date on which to return data (exclusive)
    OUTPUTS:
        spline_u:   An interpolation function that accepts array inputs ts and asteroid_id, e.g.
                    u = spline_elt(ts, element_id)
    """
    # Pad dates
    mjd0, mjd1 = pad_dates(mjd0=mjd0, mjd1=mjd1)
    # Get the orbital elements for these asteroids; these are the inputs to build the spline
    df_in = load_ast_dir(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Delegate to make_spline_df
    cols_spline = ['ux', 'uy', 'uz', 'LightTime']
    time_col = 'tObs'
    id_col = 'AsteroidID'
    spline_u_np = make_spline_df(df=df_in, cols_spline=cols_spline, time_col=time_col, id_col=id_col)

    # Wrap up a function that returns the direction vector with named arguments
    def spline_u(ts: np.ndarray, asteroid_id: np.ndarray):
        """Splined direction as a function of time and asteroid_id"""
        # Evaluate the spline
        spline_out = spline_u_np(ts, asteroid_id)
        # Unpack the direction and light time
        u = spline_out[:, 0:3]
        light_time = spline_out[:, 3]
        return u, light_time

    # Return the function that calculates direction directly
    return spline_u
