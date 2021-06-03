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
from scipy.interpolate import CubicSpline, RectBivariateSpline

# Astronomy
import astropy
from astropy.units import au, day

# Local imports
from planets_interp import get_earth_elts, get_sun_vectors
from asteroid_data import load_ast_data
from orbital_element import unpack_elt_df, elt2vec, elt2pos, anomaly_f2M
from orbital_element_test import report_test

# Types
from typing import List, Tuple, Optional

# ********************************************************************************************************************* 
# Speed of light; express this in AU / day
light_speed_au_day = astropy.constants.c.to(au / day).value

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
                   id_col: Optional[str]=None, time_col: Optional[str]=None) -> pd.DataFrame:
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

    # Data to be splined: y axis is the AsteroidID; this isn't *really* splined, it just uniquely identifies a row
    y = df[id_col].values[0:N_row:N_t_in]

    # Set splining order: cubic on x axis (time), linear on dummy AsteroidID axis
    kx=3
    ky=1

    # Construct a splining function for each splined columns
    spline_funcs = []
    for i, col in enumerate(cols_spline):
        zi = df[col].values.reshape((N_obj, N_t_in)).T
        spline_func_i = RectBivariateSpline(x=x, y=y, z=zi, kx=3, ky=1)
        spline_funcs.append(spline_func_i)

    # Wrap the column splines into a single spline function that will output one array of shape (sz, k)
    def spline_func(x: np.ndarray, y: np.ndarray):
        z = np.zeros((x.size, k))
        for i in range(k):
            # z[:, i] = spline_funcs[i](x, y).flatten()
            z[:, i] = spline_funcs[i].ev(x, y)
        return z

    # Return the assembled spline function
    return spline_func

# # ********************************************************************************************************************* 
# def make_spline_df_v1(df: pd.DataFrame, cols_spline: List[str],
#                    id_col: Optional[str]=None, time_col: Optional[str]=None) -> pd.DataFrame:
#     """
#     Build a splining interpolator from a DataFrame
#     INPUTS:
#         df:             DataFrame with key column, time column, and data columns to be splined
#         cols_spline:    List of column names to be splined
#         id_col:         Name of the column identifying entities (e.g. AsteroidID)
#         time_col:       Name of the column with the spline times
#     OUTPUTS:
#         spline_func:    New data frame with key columns and splined data columns
#     """
#     # If the id column was not provided, default to the first column
#     id_col = id_col or df.columns[0]
#     # If the time column was not provided, default to the second column
#     time_col = time_col or df.columns[1]

#     # Get the size of the input DataFrame
#     N_obj, N_t_in = get_df_shape(df=df, id_col=id_col)

#     # Data to be splined: x axis is time; column name in the asteroids DataFrame is mjd
#     x_spline = df[time_col].values[0:N_t_in]

#     # the indexing of y_spline is (ast_num, time_step, data_dim) with shape e.g. (16, 3653, 16)
#     y_spline = df[cols_spline].values.reshape(N_obj, N_t_in, -1)

#     # splining function for the splined columns
#     spline_func = CubicSpline(x=x_spline, y=y_spline, axis=1)

#     return spline_func

# ********************************************************************************************************************* 
def spline_data(df: pd.DataFrame, ts: np.ndarray, ids: np.ndarray, cols_spline: List[str], 
                id_col: Optional[str]=None, time_col: Optional[str]=None) -> pd.DataFrame:
    """
    Build a splining interpolator from a DataFrame
    INPUTS:
        df:             DataFrame with key column, time column, and data columns to be splined
        ts:             Array of times at which splined output is desired
        ids:            Array of ids for which splined output is desired; must match length of ts
        cols_spline:    List of column names to be splined
        id_col:         Name of the column identifying entities (e.g. AsteroidID)
        time_col:       Name of the column with the spline times
    OUTPUTS:
        df_out:         New data frame with key columns and splined data columns
    """
    # If the id column was not provided, default to the first column
    id_col = id_col or df.columns[0]
    # If the time column was not provided, default to the second column
    time_col = time_col or df.columns[1]

    # Delegate to make_spline_df to build an interpolating spline
    spline_func = make_spline_df(df=df, cols_spline=cols_spline, id_col=id_col, time_col=time_col)

    # Splined output
    spline_data = spline_func(x=ts, y=ids)

    # Key fields in the output (splined) DataFrame
    # ts_out = np.tile(ts, N_obj)
    out_dict_keys = {
        id_col: ids,
        time_col: ts
    }

    # The splined asteroid state vectos
    out_dict_splined = {col:spline_data[:,k] for k, col in enumerate(cols_spline)}
    
    # Wrap the output arrays into one dictionary and DataFrame
    out_dict = dict(**out_dict_keys, **out_dict_splined)
    df_out = pd.DataFrame(out_dict)

    return df_out

# ********************************************************************************************************************* 
# Specialize splines for asteroids - either directly spline state vectors or spline elements, then transform them
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def spline_ast_data(df_ast: pd.DataFrame, ts: np.ndarray, ids: np.ndarray, cols_spline: List[str]) -> pd.DataFrame:
    """
    Spline the integrated state vectors of asteroids and earth at the desired times
    INPUTS:
        df_ast:         DataFrame with asteroid position and velocity
        ts:             Array of times at which splined output is desired
        ids:            Array of AsteroidIDs for which splined output is desired
        cols_spline:    List of column names to be splined
    OUTPUTS:
        df_out:         New data frame with key columns and splined data columns
    """

    # Delegate to spline_data
    id_col = 'AsteroidID'
    time_col = 'mjd'
    df_out = spline_data(df=df_ast, ts=ts, ids=ids, cols_spline=cols_spline, id_col=id_col, time_col=time_col)
    return df_out

# ********************************************************************************************************************* 
def spline_ast_vec(vec: pd.DataFrame, ts: np.ndarray, ids: np.ndarray) -> pd.DataFrame:
    """
    Spline the integrated state vectors of asteroids and earth at the desired times
    INPUTS:
        vec:       DataFrame with asteroid state vectors
        ts:        Array of times at which splined output is desired
        ids:       Array of AsteroidIDs for which splined output is desired
    OUTPUTS:
        df_out:    Position & velocity of asteroids in barycentric frame
    """
    # The columns to spline - six cartesian coordinates
    cols_spline = ['qx', 'qy', 'qz', 'vx', 'vy', 'vz']

    # Spline these output columns
    df_out = spline_ast_data(df_ast=vec, ts=ts, ids=ids, cols_spline=cols_spline)

    return df_out

# ********************************************************************************************************************* 
def spline_ast_elt(elt: pd.DataFrame, ts: np.ndarray, ids: np.ndarray) -> pd.DataFrame:
    """
    Spline the integrated orbital elements of asteroids and earth at the desired times.
    Return equivalent state vectors.
    INPUTS:
        elt:      DataFrame with asteroid orbital elements
        ts:        Array of times at which splined output is desired
        ids:       Array of AsteroidIDs for which splined output is desired
    OUTPUTS:
        df_out:    Position & velocity of asteroids in barycentric frame
    """
    # Compute cosine and sine of f for splining
    f_in = elt.f.values
    elt['fx'] = np.cos(f_in)
    elt['fy'] = np.sin(f_in)

    # The columns to spline
    cols_spline = ['a', 'e', 'inc', 'Omega', 'omega', 'fx', 'fy']

    # Spline these output columns
    elt_out = spline_ast_data(df_ast=elt, ts=ts, ids=ids, cols_spline=cols_spline)
    # Compute the splined f using atan2
    elt_out['f'] = np.arctan2(elt_out.fy, elt_out.fx)
    # Compute the splined M from f
    elt_out['M'] = anomaly_f2M(f=elt_out.f, e=elt_out.e)

    # Array of output times including repeated times for each asteroid
    ts_out = elt_out.mjd.values

    # Unpack orbital elements
    a, e, inc, Omega, omega, f = unpack_elt_df(elt_out)

    # Compute q and v in the heliocentric frame
    q_hel, v_hel = elt2vec(a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)

    # Position and velocity of the sun
    q_sun, v_sun = get_sun_vectors(ts=ts_out)

    # The recovered position and velocity of the asteroid in the barycentric frame
    q: np.ndarray = q_hel + q_sun
    v: np.ndarray = v_hel + v_sun

    # Build output state vectors; the two key columns in the output are AsteroidID and mjd
    key_cols = elt_out.columns[0:2]
    # Cartesian coordinate columns
    cols_q = ['qx', 'qy', 'qz']
    cols_v = ['vx', 'vy', 'vz']
    # Need to copy from elts_out to avoid dreaded "setting with copy" warning in Pandas
    vec_out = elt_out[key_cols].copy()
    vec_out[cols_q] = q
    vec_out[cols_v] = v

    return vec_out

# ********************************************************************************************************************* 
def test_ast_spline_elt():
    """Test splining of asteroid orbital elements"""
    # Get test orbital elements from the first 10 asteroids
    df_ast = load_ast_data(n0=1, n1=11)

    # Choose a set of interpolation times that exactly match the ORIGINAL sample times
    # Number of asteroids
    N_ast = np.unique(df_ast.AsteroidID.values).size
    # Number of times per asteroid
    N_t = df_ast.shape[0] // N_ast
    
    # Array of times for output
    ts = df_ast.mjd.values
    # Array of AsteroidIDs for output
    ids = df_ast.AsteroidID.values

    # Down-sample from every 4 days to every 8 days
    interval = 24*60*8
    mask = (df_ast.TimeID % interval == 0)
    df_in = df_ast[mask].copy()

    # Spline vectors directly
    df_vec = spline_ast_vec(vec=df_in, ts=ts, ids=ids)

    # Spline vectors via orbital elements
    df_elt = spline_ast_elt(elt=df_in, ts=ts, ids=ids)

    # Unpack the position from (1) original data (2) interpolated via vectors (3) interpolated via elements
    cols_q = ['qx', 'qy', 'qz']
    q0 = df_ast[cols_q].values
    q1 = df_vec[cols_q].values
    q2 = df_elt[cols_q].values

    # Unpack the velocity from same 3 sources
    # cols_v = ['vx', 'vy', 'vz']
    # v0 = df_elt[cols_v].values
    # v1 = df_vec[cols_v].values
    # v1 = df_elt[cols_v].values

    # Reconstruction error
    dq1: np.ndarray = q1 - q0
    dq2: np.ndarray = q2 - q0
    # dv1: np.ndarray = v1 - v0
    # dv1: np.ndarray = v2 - v0
    err_q1: np.ndarray = np.sqrt(np.sum(np.square(dq1), axis=-1))
    err_q2: np.ndarray = np.sqrt(np.sum(np.square(dq2), axis=-1))   
    # err_v1: np.ndarray = np.sqrt(np.sum(np.square(dv1), axis=-1))
    # err_v2: np.ndarray = np.sqrt(np.sum(np.square(dv2), axis=-1))

    # Report the results
    print(f'Splined position using orbital elements; first 10 asteroids.')
    print(f'Original data interval is 4 days.  Spline built using interval of 8 days.')
    is_ok_q1 = report_test(err=err_q1, test_name='ast_spline_elt (position q, spline on vectors)', thresh=1.0E-4)
    is_ok_q2 = report_test(err=err_q2, test_name='ast_spline_elt (position q, spline on elements)', thresh=1.0E-4)
    return (is_ok_q1 and is_ok_q2)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    print(f'Running test suite on asteroid_spline.py...')
    test_ast_spline_elt()
