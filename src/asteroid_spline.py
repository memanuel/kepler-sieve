"""
Load integrated asteroid trajectories as Pandas DataFrames.
Add earth or sun position / state vectors to asteroid DataFrame.

Michael S. Emanuel
Sat Sep 21 10:38:38 2019
"""

# Core
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

# Astronomy
import astropy
from astropy.units import au, day

# Local imports
from planets_interp import get_earth_elts, get_sun_vectors
from asteroid_data import load_ast_data
from orbital_element import unpack_elt_df, elt2vec, elt2pos, anomaly_M2f
from orbital_element_test import report_test

# Types
from typing import List

# ********************************************************************************************************************* 
# Speed of light; express this in AU / day
light_speed_au_day = astropy.constants.c.to(au / day).value

# ********************************************************************************************************************* 
def spline_ast_data(df_ast: pd.DataFrame, ts: np.ndarray, cols_spline: List[str]) -> pd.DataFrame:
    """
    Spline the integrated state vectors of asteroids and earth at the desired times
    INPUTS:
        df_ast:         DataFrame with asteroid position and velocity
        ts:             Array of times at which splined output is desired
        cols_spline:    List of column names to be splined
    OUTPUTS:
        df_out:         New data frame with key columns and splined data columns
    """

    # Distinct asteroids or elements; get the name of the relevant column and array of distinct values
    id_col = 'AsteroidID' if 'AsteroidID' in df_ast.columns else 'ElementID'
    id_val_in = df_ast[id_col].values
    id_val_unq = np.unique(id_val_in)

    # Number of asteroids or elements
    N_ast = id_val_unq.size
    # Number of times in splined output
    N_t_out = ts.shape[0]

    # Number of rows in input and splined output
    N_row_in = df_ast.shape[0]
    N_row_out = N_ast * N_t_out

    # Calculated number of times in input
    N_t_in = N_row_in // N_ast

    # Data to be splined: x axis is time; column name in the asteroids DataFrame is mjd
    x_spline = df_ast.mjd.values[0:N_t_in]

    # the indexing of y_spline is (ast_num, time_step, data_dim) with shape e.g. (16, 3653, 16)
    y_spline = df_ast[cols_spline].values.reshape(N_ast, N_t_in, -1)

    # splining function for the splined columns
    spline_func = CubicSpline(x=x_spline, y=y_spline, axis=1)

    # splined output
    spline_data = spline_func(ts)

    # ID column (asteroid_num or element_id) corresponding to splined output
    id_val = np.repeat(id_val_unq, N_t_out)

    # Key fields in the output (splined) DataFrame
    mjd_out = np.tile(ts, N_ast)
    out_dict_keys = {
        id_col: id_val,
        'mjd': mjd_out
    }

    # The splined asteroid state vectos
    out_dict_splined = {col:spline_data[:,:,k].reshape(N_row_out) for k, col in enumerate(cols_spline)}
    
    # Wrap the output arrays into one dictionary and DataFrame
    out_dict = dict(**out_dict_keys, **out_dict_splined)
    df_out = pd.DataFrame(out_dict)

    return df_out

# ********************************************************************************************************************* 
def spline_ast_vec(vec: pd.DataFrame, ts: np.ndarray) -> pd.DataFrame:
    """
    Spline the integrated state vectors of asteroids and earth at the desired times
    INPUTS:
        vec:       DataFrame with asteroid state vectors
        ts:        Array of times at which splined output is desired
    OUTPUTS:
        df_out:    Position & velocity of asteroids in barycentric frame
    """
    # The columns to spline - six cartesian coordinates
    cols_spline = ['qx', 'qy', 'qz', 'vx', 'vy', 'vz']

    # Spline these output columns
    df_out = spline_ast_data(df_ast=vec, ts=ts, cols_spline=cols_spline)

    return df_out

# ********************************************************************************************************************* 
def spline_ast_elt(elt: pd.DataFrame, ts: np.ndarray) -> pd.DataFrame:
    """
    Spline the integrated orbital elements of asteroids and earth at the desired times.
    Return equivalent state vectors.
    INPUTS:
        elt:      DataFrame with asteroid orbital elements
        ts:        Array of times at which splined output is desired
    OUTPUTS:
        df_out:    Position & velocity of asteroids in barycentric frame
    """
    # The columns to spline
    cols_spline = ['a', 'e', 'inc', 'Omega', 'omega', 'Mx', 'My']

    # Compute cosine and sine of M for splining
    M_in = elt.M.values
    elt['Mx'] = np.cos(M_in)
    elt['My'] = np.sin(M_in)

    # Spline these output columns
    elt_out = spline_ast_data(df_ast=elt, ts=ts, cols_spline=cols_spline)
    # Compute the splined M using atan2
    elt_out['M'] = np.arctan2(elt_out.My, elt_out.Mx)
    # Compute the splined f from M
    elt_out['f'] = anomaly_M2f(M=elt_out.M, e=elt_out.e)

    # Array of output times including repeated times for each asteroid
    ts_out = elt_out.mjd.values

    # Unpack orbital elements
    a, e, inc, Omega, omega, f = unpack_elt_df(elt_out)

    # Compute q and v in the heliocentric frame
    q_hel, v_hel = elt2vec(a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)

    # Position and velocity of the sun
    q_sun, v_sun = get_sun_vectors(ts=ts_out)

    # The recovered position and velocity of the asteroid in the barycentric frame
    q: np.array = q_hel + q_sun
    v: np.array = v_hel + v_sun

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
    # Array of distinct times
    ts = df_ast.mjd.values[0:N_t]

    # Down-sample from every 4 days to every 8 days
    interval = 24*60*8
    mask = (df_ast.TimeID % interval == 0)
    df_in = df_ast[mask].copy()

    # Spline vectors directly
    df_vec = spline_ast_vec(vec=df_in, ts=ts)

    # Spline vectors via orbital elements
    df_elt = spline_ast_elt(elt=df_in, ts=ts)

    # Unpack the position from (1) original data (2) interpolated via vectors (3) interpolated via elements
    cols_q = ['qx', 'qy', 'qz']
    q0 = df_ast[cols_q].values
    q1 = df_vec[cols_q].values
    q2 = df_elt[cols_q].values

    # Unpack the velocity from same 3 sources
    cols_v = ['vx', 'vy', 'vz']
    v0 = df_elt[cols_v].values
    v1 = df_vec[cols_v].values
    v1 = df_elt[cols_v].values

    # Reconstruction error
    dq1: np.array = q1 - q0
    dq2: np.array = q2 - q0
    # dv1: np.array = v1 - v0
    # dv1: np.array = v2 - v0
    err_q1: np.array = np.sqrt(np.sum(np.square(dq1), axis=-1))
    err_q2: np.array = np.sqrt(np.sum(np.square(dq2), axis=-1))   
    # err_v1: np.array = np.sqrt(np.sum(np.square(dv), axis=-1))

    # Report the results
    print(f'Splined position using orbital elements; first 10 asteroids.')
    print(f'Original data interval is 4 days.  Spline built using interval of 8 days.')
    is_ok_q1 = report_test(err=err_q1, test_name='ast_spline_elt (position q, spline on vectors)', thresh=1.0E-4)
    is_ok_q2 = report_test(err=err_q2, test_name='ast_spline_elt (position q, spline on elements)', thresh=1.0E-4)
    # is_ok_v = report_test(err=err_v, test_name='ast_spline_elt (velocity v)', thresh=1.0E-6)
    return (is_ok_q1 and is_ok_q2)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    print(f'Running test suite on asteroid_spline.py...')
    test_ast_spline_elt()
