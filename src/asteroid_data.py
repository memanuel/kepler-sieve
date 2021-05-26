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

# Utility
# from tqdm.auto import tqdm

# Local imports
from utils import range_inc
from db_utils import sp2df
from planets_interp import get_earth_pos, get_earth_vectors, get_sun_vectors
from ra_dec import qv2dir, dir2radec, calc_topos, astrometric_dir, direction_diff

# Type names
from typing import Tuple

# ********************************************************************************************************************* 
# Speed of light; express this in AU / day
light_speed_au_day = astropy.constants.c.to(au / day).value

# Last date of calculated data; filtering for mjd <= this date will include everything
mjd_last: np.float64 = np.float64(2.0**20)

# ********************************************************************************************************************* 
def load_ast_vectors(n0: int, n1: int, mjd0: float=0.0, mjd1: float=mjd_last) -> pd.DataFrame:
    """
    Load the MSE asteroid integrations for this range of asteroids; returns just state vectors.
    INPUTS:
        n0:  First asteroid to load, e.g. 0
        n1:  Last asteroid to load, (exclusive) e.g. 1000
        mjd0:  Start modfified julian date used to filter output.
        mjd1:  Last modified julian date used to filter output.
               Default for mjd0 and mjd1 is None; then return all available time steps
    OUTPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame
    """

    # Get asteroid data (state vectors and orbital elements) from database
    sp_name = 'KS.GetAsteroidVectors'
    params = {
        'n0': n0,
        'n1': n1,
        'mjd0': mjd0,
        'mjd1': mjd1,
    }
    df_ast = sp2df(sp_name, params)

    return df_ast

# ********************************************************************************************************************* 
def load_ast_elements(n0: int, n1: int, mjd0: float=0.0, mjd1: float=mjd_last) -> pd.DataFrame:
    """
    Load the MSE asteroid integrations for this range of asteroids; returns just orbital elements.
    INPUTS:
        n0:  First asteroid to load, e.g. 0
        n1:  Last asteroid to load, (exclusive) e.g. 1000
        mjd0:  Start modfified julian date used to filter output.
        mjd1:  Last modified julian date used to filter output.
               Default for mjd0 and mjd1 is None; then return all available time steps
    OUTPUTS:
        df_ast:   Heliocentric orbital elements
    """

    # Get asteroid data (state vectors and orbital elements) from database
    sp_name = 'KS.GetAsteroidElements'
    params = {
        'n0': n0,
        'n1': n1,
        'mjd0': mjd0,
        'mjd1': mjd1,
    }
    df_ast = sp2df(sp_name, params)

    return df_ast

# ********************************************************************************************************************* 
def load_ast_data(n0: int, n1: int, mjd0: float=0.0, mjd1: float=mjd_last) -> \
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load the MSE asteroid integrations for this range of asteroids; returns both state vectors and orbital elements.
    INPUTS:
        n0:  First asteroid to load, e.g. 0
        n1:  Last asteroid to load, (exclusive) e.g. 1000
        mjd0:  Start modfified julian date used to filter output.
        mjd1:  Last modified julian date used to filter output.
               Default for mjd0 and mjd1 is None; then return all available time steps
    OUTPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
    """

    # Get asteroid data (state vectors and orbital elements) from database
    sp_name = 'KS.GetAsteroidData'
    params = {
        'n0': n0,
        'n1': n1,
        'mjd0': mjd0,
        'mjd1': mjd1,
    }
    df_ast = sp2df(sp_name, params)

    return df_ast

# ********************************************************************************************************************* 
def ast_add_earth_pos(df_ast: pd.DataFrame) -> None:
    """
    Add the splined earth position vectors to an asteroids DataFrame
    INPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
    OUTPUTS:
        None.  Modifies df_ast in place
    """
    # Spline earth and sun vectors
    ts = df_ast.mjd.values
    q_earth = get_earth_pos(ts)

    # Add new colums to DataFrame
    cols = ['earth_qx', 'earth_qy', 'earth_qz']
    df_ast[cols] = q_earth

# ********************************************************************************************************************* 
def ast_add_earth_vectors(df_ast: pd.DataFrame) -> None:
    """
    Add the splined earth position and velocity vectors to an asteroids DataFrame
    INPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
    OUTPUTS:
        None.  Modifies df_ast in place
    """
    # Spline earth and sun vectors
    ts = df_ast.mjd.values
    q_earth, v_earth = get_earth_vectors(ts)

    # Stack vectors
    qv_earth = np.hstack([q_earth, v_earth])

    # Add new colums to DataFrame
    cols_earth = ['earth_qx', 'earth_qy', 'earth_qz', 'earth_vx', 'earth_vy', 'earth_vz']
    df_ast[cols_earth] = qv_earth

# ********************************************************************************************************************* 
def ast_add_sun_vectors(df_ast: pd.DataFrame) -> None:
    """
    Add the splined sun position and velocity vectors to an asteroids DataFrame
    INPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
    OUTPUTS:
        None.  Modifies df_ast in place
    """
    # Spline earth and sun vectors
    ts = df_ast.mjd.values
    q_sun, v_sun= get_sun_vectors(ts)

    # Stack vectors
    qv_sun= np.hstack([q_sun, v_sun])

    # Add new colums to DataFrame
    cols_sun = ['sun_qx', 'sun_qy', 'sun_qz', 'sun_vx', 'sun_vy', 'sun_vz']
    df_ast[cols_sun] = qv_sun

# # ********************************************************************************************************************* 
# def spline_ast_vec_df(df_ast: pd.DataFrame, mjd: np.ndarray,
#                       include_elts: bool = True, progbar: bool = True) \
#                       -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
#     """
#     Load the MSE asteroid integrations for this range of asteroids.
#     INPUTS:
#         df_ast:    DataFrame with asteroid position and velocity
#         mjd:       Array of times at which splined output is desired
#         include_elts: Whether to include orbital elements for the asteroids and earth
#         progbar:   Whether to print a progress bar to the console
#     OUTPUTS:
#         df_ast_out:   Position & velocity of asteroids in barycentric frame
#         df_earth_out: Position & velocity of earth in barycentric frame
#         df_sun_out:   Position & velocity of sun in barycentric frame
#     """
#     # Time key from mjd at spline points
#     # time_key = np.int32(np.round(mjd*24*3600))

#     # Distinct asteroids or elements; get the name of the relevant column and array of distinct values
#     id_col = 'asteroid_num' if 'asteroid_num' in df_ast.columns else 'element_id'
#     id_val_unq = np.unique(df_ast[id_col].values)

#     # Number of asteroids or elements
#     N_ast = id_val_unq.size
#     # Number of times in input and splined output
#     N_t_in = df_earth.mjd.size
#     N_t_out = mjd.size
#     # Number of rows in input and splined output
#     N_row_in = N_ast * N_t_in
#     N_row_out = N_ast * N_t_out

#     # Data to be splined: x axis is time
#     x_spline = df_earth.mjd

#     # Desired output columns to be splined for asteroids, earth and sun
#     cols_cart = ['qx', 'qy', 'qz', 'vx', 'vy', 'vz']
#     cols_elt = ['a', 'e', 'inc', 'Omega', 'omega', 'f']
#     cols_ast = cols_cart + cols_elt if include_elts else cols_cart
#     cols_earth = cols_cart + cols_elt if include_elts else cols_cart
#     cols_sun = cols_cart

#     # the indexing of y_spline_ast is (ast_num, time_step, data_dim) with shape e.g. (16, 3653, 6)
#     y_spline_ast = df_ast[cols_ast].values.reshape(N_ast, N_t_in, -1)

#     # splining functions for asteroids, earth, and sun
#     spline_func_ast = CubicSpline(x=x_spline, y=y_spline_ast, axis=1)
#     spline_func_earth = CubicSpline(x=x_spline, y=y_spline_earth, axis=0)
#     spline_func_sun = CubicSpline(x=x_spline, y=y_spline_sun, axis=0)

#     # splined output for asteroids, earth and sun
#     spline_data_ast = spline_func_ast(mjd)
#     spline_data_earth = spline_func_earth(mjd)
#     spline_data_sun = spline_func_sun(mjd)

#     # ID column (asteroid_num or element_id) corresponding to splined output
#     id_val = np.repeat(id_val_unq, N_t_out)

#     # asteroid DataFrame
#     ast_dict_keys = {
#         id_col: id_val,
#         'mjd': np.tile(mjd, N_ast),
#         # 'time_key' : np.tile(time_key, N_ast)
#     }
#     ast_dict_data = {col:spline_data_ast[:,:,k].reshape(N_row_out) for k, col in enumerate(cols_ast)}
#     ast_dict = dict(**ast_dict_keys, **ast_dict_data)
#     df_ast_out = pd.DataFrame(ast_dict)

#     # earth DataFrame
#     earth_dict_keys = {
#         'mjd': mjd,
#         # 'time_key': time_key
#     }
#     earth_dict_data = {col:spline_data_earth[:,k] for k, col in enumerate(cols_earth)}
#     earth_dict = dict(**earth_dict_keys, **earth_dict_data)
#     df_earth_out = pd.DataFrame(earth_dict)

#     # sun DataFrame
#     sun_dict_keys = earth_dict_keys
#     sun_dict_data = {col:spline_data_sun[:,k] for k, col in enumerate(cols_sun)}
#     sun_dict = dict(**sun_dict_keys, **sun_dict_data)
#     df_sun_out = pd.DataFrame(sun_dict)
    
#     return df_ast_out, df_earth_out, df_sun_out
