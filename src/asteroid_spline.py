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
from planets_interp import get_earth_pos, get_earth_vectors, get_sun_vectors
# from ra_dec import qv2dir, dir2radec, calc_topos, astrometric_dir, direction_diff

# ********************************************************************************************************************* 
# Speed of light; express this in AU / day
light_speed_au_day = astropy.constants.c.to(au / day).value

# ********************************************************************************************************************* 
def spline_ast_vec(df_ast: pd.DataFrame, mjd: np.ndarray) -> pd.DataFrame:
    """
    Spline the integrated state vectors of asteroids and earth at the desired times
    INPUTS:
        df_ast:    DataFrame with asteroid position and velocity
        mjd:       Array of times at which splined output is desired
    OUTPUTS:
        df_out:    Position & velocity of asteroids in barycentric frame
    """
    # distinct asteroids or elements; get the name of the relevant column and array of distinct values
    id_col = 'AsteroidID' if 'AsteroidID' in df_ast.columns else 'ElementID'
    id_val_in = df_ast[id_col].values
    id_val_unq = np.unique(id_val_in)

    # Number of asteroids or elements
    N_ast = id_val_unq.size
    # Number of times in splined output
    N_t_out = mjd.shape[0]

    # Number of rows in input and splined output
    N_row_in = df_ast.shape[0]
    N_row_out = N_ast * N_t_out

    # Calculated number of times in input
    N_t_in = N_row_in // N_ast

    # Data to be splined: x axis is time
    x_spline = df_ast.mjd.values[0:N_t_in]

    # Collections of columns
    cols_cart = ['qx', 'qy', 'qz', 'vx', 'vy', 'vz']
    # cols_elt = ['a', 'e', 'inc', 'Omega', 'omega', 'f']
    cols_earth = ['earth_qx', 'earth_qy', 'earth_qz']
    # cols_sun = ['sun_qx', 'sun_qy', 'sun_qz']
    # The columns to be splined
    cols_spline = cols_cart

    # the indexing of y_spline is (ast_num, time_step, data_dim) with shape e.g. (16, 3653, 16)
    y_spline = df_ast[cols_spline].values.reshape(N_ast, N_t_in, -1)

    # splining function for the splined columns
    spline_func = CubicSpline(x=x_spline, y=y_spline, axis=1)

    # splined output
    spline_data = spline_func(mjd)

    # ID column (asteroid_num or element_id) corresponding to splined output
    id_val = np.repeat(id_val_unq, N_t_out)

    # Key fields in the output (splined) DataFrame
    mjd_out = np.tile(mjd, N_ast)
    out_dict_keys = {
        id_col: id_val,
        'mjd': mjd_out
    }

    # The splined asteroid state vectos
    out_dict_splined = {col:spline_data[:,:,k].reshape(N_row_out) for k, col in enumerate(cols_spline)}
    
    # Wrap the output arrays into one dictionary and DataFrame
    out_dict = dict(**out_dict_keys, **out_dict_splined)
    df_out = pd.DataFrame(out_dict)

    # Add the splined earth position to df_out
    q_earth = get_earth_pos(mjd_out)
    df_out[cols_earth] = q_earth
    
    return df_out
