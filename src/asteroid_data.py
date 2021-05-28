"""
Load integrated asteroid trajectories as Pandas DataFrames.
Add earth or sun position / state vectors to asteroid DataFrame.

Michael S. Emanuel
Sat Sep 21 10:38:38 2019
"""

# Core
import numpy as np
import pandas as pd

# Utility
# from tqdm.auto import tqdm

# Local imports
from db_utils import sp2df
from planets_interp import get_earth_pos, get_earth_vectors, get_earth_elt, get_sun_vectors

# Type names
from typing import Tuple

# ********************************************************************************************************************* 
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
def ast_add_earth_elt(df_ast: pd.DataFrame) -> None:
    """
    Add the splined earth orbital elements to an asteroids DataFrame
    INPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
    OUTPUTS:
        None.  Modifies df_ast in place
    """
    # Spline earth and sun vectors
    ts = df_ast.mjd.values
    elt_earth = get_earth_elt(ts)

    # Add new colums to DataFrame
    cols = ['earth_a', 'earth_e', 'earth_inc', 'earth_Omega', 'earth_omega', 'earth_f', 'earth_M']
    df_ast[cols] = elt_earth

# # ********************************************************************************************************************* 
# def ast_add_earth_vectors(df_ast: pd.DataFrame) -> None:
#     """
#     Add the splined earth position and velocity vectors to an asteroids DataFrame
#     INPUTS:
#         df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
#     OUTPUTS:
#         None.  Modifies df_ast in place
#     """
#     # Spline earth and sun vectors
#     ts = df_ast.mjd.values
#     q_earth, v_earth = get_earth_vectors(ts)

#     # Stack vectors
#     qv_earth = np.hstack([q_earth, v_earth])

#     # Add new colums to DataFrame
#     cols_earth = ['earth_qx', 'earth_qy', 'earth_qz', 'earth_vx', 'earth_vy', 'earth_vz']
#     df_ast[cols_earth] = qv_earth

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
