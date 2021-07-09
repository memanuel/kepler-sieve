"""
Load integrated asteroid trajectories as Pandas DataFrames.
Add earth or sun position / state vectors to asteroid DataFrame.

Functions in this module:
load_ast_vectors(n0, n1, mjd0, mjd1)
load_ast_elements(n0, n1, mjd0, mjd1)
load_ast_data(n0, n1, mjd0, mjd1)
ast_add_earth_pos(df_ast)
ast_add_earth_elts(df_ast)
ast_add_sun_vectors(df_ast)

Michael S. Emanuel
Sat Sep 21 10:38:38 2019
"""

# Core
import numpy as np
import pandas as pd

# Local imports
from db_utils import sp2df
from planets_interp import get_earth_pos, get_earth_elts, get_sun_vectors

# Type names
from typing import Tuple

# ********************************************************************************************************************* 
# Last date of calculated data; filtering for mjd <= this date will include everything
# The *real* last date is mjd=63000.  This is meant to be a future proof dummy "big" date (DB won't accept infinity).
mjd_last: np.float64 = np.float64(2.0**20)

# ********************************************************************************************************************* 
# Catalog of known asteroids
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def get_asteroid_ids() -> np.ndarray:
    """
    Return array asteroid_id with the AsteroidID field of all known asteroids
    """
    ast: pd.DataFrame = sp2df(sp_name='KS.GetAsteroidIDs')
    return ast.AsteroidID.values

# ********************************************************************************************************************* 
def get_asteroids(key_to_body_id: bool=False) -> pd.DataFrame:
    """
    Return list of known asteroid names and IDs
    INPUTS:
        key_to_body_id: When true, reindex this DataFrame to BodyID.  Default (false) keys to AsteroidID.
    """
    ast: pd.DataFrame = sp2df(sp_name='KS.GetAsteroids')
    if key_to_body_id:
        ast.set_index(keys='BodyID', drop=False, inplace=True)
    else:
        ast.set_index(keys='AsteroidID', drop=False, inplace=True)
    return ast

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
def load_ast_pos(n0: int, n1: int, mjd0: float=0.0, mjd1: float=mjd_last) -> pd.DataFrame:
    """
    Load the MSE asteroid integrations for this range of asteroids; returns just position.
    INPUTS:
        n0:  First asteroid to load, e.g. 0
        n1:  Last asteroid to load, (exclusive) e.g. 1000
        mjd0:  Start modfified julian date used to filter output.
        mjd1:  Last modified julian date used to filter output.
               Default for mjd0 and mjd1 is None; then return all available time steps
    OUTPUTS:
        df_ast:   Position of asteroids in barycentric frame
    """

    # Get asteroid data (state vectors and orbital elements) from database
    sp_name = 'KS.GetAsteroidPositions'
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
def load_ast_data(n0: int, n1: int, mjd0: float=0.0, mjd1: float=mjd_last) -> pd.DataFrame:
    """
    Load the MSE asteroid integrations for this range of asteroids; returns both state vectors and orbital elements.
    INPUTS:
        n0:     First asteroid to load, e.g. 0
        n1:     Last asteroid to load, (exclusive) e.g. 1000
        mjd0:   Start modfified julian date used to filter output.
        mjd1:   Last modified julian date used to filter output.
                Default for mjd0 and mjd1 is None; then return all available time steps
    OUTPUTS:
        df_ast: Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
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
def load_ast_dir(n0: int, n1: int, mjd0: float=0.0, mjd1: float=mjd_last) -> pd.DataFrame:
    """
    Load the directions from asteroids to earth geocenter from MSE integration.
    The times are tObs (observer time), NOT the time light left the asteroid!
    INPUTS:
        n0:     First asteroid to load, e.g. 0
        n1:     Last asteroid to load, (exclusive) e.g. 1000
        mjd0:   Start modfified julian date used to filter output.
        mjd1:   Last modified julian date used to filter output.
                Default for mjd0 and mjd1 is None; then return all available time steps
    OUTPUTS:
        df_ast: Direction from earth geocenter and light time
    """

    # Get asteroid directions in this time range
    sp_name = 'KS.GetAsteroidDirection'
    params = {
        'n0': n0,
        'n1': n1,
        'mjd0': mjd0,
        'mjd1': mjd1
    }
    df_ast = sp2df(sp_name=sp_name, params=params)

    return df_ast

## ********************************************************************************************************************* 
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
def ast_add_earth_elts(df_ast: pd.DataFrame) -> None:
    """
    Add the splined earth orbital elements to an asteroids DataFrame
    INPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
    OUTPUTS:
        None.  Modifies df_ast in place
    """
    # Spline earth and sun vectors
    ts = df_ast.mjd.values
    elt_earth = get_earth_elts(ts)

    # Add new colums to DataFrame
    cols = ['earth_a', 'earth_e', 'earth_inc', 'earth_Omega', 'earth_omega', 'earth_f', 'earth_M']
    df_ast[cols] = elt_earth

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
    qv_sun = np.hstack([q_sun, v_sun])

    # Add new colums to DataFrame
    cols_sun = ['sun_qx', 'sun_qy', 'sun_qz', 'sun_vx', 'sun_vy', 'sun_vz']
    df_ast[cols_sun] = qv_sun