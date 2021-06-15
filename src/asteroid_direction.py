"""
Direction from known asteroids to Earth center implied by integrated trajectories.
Example calls:
$ python asteroid_direction.py 0 1000

Functions in this module:

light_time_iter(df, t_ast, q_ast)
calc_dir_ast2obs(n0, n1)
light_time_error(df)
insert_dir_ast2obs(df)
main()

Michael S. Emanuel
2021-06-01
"""

# Core
import numpy as np
import pandas as pd

# Astronomy
import astropy
from astropy.units import au, day, minute

# Commandline arguments
import argparse

# Utility
from tqdm.auto import tqdm

# Typing
from typing import Callable, Optional
spline_type_ast = Callable[[np.ndarray, np.ndarray], np.ndarray]
spline_type_obs = Callable[[np.ndarray], np.ndarray]

# Local imports
from planets_interp import get_earth_pos
from asteroid_spline import make_spline_ast_vec, make_spline_ast_pos, make_spline_ast_pos_direct
from astro_utils import infer_shape
from topos import calc_topos
from db_utils import sp2df, df2db
from utils import arange_inc

# ********************************************************************************************************************* 
# Speed of light; express this in AU / min
c = astropy.constants.c.to(au / minute).value
# Speed of light in AU/day
c_au_day = astropy.constants.c.to(au / day).value

# Number of minutes in one day
mpd: float = 1440.0
# Number of days in one minute
dpm: float = 1.0 / mpd

# Column names used in calculations
cols_q_obs = ['qObs_x', 'qObs_y', 'qObs_z']
cols_q_ast = ['qAst_x', 'qAst_y', 'qAst_z']
cols_dir = ['ux', 'uy', 'uz']

# DataFrame of asteroids
ast = sp2df(sp_name='KS.GetAsteroids')

# *************************************************************************************************
# Utility functions - calculate distance and direction between two positions
# *************************************************************************************************

# ********************************************************************************************************************* 
def calc_distance(q0: np.ndarray, q1:np.ndarray):
    """
    Distance between two position arrays.
    INPUTS:
        q0: First array of positions; shape (N, 3,)
        q1: Second array of positions; shape (N, 3,)
    OUTPUTS:
        r:  Array of distances from q0 to q1; shape (N,)
    """
    dq = q1 - q0
    r = np.sqrt(np.sum(np.square(dq), axis=-1))
    return r

# ********************************************************************************************************************* 
def calc_direction(q0: np.ndarray, q1:np.ndarray):
    """
    Direction between two position arrays.
    INPUTS:
        q0: First array of positions; shape (N, 3,)
        q1: Second array of positions; shape (N, 3,)
    OUTPUTS:
        u:  Array of directions from q0 to q1; shape (N, 3,)
    """
    dq = q1 - q0
    r = np.sqrt(np.sum(np.square(dq), axis=-1, keepdims=True))
    u = dq / r
    return u

# ********************************************************************************************************************* 
def calc_direction_dist(q0: np.ndarray, q1:np.ndarray):
    """
    Direction and distance between two position arrays.
    INPUTS:
        q0: First array of positions; shape (N, 3,)
        q1: Second array of positions; shape (N, 3,)
    OUTPUTS:
        u:  Array of directions from q0 to q1; shape (N, 3,)
        r:  Array of distances from q0 to q1; shape (N,)    
    """
    dq = q1 - q0
    r = np.sqrt(np.sum(np.square(dq), axis=-1, keepdims=True))
    u = dq / r
    return u, r.squeeze()

# *************************************************************************************************
# Calculate direction in the linear model
# *************************************************************************************************

# ********************************************************************************************************************* 
def calc_dir_linear(q_tgt: np.ndarray, v_tgt: np.ndarray, q_obs: np.ndarray):
    """
    Compute the astrometric direction from earth to a space body as a unit displacement vector
    u = (ux, uy, uz) in the ecliptic plane.
    Uses a simple linear model where the target is assumed to be moving at constant velocity << c (speed of light).
    INPUTS:
        q_tgt:  position of target body in ecliptic coordinate frame; array in AU
        v_tgt:  velocity of target body in ecliptic coordinate frame; array in AU/day
        q_obs:  position of observer in ecliptic coordinate frame; array in AU
                typically this will be computed as q_earth + dq_topos
    RETURNS:
        u:      An array [ux, uy, uz] on the unit sphere in the ecliptic frame
        delta:  The distance from observer to target in AU.
    """
    # First iteration of distance from object to observer
    r = calc_distance(q0=q_obs, q1=q_tgt)
    # Light time in DAYS (usually do this in minutes, but here want days so can multiply by v_tgt)
    light_time = (r / c_au_day).reshape((-1, 1))

    # Adjusted relative position, accounting for light time
    dq = v_tgt * light_time
    # Position of object at time light left using linear model
    q_tgt_adj = q_tgt - dq

    # Calculate direction including down-leg light time aberration at constant velocity
    u, delta = calc_direction_dist(q0=q_obs, q1=q_tgt_adj)
    return u, delta

# ********************************************************************************************************************* 
def calc_dir_linear_topos(q_tgt: np.ndarray, v_tgt: np.ndarray, t_obs: np.ndarray, 
                          site_name: str = 'geocenter') -> np.ndarray:
    """
    Compute the direction of displacement from earth to a space body as a unit displacement vector
    u = (ux, uy, uz) in the ecliptic plane.
    Linear model including topos adjustment for position of observatory on Earth.
    INPUTS:
        q_body:         position of target body in BME; array in AU
        v_body:         velocity of target body in BME; array in AU / day
        t_obs:          observation time as a modified julian dates
        site_name:      String describing geolocation of the observatory, e.g. 'geocenter' or 'palomar'
    RETURNS:
        u:              an array [ux, uy, uz] on the unit sphere in the ecliptic frame
        delta:          the distance from observer to target in AU.
    """
    # Calculate Earth position at observation times
    q_earth = get_earth_pos(ts=t_obs) 
    # Compute the correction due to the observatory of t_obs and site_name
    dq_topos, _ = calc_topos(t_obs=t_obs, site_name=site_name)
    # Reshape dq_topos to match q_earth
    data_axis, space_axis, shape = infer_shape(q_earth)
    dq_topos = dq_topos.reshape(shape)
    # Position of the observer in space
    q_obs = q_earth + dq_topos

    # Calculate astrometric direction in the linear model
    u, delta = calc_dir_linear(q_tgt=q_tgt, v_tgt=v_tgt, q_obs=q_obs)

    return u, delta
    
# *************************************************************************************************
# Calculate direction with iterative method and an asteroid position spline
# *************************************************************************************************

# *************************************************************************************************
def light_time_iter(spline_q_ast: spline_type_ast, q_obs: np.ndarray, t_obs: np.ndarray, 
                    asteroid_id: np.ndarray, light_time: np.ndarray):
    """
    One iteration of refinement of light time calculation specialized for asteroids.
    INPUTS:
        spline_q_ast:   Spline function that returns the position of asteroids vs. time.
                        Inputs to this spline are x=t_ast, y=asteroid_id.
        q_obs:          Position of observer at observation times.
        t_obs:          Array of observation times
        asteroid_id:    Array of asteroid IDs
        light_time:     Current estimate of light time
    OUTPUTS:
        light_time:     New estimate of light time
    """
    # Calculate t_ast from t_obs and current estimate of light_time
    t_ast = t_obs - (light_time*dpm)
    # Calculate new asteroid positions at the revised times
    q_ast = spline_q_ast(t_ast, asteroid_id)
    # Compute distance from observer to asteroid
    r = calc_distance(q0=q_obs, q1=q_ast)
    # Compute light time
    return r/c

# ********************************************************************************************************************* 
def calc_dir_spline(spline_q_ast: spline_type_ast, q_obs: np.ndarray, t_obs: np.ndarray, 
                    asteroid_id: np.ndarray, light_time: np.ndarray, iters:int=2):
    """
    One iteration of refinement of light time calculation specialized for asteroids.
    INPUTS:
        spline_q_ast:   Spline function that returns the position of asteroids vs. time.
                        Inputs to this spline are x=t_ast, y=asteroid_id.
        q_obs:          Position of observer at observation times.
        t_obs:          Array of observation times
        asteroid_id:    Array of asteroid IDs
        light_time:     Initial estimate of light time
    OUTPUTS:
        u
        light_time:     Final estimate of light time
    """
    # Iterate to find the light time
    for i in range(iters):
        light_time = light_time_iter(spline_q_ast=spline_q_ast, q_obs=q_obs, t_obs=t_obs, 
                                     asteroid_id=asteroid_id, light_time=light_time)

    # Calculate t_ast from t_obs and final estimate of light time
    t_ast = t_obs - (light_time*dpm)
    # Calculate asteroid positions at the asteroid times (when light leaves)
    q_ast = spline_q_ast(t_ast, asteroid_id)
    # Calculate direction and distance
    u, delta = calc_direction_dist(q0=q_obs, q1=q_ast)
    return u, delta

# *************************************************************************************************
# End to end calculation of direction from an observatory to a known asteroid.
# *************************************************************************************************
def calc_dir_ast2obs(t_obs: np.ndarray, asteroid_id: np.ndarray, q_obs: np.ndarray):
    """
    Calculate direction for asteroids in the given range.
    Use linear method, with position and velocity of asteroid calculated via splined orbital elements.
    The time when photons leave the asteroid is fixed on a regular schedule.
    The arrival time includes the light time.
    INPUTS:
        t_obs:          Array of detection times
        asteroid_id:    Array of asteroid IDs
        q_obs:          Array of positions of the observatory in the BME frame
    RETURNS:
        df:             DataFrame with asteroid direction and light time.
                        Columns: AsteroidID, tObs, ux, uy, uz, LightTime
    """
    # Get range of asteroids
    n0: int = np.min(asteroid_id)
    n1: int = np.max(asteroid_id)+1     # Need to add 1 here b/c n1 is EXCLUSIVE in make_spline_ast_pos()
    # Get range of MJD
    mjd0: int = np.min(t_obs)
    mjd1: int = np.max(t_obs)

    # Build spline of asteroid posistion that supports this range of asteroids and dates
    spline_vec_ast = make_spline_ast_vec(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)
    # Position of asteroid at t_obs    
    q_ast, v_ast = spline_vec_ast(ts=t_obs, asteroid_id=asteroid_id)
    # Delegate to calc_dir_linear
    u, delta = calc_dir_linear(q_tgt=q_ast, v_tgt=v_ast, q_obs=q_obs)
    # Light time from delta
    light_time = delta / c

    # Save results to DataFrame
    df = pd.DataFrame()
    # Key columns
    df['AsteroidID'] = asteroid_id
    df['tObs'] = t_obs
    # Direction and light time
    df[cols_dir] = u
    df['LightTime'] = light_time

    return df

# ********************************************************************************************************************* 
def calc_dir_ast2obs_spline(t_obs: np.ndarray, asteroid_id: np.ndarray, q_obs: np.ndarray, iters: int=2):
    """
    Calculate direction for asteroids in the given range.
    Use iterated spline method. (This is slower but results appear to be basically the same.)
    The time when photons leave the asteroid is fixed on a regular schedule.
    The arrival time includes the light time
    INPUTS:
        t_obs:          Array of detection times
        asteroid_id:    Array of asteroid IDs
        q_obs:          Array of positions of the observatory in the BME frame
    RETURNS:
        df:             DataFrame with asteroid direction and light time.
                        Columns: AsteroidID, tObs, ux, uy, uz, LightTime
    """
    # Get range of asteroids
    n0: int = np.min(asteroid_id)
    n1: int = np.max(asteroid_id)+1     # Need to add 1 here b/c n1 is EXCLUSIVE in make_spline_ast_pos()
    # Get range of MJD
    mjd0: int = np.min(t_obs)
    mjd1: int = np.max(t_obs)

    # Build spline of asteroid posistion that supports this range of asteroids and dates
    spline_q_ast = make_spline_ast_pos(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Create separate array for t_ast
    t_ast = t_obs.copy()

    # Position of asteroid at t_obs    
    q_ast = spline_q_ast(ts=t_ast, asteroid_id=asteroid_id)
    # Compute initial light time
    r = calc_distance(q0=q_obs, q1=q_ast)
    light_time = r/c

    # Delegate to calc_dir_spline
    u, delta = calc_dir_spline(spline_q_ast=spline_q_ast, q_obs=q_obs, t_obs=t_obs, 
                               asteroid_id=asteroid_id, light_time=light_time, iters=iters)
    # Update light_time from delta
    light_time = delta / c

    # Save results to DataFrame
    df = pd.DataFrame()
    # Key columns
    df['AsteroidID'] = asteroid_id
    df['tObs'] = t_obs
    # Direction and light time
    df[cols_dir] = u
    df['LightTime'] = light_time

    return df

# *************************************************************************************************
# Batch process directions to a block of asteroids at regular intervals
# *************************************************************************************************

# ********************************************************************************************************************* 
def prep_ast_block(n0: int, n1: int, mjd0: int, mjd1: int, interval_min: int):
    """
    Prepare arrays used in calc_dir_ast_block.
    N_ast = number of asteroids
    N_t   = number of times
    N_row = N_ast x N_t (number of rows)
    INPUTS:
        n0:             First asteroid to process; inclusive.
        n1:             Last asteroid to process; exclusive.
        mjd0:           First date to process
        mjd1:           Last date to process
        interval_min:   Step between times in MINUTES
    OUTPUTS:
        t_obs:          Array of observation times; shape (N_row,)
        asteroid_id:    Array of Asteroid IDs; shape (N_row,)
    """
    # Start, end and step for TimeID
    TimeID_0 = mjd0 * mpd
    TimeID_1 = mjd1 * mpd
    # Array of TimeID
    TimeIDs_one = arange_inc(TimeID_0, TimeID_1, interval_min).astype(np.int64)

    # Get array of asteroid IDs between n0 and n1
    mask = (n0 <= ast.AsteroidID) & (ast.AsteroidID < n1)
    asteroid_id_one = ast.AsteroidID[mask].values

    # The number of observation times and asteroids
    N_t = TimeIDs_one.shape[0]
    N_ast = asteroid_id_one.shape[0]

    # Tile TimeIDs and calculate t_obs from TimeIDs
    TimeIDs = np.tile(TimeIDs_one, N_ast)
    t_obs = TimeIDs * dpm
    # Repeat the asteroid IDs
    asteroid_id = np.repeat(asteroid_id_one, N_t)

    return t_obs, asteroid_id
    
# ********************************************************************************************************************* 
def calc_dir_ast_block(n0: int, n1: int, mjd0: int, mjd1: int, interval_min: int):
    """
    Calculate direction from Earth geocenter to a block of asteroids
    INPUTS:
        n0:             First asteroid to process; inclusive.
        n1:             Last asteroid to process; exclusive.
        mjd0:           First date to process
        mjd1:           Last date to process
        interval_min:   Step between times in MINUTES
    OUTPUTS:
        df:             DataFrame of asteroid directions
    """
    # Build arrays of observation times and asteroid IDs
    t_obs, asteroid_id = prep_ast_block(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1, interval_min=interval_min)

    # Calculate observer position; using geocentric observer so there is no topos adjustment
    q_obs = get_earth_pos(ts=t_obs)

    # Delegate to calc_dir_ast2obs()
    df = calc_dir_ast2obs(t_obs=t_obs, asteroid_id=asteroid_id, q_obs=q_obs)
    # Add the supplemental column TimeID
    time_id = np.round(t_obs*mpd).astype(np.int64)
    df['TimeID'] = time_id
    return df

# *************************************************************************************************
# Populate DB table KS.AsteroidDirection with directions from Earth geocenter to known asteroids.
# *************************************************************************************************

# ********************************************************************************************************************* 
def insert_dir_ast2obs(df: pd.DataFrame):
    """Insert asteroid direction calculations to database"""
    # Rename column tAst back to mjd to match DB schema
    df.rename(columns={'tAst':'mjd'}, inplace=True)

    # Arguments to df2db
    schema = 'KS'
    table = 'AsteroidDirection'
    columns = ['AsteroidID', 'TimeID', 'tObs', 'ux', 'uy', 'uz', 'LightTime']
    chunksize = 2**19
    verbose = False
    progbar = False

    # Dispatch to df2db
    df2db(df=df, schema=schema, table=table, columns=columns, chunksize=chunksize, verbose=verbose, progbar=progbar)

# *************************************************************************************************
# Console program populates DB table KS.AsteroidDirection
# *************************************************************************************************

# ********************************************************************************************************************* 
def main():
    """Calculate the direction and light time of the selected batch of asteroid"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description='Calculated direction from known asteroids to Earth center '
    'implied by rebound integration.  Populates DB table KS.AsteroidDirection.')
    parser.add_argument('n0', nargs='?', metavar='n0', type=int, default=0,
                        help='the first asteroid number to process')
    parser.add_argument('n_ast', nargs='?', metavar='B', type=int, default=1000,
                        help='the number of asteroids to process in this batch'),
    
    # Unpack command line arguments
    args = parser.parse_args()
    
    # Block of asteroids to integrate
    n0: int = args.n0
    n1: int = n0 + args.n_ast

    # Set the time range by policy to match AsteroidVectors and AsteroidElements
    mjd0: int = 48000
    mjd1: int = 63000
    interval_min: int = int(4*mpd)

    # Report arguments
    print(f'Processing asteroid directions for asteroid number in range [{n0}, {n1})...')
    # Set the batch size
    b: int = 100
    k0: int = n0 // b
    k1: int = max(n1 // b, k0+1)
    for k in tqdm(range(k0, k1)):
        # Start and end of this batch
        n0_i = k*b
        n1_i = min(n0_i + b, n1)
        # Calculate the direction and light time
        df = calc_dir_ast_block(n0=n0_i, n1=n1_i, mjd0=mjd0, mjd1=mjd1, interval_min=interval_min)
        # Insert results to database
        insert_dir_ast2obs(df=df)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
