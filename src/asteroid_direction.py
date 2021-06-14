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
from asteroid_spline import make_spline_ast_pos
from astro_utils import infer_shape
from topos import calc_topos
from db_utils import df2db

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
        q_tgt:  position of target body in ecliptic coordinate frame; either with astropy units or an array in AU
        v_tgt:  velocity of target body in ecliptic coordinate frame; either with astropy units or an array in AU/day
        q_obs:  position of observer in ecliptic coordinate frame; either with astropy units or an array in AU
                typically this will be computed as q_earth + dq_topos
    RETURNS:
        u:      An array [ux, uy, uz] on the unit sphere in the ecliptic frame
        delta:  The distance from observer to target in AU.
    """
    # Convert q_tgt and q_obs to be in units of AU
    if isinstance(q_tgt, astropy.units.quantity.Quantity):
        q_tgt = q_tgt.to(au).values
    if isinstance(q_obs, astropy.units.quantity.Quantity):
        q_obs = q_obs.to(au).values
    # Convert v_tgt to be in units of AU/day
    if isinstance(v_tgt, astropy.units.quantity.Quantity):
        v_tgt = v_tgt.to(au/day).values

    # At the end of this block, q_tgt, q_obs, and v_tgt are plain old numpy arrays
    # q_tgt and q_obs are in AU; v_tgt is in AU/day

    # First iteration of distance from object to observer
    r = calc_distance(q0=q_obs, q1=q_tgt)
    # Light time in days
    light_time = (r / c_au_day).reshape((-1, 1))

    # Adjusted relative position, accounting for light time
    dq = v_tgt * light_time
    # Position of object at time light left using linear model
    q_tgt_adj = q_tgt - dq

    # Calculate direction including down-leg light time aberration at constant velocity
    u, delta = calc_direction_dist(q0=q_obs, q1=q_tgt_adj)
    return u, delta

# ********************************************************************************************************************* 
def calc_dir_linear_topos(q_tgt: np.ndarray, v_tgt: np.ndarray, q_earth: np.ndarray, 
                          obstime_mjd: Optional[np.ndarray] = None, site_name: str = 'geocenter') -> np.ndarray:
    """
    Compute the direction of displacement from earth to a space body as a unit displacement vector
    u = (ux, uy, uz) in the ecliptic plane.
    Linear model including topos adjustment for position of observatory on Earth.
    INPUTS:
        q_body:         position of target body in BME; either with astropy units or an array in AU
        v_body:         velocity of target body in BME; passed with units (default AU / day)
        q_earth:        position of earth in BME; either with astropy units or an array in AU
        obstime_mjd:    observation time as a modified julian date; only required if passing obsgeoloc 
        site_name:      String describing geolocation of the observatory, e.g. 'geocenter' or 'palomar'
    RETURNS:
        u:              an array [ux, uy, uz] on the unit sphere in the ecliptic frame
        delta:          the distance from observer to target in AU.
    EXAMPLE:
        u, delta = calc_dir_linear_topos(q_body=np.array([-0.328365, 1.570624, 0.040733])*au, 
                   v_body=np.array([-0.013177, -0.001673, 0.000288])*au/day,
                   q_earth=np.array([-0.813785, -0.586761, -0.000003])*au,
                   obsgeoloc=[-2410346.78217658, -4758666.82504051, 3487942.97502457] * meter)
    """
    # Compute the correction due to the observatory of obstime_mjd and site_name are passed
    dq_topos, _ = calc_topos(obstime_mjd=obstime_mjd, site_name=site_name)

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
def calc_dir_ast2obs(t_obs: np.ndarray, asteroid_id: np.ndarray, q_obs: np.ndarray, iters: int=2):
    """
    Calculate direction for asteroids in the given range.
    The time when photons leave the asteroid is fixed on a regular schedule.
    The arrival time includes the light time
    INPUTS:
        t_obs:          Array of detection times
        asteroid_id:    Array of asteroid IDs
        q_obs:          Array of positions of the observatory in the BME frame
    RETURNS:
        df:     DataFrame with asteroid direction and light time.
                Columns: AsteroidID, TimeID, tAst, qAst_x, qAst_y, qAst_z, LightTime, tObs, qObs_x, qObs_y, qObs_z
    """
    # Get range of asteroids
    n0: int = np.min(asteroid_id)
    n1: int = np.max(asteroid_id)
    # Get range of MJD
    pad: int = 32
    mjd0: int = np.min(t_obs) - pad
    mjd1: int = np.max(t_obs) + pad

    # Build spline of asteroid posistion that supports this range of asteroids and dates
    spline_q_ast = make_spline_ast_pos(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Create seperate array for t_ast
    t_ast = t_obs.copy()

    # # Position of asteroid at t_obs
    # q_ast = spline_q_ast(ts=t_ast, asteroid_id=asteroid_id)
    # # Compute initial light time
    # r = calc_distance(q0=q_obs, q1=q_ast)
    # light_time = r/c
    light_time = np.zeros_like(t_ast)

    # Delegate to calc_dir_spline
    u, delta = calc_dir_spline(spline_q_ast=spline_q_ast, q_obs=q_obs, t_obs=t_obs, 
                               asteroid_id=asteroid_id, light_time=light_time, iters=iters)

    # Save results to DataFrame
    df = pd.DataFrame()
    # Key columns
    df['AsteroidID'] = asteroid_id
    df['tObs'] = t_obs
    # Direction and light time
    df[cols_dir] = u
    df['LightTime'] = light_time

    # TESTING
    t_ast = t_obs - light_time / mpd
    df['tAst'] = t_ast
    df[cols_q_obs] = q_obs
    df[cols_q_ast] = spline_q_ast(t_ast)

    return df

# ********************************************************************************************************************* 
def calc_dir_ast2obs_topos(ts: np.ndarray, asteroid_id: np.ndarray, site_name: str):
    """Calculate direction of asteroids to a named observatory"""
    # The splined Earth position at the observer times
    q_earth = get_earth_pos(ts=ts)
    # Compute the correction due to topos
    dq_topos, _ = calc_topos(obstime_mjd=ts, site_name=site_name)
    data_axis, space_axis, shape = infer_shape(q_earth)
    dq_topos = dq_topos.reshape(shape)
    # Position of the observer in BME frame
    q_obs = q_earth + dq_topos
    # Delegate to calc_dir_ast2obs()
    pass

# ********************************************************************************************************************* 
def light_time_error(df: pd.DataFrame):
    """Calculate the light time error in minutes and print it to screen"""
    err = df.LightTime - (df.tObs - df.tAst)*mpd
    mean_err = np.mean(np.abs(err))
    max_err = np.max(np.abs(err))
    print(f'Light time error (in minutes):')
    print(f'mean error: {mean_err:5.3e}')
    print(f'max error : {max_err:5.3e}')

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
    table = 'AsteroidDirections'
    columns = ['AsteroidID', 'TimeID', 'tObs', 'ux', 'uy', 'uz', 'LightTime']
    chunksize = 2**19
    verbose = False
    progbar = False

    # Dispatch to df2db
    df2db(df=df, schema=schema, table=table, columns=columns, chunksize=chunksize, verbose=verbose, progbar=progbar)

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

    # Report arguments
    print(f'Processing asteroid directions for asteroid number {n0} <= AsteroidID < {n1}...')
    # Set the batch size
    b: int = 200
    k0: int = n0 // b
    k1: int = n1 // b
    for k in tqdm(range(k0, k1)):
        # Start and end of this batch
        n0_i = k*b
        n1_i = n0_i + b
        # Calculate the direction and light time
        df = calc_dir_ast2obs(n0=n0_i, n1=n1_i)
        # Insert results to database
        insert_dir_ast2obs(df=df)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
