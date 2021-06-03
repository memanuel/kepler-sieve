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
from astropy.units import au, minute

# Commandline arguments
import argparse

# Utility
from tqdm.auto import tqdm

# Typing
from typing import Callable

# Local imports
from asteroid_data import load_ast_pos
from planets_interp import get_earth_pos
from asteroid_spline import make_spline_df
from db_utils import df2db

# ********************************************************************************************************************* 
# Speed of light; express this in AU / min
c = astropy.constants.c.to(au / minute).value

# Column names used in calculations
cols_q_ast = ['qAst_x', 'qAst_y', 'qAst_z']
cols_q_obs = ['qObs_x', 'qObs_y', 'qObs_z']

# ********************************************************************************************************************* 
def light_time_iter(df: pd.DataFrame, spline_q_ast: Callable[[np.ndarray, np.ndarray], np.ndarray]):
    """
    One iteration of refinement of light time calculation; modifies DataFrame in place.
    INPUTS:
        df:             DataFrame that includes columns:
                        AsteroidID, tAst, qAst_x, qAst_y, qAst_z, tObs, qObs_x, qObs_y, qObs_z, LightTime
        spline_q_ast:   Spline function that returns the position of asteroids vs. time.
    OUTPUTS:
        None.           Modifies df in place.
    """
    # Unpack vectors from DataFrame
    asteroid_id = df.AsteroidID.values
    t_obs = df.tObs.values
    q_obs = df[cols_q_obs].values
    light_time = df.LightTime.values

    # Calculate t_ast from t_obs and current estimate of light_time
    t_ast = t_obs - light_time/1440.0             # 1440 is number of minutes in one day
    # Save revised time light leaves asteroid to DataFrame
    df['tAst'] = t_ast

    # Calculate new asteroid positions at the revised times and save to DataFrame
    q_ast = spline_q_ast(x=t_ast, y=asteroid_id)
    df[cols_q_ast] = q_ast

    # Compute position difference and distance from asteroid to earth
    dq = q_ast-q_obs
    r = np.sqrt(np.sum(dq*dq, axis=-1))

    # Compute light time and update it on DataFrame
    light_time = r/c
    df['LightTime'] = light_time

# ********************************************************************************************************************* 
def calc_dir_ast2obs(n0: int, n1: int):
    """
    Calculate light time and direction for asteroids in the given range.
    The time when photons leave the asteroid is fixed on a regular schedule.
    The arrival time includes the light time
    INPUTS:
        n0: Fist asteroid number to include; inclusive.
        n1: Last asteroid number to include; exclusive.
    RETURNS:
        df: DataFrame with asteroid direction and light time.
            Columns: AsteroidID, TimeID, tAst, qAst_x, qAst_y, qAst_z, LightTime, tObs, qObs_x, qObs_y, qObs_z
    """
    # Load asteroid positions
    ast_pos = load_ast_pos(n0=n0, n1=n1)

    # Reorder columns to put AsteroidID first
    cols = ['AsteroidID', 'TimeID'] + list(ast_pos.columns[2:])
    ast_pos = ast_pos[cols]

    # Rename columns
    col_tbl = {
        'mjd': 'tAst',
        'qx': 'qAst_x',
        'qy': 'qAst_y',
        'qz': 'qAst_z',
    }
    ast_pos.rename(columns=col_tbl, inplace=True)

    # Copy asteroid positions into DataFrame used for directions
    df = ast_pos.copy()

    # The time the light leaves the asteroid and arrives at the observer
    t_ast = df.tAst.values
    t_obs = t_ast
    # Save t_obs to DataFrame
    df['tObs'] = t_obs

    # Spline earth vectors and update qObs on DataFrame; this will not change
    q_obs = get_earth_pos(ts=t_obs)
    df[cols_q_obs] = q_obs

    # Extract asteroid position vectors
    q_ast = ast_pos[cols_q_ast].values

    # Compute position difference and distance from asteroid to earth
    dq = q_ast-q_obs
    r = np.sqrt(np.sum(dq*dq, axis=-1))

    # Compute light time and update it on DataFrame
    light_time = r/c
    df['LightTime'] = light_time

    # Build asteroid position spline
    id_col = 'AsteroidID'
    time_col = 'tAst'
    spline_q_ast = make_spline_df(df=ast_pos, cols_spline=cols_q_ast, id_col=id_col, time_col=time_col)

    # Experiments show that 3 iterations is sufficient to achieve full convergence
    # Error in light_time (number of minutes) around 5E-9 at this point
    for _ in range(4):
        light_time_iter(df=df, spline_q_ast=spline_q_ast)

    # Compute position difference and distance from asteroid to earth
    q_obs = df[cols_q_obs].values
    dq = q_ast - q_obs
    r = np.sqrt(np.sum(np.square(dq), axis=-1, keepdims=True))

    # Calculate the direction and save it to DataFrame
    u = dq/r
    cols_dir = ['ux', 'uy', 'uz']
    df[cols_dir] = u

    return df

# ********************************************************************************************************************* 
def light_time_error(df: pd.DataFrame):
    """Calculate the light time error in minutes and print it to screen"""
    err = df.LightTime - (df.tObs - df.tAst)*1440.0
    mean_err = np.mean(np.abs(err))
    max_err = np.max(np.abs(err))
    print(f'Light time error (in minutes):')
    print(f'mean error: {mean_err:5.3e}')
    print(f'max error : {max_err:5.3e}')

# ********************************************************************************************************************* 
def insert_dir_ast2obs(df: pd.DataFrame):
    """Insert asteroid direction calculations """
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
    'implied by rebound integration.  Populates DB table KS.AsteroidDirections.')
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
