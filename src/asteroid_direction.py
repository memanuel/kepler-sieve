"""
Direction from known asteroids to Earth center implied by integrated trajectories.
Example calls:
$ python asteroid_direction.py 0 1000

Functions in this module:
calc_light_time(n0, n1)
light_time_iter(df, t_ast, q_ast)
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

# Local imports
from asteroid_data import load_ast_vectors
from planets_interp import get_earth_pos
from db_utils import df2db

# ********************************************************************************************************************* 
# Speed of light; express this in AU / min
c = astropy.constants.c.to(au / minute).value

# ********************************************************************************************************************* 
def light_time_iter(df: pd.DataFrame, t_ast: np.ndarray, q_ast: np.ndarray):
    """
    One iteration of refinement of light time calculation; modifies DataFrame in place.
    INPUTS:
        df:     DataFrame assembled in function calc_light_time.
        t_ast:  Array of times when photons leave asteroid.
        q_ast:  Position of asteroids at these times.
    OUTPUTS:
        None.  Modifies df in place.
    """
    # Column names
    cols_q_obs = ['qObs_x', 'qObs_y', 'qObs_z']

    # Calculate t_obs from t_ast and current estimate of light_time
    light_time = df['LightTime'].values
    t_obs = t_ast + light_time / 1440.0             # 1440 is number of minutes in one day
    df['tObs'] = t_obs

    # Spline earth vectors and update qObs on DataFrame
    q_obs = get_earth_pos(ts=t_obs)
    df[cols_q_obs] = q_obs

    # Compute position difference and distance from asteroid to earth
    dq = q_ast - q_obs
    r = np.sqrt(np.sum(dq*dq, axis=1))

    # Compute light time and update it on DataFrame
    light_time = r / c
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
    # Load asteroid state vectors
    df = load_ast_vectors(n0=n0, n1=n1)

    # Reorder columns to put AsteroidID first
    cols = ['AsteroidID', 'TimeID'] + list(df.columns[2:])
    df = df[cols]

    # Sort by AsteroidID then TimeID
    df.sort_values(by=['AsteroidID', 'TimeID'], inplace=True)

    # Rename columns
    col_tbl = {
        'mjd': 'tAst',
        'qx': 'qAst_x',
        'qy': 'qAst_y',
        'qz': 'qAst_z',
        'vx': 'vAst_x',
        'vy': 'vAst_y',
        'vz': 'vAst_z',
    }
    df.rename(columns=col_tbl, inplace=True)

    # Column names used in calculations
    cols_q_ast = ['qAst_x', 'qAst_y', 'qAst_z']
    cols_q_obs = ['qObs_x', 'qObs_y', 'qObs_z']

    # The time the light leaves the asteroid
    t_ast = df.tAst.values

    # Extract asteroid position vectors
    q_ast = df[cols_q_ast].values

    # Add light time column; Initial guess is 0
    df['LightTime'] = 0.0

    # Initial guess is that observation time equals asteroid time
    df['tObs'] = t_ast

    # Add new colums to DataFrame for observer position; put in dummy value of 0 for iteration 0
    df[cols_q_obs] = 0.0

    # Experiments show that 4 iterations is sufficient to achieve full convergence
    # Error in light_time (number of minutes) around 5E-9 at this point
    for _ in range(4):
        light_time_iter(df=df, t_ast=t_ast, q_ast=q_ast)

    # Compute position difference and distance from asteroid to earth
    q_obs = df[cols_q_obs].values
    dq = q_ast - q_obs
    r = np.sqrt(np.sum(dq*dq, axis=1)).reshape((-1, 1))

    # Calculate the direction and save it to DataFrame
    u = dq / r
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
    columns = ['AsteroidID', 'TimeID', 'mjd', 'ux', 'uy', 'uz', 'LightTime']
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
    # parser.add_argument('n1', nargs='?', metavar='n1', type=int, default=1000,
    #                     help='the first asteroid number to process')
    
    # Unpack command line arguments
    args = parser.parse_args()
    
    # Block of asteroids to integrate
    n0: int = args.n0
    n1: int = n0 + args.n_ast
    # n1: int = args.n1

    # Report arguments
    print(f'Processing asteroid directions for asteroid number {n0} <= AsteroidID < {n1}...')
    # Set the batch size
    b: int = 200
    for k0 in tqdm(range(n0, n1, b)):
        k1: int = min(k0 + b, n1)
        # Calculate the direction and light time
        df = calc_dir_ast2obs(n0=k0, n1=k1)
        # Insert results to database
        insert_dir_ast2obs(df=df)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
