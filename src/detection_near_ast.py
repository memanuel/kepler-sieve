"""
Calculate set of asteroid detections near to known asteroid trajectories.
Example calls:
$ python detection_near_ast.py 0 1000

Functions in this module:
main()

Michael S. Emanuel
2021-06-02
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
from asteroid_data import get_asteroid_ids
from asteroid_spline import make_spline_ast_dir
from asteroid_direction import calc_distance
from astro_utils import sec2dist
from db_utils import sp2df, df2db
from utils import print_stars

# ********************************************************************************************************************* 
# Speed of light; express this in AU / min
c = astropy.constants.c.to(au / minute).value

# Margin for preliminary search
arcmin_margin: float = 15.0

# Column groups
cols_u_obs = ['uObs_x', 'uObs_y', 'uObs_z']
cols_u_ast = ['uAst_x', 'uAst_y', 'uAst_z']

# ********************************************************************************************************************* 
def calc_det_near_ast(n0: int, n1: int, arcsec_max: float):
    """Calculate asteroid detections near a batch of known asteroids """
    # Convert threshold from arc seconds to distance
    s_max = sec2dist(arcsec_max)

    # Get batch of candidates
    sp_name = 'KS.GetDetectionNearAstCand'
    params = {
        'AsteroidID_0': n0,
        'AsteroidID_1': n1,
    }
    df = sp2df(sp_name=sp_name, params=params)

    # Extract arrays of observation times and asteroid_id
    t_obs: np.ndarray = df.tObs.values
    asteroid_id: np.ndarray = df.AsteroidID.values

    # Extract array of detection directions
    u_obs = df[cols_u_obs].values

    # Get date range
    mjd0: float = np.min(t_obs)
    mjd1: float = np.max(t_obs)

    # Build spline of asteroid direction
    spline_u = make_spline_ast_dir(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Calculate asteroid directions and light time from spline; add light time to DataFrame
    u_ast, light_time = spline_u(t_obs, asteroid_id)
    # df[cols_u_ast] = u_ast
    df['LightTime'] = light_time

    # Calculate distance between the detection and the asteroid and add it to DataFrame
    s = calc_distance(u_obs, u_ast)
    df['s'] = s

    # Mask down to only rows within the maximum distance
    mask = (s < s_max)
    df = df[mask].reset_index(drop=True)
    return df
    
# ********************************************************************************************************************* 
def insert_det_near_ast(df: pd.DataFrame):
    """Insert asteroid direction calculations """
    # Arguments to df2db
    schema = 'KS'
    table = 'DetectionNearAsteroid'
    columns = ['DetectionID', 'AsteroidID', 's', 'LightTime']
    chunksize = 2**19
    verbose = False
    progbar = False

    # Dispatch to df2db
    df2db(df=df, schema=schema, table=table, columns=columns, chunksize=chunksize, verbose=verbose, progbar=progbar)
    
# ********************************************************************************************************************* 
def main():
    """Calculate the detection near asteroids in the selected batch of asteroids"""
    # Process command line arguments
    parser = argparse.ArgumentParser(description='Populates DB table KS.AsteroidSkypatch.')
    parser.add_argument('jn', nargs='?', metavar='jn', type=int, 
                        help='the job number; job jn processes asteroids in block [jn*sz, (jn+1)*sz)')
    parser.add_argument('sz', nargs='?', metavar='sz', type=int, default=25000,
                        help='the number of asteroids to process in this job'),
    parser.add_argument('arcsec_max', nargs='?', metavar='B', type=float, default=60.0,
                        help='the maximum distance in arc seconds between the observation and asteroid direction.'),

    # Unpack command line arguments
    args = parser.parse_args()
    jn: int = args.jn
    sz: int = args.sz
    arcsec_max: float = args.arcsec_max

    # Get array of all known asteroid_ids
    asteroid_id = get_asteroid_ids()
    ast_count = asteroid_id.shape[0]
    # Append dummy entry to end of asteroid_id to avoid overflowing on the last asteroid
    ast_id_dummy = np.max(asteroid_id)+1
    asteroid_id = np.append(asteroid_id, [ast_id_dummy])
    # Get start and end of block of sz asteroids for this job
    i0_job = sz*jn
    i1_job = min(i0_job+sz, ast_count)
    # Handle edge case where jn*sz is too big
    if i0_job > ast_count:
        print(f'There are only {ast_count} asteroids and this job starts at {i0_job}. Quitting early!')
        quit()
    # Get asteroid_id for start and end of job; used only in status message, not loop control
    n0_job = asteroid_id[i0_job]
    n1_job = asteroid_id[i1_job]

    # Report arguments
    print(f'Processing detections near asteroids for job number {jn} and job size {sz}...')
    print(f'Range of asteroid index: [{i0_job}, {i1_job})')
    print(f'Range of AsteroidID:     [{n0_job}, {n1_job})')
    print(f'Threshold distance:     {arcsec_max:5.1f} arc seconds')

    # Set the batch size
    b: int = 10
    # Loop over asteroids in batches of size b
    for i0 in tqdm(range(i0_job, i1_job, b)):
        # Ending asteroid index of this batch
        i1: int = min(i0 + b, i1_job)
        # Asteroid numbers in this batch
        n0: int = asteroid_id[i0]
        n1: int = asteroid_id[i1]
        # Calculate the asteroid skypatch IDs
        df = calc_det_near_ast(n0=n0, n1=n1, arcsec_max=arcsec_max)
        # Insert results to database
        insert_det_near_ast(df=df)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
