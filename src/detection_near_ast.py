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
from asteroid_spline import make_spline_ast_vec
from asteroid_direction import calc_distance, calc_dir_linear
from astro_utils import sec2dist
from db_utils import sp2df, df2db

# ********************************************************************************************************************* 
# Speed of light; express this in AU / min
c = astropy.constants.c.to(au / minute).value

# Margin for preliminary search
arcmin_margin: float = 15.0

# Column groups
cols_q_obs = ['qObs_x', 'qObs_y', 'qObs_z']
cols_q_ast = ['qAst_x', 'qAst_y', 'qAst_z']
cols_v_ast = ['vAst_x', 'vAst_y', 'vAst_z']
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
        'AsteroidID_1': n1
    }
    df = sp2df(sp_name=sp_name, params=params)

    # Handle corner case where t_obs has no entries
    if df.shape[0]==0:
        return df
    
    # Extract arrays from candidates DataFrame
    asteroid_id: np.ndarray = df.AsteroidID.values
    t_obs: np.ndarray = df.tObs.values
    q_obs: np.ndarray = df[cols_q_obs].values
    u_obs = df[cols_u_obs].values

    # Get date range of these observation times; it is padded downstream by make_spline_ast_vec, don't pad again
    mjd0: float = np.min(t_obs)
    mjd1: float = np.max(t_obs)

    # Build spline of asteroid vectors
    spline_ast_vec = make_spline_ast_vec(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Calculate asteroid directions and light time from spline and add to DataFrame
    # Note: while it's tempting to use spline_ast_dir, this splines the direction from geocenter only;
    # it doesn't handle arbitrary observatory locations on Earth.
    # That approach would be close, but introduces errors on the order of 2-3 arc seconds.
    q_ast, v_ast = spline_ast_vec(ts=t_obs, asteroid_id=asteroid_id)

    # Delegate to calc_dir_linear
    u_ast, delta = calc_dir_linear(q_tgt=q_ast, v_tgt=v_ast, q_obs=q_obs)

    # Light time from delta
    light_time = delta / c

    # Calculate distance between the detection and the asteroid
    s = calc_distance(u_obs, u_ast)

    # Save the distance and light time to the DataFrame; only columns we want to write to DB later
    df['s'] = s
    df['LightTime'] = light_time

    # Mask down to only rows withing the maximum distance
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
    parser.add_argument('arcsec_max', nargs='?', metavar='B', type=float, default=10.0,
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
    b: int = 20
    # Loop over asteroids in batches of size b
    for i0 in tqdm(range(i0_job, i1_job, b)):
        # Ending asteroid index of this batch
        i1: int = min(i0 + b, i1_job)
        # Job number for this batch; based on sz=25000 used in AsteroidSkyPatch job
        jn_asp: int = i0 // 25000
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
