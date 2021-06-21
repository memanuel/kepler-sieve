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
from db_utils import sp2df, df2db
from utils import print_stars

# ********************************************************************************************************************* 
# Speed of light; express this in AU / min
c = astropy.constants.c.to(au / minute).value

# Margin for preliminary search
arcmin_margin: float = 15.0

# Column groups - detections
cols_u_obs = ['uObs_x', 'uObs_y', 'uObs_z']
cols_q_obs = ['qObs_x', 'qObs_y', 'qObs_z']

# Column groups - asteroid
cols_u_ast = ['uAst_x', 'uAst_y', 'uAst_z']
cols_q_ast = ['qAst_x', 'qAst_y', 'qAst_z']
cols_v_ast = ['vAst_x', 'vAst_y', 'vAst_z']
cols_spline_ast = cols_u_ast + ['LightTime']

# Get detection times
dt = sp2df(sp_name='KS.GetDetectionTimes', params=dict())

# Get range of asteroids
anr = sp2df(sp_name='KS.GetAsteroidNumberRange')
n0_num = anr.loc[0, 'AsteroidID_0']
n1_num = anr.loc[0, 'AsteroidID_1']
n0_unum = anr.loc[1, 'AsteroidID_0']
n1_unum = anr.loc[1, 'AsteroidID_1']

# # TEST
# n0_num = 0
# n1_num = 1000
# n0_unum = 1000000
# n1_unum = 1001000

# ********************************************************************************************************************* 
def calc_det_near_ast(n0: int, n1: int, arcmin_max: float):
    """Calculate asteroid detections near a batch of known asteroids """
    # Get DataFrame of candidate interactions

    # Calculate the the asteroid directions

    
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
    parser.add_argument('arcmin_max', nargs='?', metavar='B', type=float, default=15.0,
                        help='the maximum distance in arcminutes between the observation and asteroid direction.'),

    # Unpack command line arguments
    args = parser.parse_args()
    jn: int = args.jn
    sz: int = args.sz
    arcmin_max: float = args.arcmin_max

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
    print(f'Threshold distance for a detection to be near: {arcmin_max:6.1f} arc minutes')

    # Set the batch size
    b: int = 100
    # Loop over asteroids in batches of size b
    for i0 in tqdm(range(i0_job, i1_job, b)):
        # Ending asteroid index of this batch
        i1: int = min(i0 + b, i1_job)
        # Asteroid numbers in this batch
        n0: int = asteroid_id[i0]
        n1: int = asteroid_id[i1]
        # Calculate the asteroid skypatch IDs
        df = insert_det_near_ast(n0=n0, n1=n1, arcmin_max=arcmin_max)
        # Insert results to database
        insert_det_near_ast(df=df, jn=jn)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
