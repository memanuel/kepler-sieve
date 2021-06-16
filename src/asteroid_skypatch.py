"""
Calculate sky patch occupied by known asteroids over time and populate DB table KS.AsteroidSkyPatch.
Example calls:
$ python asteroid_skypatch.py 0 1000

Functions in this module:

main()

Michael S. Emanuel
2021-06-15
"""

# Core
import numpy as np
import pandas as pd

# Commandline arguments
import argparse

# Utility
from tqdm.auto import tqdm

# Local
from asteroid_data import get_asteroid_ids
from asteroid_spline import make_spline_ast_dir
from asteroid_direction import prep_ast_block
from sky_patch import dir2SkyPatchID
from db_utils import sp2df, df2db

# ********************************************************************************************************************* 
# SkyPatch grid size
N_sp = 1024

# Minutes in one day
mpd: int = 1440

# ********************************************************************************************************************* 
# Calculate the skypatch of known asteroids in contiguous segments.
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def calc_ast_skypatch(n0: int, n1: int, mjd0: int, mjd1: int, interval_min: int) -> pd.DataFrame:
    """
    Calculate the skypatch over time for a block of known asteroids.
    Output is a DataFrame matching the specifications of DB table KS.AsteroidSkyPatch.
    INPUTS:
        n0:             The first asteroid to process (inclusive)
        n1:             The last asteroid to process (inclusive)
        mjd0:           The first time to process as an MJD; (inclusive)
        mjd1:           The last time to process as an MJD; (inclusive)
        interval_min:   Interval in minutes for the initial schedule; resolution of the output
    OUTPUTS:
        df:             DataFrame with the sky patch of this asteroid.  Columns include
                        AsteroidID, Segment, SkyPatchID, TimeID_0, TimeID_1
    """
    # Prepare asteroid block
    t_obs, asteroid_id = prep_ast_block(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1, interval_min=interval_min)
    # The TimeID
    time_id = np.round(t_obs * mpd).astype(np.int64)

    # Build spline of asteroid direction
    spline_u = make_spline_ast_dir(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)
    # Calculate asteroid directions from position and velocity
    u_ast, _ = spline_u(t_obs, asteroid_id)
    # Calculate the SkyPatchID
    sky_patch_id = dir2SkyPatchID(dir=u_ast, N=N_sp)

    # Calculate the asteroid number in this batch from asteroid_id
    asteroid_id_unq = np.unique(asteroid_id)
    ast_bn = np.searchsorted(asteroid_id_unq, asteroid_id)

    # Number of rows
    N_row = asteroid_id.shape[0]
    # Rows with a new asteroid
    is_new_ast = np.ones(N_row, dtype=np.bool)
    is_new_ast[1:N_row] = (asteroid_id[1:N_row] != asteroid_id[0:N_row-1])
    # Rows with a new sky patch entry
    is_new_skypatch = np.ones(N_row, dtype=np.bool)
    is_new_skypatch[1:N_row] = (sky_patch_id[1:N_row] != sky_patch_id[0:N_row-1])
    # Mask for first pass; these are rows that begin a new segment
    mask = (is_new_ast | is_new_skypatch)

    # The segment number - counter of distinct segments for each asteroid
    sn_cum = np.cumsum(mask)
    sn_base = sn_cum[is_new_ast]
    sn = sn_cum - sn_base[ast_bn]

    # Wrap into DataFrame
    df_tbl = {
        'AsteroidID': asteroid_id,
        'Segment': sn,
        'TimeID': time_id,
        'SkyPatchID': sky_patch_id,
    }
    df_all = pd.DataFrame(df_tbl)

    # Create DataFrame for just the rows that start segments (output rows)
    cols = ['AsteroidID', 'Segment', 'SkyPatchID', 'TimeID']
    df = df_all.loc[mask, cols].copy().reset_index(drop=True)

    # Index of next entry on the output DataFrame
    idx0 = df.index.values
    idx1 = idx0+1
    idx1[-1] = idx0[0]

    # Populate TimeID_0 and TimeID_1
    df['TimeID_0'] = df.TimeID[idx0].values
    df['TimeID_1'] = df.TimeID[idx1].values
    df.drop('TimeID', axis='columns', inplace=True)

    # Flag for last segment on each asteroid
    is_last_segment = (df.AsteroidID[idx0].values != df.AsteroidID[idx1].values)
    # Set TimeID_1 on last segment of each asteroid 200 minutes ahead of TimeID_0 (average found on first 100 asteroids)
    df.loc[is_last_segment, 'TimeID_1'] = df.loc[is_last_segment, 'TimeID_0'] + 200

    return df

# ********************************************************************************************************************* 
# Save asteroid skypatch calculations to DB table KS.AsteroidSkyPatch.
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def insert_ast_skypatch(df: pd.DataFrame, jn: int) -> None:
    """
    Insert asteroid direction calculations to database
    INPUTS:
        df:     The DataFrame of rows to be inserted
        jn:     The job number; supports parallel processing.
    """
    # Arguments to df2db
    schema = 'KS'
    table = f'AsteroidSkyPatch_Stage_{jn:02d}'
    columns = ['AsteroidID', 'Segment', 'SkyPatchID', 'TimeID_0', 'TimeID_1']
    chunksize = 2**19
    verbose = False
    progbar = False

    # Dispatch to df2db
    df2db(df=df, schema=schema, table=table, columns=columns, chunksize=chunksize, verbose=verbose, progbar=progbar)

# *************************************************************************************************
# Console program populates DB table KS.AsteroidSkyPatch
# *************************************************************************************************

# ********************************************************************************************************************* 
def main():

    # Process command line arguments
    parser = argparse.ArgumentParser(description='Populates DB table KS.AsteroidSkypatch.')
    parser.add_argument('jn', nargs='?', metavar='jn', type=int, 
                        help='the job number; job jn processes asteroids in block [jn*sz, (jn+1)*sz)')
    parser.add_argument('sz', nargs='?', metavar='sz', type=int, default=25000,
                        help='the number of asteroids to process in this job'),
    
    # Unpack command line arguments
    args = parser.parse_args()
    jn: int = args.jn
    sz: int = args.sz

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

    # Set the time range to approximately match available data on detections
    # Exact date range on 15Jun2021 is MJD 58270-59295
    # Round this to MJD 58000-60000
    mjd0: int = 58000
    mjd1: int = 60000
    # Use a resolution of every 15 minutes.
    # The average segment length on first 100 asteroids is about 200 minutes.
    interval_min: int = 15

    # Report arguments
    print(f'Processing asteroid sky patch for job number {jn} and job size {sz}...')
    print(f'Range of asteroid index: [{i0_job}, {i1_job})')
    print(f'Range of AsteroidID:     [{n0_job}, {n1_job})')

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
        df = calc_ast_skypatch(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1, interval_min=interval_min)
        # Insert results to database
        insert_ast_skypatch(df=df, jn=jn)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
