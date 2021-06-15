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
from asteroid_spline import make_spline_ast_dir
from asteroid_direction import prep_ast_block
from sky_patch import dir2SkyPatchID
from db_utils import sp2df, df2db
from astro_utils import dist2deg

# Minutes in one day
mpd: int = 1440

# ********************************************************************************************************************* 
def calc_ast_skypatch(n0: int, n1: int, mjd0: int, mjd1: int, interval_min: int):
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
    u_ast, light_time = spline_u(t_obs, asteroid_id)
    # Calculate the SkyPatchID
    N_sp = 1024
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
def test_skypatch_dir(name1: str, name2: str, u1: np.ndarray, u2: np.ndarray, verbose: bool=False) -> float:
    """
    Report 
    INPUTS:
        name1: Descriptive name of the first source, e.g. 'spline'
        name2: Descriptive name of the second source, e.g. 'skypatch'
        u1:    Array of directions from source 1; shape Nx3
        u2:    Array of directions from source 2; shape Nx3
        verbose: Whether to report results to console
    """
    # Difference in unit directions
    du = u2 - u1
    du_norm = np.sqrt(np.sum(np.square(du), axis=-1))
    du_deg = dist2deg(du_norm)
    du_sec = du_deg * 3600.0

    # Calculate mean, median and max difference in degrees
    du_mean = np.mean(du_sec)
    du_median = np.median(du_sec)
    du_max = np.max(du_sec)

    if verbose:
        print(f'Angle Difference: {name2} vs. {name1} in arc seconds')
        print(f'*Mean  : {du_mean:9.1f}*')
        print(f' Median: {du_median:9.1f}')
        print(f' Max   : {du_max:9.1f}')

    # Return the difference of direction vectors in seconds of arc
    return du_mean
    
# ********************************************************************************************************************* 
def main():

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

    # Set the time range by policy to approximately match available detections
    mjd0: int = 58000
    mjd1: int = 60000
    interval_min: int = 5

    # Report arguments
    print(f'Processing asteroid sky patch for asteroid number in range [{n0}, {n1})...')
    # Set the batch size
    b: int = 100
    k0: int = n0 // b
    k1: int = max(n1 // b, k0+1)
    for k in tqdm(range(k0, k1)):
        # Start and end of this batch
        n0_i = k*b
        n1_i = min(n0_i + b, n1)
        # Calculate the direction and light time
        df = 0
        # Insert results to database
        
