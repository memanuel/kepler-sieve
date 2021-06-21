"""
Calculate set of asteroid detections near to known asteroid trajectories.
Example calls:
$ python detection_near_ast.py 0 1000

Functions in this module:
get_data_detections(d0, d1)
get_data_ast_dir(n0, n1)
asteroid_batch_prelim(det, n0, n1, arcmin_max)
light_time_adj(df, verbose)
asteroid_batch_vec(dna, n0, n1, arcmin_max)
asteroid_batch_elt(dna, n0, n1, arcmin_max)
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

# Typing
from typing import Callable

# Local imports
from asteroid_data import load_ast_elements, load_ast_vectors, load_ast_dir
from asteroid_spline import spline_ast_data, make_spline_df
from astro_utils import deg2dist
from orbital_element import elt2pos, anomaly_M2f
from planets_interp import get_sun_pos
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

# Calculate range of times required for splines
pad = 16.0
mjd0 = np.floor(np.min(dt.mjd)) - pad
mjd1 = np.ceil(np.max(dt.mjd)) + pad

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
def get_data_detections(d0: int, d1: int):
    """
    Get DataFrame of asteroid detections
    INPUTS:
        d0: First DetectionID to process; inclusive.
        d1: Last DetectionID to process; exclusive.
    OUTPUTS:
        det: DataFrame of detections; payload includes DetectionID, tObs, and uObs.
    """

    # Get block of asteroid detections (directions only)
    sp_name = 'KS.GetDetectionDirections'
    params = {
        'd0': d0,
        'd1': d1,
    }
    det = sp2df(sp_name=sp_name, params=params)

    # Rename columns in detections DataFrame
    col_tbl_obs = {
        'mjd': 'tObs',
        'ux': 'uObs_x',
        'uy': 'uObs_y',
        'uz': 'uObs_z',
    }
    det.rename(columns=col_tbl_obs, inplace=True)

    # Set index to be DetectionID
    det.set_index(keys='DetectionID', drop=False, inplace=True)

    return det

# ********************************************************************************************************************* 
def get_data_ast_dir(n0: int, n1: int):
    """
    Get DataFrame of asteroid directions
    INPUTS:
        n0:     First AsteroidID to process; inclusive.
        n1:     Last AsteroidID to process; exclusive.
    OUTPUTS:
        ast_in: DataFrame including asteroid directions and light time for selected asteroids and date range.
    """
    # Get asteroid directions in this time range
    ast_in = load_ast_dir(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Rename columns to disambiguate asteroid direction from observer direction
    col_tbl_ast = {
        'ux': 'uAst_x',
        'uy': 'uAst_y',
        'uz': 'uAst_z',
    }
    ast_in.rename(columns=col_tbl_ast, inplace=True)

    return ast_in

# ********************************************************************************************************************* 
def make_ast_num_df(df: pd.DataFrame):
    """Build a DataFrame mapping an AsteroidID to an ast_num field, n, that counts asteroids in this batch"""

    # Number of distinct asteroids in the input data
    ast_unq = np.unique(df.AsteroidID.values)
    N_ast = ast_unq.size
    
    # Build data frame with index = AsteroidID, payload n (asteroid number in this batch)
    dfa_tbl = {
        'AsteroidID': ast_unq,
        'n': np.arange(N_ast, dtype=np.int)
    }
    dfa = pd.DataFrame(dfa_tbl)    
    dfa.set_index(keys='AsteroidID', drop=False, inplace=True)
    return dfa
    
# ********************************************************************************************************************* 
def asteroid_batch_prelim(det: pd.DataFrame, n0: int, n1: int, arcmin_max: float):
    """
    Process a batch of near asteroid calculations.
    Preliminary calculation using a spline of asteroid directions.
    Spline directly on components of position; ignore light time and set t_ast = t_obs.
    Sufficient to identify which interactions (DetectionID, AsteroidID) are close together.
    INPUTS:
        det:            DataFrame of asteroid detections.
        n0:             First asteroid to process, inclusive.
        n1:             Last asteroid to process, exclusive.
        arcmin_max:     Largest angular distance between a detection and an asteroid direction to include.
    OUTPUTS:
        dna:            DataFrame of detections near asteroids. Columns include:
        AsteroidID, DetectionID, tObs, uObs_x, uObs_y, uObs_z, qObs_x, qObs_y, qObs_z, uAst_x, uAst_y, uAst_z, s
    """

    # Get asteroid directions
    ast_dir = get_data_ast_dir(n0=n0, n1=n1)

    # Array of distinct asteroids
    asteroid_id_unq = np.unique(ast_dir.AsteroidID.values)

    # Shape of data
    N_det = det.shape[0]
    N_ast = asteroid_id_unq.shape[0]

    # Build spline for asteroid direction and light time
    spline_u_ast = make_spline_df(df=ast_dir, cols_spline=cols_u_ast, time_col='tObs', id_col='AsteroidID')
    spline_LT = make_spline_df(df=ast_dir, cols_spline=['LightTime'], time_col='tObs', id_col='AsteroidID')

    # First block  of interactions: tile the detections DataFrame N_ast times
    da = pd.concat([det] * N_ast)
    da.reset_index(drop=True, inplace=True)

    # Add column to dna for AsteroidID
    da.insert(1, 'AsteroidID', np.tile(asteroid_id_unq, N_det))

    # Extract arrays of DetectionID and AsteroidID
    # detection_id = da.DetectionID.values
    asteroid_id = da.AsteroidID.values

    # Extract array of tObs and uObs
    t_obs = da.tObs.values
    u_obs = da[cols_u_obs].values

    # Spline light time
    light_time = spline_LT(x=t_obs, y=asteroid_id).flatten()
    # First estimate of asteroid time
    t_ast = t_obs - light_time / 1440.0    

    # Spline asteroid directions at the observation times
    u_ast = spline_u_ast(x=t_obs, y=asteroid_id)
    # Add columns for tAst and uAst
    da['tAst'] = t_ast
    da[cols_u_ast] = u_ast

    # Calculate cartesian distance and save to DataFrame
    s = np.sqrt(np.sum(np.square(u_ast - u_obs), axis=-1))
    da['s'] = s
    # Save LightTime to DataFrame
    da['LightTime'] = light_time

    # Convert angular threshold to Cartesian distance
    # Add a cushion for preliminary calculation
    arcmin_max_prelim = arcmin_max + arcmin_margin
    deg_max: float = arcmin_max_prelim / 60.0
    s_max: float = deg2dist(deg_max)

    # Build the mask for close interactions
    mask = (s < s_max)
    dna = da[mask].copy()

    # Reset the index so row_number works as expected downstream
    dna.reset_index(drop=True, inplace=True)

    return dna

# ********************************************************************************************************************* 
def light_time_adj(df: pd.DataFrame, verbose: bool):
    """
    Calculate adjustment required to light time
    INPUTS:
        df:     DataFrame including tObs, qObs, tAst, qAst
    OUTPUTS:
        t_adj:  Adjustment required to match on light time.
    """
    # The light time used to build this frame
    light_time_used = (df.tObs.values - df.tAst.values) * 1440.0

    # The time of flight for light
    q_ast = df[cols_q_ast].values
    q_obs = df[cols_q_obs].values
    dq = q_obs - q_ast
    r = np.sqrt(np.sum(np.square(dq), axis=-1))
    light_time_true = r / c

    # The adjustment required; this should be added to tAst, i.e. df.tAst += t_adj
    # Note that adding to tAst will reduce the calculated light time in the next iteration,
    # hence the sign flip in the formula below.
    t_adj = (light_time_used - light_time_true) / 1440.0
    
    # Report error in minutes if requested
    if verbose:
        mean_light_time_err = np.mean(np.abs(t_adj))*1440.0
        print(f'Mean light time error: {mean_light_time_err:5.2e} minutes')

    return t_adj
    
# ********************************************************************************************************************* 
def asteroid_batch_finish(dna: pd.DataFrame, q_ast_func: Callable, arcmin_max: float):
    """
    Complete processing of a batch of detections near asteroids.
    Same whether batch is computed by interpolating state vectors or orbital elements.
    INPUTS:
        dna:            DataFrame of detections near asteroids with preliminary calculations.
        q_ast_func:     Function taking 
    """
    # Unpack arrays
    t_ast = dna.tAst.values
    q_obs = dna[cols_q_obs].values
    u_obs = dna[cols_u_obs].values

    # Splined asteroid position
    q_ast =  q_ast_func(t_ast)
    dna[cols_q_ast] = q_ast

    # Sharpen the asteroid time
    t_adj = light_time_adj(df=dna, verbose=False)
    dna.tAst += t_adj

    # Update the asteroid position at this new time
    q_ast =  q_ast_func(t_ast)
    dna[cols_q_ast] = q_ast

    # Displacement and distance from asteroid to observer
    dq = q_ast - q_obs
    r = np.sqrt(np.sum(np.square(dq), axis=-1, keepdims=True))

    # Save the improved light time
    light_time = r / c
    dna.LightTime = light_time

    # Direction from asteroid to observer
    u_ast = dq / r
    dna[cols_u_ast] = u_ast

    # Calculate distance between the two directions, s, and save it to DataFrame
    s = np.sqrt(np.sum(np.square(u_obs - u_ast), axis=-1))
    dna.s = s

    # Mask down to rows below the threshold
    deg_max: float = arcmin_max / 60.0
    s_max: float = deg2dist(deg_max)
    mask = (dna.s < s_max)
    dna = dna[mask]

    # Reset the index
    dna.reset_index(drop=True, inplace=True)

# ********************************************************************************************************************* 
def asteroid_batch_vec(dna: pd.DataFrame, n0: int, n1: int, arcmin_max: float):
    """
    Sharpen preliminary asteroid direction calculations by splining on asteroid state vectors.
    INPUTS:
        dna:            Preliminary calculation of asteroid detections near known asteroids.
        n0:             First asteroid to process; inclusive.
        n1:             Last asteroid to process; exclusive.
        arcmin_max:     Largest angular distance between a detection and an asteroid direction to include.
    OUTPUTS:
        None.  Modifies dna in place.
    """

    # Unpack DataFrame to numpy arrays
    asteroid_id = dna.AsteroidID.values
    t_obs = dna.tObs.values
    t_ast = dna.tAst.values
    light_time = dna.LightTime.values

    # Refine estimate of t_ast and write to DataFrame
    t_ast = t_obs - light_time / 1440.0
    dna['tAst'] = t_ast

    # Get the asteroid vectors
    ast_vec = load_ast_vectors(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Build asteroid vectors spline
    cols_spline_q_ast = ['qx', 'qy', 'qz']
    spline_q_ast = make_spline_df(df=ast_vec, id_col='AsteroidID', time_col='mjd', cols_spline=cols_spline_q_ast)

    # Subfunction: get position from splined orbital elements
    def q_ast_func(t_ast: np.ndarray):
        q_ast = spline_q_ast(x=t_ast, y=asteroid_id)
        return q_ast

    # Complete the post-processing step
    asteroid_batch_finish(dna=dna, q_ast_func=q_ast_func, arcmin_max=arcmin_max)

# ********************************************************************************************************************* 
def asteroid_batch_elt(dna: pd.DataFrame, n0: int, n1: int, arcmin_max: float):
    """
    Sharpen asteroid direction calculations by splining on asteroid orbital elements.
    INPUTS:
        dna:            Preliminary calculation of asteroid detections near known asteroids.
        n0:             First asteroid to process; inclusive.
        n1:             Last asteroid to process; exclusive.
        arcmin_max:     Largest angular distance between a detection and an asteroid direction to include.
    OUTPUTS:
        None.  Modifies dna in place.
    """

    # Unpack DataFrame to numpy arrays
    asteroid_id = dna.AsteroidID.values
    t_obs = dna.tObs.values
    t_ast = dna.tAst.values
    light_time = dna.LightTime.values

    # Refine estimate of t_ast and write to DataFrame
    t_ast = t_obs - light_time / 1440.0
    dna['tAst'] = t_ast

    # Get the asteroid vectors
    ast_elt = load_ast_elements(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Compute cosine and sine of f for splining
    f_in = ast_elt.f.values
    ast_elt['fx'] = np.cos(f_in)
    ast_elt['fy'] = np.sin(f_in)

    # Build asteroid elements spline
    cols_spline_elt = ['a', 'e', 'inc', 'Omega', 'omega', 'M']
    spline_elt = make_spline_df(df=ast_elt, id_col='AsteroidID', time_col='mjd', cols_spline=cols_spline_elt)

    # Subfunction: get position from splined orbital elements
    def q_ast_func(t_ast: np.ndarray):
        # The splined orbital elements for the selected asteroid and t_ast
        # This is an array of shape (N_pair, 7) with columns matching cols_spline_elt
        elt = spline_elt(x=t_ast, y=asteroid_id)

        # Unpack the elements from the splined array
        a = elt[:, 0]
        e = elt[:, 1]
        inc = elt[:, 2]
        Omega = elt[:, 3]
        omega = elt[:, 4]
        M = elt[:, 5]
        # Compute f from M and e
        f = anomaly_M2f(M=M, e=e)

        # Compute q in the heliocentric frame
        q_hel = elt2pos(a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)
        # Position and velocity of the sun
        q_sun = get_sun_pos(ts=t_ast)
        # Position of the selected asteroid in barycentric frame is sum
        q_ast = q_hel + q_sun
        return q_ast

    # Complete the post-processing step
    asteroid_batch_finish(dna=dna, q_ast_func=q_ast_func, arcmin_max=arcmin_max)

# ********************************************************************************************************************* 
def insert_det_near_ast(dna: pd.DataFrame):
    """Insert asteroid direction calculations """
    # Rename column tAst back to mjd to match DB schema
    dna.rename(columns={'tAst':'mjd'}, inplace=True)

    # Arguments to df2db
    schema = 'KS'
    table = 'DetectionNearAsteroid'
    columns = ['DetectionID', 'AsteroidID', 's', 'LightTime']
    chunksize = 2**19
    verbose = False
    progbar = False

    # Dispatch to df2db
    df2db(df=dna, schema=schema, table=table, columns=columns, chunksize=chunksize, verbose=verbose, progbar=progbar)
    
# ********************************************************************************************************************* 
def process_batch(det: pd.DataFrame, n0: int, n1: int, arcmin_max: float):
    """Process one batch -- a single block of asteroids"""
    # Run preliminary search
    dna = asteroid_batch_prelim(det=det, n0=n0, n1=n1, arcmin_max=arcmin_max)
    # Refine using vectors (same answer as orbital elements, simpler and slightly faster)
    asteroid_batch_vec(dna=dna, n0=n0, n1=n1, arcmin_max=arcmin_max)
    # Insert into the database
    insert_det_near_ast(dna=dna)
    # Return the number of rows inserted
    return dna.shape[0]

# ********************************************************************************************************************* 
def process_all_asteroids(det: pd.DataFrame, arcmin_max: float, b: int):
    """Process a set of detections against all known asteroids"""
    # Count number of rows
    row_count: int = 0
    
    # Loop over numbered asteroids
    print(f'Processing numbered asteroids with AsteroidID between {n0_num} and {n1_num}')
    k0: int = n0_num // b
    k1: int = n1_num // b
    for k in tqdm(range(k0, k1)):
        n0: int = b*(k+0)
        n1: int = b*(k+1)
        row_count += process_batch(det=det, n0=n0, n1=n1, arcmin_max=arcmin_max)

    # Loop over unnumbered asteroids
    print(f'Processing numbered asteroids with AsteroidID between {n0_unum} and {n1_unum}')
    k0: int = n0_unum // b
    k1: int = n1_unum // b
    for k in tqdm(range(k0, k1)):
        n0: int = b*(k+0)
        n1: int = b*(k+1)
        row_count += process_batch(det=det, n0=n0, n1=n1, arcmin_max=arcmin_max)

    return row_count

# ********************************************************************************************************************* 
def main():
    """Calculate the direction and light time of the selected batch of asteroid"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description=
    'Calculate pairs (DetectionID, AsteroidID) where an asteroid detection '
    'is near to the direction of a known asteroid.')
    parser.add_argument('d0', nargs='?', metavar='d0', type=int, default=0,
                        help='the first DetectionID to process')
    parser.add_argument('n_det', nargs='?', metavar='B', type=int, default=5000000,
                        help='the number of detections to process in this batch'),
    parser.add_argument('arcmin_max', nargs='?', metavar='B', type=float, default=60.0,
                        help='the maximum distance in arcminutes between the observation and asteroid direction.'),

    # Unpack command line arguments
    args = parser.parse_args()
    d0: int = args.d0
    d1: int = d0 + args.n_det
    arcmin_max: float = args.arcmin_max

    # Report arguments
    print(f'Processing asteroid directions for DetectionID in between {d0} and {d1} '
          f'with s < {arcmin_max} arc minutes...')

    # Set batch size for asteroids and detecions
    b_det: int = 200000
    b_ast: int = 1000

    # Loop over detection batches
    k0: int = d0 // b_det
    k1: int = d1 // b_det
    for k in range(k0, k1):
        # Range of detections
        d0_k = k*b_det
        d1_k = d0_k + b_det
        # Status
        print()
        print_stars()
        print(f'Processing detection IDs between {d0} and {d1}...')
        # Get asteroid detections in the input range
        det = get_data_detections(d0=d0_k, d1=d1_k)
        # Process all asteroids 
        row_count = process_all_asteroids(det=det, arcmin_max=arcmin_max, b=b_ast)
        # Report total number of rows inserted
        print(f'Inserted {row_count} rows of detections near a known asteroid.')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
