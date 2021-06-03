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
# from scipy.interpolate import CubicSpline

# Astronomy
import astropy
from astropy.units import au, minute

# Commandline arguments
import argparse

# Utility
from tqdm.auto import tqdm

# Local imports
from asteroid_data import load_ast_elements, load_ast_vectors
from asteroid_spline import spline_ast_data, make_spline_df
from astro_utils import deg2dist
from orbital_element import unpack_elt_df, elt2pos
from planets_interp import get_sun_pos
from db_utils import sp2df

# ********************************************************************************************************************* 
# Speed of light; express this in AU / min
c = astropy.constants.c.to(au / minute).value

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
        n0: First AsteroidID to process; inclusive.
        n1: Last AsteroidID to process; exclusive.
    OUTPUTS:

    """
    # Get asteroid directions in this time range
    sp_name = 'KS.GetAsteroidDirections'
    params = {
        'n0': n0,
        'n1': n1,
        'mjd0': mjd0,
        'mjd1': mjd1
    }
    ast_in = sp2df(sp_name=sp_name, params=params)

    # Rename columns
    col_tbl_ast = {
        'ux': 'uAst_x',
        'uy': 'uAst_y',
        'uz': 'uAst_z',
    }
    ast_in.rename(columns=col_tbl_ast, inplace=True)

    return ast_in

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
    ast_in = get_data_ast_dir(n0=n0, n1=n1)

    # Number of distinct asteroids in the input data
    N_ast = np.unique(ast_in.AsteroidID.values).size

    # Extract array of detection times
    t_obs = det.tObs.values

    # Data shape
    N_det = t_obs.shape[0]
    # print('Shape of input data:')
    # print(f'N_ast = {N_ast}')
    # print(f'N_det = {N_det}')

    # Spline asteroid directions to match these times
    ast = spline_ast_data(df_ast=ast_in, ts=t_obs, cols_spline=cols_spline_ast)

    # Add column for the DetectionID, 
    ast['DetectionID'] = np.tile(det.DetectionID.values, N_ast)

    # Rename column from mjd to tAst
    ast.rename(columns={'mjd':'tAst'}, inplace=True)

    # Extract detection_id and asteroid_id
    detection_id = ast.DetectionID.values
    light_time = ast.LightTime.values

    # Direction of detections
    u_obs = det[cols_u_obs].values.reshape((1, N_det, 3))

    # Direction of asteroids at detection times
    u_ast = ast[cols_u_ast].values.reshape((N_ast, N_det, 3))

    # The squared cartesian distance
    s2 = np.sum(np.square(u_ast - u_obs), axis=2)

    # Convert angular threshold to Cartesian distance
    # Add a cushion of 5 arc minutes for preliminary calculation
    arcmin_max_prelim = arcmin_max + 5.0
    deg_max: float = arcmin_max_prelim / 60.0
    s_max: float = deg2dist(deg_max)
    s2_max = np.square(s_max)

    # Build the mask for close interactions
    mask_2d = (s2 < s2_max)
    mask = mask_2d.flatten()

    # Arrays of near interactions only
    detection_id_near = detection_id[mask]
    light_time_near = light_time[mask]

    # Look up the distance on the 2D table
    s_near = np.sqrt(s2[mask_2d])

    # Look up detection data at near interactions only
    t_obs_near = det.loc[detection_id_near, 'tObs'].values
    u_obs_near = det.loc[detection_id_near, cols_u_obs].values
    q_obs_near = det.loc[detection_id_near, cols_q_obs].values

    # Copy splined asteroids that are near detections
    dna = ast[mask].copy()

    # Add columns with data from the detections
    dna['tObs'] = t_obs_near
    dna[cols_u_obs] = u_obs_near
    dna[cols_q_obs] = q_obs_near
    dna['LightTime'] = light_time_near
    dna['s'] = s_near

    # Reorder columns
    cols_dna_obs = ['AsteroidID', 'DetectionID', 'tObs'] + cols_u_obs + cols_q_obs
    cols_dna_ast = ['tAst', 'LightTime'] + cols_u_ast + ['s'] 
    cols_dna =  cols_dna_obs + cols_dna_ast
    dna = dna.reindex(columns=cols_dna)

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
    # The light time used in this frame
    light_time_used = (df.tObs.values - df.tAst.values) * 1440.0

    # The time of flight for light
    q_ast = df[cols_q_ast].values
    q_obs = df[cols_q_obs].values
    dq = q_obs - q_ast
    r = np.sqrt(np.sum(np.square(dq), axis=-1))
    light_time_true = r / c

    # The adjustment required; this should be added to tObs, i.e. df.tObs += t_adj
    t_adj = (light_time_used - light_time_true) / 1440.0
    
    # Report error in minutes if requested
    if verbose:
        mean_light_time_err = np.mean(np.abs(t_adj))*1440.0
        print(f'Mean light time error: {mean_light_time_err:5.2e} minutes')

    return t_adj
    
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
    u_obs = dna[cols_u_obs].values
    q_obs = dna[cols_q_obs].values
    t_ast = dna.tAst.values
    light_time = dna.LightTime.values

    # Refine estimate of t_ast and write to DataFrame
    t_ast = t_obs - light_time / 1440.0
    dna['tAst'] = t_ast

    # Calculate the asteroid number in this batch, n
    n = asteroid_id - n0

    # The row number - for indexing into splines
    row_num = dna.index.values

    # Get the asteroid vectors
    ast_vec = load_ast_vectors(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Build asteroid vectors spline
    cols_spline_q_ast = ['qx', 'qy', 'qz']
    spline_q_ast = make_spline_df(df=ast_vec, id_col='AsteroidID', time_col='mjd', cols_spline=cols_spline_q_ast)

    # Splined asteroid position
    q_ast =  spline_q_ast(t_ast)[n, row_num]
    dna[cols_q_ast] = q_ast

    # Sharpen the asteroid time
    t_adj = light_time_adj(df=dna, verbose=False)
    dna.tAst += t_adj

    # Save the improved light time
    light_time = (dna.tObs - dna.tAst)*1440.0
    dna.LightTime = light_time

    # Update the asteroid position at this new time
    q_ast =  spline_q_ast(t_ast)[n, row_num]
    dna[cols_q_ast] = q_ast

    # Displacement and distance from asteroid to observer
    dq = q_ast - q_obs
    r = np.sqrt(np.sum(np.square(dq), axis=-1, keepdims=True))

    # Direction from asteroid to observer
    u_ast = dq / r
    dna[cols_u_ast] = u_ast

    # Distance between two directions
    s = np.sqrt(np.sum(np.square(u_obs - u_ast), axis=-1))
    dna.s = s

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
    u_obs = dna[cols_u_obs].values
    q_obs = dna[cols_q_obs].values
    t_ast = dna.tAst.values
    light_time = dna.LightTime.values

    # Refine estimate of t_ast and write to DataFrame
    t_ast = t_obs - light_time / 1440.0
    dna['tAst'] = t_ast

    # Calculate the asteroid number in this batch, n
    n = asteroid_id - n0

    # The row number - for indexing into splines
    row_num = dna.index.values

    # Get the asteroid vectors
    ast_elt = load_ast_elements(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Compute cosine and sine of f for splining
    f_in = ast_elt.f.values
    ast_elt['fx'] = np.cos(f_in)
    ast_elt['fy'] = np.sin(f_in)

    # Build asteroid elements spline
    cols_spline_elt = ['a', 'e', 'inc', 'Omega', 'omega', 'fx', 'fy']
    spline_elt = make_spline_df(df=ast_elt, id_col='AsteroidID', time_col='mjd', cols_spline=cols_spline_elt)

    # Subfunction: get position from splined orbital elements
    def q_ast_func(t_ast: np.ndarray):
        # The splined orbital elements for the selected asteroid and t_ast
        # This is an array of shape (N_pair, 7) with columns matching cols_spline_elt
        elt =  spline_elt(t_ast)[n, row_num]
        # Unpack the elements from the splined array
        a = elt[:, 0]
        e = elt[:, 1]
        inc = elt[:, 2]
        Omega = elt[:, 3]
        omega = elt[:, 4]
        fx = elt[:, 5]
        fy = elt[:, 6]
        # Compute the splined f using atan2
        f = np.arctan2(fy, fx)
        # Compute q in the heliocentric frame
        q_hel = elt2pos(a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)
        # Position and velocity of the sun
        q_sun = get_sun_pos(ts=t_ast)
        # Position of the selected asteroid in barycentric frame is sum
        q_ast = q_hel + q_sun
        return q_ast

    # Splined asteroid position
    q_ast =  q_ast_func(t_ast)
    dna[cols_q_ast] = q_ast

    # Sharpen the asteroid time
    t_adj = light_time_adj(df=dna, verbose=False)
    dna.tAst += t_adj

    # Save the improved light time
    light_time = (dna.tObs - dna.tAst)*1440.0
    dna.LightTime = light_time

    # Update the asteroid position at this new time
    q_ast =  q_ast_func(t_ast)
    dna[cols_q_ast] = q_ast

    # Displacement and distance from asteroid to observer
    dq = q_ast - q_obs
    r = np.sqrt(np.sum(np.square(dq), axis=-1, keepdims=True))

    # Direction from asteroid to observer
    u_ast = dq / r
    dna[cols_u_ast] = u_ast

    # Distance between two directions
    s = np.sqrt(np.sum(np.square(u_obs - u_ast), axis=-1))
    dna.s = s

# ********************************************************************************************************************* 
def main():
    """Calculate the direction and light time of the selected batch of asteroid"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description=
    'Calculate pairs (DetectionID, AsteroidID) where an asteroid detection '
    'is near to the direction of a known asteroid.')
    parser.add_argument('d0', nargs='?', metavar='d0', type=int, default=0,
                        help='the first DetectionID to process')
    parser.add_argument('n_det', nargs='?', metavar='B', type=int, default=10000,
                        help='the number of detections to process in this batch'),
    parser.add_argument('arcmin_max', nargs='?', metavar='B', type=float, default=5.0,
                        help='the maximum distance in arcminutes between the observation and asteroid direction.'),
    
    # Unpack command line arguments
    args = parser.parse_args()
    d0: int = args.d0
    d1: int = d0 + args.n_det
    arcmin_max: float = args.arcmin_max

    # Report arguments
    print(f'Processing asteroid directions for DetectionID {d0} <= DetectionID < {d1} '
           'with s < {arcmin_max} arc minutes...')

    # Get asteroid detections in the input range
    det = get_data_detections(d0=d0, d1=d1)

    # Set the batch size
    b: int = 100
    # Loop over asteroid batches
    n0_loop: int = 0 
    n1_loop: int = 1000
    for n0 in tqdm(range(n0_loop, n1_loop, b)):
        n1: int = min(n0 + b, n1_loop)
        dna = asteroid_batch_prelim(det=det, n0=n0, n1=n1, arcmin_max=arcmin_max)
        asteroid_batch_vec(dna=dna, n0=n0, n1=n1, arcmin_max=arcmin_max)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
