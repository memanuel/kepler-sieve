"""
Harvard IACS Masters Thesis
ZTF Data
Calculate distances between ZTF observations and candidate orbital elements.
Load DataFrame of ZTF observations and their interactions candidate elements within a threshold.

Michael S. Emanuel
27-Feb-2020
"""

# Libraries for getting Alerce data out of ZTF2 database
import json
import psycopg2
from alerce.api import AlerceAPI

# Standard libraries
import numpy as np
import pandas as pd

# Astronomy related
from astropy.units import deg
from scipy.interpolate import CubicSpline

# Utility
import os
import datetime
from datetime import date
from tqdm.auto import tqdm

# MSE imports
from utils import range_inc
from astro_utils import date_to_mjd, deg2dist, dist2deg, dist2sec
from ra_dec import radec2dir
from asteroid_dataframe import calc_ast_data, spline_ast_vec_df, calc_ast_dir, spline_ast_vec_dir
from ztf_data import load_ztf_det_all
from ztf_ast import load_ztf_nearest_ast

# Typing
from typing import Optional, Dict

# ********************************************************************************************************************* 
def ztf_elt_hash(elts: pd.DataFrame, ztf: pd.DataFrame=None, thresh_deg: float = 1.0, near_ast: bool = False):
    """
    Load or generate a ZTF batch with all ZTF observations within a threshold of the given elements
    INPUTS:
        elts:       Dataframe including element_id; 6 orbital elements
        ztf:        Dataframe of ZTF observations used
        thresh_deg: Threshold in degrees; only observations this close to elements are returned
        near_ast:   Whether to include data on the neareast asteroid
    OUTPUTS:
        hash_id:    Unique ID for these inputs
    """
    # Columns of the Dataframe to hash
    cols_to_hash_elt = ['a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch']
    # Tuple of int64; one per orbital element candidate
    hash_elts = tuple((pd.util.hash_pandas_object(elts[cols_to_hash_elt])).values)
    # Hash of the ZTF if it was input
    if ztf is not None:
        cols_to_hash_ztf = ['mjd',]
        hash_ztf = tuple((pd.util.hash_pandas_object(ztf[cols_to_hash_ztf])).values)
        near_ast = True
    else:
        hash_ztf = tuple()
    # Combine the element hash tuple with the threshold and near_ast flags
    thresh_int = int(thresh_deg*2**48)
    hash_id = abs(hash(hash_elts + hash_ztf + (thresh_int, near_ast,)))

    return hash_id

# ********************************************************************************************************************* 
def make_ztf_near_elt(ztf: pd.DataFrame, 
                      df_dir: pd.DataFrame, 
                      df_ast: pd.DataFrame,
                      thresh_deg: float, 
                      progbar: bool=False) \
                      -> Dict[np.int32, pd.DataFrame]:
    """
    Assemble Python dict of DataFrames with ZTF observations near a batch of orbital elements.
    INPUTS:
        ztf:        DataFrame of candidate ZTF observations
        df_dir:     DataFrame of splined directions for elements at the unique observation times
        thresh_deg: Threshold for a close observation in degrees
        progbar:    Whether to display a progress bar
    """
    # Get unique element IDs
    element_ids = df_dir.element_id.values
    element_ids_unq = np.unique(element_ids)

    # Column collections used on ztf
    cols_catalog = ['ObjectID', 'CandidateID', 'TimeStampID', 'mjd']
    cols_radec = ['ra', 'dec']
    cols_dir = ['ux', 'uy', 'uz']
    cols_elt_qv = ['qx', 'qy', 'qz', 'vx', 'vy', 'vz']
    cols_elt_dir = ['elt_ux', 'elt_uy', 'elt_uz']
    # Add nearest asteroid data to columns if available
    cols_nearest_ast_all = ['nearest_ast_num', 'nearest_ast_dist', 
                            'ast_ra', 'ast_dec', 'ast_ux', 'ast_uy', 'ast_uz']
    cols_nearest_ast = [col for col in cols_nearest_ast_all if col in ztf.columns]
    # All the columns pertaining to ZTF and nearest_ast
    # The calculated information with the position, velocity and direction is separate
    cols_ztf = cols_catalog + cols_radec + cols_dir + cols_nearest_ast
    # Flag indicating whether we have nearest asteroid data
    is_near_ast: bool = 'nearest_ast_num' in ztf.columns

    # Extract TimeStampID as (M,) array; name it row_num to emphasize that we use it to index into u_elt
    row_num = ztf.TimeStampID.values

    # Extract directions of the ZTF observations as an Mx3 array
    u_ztf = ztf[cols_dir].values

    # Threshold as Cartesian distance
    thresh_s = deg2dist(thresh_deg)    
    # Threshold value of z = 1- s^2 / 2
    thresh_z = 1.0 - thresh_s**2 / 2.0

    # Theshold for hits is 2.0 arc seconds (completely separate from threshold above, for inclusion)
    thresh_hit_sec = 2.0

    # Dictionary of DataFrames that are close
    ztf_tbl = dict()

    # Iterate over distinct element IDs
    iterates = tqdm(element_ids_unq) if progbar else element_ids_unq
    for element_id in iterates:
        # Projected directions with this element id
        mask_elt = (df_dir.element_id == element_id)
        # Position and velocity of this candidate element
        qv_elt = df_ast.loc[mask_elt, cols_elt_qv].values
        # Directions of this candidate element at the unique time stamps
        u_elt = df_dir.loc[mask_elt, cols_dir].values
        r_elt = df_dir.loc[mask_elt, 'delta'].values
        # Difference bewteen ztf direction and this element's direction
        s_all = np.linalg.norm(u_ztf - u_elt[row_num], axis=1)
        # Which rows are within the threshold distance?
        mask_close = (s_all < thresh_s)
        # Distances for the close rows only
        s = s_all[mask_close]
        # Row numbers corresponding to close observations
        row_num_close = row_num[mask_close]
        # Copy this slice
        ztf_i = ztf.loc[mask_close, cols_ztf].copy()

        # Insert columns with the ztf_id and element_id
        ztf_i.insert(loc=0, column='ztf_id', value=ztf_i.index.values)
        ztf_i.insert(loc=1, column='element_id', value=element_id)

        # Insert columns with the predicted position and velocity
        for col in cols_elt_qv:
            ztf_i.insert(loc=ztf_i.columns.size, column=col, value=0.0)
        ztf_i[cols_elt_qv] = qv_elt[row_num_close]

        # Insert columns with the predicted directions of the elements
        for col in cols_elt_dir:
            ztf_i.insert(loc=ztf_i.columns.size, column=col, value=0.0)
        ztf_i[cols_elt_dir] = u_elt[row_num_close]
        ztf_i['elt_r'] = r_elt[row_num_close]

        # Save the distance in a column named s
        ztf_i['s'] = s
        # Convert distance between element and observation to degrees
        # dist_deg = dist2deg(dist_i[mask_close])
        s_sec = dist2sec(s)
        # Save distance in arc seconds
        ztf_i['s_sec'] = s_sec
        # Save the quantity z = 1 - s^2 / 2
        ztf_i['z'] = 1.0 - s**2 / 2.0
        # Save the quantity v = (1 - z) / (1 - thresh_z); for random observations, v will uniform on [0, 1]
        # More numerically accurate calculation simplifies 1-z = s^s so (1-z)/(1-thresh_z) = (s/thresh_s)^2
        # ztf_i['v'] = (1.0 - ztf_i.z) / (1.0 - thresh_z)
        ztf_i['v'] = (s/thresh_s)**2
        # Compute score function with given resolution
        # score_arg = -0.5 * (dist_deg / R_deg)**2
        # ztf_i['score'] = np.exp(score_arg)
        # Flag for hits; an entry is a hit if these elements are within threshold of the ZTF observation
        ztf_i['is_hit'] = (ztf_i.s_sec < thresh_hit_sec)
        # Flag for matches; a match is when the element_id matches the nearest asteroid
        # ztf_i['is_match'] = False
        if is_near_ast:
            ztf_i['is_match'] = (ztf_i.nearest_ast_num == ztf_i.element_id)
        # Save ztf_i to the table (dict) of elements
        ztf_tbl[element_id] = ztf_i

    # Combine rows into one big dataframe
    ztf_batch = pd.concat(ztf_tbl.values())

    # Reset the index so it does not have duplicate indices
    ztf_batch.reset_index(inplace=True, drop=True)

    return ztf_batch

# ********************************************************************************************************************* 
def make_ztf_batch(elts: pd.DataFrame, ztf: pd.DataFrame=None, thresh_deg: float = 1.0, near_ast: bool=False):
    """
    Generate a ZTF batch with all ZTF observations within a threshold of the given elements
    INPUTS:
        elts:       Dataframe including element_id; 6 orbital elements
        ztf:        DataFrame of ZTF observations, possibly with asteroid data.
                    Default of None means to use all of the saved ZTF observations
        thresh_deg: Threshold in degrees; only observations this close to elements are returned
        near_ast:   Whether to include data on the neareast asteroid
    """
    # Load all ZTF observations; include nearest asteroid data if requested
    if ztf is None and near_ast:
        # ZTF detections w/ nearest ast; need to regenerate unique times
        ztf = load_ztf_nearest_ast() 
        mjd = np.unique(ztf.mjd)
    elif ztf is None and not near_ast:
        # just the ZTF detections w/o nearest asteroid data
        ztf, mjd = load_ztf_det_all()
    else:
        # Manual treatment of the dates and time stamps if we used a custom ZTF frame
        # Generate the TimeStampID column.  This is integer key counting the unique times of observations.
        mjd, time_stamp_id = np.unique(ztf.mjd.values, return_inverse=True)
        # Convert time_stamp_id to 32 bit integer
        time_stamp_id = time_stamp_id.astype(np.int32)
        # Update TimeStampID on the DataFrame
        ztf['TimeStampID'] = time_stamp_id

    # element_id is the unique identifier for each orbital element considered
    element_id = elts.element_id.values

    # Compute mjd0 and mjd1 from mjd_unq
    mjd0 = np.floor(np.min(mjd))
    mjd1 = np.ceil(np.max(mjd))

    # Calculate positions in this date range, sampled daily
    df_ast_daily, df_earth_daily, df_sun_daily = calc_ast_data(elts=elts, mjd0=mjd0, mjd1=mjd1, element_id=element_id)

    # Spline positions at ztf times
    df_ast, df_earth, df_sun = spline_ast_vec_df(df_ast=df_ast_daily, df_earth=df_earth_daily, df_sun=df_sun_daily, 
                                                 mjd=mjd, include_elts=False)
    # Direction from palomar
    df_dir = calc_ast_dir(df_ast=df_ast, df_earth=df_earth, site_name='palomar')

    # Calculate subset of ZTF data within threshold of this batch
    progbar = True
    ztf_batch = make_ztf_near_elt(ztf=ztf, df_dir=df_dir, df_ast=df_ast,
                                  thresh_deg=thresh_deg, progbar=progbar)

    return ztf_batch

# ********************************************************************************************************************* 
def load_ztf_batch(elts: pd.DataFrame,
                   ztf: pd.DataFrame=None,
                   thresh_deg: float = 1.0, 
                   near_ast: bool = False,
                   regenerate: bool = False):
    """
    Load or generate a ZTF batch with all ZTF observations within a threshold of the given elements
    INPUTS:
        elts:       Dataframe including element_id; 6 orbital elements
        thresh_deg: Threshold in degrees; only observations this close to elements are returned
        near_ast:   Whether to include data on the neareast asteroid
    OUTPUTS:
        ztf_elt:    Dataframe including ZTF observation (mjd, ra, dec, ...)
                    elt_ux, elt_uy, elt_uz, 
    """
    # Get hash of arguments
    hash_id = ztf_elt_hash(elts=elts, ztf=ztf, thresh_deg=thresh_deg, near_ast=near_ast)

    # Name of file
    file_path = f'../data/ztf_elt/ztf_elt_{hash_id}.h5'

    # Try to load file if available
    try:
        if regenerate:
            raise FileNotFoundError
        ztf_elt = pd.read_hdf(file_path)
    # Generate it on the fly if it's not available
    except FileNotFoundError:
        ztf_elt = make_ztf_batch(elts=elts, ztf=ztf, thresh_deg=thresh_deg, near_ast=near_ast)
        ztf_elt.to_hdf(file_path, key='ztf_elt', mode='w')
    
    return ztf_elt

# ********************************************************************************************************************* 
def ztf_score_by_elt(ztf_elt: pd.DataFrame):
    """Calculate DataFrame with summary score for elements"""
    # Score by element; use log_v as a proxy.  This has E[log(v)] = 0, Var[log(v)] = 1 b/c V ~ Unif[0, 1]
    score_func = lambda x: -1.0 - np.log(x)
    # ztf_elt['score'] = 1.0 - np.log(ztf_elt.v)
    # score_by_elt = ztf_elt['v'].apply(np.log).groupby(ztf_elt.element_id).agg(['sum', 'count'])
    score_by_elt = ztf_elt['v'].apply(score_func).groupby(ztf_elt.element_id).agg(['sum', 'count'])
    score_by_elt.rename(columns={'sum': 'score_sum', 'count': 'num_obs'}, inplace=True)
    score_by_elt['t_score'] = score_by_elt['score_sum'] / np.sqrt(score_by_elt['num_obs'])    
    
    return score_by_elt

# ********************************************************************************************************************* 
def ztf_elt_summary(ztf_elt: pd.DataFrame, score_by_elt, elt_name: str):
    """Report summary attributes of a ztf_elt dataframe"""
    # Score by element; use log_v as a proxy.  This has E[log(v)] = 0, Var[log(v)] = 1 b/c V ~ Unif[0, 1]
    score_func = lambda x: -1.0 - np.log(x)
    score = score_func(ztf_elt.v)

    # Calculate summary statistics
    num_obs = ztf_elt.shape[0]
    batch_size = np.unique(ztf_elt.element_id).size
    obs_per_batch = num_obs / batch_size
    num_hits = np.sum(ztf_elt.is_hit)
    hits_per_batch = num_hits / batch_size
    hit_rate = np.mean(ztf_elt.is_hit)

    # Summary statistics after grouping by element_id
    mean_score_sum = np.mean(np.mean(score_by_elt.score_sum))
    mean_t_score= np.mean(np.mean(score_by_elt.t_score))

    # Report results
    print(f'ZTF Element Dataframe {elt_name}:')
    print(f'                  Total     (Per Batch)')
    print(f'Observations   : {num_obs:8d}   ({obs_per_batch:9.0f})')
    # print(f'Hits           : {num_hits:8d}   ({hits_per_batch:9.2f})')
    
    print(f'\nSummarize score = sum(-1.0 - log(v)) by batch.  (Mean=0, Variance=num_obs)')
    print(f'Mean score     :  {mean_score_sum:9.2f}')
    print(f'Sqrt(batch_obs):  {np.sqrt(obs_per_batch):9.2f}')
    print(f'Mean t_score   :  {mean_t_score:9.2f}')
