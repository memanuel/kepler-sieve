"""
Harvard IACS Masters Thesis
asteroid_search.py: Search for orbital elements of asteroids given observational data.

Michael S. Emanuel
Thu Oct 17 15:24:10 2019
"""

# Core
import numpy as np
import pandas as pd

# Tensorflow / ML
import tensorflow as tf

# Utility
import time
import datetime
import argparse
import os
import warnings

# ********************************************************************************************************************* 
# Tensorflow config
from tf_utils import gpu_grow_memory, get_gpu_device, tf_quiet
# Need to do this before importing local modules with TensorFlow variables / constants stored on GPU

# Run TF quietly
tf_quiet()

# Silence irrelevant warnings
warnings.filterwarnings('ignore')

# Configure TF GPU growth
# gpu_grow_memory(verbose=True)

# ********************************************************************************************************************* 
# MSE Imports
from asteroid_element import load_ast_elt
from candidate_element import asteroid_elts, perturb_elts, random_elts, elts_add_mixture_params, elts_add_H
from random_elements import load_best_random_elts, make_ztf_ast
from ztf_ast import load_ztf_nearest_ast, calc_hit_freq
from ztf_element import load_ztf_batch, make_ztf_batch, ztf_score_by_elt, ztf_elt_summary
from asteroid_model import AsteroidPosition, AsteroidDirection, make_model_ast_pos
from asteroid_search_layers import CandidateElements, MixtureParameters, TrajectoryScore
from asteroid_search_model import AsteroidSearchModel, save_dir
from asteroid_search_report import traj_diff
from nearest_asteroid import nearest_ast_elt_cart, nearest_ast_elt_cov, elt_q_norm
from element_eda import score_by_elt
from asteroid_dataframe import calc_ast_data, spline_ast_vec_df
from astro_utils import deg2dist, dist2deg, dist2sec
from utils import print_header
from tf_utils import tf_quiet, gpu_grow_memory, get_gpu_device

# Typing
from typing import Dict

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
# Constants
dtype = tf.float32
dtype_np = np.float32
space_dims = 3

# Orbital elements of known asteroids
ast_elt = load_ast_elt()

# ********************************************************************************************************************* 
def file_name_model(seed: int, known_ast: bool, batch_size_init: int, batch_size: int, thresh_deg: float):
    """File name for saved model"""
    known_token = 'known' if known_ast else 'unknown'
    thresh_sec = int(thresh_deg * 3600)
    file_name = f'candidate_elt_{known_token}_seed_{seed:03d}_size_{batch_size}_of_{batch_size_init}_thresh_{thresh_sec}.h5'
    return file_name

# ********************************************************************************************************************* 
def file_path_fitted_elts(known_ast: bool):
    """
    File name for global fitted elements file
    INPUTS:
        known_ast: Whether elements are fitted against known or unknown asteroid observations.    
    """
    known_type = 'known' if known_ast else 'unknown'
    file_name = f'fitted_elts_{known_type}.h5'
    file_path = os.path.join(save_dir, file_name)
    return file_path

# ********************************************************************************************************************* 
def file_path_ztf_hits(known_ast: bool):
    """
    File name for global ZTF hits file
    INPUTS:
        known_ast: Whether hits are against known or unknown asteroid observations.    
    """
    known_type = 'known' if known_ast else 'unknown'
    file_name = f'ztf_hits_{known_type}.h5'
    file_path = os.path.join(save_dir, file_name)
    return file_path

# ********************************************************************************************************************* 
def load_fitted_elts(known_ast: bool, display:bool=True, min_hits: int=0):
    """
    Load the fitted orbital elements
    INPUTS:
        known_ast: Whether elements are fitted against known or unknown asteroid observations.
        display:   When true, only show the columns formatted for user friendly display mode
        min_hits:  Minimum number of hits
    """

    # Location of the file
    file_path = file_path_fitted_elts(known_ast=known_ast)

    # Load the fitted elements if available, or return an empty DataFrame if not
    fitted_elts: pd.DataFrame
    try:
        fitted_elts = pd.read_hdf(file_path)
    except:
        fitted_elts = pd.DataFrame()

    # If display was requested, return only columns for human consumers
    if display and fitted_elts.shape[0] > 0:
        cols_display = [
            'element_id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch',
            'num_hits', 'R_sec', 'thresh_sec', 'log_like', 'hits', 'num_rows_close', 'timestamp']
        fitted_elts = fitted_elts[cols_display]
        # Deduplicate
        fitted_elts = fitted_elts.loc[~fitted_elts.index.duplicated(keep='last')]
        # Sort to show the best fitted elements first
        fitted_elts.sort_values(by=['hits', 'log_like'], ascending=False, inplace=True)
    
    # If min_hits was passed, limit output to only elements with the required number of hits
    if min_hits and fitted_elts.shape[0] > 0:
        # Is this element_id successfully fit with?
        is_good_fit = (min_hits <= np.ceil(fitted_elts.hits))
        # Apply this filter
        fitted_elts = fitted_elts[is_good_fit]

    return fitted_elts

# ********************************************************************************************************************* 
def append_fitted_elt(fitted_elt: pd.DataFrame, known_ast: bool):
    """
    Append the current batch of fitted elements to the global store
    INPUTS:
        fitted_elt:  The new batch of fitted elements to be added
        known_ast:    Whether elements are fitted against known or unknown asteroid observations.    
    """

    # Location of the file
    file_path = file_path_fitted_elts(known_ast=known_ast)
    
    # Load the file
    fitted_elts = load_fitted_elts(known_ast=known_ast, display=False)
    
    # Append the new fitted elements
    fitted_elts = fitted_elts.append(fitted_elt)
    
    # Save the modified file
    fitted_elts.to_hdf(file_path, key='fitted_elts')
    
    # Save the display version
    fitted_elts_disp = load_fitted_elts(known_ast=known_ast, display=True)
    file_path_csv = file_path.replace('.h5', '.csv')
    fitted_elts_disp.to_csv(file_path_csv)

# ********************************************************************************************************************* 
def load_ztf_hits(known_ast: bool, display: bool=True, min_hits: int=0):
    """
    Load ZTF hits
    INPUTS:
        known_ast:  Whether hits are for known or unknown asteroid observations
        display:    Whether to return output in display mode or all columns
        min_hits:   Minimum number of hits
    """

    # Location of the file
    file_path = file_path_ztf_hits(known_ast=known_ast)
    
    # Load the fitted elements if available, or return an empty DataFrame if not
    ztf_hits: pd.DataFrame
    try:
        ztf_hits = pd.read_hdf(file_path)
    except:
        ztf_hits = pd.DataFrame()

    # If display was requested, return only columns for human consumers
    if display and ztf_hits.shape[0] > 0:
        cols_display = [
            'element_id', 'ztf_id', 'ObjectID', 'CandidateID', 
            'mjd', 'ra', 'dec', 'mag_app', 'ux', 'uy', 'uz',
            'elt_ux', 'elt_uy', 'elt_uz', 's_sec',
            'timestamp']
        ztf_hits = ztf_hits[cols_display]
        # Deduplicate
        ztf_hits = ztf_hits.loc[~ztf_hits.index.duplicated(keep='last')]
        # Don't need to sort; index is already sorted (also it fails b/c of horrible Pandas treatment of indexes)

    # If good was specified, limit to only those elements that were well fit
    if min_hits:
        fitted_elts_good = load_fitted_elts(known_ast=known_ast, display=False, min_hits=min_hits)
        if fitted_elts_good.shape[0] > 0:
            element_id_good = np.unique(fitted_elts_good.element_id)
            # Filter ztf_hits on the first level of the multi-index
            ztf_hits = ztf_hits.loc[(element_id_good),:]
        else:
            ztf_hits = pd.DataFrame()
    
    return ztf_hits

# ********************************************************************************************************************* 
def append_ztf_hit(ztf_hit: pd.DataFrame, known_ast: bool):
    """
    Append the current batch of ztf hits to the global store
    INPUTS:
        ztf_hit:   The new batch of ZTF hits
        known_ast: Whether hits are against known or unknown asteroid observations.    
    """

    # Location of the file
    file_path = file_path_ztf_hits(known_ast=known_ast)
    
    # Load the file
    ztf_hits = load_ztf_hits(known_ast=known_ast, display=False)
    
    # Append the new fitted elements
    ztf_hits = ztf_hits.append(ztf_hit)
    
    # Save the modified file
    ztf_hits.to_hdf(file_path, key='ztf_hits')

    # Save the display version
    ztf_hits_disp = load_ztf_hits(known_ast=known_ast, display=True)
    file_path_csv = file_path.replace('.h5', '.csv')
    ztf_hits_disp.to_csv(file_path_csv)

# ********************************************************************************************************************* 
def remove_prior_hits(ztf_elt, elts, min_hits: int=10):
    """
    Filter out rows in ztf_elt that correspond to hits in ztf_hits for prior elements to those in elts.
    INPUTS:
        ztf_elt:    ZTF observations near the candidate elements
        elts:       The batch of candidate elements
        min_hits:   Minimum number of hits an element has 
    OUTPUTS:
        ztf_elt:    Filtered version of ztf_elt with prior hits removed
    """
    # First element_id in the candidate elements
    element_id_min: np.int32 = np.min(elts.element_id)

    # Load current ztf hits (good only)
    ztf_hits: pd.DataFrame = load_ztf_hits(known_ast=known_ast, display=False, min_hits=min_hits)
    
    # Check for edge case that either fitted_elts or ztf_hits is empty
    if ztf_hits.shape[0] == 0:
        return ztf_elt

    # Which hits are for prior elements to the current one?
    is_prior_elt: pd.Series = (ztf_hits.element_id < element_id_min)

    # Unique ZTF IDs appearing in hits for prior elements
    ztf_id_unq: np.ndarray = np.unique(ztf_hits[is_prior_elt].ztf_id)

    # Is this observation a duplicate of one the existing hits?
    is_dup: pd.Series = ztf_elt.ztf_id.isin(ztf_id_unq)

    # Count the number of duplicates
    dup_count: np.int32 = np.sum(is_dup)
    # Status
    if dup_count > 0:
        print(f'Removing {dup_count} duplicate hits from ztf_elt for elements with at least {min_hits} hits.')

    # Return ztf_elt filtered to remove duplicates
    return ztf_elt[~is_dup]

# ********************************************************************************************************************* 
def train_one_batch(i: int, random_seed: int, ztf_ast: pd.DataFrame):
    """Train one batch of candidate elements"""
    # Arguments to make_ztf_batch
    near_ast = False
    regenerate = False

    # Mixture parameters
    num_hits: int = 10
    R_deg: float = 0.5

    # Observatory for ZTF data is Palomar Mountain
    site_name = 'palomar'

    # Training parameters
    learning_rate = 2.0**-12
    clipnorm = 1.0

    # File name
    file_name = file_name_model(seed=random_seed, known_ast=known_ast, 
                                batch_size_init=batch_size_init, batch_size=batch_size, 
                                thresh_deg=thresh_deg)
    file_path = os.path.join(save_dir, file_name)

    # If this file already exists, quit early
    if os.path.isfile(file_path):
        print(f'Found file for trained model {file_path}, moving on to next one.')
        return

    # Status
    time_str = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'\nItem {i:3}, random_seed {random_seed}. {time_str}')
    
    # Load the best random elements (or recompute if unavailable)
    elts = load_best_random_elts(
        random_seed=random_seed, known_ast=known_ast,
        batch_size_init=batch_size_init, batch_size=batch_size, thresh_deg=thresh_deg,
        ztf_ast=ztf_ast)

    # Load ztf_elt data
    ztf_elt = load_ztf_batch(elts=elts, ztf=ztf_ast, thresh_deg=thresh_deg, 
                             near_ast=near_ast, regenerate=regenerate)

    # Remove rows corresponding to prior hits
    ztf_elt = remove_prior_hits(ztf_elt=ztf_elt, elts=elts)

    # Add mixture parameters to candidate elements
    elts_add_mixture_params(elts=elts, num_hits=num_hits, R_deg=R_deg, thresh_deg=thresh_deg)

    # Add brightness parameter H
    elts_add_H(elts=elts)

    # Build asteroid search model
    model = AsteroidSearchModel(
                    elts=elts, ztf_elt=ztf_elt, 
                    site_name=site_name, thresh_deg=thresh_deg, 
                    learning_rate=learning_rate, clipnorm=clipnorm,
                    name='model', file_name=file_name)

    # Load the trained model if possible
    model.load()

    # Report before training starts
    model.report()

    # If the model has already been trained, quit early
    if model.current_batch > 10000:
        return

    # Train the model    
    model.sieve()

    # Save the trained model
    model.save_state()

    # Get a common timestamp for both entries (elements and hits)
    timestamp = pd.to_datetime('today')

    # Extract the fitted elements
    fitted_elt = model.candidates_df()
    # Reindex fitted_elt by element_id column so data can be merged
    fitted_elt.set_index(keys='element_id', drop=False, inplace=True)
    # Add a timestamp column
    fitted_elt['timestamp'] = timestamp
    # Append fitted elements to global store
    append_fitted_elt(fitted_elt=fitted_elt, known_ast=known_ast)

    # Generate ZTF hits
    print(f'\nGenerating ZTF hits:')
    ztf_hit = model.calc_ztf_hits()
    # Reindex ztf_hit by element_id column so data can be merged
    ztf_hit.set_index(keys=['element_id', 'ztf_id'], drop=False, inplace=True)
    # Add a timestamp column
    ztf_hit['timestamp'] = timestamp
    # Append ZTF hits to global store
    append_ztf_hit(ztf_hit=ztf_hit, known_ast=known_ast)

# ********************************************************************************************************************* 
def main(seed0: int, seed1: int, stride: int, 
         known_ast: bool, R_deg: float, thresh_deg: float, 
         batch_size_init: int, batch_size: int):
    """
    Main program body
    INPUTS:
        seed0:              The first random seed, inclusive
        seed1:              The last random seed, exclusive
        stride:             Stride for stepping through the random seeds
        known_ast:          True: include only hits (<2 arc sec) vs. a known ast; False: only non-hits
        R_deg:              Resolution in degrees, e.g. 0.5 degrees
        thresh_deg:         Threshold for observations to be included, e.g. 2.0 degrees
        batch_size_init:    Batch size for large initial batch, e.g. 1024
        batch_size:         Batch size for small final batch, e.g. 64 best out of 1024
    """
    
    # Build the ztf_ast DataFrame
    known_type = 'known' if known_ast else 'unknown'
    print(f'Building ZTF DataFrame: {known_type} asteroids; thresh = {thresh_deg}\n')
    ztf_ast = make_ztf_ast(known_ast=known_ast)
    print(f'Loaded ztf_ast with {ztf_ast.shape[0]} rows.')

    # Iterate over random seeds in the specified range
    random_seeds = list(range(seed0, seed1, stride))
    for i, random_seed in enumerate(random_seeds):
        train_one_batch(i=i, random_seed=random_seed, ztf_ast=ztf_ast)

# ********************************************************************************************************************* 
if __name__ == '__main__':    
    # Process command line arguments

    parser = argparse.ArgumentParser(description='Search for asteroids.')
    parser.add_argument('-seed0', nargs='?', metavar='s0', type=int, default=-1,
                        help='First random seed to process, inclusive.')
    parser.add_argument('-seed1', nargs='?', metavar='s1', type=int, default=1024,
                        help='Last random seed to process, exclusive.')
    parser.add_argument('-stride', nargs='?', metavar='str', type=int, default=4,
                        help='Stride for stepping through seeeds.')
    known_ast_help_msg= 'when true, match ZTF observations that match a known asteroid to 2.0 arc seconds.\n' \
                        'when false, match ZTF observations that do not match a known asteroid.'
    parser.add_argument('-known_ast', default=False, action='store_true',
                        help=known_ast_help_msg)
    # parser.add_argument('-gpu_num', nargs='?', metavar='gpu_num', type=int, default=0,
    #                     help='the GPU to use; defaults to 0.')
    parser.add_argument('-batch_size_init', nargs='?', metavar='bsi', type=int, default=1024,
                        help='Large batch size for initial pass.')
    parser.add_argument('-batch_size', nargs='?', metavar='bs', type=int, default=64,
                        help='Small batch size; final result, after filtering.')
    parser.add_argument('-R_deg', nargs='?', metavar='R_deg', type=float, default=0.5,
                        help='Initial resolution in degrees.')
    parser.add_argument('-thresh_deg', nargs='?', metavar='thresh_deg', type=float, default=2.0,
                        help='Threshold in degrees; max distance from elements direction to observation.')

    # Parse arguments
    args = parser.parse_args()

    # Alias arguments
    seed0 = args.seed0
    seed1 = args.seed1
    stride = args.stride
    known_ast = args.known_ast
    # gpu_num = args.gpu_num
    batch_size_init = args.batch_size_init
    batch_size = args.batch_size
    R_deg = args.R_deg
    thresh_deg = args.thresh_deg

    # Get the environment variable with cuda visible devices
    cuda_visible_devices = os.getenv('CUDA_VISIBLE_DEVICES')
    print(f'Environment variables:')
    print(f'CUDA_VISIBLE_DEVICES = {cuda_visible_devices}')

    # Convert the environment variable to a gpu_num; default to 0
    gpu_num: int =  int(cuda_visible_devices[0]) if cuda_visible_devices is not None else 0

    # If the default start of -1 was passed, use the gpu_num
    if seed0 < 0 :
        seed0 = gpu_num

    # Report inputs
    print(f'\nInputs to asteroid search:')
    print(f'seed0           = {seed0:5}.     First random seed for random elements (inclusive).')
    print(f'seed1           = {seed1:5}.     Last random seed for random elements (exclusive).')
    print(f'stride          = {stride:5}.     Stride for stepping through random seeds.')
    print(f'known_ast       = {known_ast:5}.     '
          f'Match ZTF observations < 2.0 arc sec to known asteroids (true) or >= 2.0 sec (false).')
    print(f'gpu_num         = {gpu_num:5}.     GPU on which computations are run.')
    print(f'batch_size_init = {batch_size_init:5}.     Size of initial batch of random elements.')
    print(f'batch_size      = {batch_size:5}.     Size of final batch of elements; top scoring.')
    print(f'R_deg           = {R_deg:5}.     Resolution in degrees.')
    print(f'thresh_deg      = {thresh_deg:5}.     '
          f'Maximum distance in sky between predicted direction of elements and ZTF observation.\n')

    main(seed0=seed0, seed1=seed1, stride=stride, 
            known_ast=known_ast, R_deg=R_deg, thresh_deg=thresh_deg,
            batch_size_init=batch_size_init, batch_size=batch_size)
