"""
Harvard IACS Masters Thesis
random_element_batch.py: 
Compute the ztf_elt array for batches of random orbital elements

Michael S. Emanuel
Fri Apr 10 10:24 2020
"""

# Core
import numpy as np
import pandas as pd

# Utility
import argparse
import datetime
from tqdm.auto import tqdm

# MSE imports
from asteroid_element import load_ast_elt
from ztf_ast import load_ztf_nearest_ast
from candidate_element import asteroid_elts, random_elts
from ztf_element import load_ztf_batch, ztf_score_by_elt
from astro_utils import deg2dist

# Typing
from typing import List, Tuple, Dict, Optional

# ********************************************************************************************************************* 
# Constants
# dtype = tf.float32
dtype_np = np.float32
space_dims = 3

# ********************************************************************************************************************* 
def make_ztf_ast(known_ast: bool):
    """Build the ZTF DataFrame including nearest asteroid information"""
    # Load ztf nearest asteroid data
    ztf_ast = load_ztf_nearest_ast()
 
    # Filter ztf_ast to only include hits
    hit_thresh_sec = 2.0
    hit_thresh_s = deg2dist(hit_thresh_sec / 3600.0)
    ztf_ast['is_hit'] = ztf_ast.nearest_ast_dist < hit_thresh_s
    # mask = np.logical_xor(ztf_ast.is_hit.values, not known_ast)
    mask = ztf_ast.is_hit if known_ast else ~ztf_ast.is_hit
    ztf_ast = ztf_ast[mask]

    return ztf_ast

# ********************************************************************************************************************* 
def best_elt_file_path(seed: int, known_ast: bool, batch_size_init: int, batch_size: int, thresh_deg: float):
    """File name for best elements"""
    known_token = 'hit' if known_ast else 'miss'
    thresh_sec = int(thresh_deg * 3600)
    file_path = f'../data/ztf_elt/random_elts_{known_token}_seed_{seed:03d}_size_{batch_size}_of_{batch_size_init}_thresh_{thresh_sec}.h5'
    return file_path

# ********************************************************************************************************************* 
def select_best_elts(elts_init: pd.DataFrame, ztf_elt: pd.DataFrame, batch_size: int, element_id_start: int):
    """
    Extract the best elements from a batch of candidate random elements
    """
    # Score by element on the original batch
    score_by_elt = ztf_score_by_elt(ztf_elt)
    
    # Sort score_by_elt by score descending; extract best scores
    score_by_elt.sort_values(by='score_sum', ascending=False, inplace=True)
    element_id_best = score_by_elt.index[0:batch_size].values
    num_obs = score_by_elt.num_obs[0:batch_size].values
    score = score_by_elt.score_sum[0:batch_size].values
    t_score = score_by_elt.t_score[0:batch_size].values

    # Build new elements: copy the best elements in sorted order
    elts = elts_init.loc[element_id_best].copy()
    # Reindex elts from 0; generate new element_id starting from designated element_id_start
    elts.reset_index(drop=True, inplace=True)
    element_id = np.arange(element_id_start, element_id_start+batch_size, dtype=np.int32)
    # Add additional columns to these random elements
    elts['element_id'] = element_id
    elts['num_obs'] = num_obs
    elts['score'] = score
    elts['t_score'] = t_score
    
    return elts

# ********************************************************************************************************************* 
def calc_best_random_elts(random_seed: int, known_ast: bool, 
                          batch_size_init: int=1024, batch_size: int=64, thresh_deg: float=2.0,
                          ztf_ast: Optional[pd.DataFrame]=None):
    """
    Load the best random elements if available, or recompute on the fly
    INPUTS:
        random_seed:        Random seed used to generate random orbital elements
        known_ast:          True: include only hits (<2 arc sec) vs. a known ast; False: only non-hits
        batch_size_init:    Batch size for large initial batch, e.g. 1024
        batch_size:         Batch size for small final batch, e.g. 64 best out of 1024
        thresh_deg:         Threshold for observations to be included, e.g. 2.0 degrees
    """
    # Path where this file will be saved
    file_path = best_elt_file_path(seed=random_seed, known_ast=known_ast, batch_size_init=batch_size_init, 
                                    batch_size=batch_size, thresh_deg=thresh_deg)

    # Build the ztf_ast DataFrame if it was not passed
    if ztf_ast is None:
        known_type = 'HITS' if known_ast else 'NO HITS'
        print(f'Building ZTF DataFrame: {known_type} vs. known asteroids; thresh = {thresh_deg}\n')
        ztf_ast = make_ztf_ast(known_ast=known_ast)

    # Build the initial array of random elements with element_id starting at 0
    element_id_start_init = 0
    # Additional inputs for load_ztf_batch
    near_ast = False
    regenerate = False

    # Inital large batch of random elements
    elts_init = random_elts(element_id_start=element_id_start_init, size=batch_size_init,
                            random_seed=random_seed, dtype=dtype_np)

    # Load or generate ZTF batch for the big random elements
    ztf_elt_init = load_ztf_batch(elts=elts_init, ztf=ztf_ast, thresh_deg=thresh_deg, 
                                    near_ast=near_ast, regenerate=regenerate)

    # Extract the best elements from the candidates
    element_id_start = random_seed * batch_size
    elts = select_best_elts(elts_init=elts_init, ztf_elt=ztf_elt_init, 
                            batch_size=batch_size, element_id_start=element_id_start)

    # Save the best elements        
    file_path = best_elt_file_path(seed=random_seed, known_ast=known_ast, batch_size_init=batch_size_init, 
                                    batch_size=batch_size, thresh_deg=thresh_deg)
    elts.to_hdf(file_path, key='elts', mode='w')
    print(f'Saved random elements in {file_path}.')

    # Regenerate ztf_elt for just the batch of good orbital elements
    load_ztf_batch(elts=elts, ztf=ztf_ast, thresh_deg=thresh_deg, near_ast=near_ast, regenerate=regenerate)        

    return elts

# ********************************************************************************************************************* 
def load_best_random_elts(random_seed: int, known_ast: bool, 
                          batch_size_init: int=1024, batch_size: int=64, thresh_deg: float=2.0,
                          ztf_ast: Optional[pd.DataFrame]=None):
    """
    Load the best random elements if available, or recompute on the fly
    INPUTS:
        random_seed:        Random seed used to generate random orbital elements
        known_ast:          True: include only hits (<2 arc sec) vs. a known ast; False: only non-hits
        batch_size_init:    Batch size for large initial batch, e.g. 1024
        batch_size:         Batch size for small final batch, e.g. 64 best out of 1024
        thresh_deg:         Threshold for observations to be included, e.g. 2.0 degrees
        ztf_ast:            Optional precomputed DataFrame of ZTF observations 
                            filtered for either known or unkown asteroids
    """
    # Path where this file is expected to be
    file_path = best_elt_file_path(seed=random_seed, known_ast=known_ast, batch_size_init=batch_size_init, 
                                    batch_size=batch_size, thresh_deg=thresh_deg)

    # Load the file if available
    try:
        elts = pd.read_hdf(file_path, key='elts')
        print(f'Loaded random elements in {file_path}.')
    # Otherwise regenerate it
    except FileNotFoundError:
        elts = calc_best_random_elts(random_seed=random_seed, known_ast=known_ast,
                                     batch_size_init=batch_size_init, thresh_deg=thresh_deg,
                                     ztf_ast=None)
    return elts

# ********************************************************************************************************************* 
def main(seed0: int, seed1: int, stride: int,
         batch_size_init: int, batch_size: int, 
         known_ast: bool, thresh_deg: float):
    """
    Main program body
    INPUTS:
        seed0:              The first random seed, inclusive
        seed1:              The last random seed, exclusive
        stride:             Stride for stepping through the random seeds
        batch_size_init:    Batch size for large initial batch, e.g. 1024
        batch_size:         Batch size for small final batch, e.g. 64 best out of 1024
        known_ast:          True: include only hits (<2 arc sec) vs. a known ast; False: only non-hits
        thresh_deg:         Threshold for observations to be included, e.g. 2.0 degrees
    """
    
    # Build the ztf_ast DataFrame
    known_type = 'HITS' if known_ast else 'NO HITS'
    print(f'Building ZTF DataFrame: {known_type} vs. known asteroids; thresh = {thresh_deg}\n')
    ztf_ast = make_ztf_ast(known_ast=known_ast)

    # Iterate over random seeds in the specified range
    random_seeds = list(range(seed0, seed1, stride))
    for i, random_seed in enumerate(random_seeds):
        # Status
        time_str = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f'\nItem {i:3}, random_seed {random_seed}. {time_str}')
        
        # Delegate to load_best_random_elts using precomputed ztf_ast
        load_best_random_elts(random_seed=random_seed, known_ast=known_ast,
                              batch_size_init=batch_size_init, batch_size=batch_size, thresh_deg=thresh_deg,
                              ztf_ast=ztf_ast)

# ********************************************************************************************************************* 
if __name__ == '__main__':    
    # Process command line arguments
    parser = argparse.ArgumentParser(description='Build ZTF DataFrames for best random elements.')
    parser.add_argument('-seed0', nargs='?', metavar='s0', type=int, default=0,
                        help='First random seed to process, inclusive.')
    parser.add_argument('-seed1', nargs='?', metavar='s1', type=int, default=256,
                        help='Last random seed to process, exclusive.')
    parser.add_argument('-stride', nargs='?', metavar='str', type=int, default=4,
                        help='Stride for stepping through seeeds.')
    parser.add_argument('-batch_size_init', nargs='?', metavar='bsi', type=int, default=1024,
                        help='Large batch size for initial pass.')
    parser.add_argument('-batch_size', nargs='?', metavar='bs', type=int, default=64,
                        help='Small batch size; final result, after filtering.')
    known_ast_help_msg= 'when true, match ZTF observations that match a known asteroid to 2.0 arc seconds.\n' \
                        'when false, match ZTF observations that do not match a known asteroid.'
    parser.add_argument('-known_ast', default=False, action='store_true',
                        help=known_ast_help_msg)
    parser.add_argument('-thresh_deg', nargs='?', metavar='thresh_deg', type=float, default=2.0,
                        help='Threshold in degrees; max distance from elements direction to observation.')

    # Parse arguments
    args = parser.parse_args()

    # Alias arguments
    seed0 = args.seed0
    seed1 = args.seed1
    stride = args.stride
    batch_size_init = args.batch_size_init
    batch_size = args.batch_size
    known_ast = args.known_ast    
    thresh_deg = args.thresh_deg

    # Report inputs
    print(f'\nInputs to random elements batch:')
    print(f'seed0  = {seed0}. First random seed for random elements (inclusive).')
    print(f'seed1  = {seed1}. Last random seed for random elements (exclusive).')
    print(f'stride = {stride:5}.     Stride for stepping through random seeds.')
    print(f'batch_size_init = {batch_size_init}. Size of initial batch of random elements.')
    print(f'batch_size = {batch_size}. Size of final batch of elements; top scoring.')
    print(f'known_ast = {known_ast}. Match ZTF observations < 2.0 arc sec to known asteroids (true) or >= 2.0 sec (false).')
    print(f'thresh_deg = {thresh_deg}. Maximum distance in sky between predicted direction of elements and ZTF observation.\n')

    # Main body
    main(seed0=seed0, seed1=seed1, stride=stride,
         batch_size_init=batch_size_init, batch_size=batch_size,
         known_ast=known_ast, thresh_deg=thresh_deg)