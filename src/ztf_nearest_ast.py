"""
Harvard IACS Masters Thesis
ztf_nearest.py
Search ZTF DataFrame against calculated asteroid orbits and implied direction from Palomar.
Find the nearest asteroid to each ZTF observation.

Michael S. Emanuel
11-Mar-2020
"""

# Libary imports
import numpy as np
import pandas as pd
import os
import argparse
from multiprocessing import Pool, cpu_count
from tqdm import auto

# MSE imports
from asteroid_integrate import load_data as load_orbital_elts
from ztf_data import load_ztf_det_all, ztf_load_nearest_ast, ztf_ast_file_path

# Typing
from typing import Optional

# ********************************************************************************************************************* 
def run_batch(n0: int, n1: int, progbar: bool, verbose: bool) -> None:
    """
    Run one batch of finding the nearest asteroid to ZTF data.
    INPUTS:
        n0: First asteroid to process, e.g 0. Inclusive.
        n1: Last asteroid to process, e.g. 1000.  Exclusive.
        progbar: Whether to include a progress bar on console.
        verbose: Whether to print status updates to console.
    RETURNS:
        nearest_ast: DataFrame aligned with ztf; 
                     contains columns nearest_ast_num, nearest_ast_dist
    """
    # File name and path
    file_path = ztf_ast_file_path(n0=n0, n1=n1)
    
    # Load file from disk if it already exists
    if os.path.isfile(file_path):
        if verbose:
            print(f'Found file {file_path}.')
        ztf_ast = pd.read_hdf(file_path)
    # Otherwise generate it
    else:
        # Load ZTF detections as a DataFrame (no data about nearest asteroids yet)
        ztf, mjd_unq = load_ztf_det_all(verbose=False)

        # Set threshold to 180.0 degrees; always find the closest even if it's far away
        thresh_deg = 180.0
        # Set regen flag to false; don't want to rebuild the file
        regen = False

        # Calculate nearest asteroid to ZTF DataFrame (already know it doesn't exist)
        # This will also save it to disk if it's built for first time
        ztf_ast = ztf_nearest_ast(ztf=ztf, n0=n0, n1=n1, thresh_deg=thresh_deg, regen=regen, 
                                  progbar=progbar, verbose=verbose)

# ********************************************************************************************************************* 
def process_all(block_size: int, max_blocks: Optional[int] = None) -> np.ndarray:
    """Process all the asteroids to find the one nearest to each ZTF observation"""
    # Load distinct asteroid numbers
    orb_elt = load_orbital_elts()
    ast_nums = orb_elt.Num.values.astype(np.int32)
    # Extract distinct blocks of data
    ast_blocks = np.unique(ast_nums // block_size)
    # Limit ast_blocks to max_blocks elements if it was specifiec
    if max_blocks is not None:
        ast_blocks = ast_blocks[0:max_blocks]

    # Create a list of arguments to run_batch
    run_batch_args = []
    for k, ast_block in enumerate(ast_blocks):
        n0 = ast_block * block_size
        n1 = n0 + block_size
        # progress bar only on the first batch
        progbar = (k == 0)
        # set verbose argument to false when processing all batches
        verbose = False
        # wrap the arguments as a tuple
        args = (n0, n1, progbar, verbose)
        # append args to full list of args
        run_batch_args.append(args)

    # Set number of processes based on number of available CPUs
    processes = min(np.round(cpu_count()*0.80), ast_blocks.size)

    # Create a thread pool and run the batches in parallel
    with Pool(processes=processes) as thread_pool:
        print(f'Created thread pool with {processes} threads.')
        # nearest_ast_blocks = thread_pool.starmap(func=run_batch, iterable=run_batch_args)
        thread_pool.starmap(func=run_batch, iterable=run_batch_args)

    # Status message
    if max_blocks is not None:
        ast_count = max_blocks * block_size
        msg = f'Processed {max_blocks} blocks of size {block_size} asteroids up to asteroid number {ast_count}.'
    else:
        msg = f'Processed all asteroids in blocks of size {block_size}.'
    print(msg)

    # Return list with the nearest_ast_block frame for each block processed along with the blocks
    return ast_blocks

# ********************************************************************************************************************* 
def nearest_ast_reduction(ast_blocks: np.ndarray, block_size: int, progbar: bool = False) -> np.ndarray:
    """
    Perform a reduction operation on a list of nearest_ast_block DataFrames.
    Each frame in the list is the nearest asteroid out of a block of asteroids.
    The reduction will return the closest asteroid number and its distance globally.
    INPUTS:
        ast_blocks: Array with the block numbers.  
                    Block number bn contains asteroids [bn*s, (bn+1)*s) where s is the block size.
                    e.g. when s = 1000, block 0 is [0, 1000), block 1 is [1000, 2000) etc.
        block_size: Size of each block, e.g. 1000
        progbar:    Whether to display a tqdm progress bar
    """

    # First block
    n0 = ast_blocks[0] * block_size
    n1 = n0 + block_size
    
    # load ztf_ast for first block and initialize it as the master
    ztf = ztf_load_nearest_ast(n0=n0, n1=n1)

    # Full range of asteroids
    ast_num_min = n0
    ast_num_max = ast_blocks[-1] * block_size + block_size

    # Column names of data pertaining to the nearest asteroid
    cols = ['nearest_ast_num', 'nearest_ast_dist', 'ast_ra', 'ast_dec', 'ast_ux', 'ast_uy', 'ast_uz']

    # Iterate through the remaining blocks.  Overwrite nearest asteroid data where a new minimum distance is found
    iterates = ast_blocks[1:]
    if progbar:
        iterates = tqdm(iterates)
    for block in iterates:
        # Indices for this block
        n0 = block * block_size
        n1 = n0 + block_size
        # Load block i
        ztf_i = ztf_load_nearest_ast(n0=n0, n1=n1)
        # Mask for rows that need to be updated
        mask = ztf_i.nearest_ast_dist < ztf.nearest_ast_dist
        # Overwrite data columns for nearest asteroid
        ztf.loc[mask, cols] = ztf_i.loc[mask, cols]

    # Save this file
    file_path = ztf_ast_file_path(n0=ast_num_min, n1=ast_num_max)
    ztf.to_hdf(file_path, key='ztf_ast', mode='w')
  
# ********************************************************************************************************************* 
def main():
    """Main routine for integrating the orbits of known asteroids"""
    # Process command line arguments
    parser = argparse.ArgumentParser(description='Find the nearest asteroid to ZTF data.')
    parser.add_argument('n0', nargs='?', metavar='n0', type=int, default=0,
                        help='the first asteroid number to process')
    parser.add_argument('n_ast', nargs='?', metavar='B', type=int, default=1000,
                        help='the number of asteroids to process in this batch')
    parser.add_argument('--progress', default=False, action='store_true',
                        help='display progress bar')
    parser.add_argument('--all', default=False, action='store_true',
                        help='process all the data')
    parser.add_argument('--test', default=False, action='store_true',
                        help='run in test mode')
    args = parser.parse_args()
    
    # Alias progbar
    progbar: bool = args.progress

    # If run in test mode, run tests without processing any asteroid trajectories
    if args.test:
        # Test something
        # TODO put something here
        
        # Quit early in test mode: don't want to do any calculations
        print()
        exit()

    # If process_all, use mutlprocessing to run everything in parallel
    elif args.all:
        # Inputs for processing all the asteroids
        block_size = 5000
        max_blocks = None
        # Generate ztf_nearest_ast for each block
        ast_blocks = process_all(block_size=block_size, max_blocks=max_blocks)
        # Perform the reduction to get the globally closest asteroid
        nearest_ast_reduction(ast_blocks=ast_blocks, block_size=block_size, progbar=progbar)
        # Status update
        ast_num_min = ast_blocks[0]*block_size
        ast_num_max = (ast_blocks[-1]+1)*block_size
        print(f'Completed reduction of asteroids from {ast_num_min} to {ast_num_max}.')
        exit()

    # Unpack command line arguments to run one batch
    else:
        n0: int = args.n0
        n1: int = n0 + args.n_ast        
        run_batch(n0=n0, n1=n1, progbar=progbar, verbose=True)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
    
