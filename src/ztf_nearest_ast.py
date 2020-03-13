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
from tqdm.auto import tqdm

# MSE imports
from asteroid_integrate import load_data as load_orbital_elts
from ztf_data import load_ztf_det_all, load_ztf_nearest_ast, ztf_nearest_ast, ztf_ast_file_path

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
def calc_ast_block_numbers(block_size: int, n0: int=0, n1: Optional[int] = None) -> np.ndarray:
    """
    Compute the asteroid block numbers in a range of asteroid numbers.
    INPUTS:
        block_size: Size of batches to run, e.g. 5000
        n0:         First asteroid number to process, inclusive, e.g. 0
        n1:         Last asteroid number to process, exclusive, e.g. 1,255,514.  
                    None means no upper limit
    OUTPUTS:
        ast_blocks: Numpy array of the the block numbers for asteroids in this range.
                    e.g. if asteroid numbers are in [1, 542,000) and [1,000,000, 1,256,000)
                    and block_size=1000, then blocks = [0, 1, ... 541, 1000, 1001, ... 1256]    
    """
    # Load distinct asteroid numbers
    orb_elt = load_orbital_elts()
    ast_nums = orb_elt.Num.values.astype(np.int32)
    # Filter for n0 <= ast_num
    filter_left = (n0 <= ast_nums)
    ast_nums = ast_nums[filter_left]
    # Filter for ast_num < n1
    if n1 is not None:
        filter_right = (ast_nums < n1)
        ast_nums = ast_nums[filter_right]

    # Extract distinct blocks of data
    ast_blocks = np.unique(ast_nums // block_size)
    return ast_blocks

# ********************************************************************************************************************* 
def process_ast_blocks(block_size: int, n0: int = 0, n1: Optional[int] = None) -> None:
    """
    Process multiple blocks of asteroids to find the one nearest to each ZTF observation
    INPUTS:
        block_size: Size of batches to run, e.g. 5000
        n0:         First asteroid number to process, inclusive, e.g. 0
        n1:         Last asteroid number to process, exclusive, e.g. 1,255,514
    OUTPUTS:
        None.  Saves the assembled files ztf_nearest_ast_n0_n1.h5 to disk.
    """
    # Calculate distinct blocks of asteroids from n0, n1 and block_size
    ast_blocks = calc_ast_block_numbers(block_size=block_size, n0=n0, n1=n1)

    # Create a list of arguments to run_batch
    run_batch_args = []
    for k, ast_block in enumerate(ast_blocks):
        # The start and end for each batch
        bn0 = ast_block * block_size
        bn1 = bn0 + block_size
        # progress bar only on the first batch
        progbar = (k == 0)
        # set verbose argument to false when processing all batches
        verbose = False
        # wrap the arguments as a tuple
        args = (bn0, bn1, progbar, verbose)
        # append args to full list of args
        run_batch_args.append(args)

    # Set number of processes based on number of available CPUs, and at most the number of blocks
    processes: int = int(np.round(cpu_count()*0.80))
    processes = min(processes, ast_blocks.size)

    # Create a thread pool and run the batches in parallel
    with Pool(processes=processes) as thread_pool:
        print(f'Created thread pool with {processes} threads.')
        # nearest_ast_blocks = thread_pool.starmap(func=run_batch, iterable=run_batch_args)
        thread_pool.starmap(func=run_batch, iterable=run_batch_args)

    # Status message
    print(f'Processed {n1-n0} asteroids from {n0} to {n1} in blocks of {block_size}.')

# ********************************************************************************************************************* 
def nearest_ast_reduction(ast_blocks: np.ndarray, block_size: int, progbar: bool = False) -> None:
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
    
    # Full range of asteroids
    ast_num_min = n0
    ast_num_max = ast_blocks[-1] * block_size + block_size

    # Check if this file exists; if so, quit early
    file_path = ztf_ast_file_path(n0=ast_num_min, n1=ast_num_max)
    if os.path.isfile(file_path):
        print(f'Found reduction file {file_path}.')
        return

    # load ztf_ast for first block and initialize it as the master
    ztf = load_ztf_nearest_ast(n0=n0, n1=n1)

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
        ztf_i = load_ztf_nearest_ast(n0=n0, n1=n1)
        # Mask for rows that need to be updated
        mask = ztf_i.nearest_ast_dist < ztf.nearest_ast_dist
        # Overwrite data columns for nearest asteroid
        ztf.loc[mask, cols] = ztf_i.loc[mask, cols]

    # Save reduced DataFrame to disk
    ztf.to_hdf(file_path, key='ztf_ast', mode='w')

# ********************************************************************************************************************* 
def main():
    """Main routine for integrating the orbits of known asteroids"""
    # Process command line arguments
    parser = argparse.ArgumentParser(description='Find the nearest asteroid to ZTF data.')
    parser.add_argument('-n0', nargs='?', metavar='n0', type=int, default=0,
                        help='the first asteroid number to process')
    parser.add_argument('-n1', nargs='?', metavar='n1', type=int, default=5000,
                        help='the last asteroid number to process')
    parser.add_argument('-block_size', nargs='?', metavar='block_size', type=int, default=1000,
                        help='the number of asteroids in each block')
    parser.add_argument('--all', default=False, action='store_true',
                        help='process all the data')
    parser.add_argument('--reduce', default=False, action='store_true',
                        help='perform reduction step only (consolidate nearest across multiple blocks)')
    # parser.add_argument('--test', default=False, action='store_true',
    #                     help='run in test mode')
    parser.add_argument('--progress', default=True, action='store_true',
                        help='display progress bar')

    # Parse arguments and alias inputs
    args = parser.parse_args()    
    n0: int = args.n0
    n1: int = args.n1
    block_size = args.block_size
    progbar: bool = args.progress

    # Settings to process all asteroids
    if args.all:
        n0 = 0
        n1 = None
        block_size = 1000

    # Status
    print(f'Processing asteroids from n0={n0} to n1={n1} with block_size={block_size}. progbar={progbar}.')

    # # If run in test mode, run tests without processing any ztf calculations
    # if args.test:
    #     # Test something
    #     # TODO put something here
        
    #     # Quit early in test mode: don't want to do any calculations
    #     print()
    #     exit()

    # Calculate distinct blocks of asteroids from n0, n1 and block_size
    ast_blocks = calc_ast_block_numbers(block_size=block_size, n0=n0, n1=n1)

    # Generate ztf_nearest_ast for each block; skip this if running in reduce mode
    if not args.reduce:
        process_ast_blocks(block_size=block_size, n0=n0, n1=n1)

    # Perform the reduction to get the globally closest asteroid
    nearest_ast_reduction(ast_blocks=ast_blocks, block_size=block_size, progbar=progbar)

    # Status update
    print(f'Completed reduction of asteroids from {n0} to {n1}.')
    exit()

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
