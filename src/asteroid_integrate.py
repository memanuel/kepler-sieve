"""
Trajectories for asteroids
Compute in rebound, save to database.
Example calls:
$ python asteroid_integrate.py 0 1000 --epoch 59000
$ python asteroid_integrate.py 0 1000

Michael S. Emanuel
2021-01-13
"""

# Core
import numpy as np
import pandas as pd

# Astronomy
import rebound

# File system
import os
import glob
from pathlib import Path

# Commandline arguments
import argparse
import sys

# MSE imports
from utils import print_stars
from astro_utils import mjd_to_date
from asteroid_element import make_sim_asteroids, get_asteroids
from rebound_integrate import integrate_df
from db_utils import df2csv, csvs2db, sp2df, clean_empty_dirs

# Typing
from typing import List, Tuple

# ********************************************************************************************************************* 
# Table names for state vectors and orbital elements
schema = 'KS'
table_vec = f'AsteroidVectors'
table_elt = f'AsteroidElements'

# Columns for state vectors and orbital element frames
cols_vec_df = ['TimeID', 'AsteroidID', 'mjd', 'qx', 'qy', 'qz', 'vx', 'vy', 'vz']
cols_elt_df = ['TimeID', 'AsteroidID', 'mjd', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'M']

# Mapping to rename orbital element columnd names
# This is necessary because DB column names are case insensitive, causing a name clash between 'Omega' and 'omega'
elt_col_map = {
    'Omega':'Omega_node', 
    'omega': 'omega_peri',
}

# The column names for the DB tables
cols_vec_db = cols_vec_df
cols_elt_db = [elt_col_map.get(c, None) or c for c in cols_elt_df]

# Directory where CSVs are saved
dir_csv: str = '../data/df2db'
# Get the process_id to avoid collisions on CSV file names for different threads
pid: int = os.getpid()
pid_str: str = f'pid_{pid:07d}'
# Base CSV files names
fname_csv_vec = os.path.join(dir_csv, table_vec, pid_str, f'{table_vec}.csv')
fname_csv_elt = os.path.join(dir_csv, table_elt, pid_str, f'{table_elt}.csv')

# ********************************************************************************************************************* 
def integrate_ast(sim: rebound.Simulation, mjd0: int, mjd1: int, 
                  interval: int, progbar:bool) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Integrate a simulation and return two DataFrames with state vectors and orbital elements
    INPUTS:
        sim:            Simulation for the desired collection of bodies as of the start epoch
        mjd0:           First date to process
        mjd1:           Last date to process
        interval:       Number of days between output frames saved out
        progbar:        Whether to show a progress bar
    """
    # Test whether we have any asteroids at all
    has_ast: bool = sim.asteroid_ids.size > 0
    # In the corner case where no asteroids included, don't integrate the planets to save time
    if not has_ast:
        mjd1 = mjd0

    # Status
    if progbar:        
        n0: int = sim.asteroid_ids[0] if has_ast else 0
        n1: int = sim.asteroid_ids[-1] if has_ast else 0
        n_ast: int = len(sim.asteroid_ids)
        print()
        print_stars()
        print(f'Integrating {n_ast} asteroids in {n0:07d}-{n1:07d} from {mjd0} to {mjd1}...')

    # Run the simulation and save as a DataFrame
    df = integrate_df(sim_epoch=sim, mjd0=mjd0, mjd1=mjd1, interval_p=interval, interval_q=1,
                      save_elements=True, progbar=progbar)

    # Mask down to only the asteroids
    mask = (df.BodyID > 1000000)
    df = df[mask]

    # Mapping table used to add the AsteroidID to the output frames
    ast_tbl = {
        'BodyID': sim.body_ids[sim.body_ids>1000000],
        'AsteroidID': sim.asteroid_ids,
    }
    ast = pd.DataFrame(ast_tbl)
    ast.set_index('BodyID', drop=False, inplace=True)    

    # Add the AsteroidID column to the DataFrame; search common asteroids frame for this
    df['AsteroidID'] = ast.loc[df.BodyID, 'AsteroidID'].values
    # DataFrame with the state vectors
    df_vec = df[cols_vec_df]
    # DataFrame with the orbital elements
    df_elt = df[cols_elt_df].rename(columns=elt_col_map)

    return df_vec, df_elt

# ********************************************************************************************************************* 
def report_csv_files(fnames_csv_vec, fnames_csv_elt, verbose: bool):
    """Report CSV file names"""
    if verbose:
        nf_vec = len(fnames_csv_vec)
        nf_elt = len(fnames_csv_elt)
        print(f'CSV files: {nf_vec} vectors and {nf_elt} elements:')
        # if nf_vec > 0:
        #     print(fnames_csv_vec[0])
        # if nf_elt > 0:
        #     print(fnames_csv_elt[0])

# ********************************************************************************************************************* 
def save_csvs(df_vec: pd.DataFrame, df_elt: pd.DataFrame, verbose:bool) -> [List[str], List[str]]:
    """
    Save DataFrames to CSV files
    INPUTS:
        df_vec:         DataFrame of asteroid state vectors
        df_elt:         DataFrame of asteroid orbital elements
        progbar:        Whether to show a progress bar
    """

    # Set chunksize for writing out DataFrame to database
    chunksize: int = 2**19

    # Initial status
    if verbose:
        print_stars()

    # Save a batch of CSV files for the state vectors
    if verbose:
        print(f'Saving from state vectors DataFrame to CSV files...')
    fnames_csv_vec = df2csv(df=df_vec, fname_csv=fname_csv_vec, columns=cols_vec_db, chunksize=chunksize)

    # Save a batch of CSV files for the orbital elements
    if verbose:
        print(f'Saving from orbital elements DataFrame to CSV files...')
    fnames_csv_elt = df2csv(df=df_elt, fname_csv=fname_csv_elt, columns=cols_elt_db, chunksize=chunksize)

    # Report results
    report_csv_files(fnames_csv_vec=fnames_csv_vec, fnames_csv_elt=fnames_csv_elt, verbose=verbose)

    # Return the names of the files
    return fnames_csv_vec, fnames_csv_elt

# ********************************************************************************************************************* 
def find_fnames_csv(verbose: bool) -> Tuple[List[str], List[str]]:
    """Generate a list of CSV file names for state vectors and orbital elements"""
    search_path_vec = os.path.join(dir_csv, table_vec, f'pid_*', f'{table_vec}-chunk-*.csv')
    search_path_elt = os.path.join(dir_csv, table_elt, f'pid_*', f'{table_elt}-chunk-*.csv')
    fnames_csv_vec = sorted(glob.glob(search_path_vec))
    fnames_csv_elt = sorted(glob.glob(search_path_elt))

    # Report results
    report_csv_files(fnames_csv_vec=fnames_csv_vec, fnames_csv_elt=fnames_csv_elt, verbose=verbose)

    return fnames_csv_vec, fnames_csv_elt

# ********************************************************************************************************************* 
def insert_csvs(fnames_csv_vec: List[str], fnames_csv_elt: List[str], progbar: bool) -> None:
    """
    Inserts CSV files into the database.
    INPUTS:
        fnames_csv_vec:   File names for CSVs with the state vectors
        fnames_csv_elt:   File names for CSVs with the orbital elements
        progbar:          Whether to display a progress bar
    OUTPUTS:
        None.  Loads tables into the database.
    """
    # Delegate to csvs2db
    csvs2db(schema=schema, table=table_vec, columns=cols_vec_db, fnames_csv=fnames_csv_vec, progbar=progbar)
    csvs2db(schema=schema, table=table_elt, columns=cols_elt_db, fnames_csv=fnames_csv_elt, progbar=progbar)

    # If the DB insert is successful, csvs2db deletes the CSV file. Clean up any empty directories now.
    clean_empty_dirs(table=table_vec)
    clean_empty_dirs(table=table_elt)

# ********************************************************************************************************************* 
def main():
    """Integrate the orbits of the planets and selected batch of asteroids"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description='Integrate planets and selected asteroids in Rebound '
            'with initial conditions from Horizons (planets) and asteroid orbital elements.')
    parser.add_argument('n0', nargs='?', metavar='n0', type=int, default=0,
                        help='the first asteroid number to process')
    parser.add_argument('n_ast', nargs='?', metavar='B', type=int, default=1000,
                        help='the number of asteroids to process in this batch'),
    parser.add_argument('mode', nargs='?', metavar='MODE', type=str, default='DB',
                        help='Mode of operation. Three valid choices DB, CSV, and INS. '
                        'DB: insert to DB via CSVs.'
                        'CSV: Calculate and save to CSVs.'
                        'INS: Insert CSVs from previous run into DB.'),
    parser.add_argument('--epoch', nargs='?', metavar='EP', type=int, default=59000,
                        help='epoch of the base simulation that is integrated forward and backwards, as an MJD')
    parser.add_argument('--mjd0', nargs='?', metavar='t0', type=int, default=48000, # originally 40400
                        help='epoch of the first date in the integration, as an MJD.')
    parser.add_argument('--mjd1', nargs='?', metavar='t1', type=int, default=63000, # originally 77600
                        help='epoch of the last date in the integration, as an MJD.')
    parser.add_argument('--interval', nargs='?', metavar='SPD', type=int, default=4,
                        help='the number of days between frames saved to the database')
    parser.add_argument('--run_all', const=True, default=False, action='store_const',
                        help='when true, run ALL asteroids; default is just to integrate missing ones.')
    parser.add_argument('--quiet', const=True, default=False, action='store_const',
                        help='run in quiet mode (hide progress bar')
    parser.add_argument('--dry_run', dest='dry_run', action='store_const', const=True, default=False,
                        help='Dry run: report parsed inputs then quit.')
    
    # Unpack command line arguments
    args = parser.parse_args()
    
    # Block of asteroids to integrate and epoch
    n0: int = args.n0
    n1: int = n0 + args.n_ast
    epoch: int = args.epoch

    # Operation mode
    mode: str = args.mode.upper()
    if mode not in ('DB', 'CSV', 'INS'):
        raise ValueErrror("Mode must be one of 'DB', 'CSV' or 'INS'.")
    mode_description_tbl = {
        'DB':  'Insert to database via CSVs.',
        'CSV': 'Calculate and save to CSVs to disk; must insert them later.',
        'INS': 'Insert CSVs from previous run into database.',
    }
    mode_description: str = mode_description_tbl[mode]

    # Flags
    run_all: bool = args.run_all
    missing: bool = not run_all
    verbose: bool = not args.quiet
    progbar: bool = not args.quiet
    dry_run: bool = args.dry_run

    # Date range for integration
    mjd0: int = args.mjd0
    mjd1: int = args.mjd1
    interval: int = args.interval
    # Epoch as a date for reporting
    epoch_dt = mjd_to_date(epoch)
    mjd0_dt = mjd_to_date(mjd0)
    mjd1_dt = mjd_to_date(mjd1)
    width_yrs: float = (mjd1 - mjd0) / 365.25
    times_saved: int = np.int32(np.ceil((mjd1-mjd0) / interval))

    # Integrator settings
    integrator: str = 'ias15'
    epsilon: float = 2.0**-32

    # Report arguments and integrator settings
    if verbose:
        print_stars()
        print(f'*n0             : {n0:06d}')
        print(f'*n1             : {n1:06d}')
        print(f'*epoch          : {epoch} ({epoch_dt})')
        print(f' date range mjd : {mjd0} to {mjd1}')
        print(f'*interval       : {interval}')
        print(f' times to save  : {times_saved}')
        print(f'*mode           : {mode}: {mode_description}')
        print(f'*run_all        : {run_all}')
        print(f'*dry_run        : {dry_run}')

    # Quit early if it was a dry run
    if dry_run:
        print('\n This was a dry run.  Bye!')
        sys.exit()

    # Set chunk_size for writing out DataFrame to database
    chunk_size: int = 2**19

    # Simulation with initial configuration for planets and selected asteroids; only take missing ones
    sim = make_sim_asteroids(epoch=epoch, n0=n0, n1=n1, missing=missing)
  
    # Delegate to appropriate functions depending on the mode
    # We need to integrate the asteroid orbits in either DB or CSV mode, but not in INS mode.
    # We also need to save the DataFrame to CSV in these modes.
    if mode in ('DB', 'CSV'):
        # Integrate the asteroids
        df_vec, df_elt = integrate_ast(sim=sim, mjd0=mjd0, mjd1=mjd1, interval=interval, progbar=progbar)
        # Save the DataFrames to CSV
        fnames_csv_vec, fnames_csv_elt = save_csvs(df_vec=df_vec, df_elt=df_elt, verbose=verbose)

    # If we are in insert mode, we don't have a list of file names yet.
    # Generate it by searching for matching files in the designated directory.
    if mode == 'INS':
        fnames_csv_vec, fnames_csv_elt = find_fnames_csv(verbose=verbose)

    # If we are in either DB or INS mode, we need to insert the CSV files to the database now
    if mode in ('DB', 'INS'):
        insert_csvs(fnames_csv_vec=fnames_csv_vec, fnames_csv_elt=fnames_csv_elt, progbar=progbar)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
