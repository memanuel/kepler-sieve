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

# Utility
import argparse
import sys

# MSE imports
from utils import print_stars
from astro_utils import mjd_to_date
from asteroid_element import make_sim_asteroids, get_asteroids
from rebound_integrate import integrate_df
from db_utils import df2db, csvs2db, sp2df, truncate_table

# Typing
from typing import List

# ********************************************************************************************************************* 
# Type of calculation: state vectors or orbital elements
calc_type_tbl = {
    'vec' : 'StateVectors',
    'elt' : 'OrbitalElements'
}

# Load a single copy of the Asteroid list keyed by BodyID
ast = get_asteroids(key_to_body_id=True)

# ********************************************************************************************************************* 
def process_sim(sim, n0: int, n1: int, mjd0: int, mjd1: int, steps_per_day: int,  
                truncate: bool, single_thread: bool, progbar:bool=False):
    """
    Integrate a simulation and save it to database
    INPUTS:
        sim:            Simulation for the desired collection of bodies as of the start epoch
        n0:             First asteroid number to process (inclusive)
        n1:             Last asteroid number to process (exclusive)
        mjd0:           First date to process
        mjd1:           Last date to process
        single_thread:  Whether to run in single threaded mode
        truncate:       Flag indicating whether to truncate DB tables
    """
    # Columns for state vectors and orbital element frames
    cols_vec = ['TimeID', 'AsteroidID', 'MJD', 'qx', 'qy', 'qz', 'vx', 'vy', 'vz']
    cols_elt = ['TimeID', 'AsteroidID', 'MJD', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'M']

    # Mapping to rename columns to match DB; MariaDB has case insensitive column names :(
    elt_col_map = {
        'Omega':'Omega_node', 
        'omega': 'omega_peri',
    }

    # Generate table names; add suffix with collection name for state vectors,
    schema = 'KS'
    table_name_vec = f'AsteroidVectors'
    table_name_elt = f'AsteroidElements'

    # Set chunksize for writing out DataFrame to database
    chunksize: int = 2**19
    # Verbosity for DF to DB upload
    verbose: bool = progbar

    # Status
    print()
    print_stars()
    print(f'Integrating asteroids {n0:06d}-{n1:06d} from {mjd0} to {mjd1}...')

    # Run the simulation and save as a DataFrame
    df = integrate_df(sim_epoch=sim, mjd0=mjd0, mjd1=mjd1, steps_per_day=steps_per_day, 
                      save_elements=True, progbar=progbar)

    # Mask down to only the asteroids
    mask = (df.BodyID > 1000000)
    df = df[mask]

    # Add the AsteroidID column to the DataFrame
    df['AsteroidID'] = ast.loc[df.BodyID, 'AsteroidID'].values

    # DataFrame with the state vectors
    df_vec = df[cols_vec]

    # DataFrame with the orbital elements
    df_elt = df[cols_elt].rename(columns=elt_col_map)

    # Status
    print()
    print_stars()
    print(f'Saving from DataFrame to {table_name_vec}...')

    # Insert to StateVectors_<CollectionName> DB table
    try:
        df2db(df=df_vec, schema=schema, table=table_name_vec, truncate=truncate, 
              chunksize=chunksize, single_thread=single_thread, verbose=verbose)
    except:
        print("Problem with DB insertion... Continuing to save orbital elements.")

    # Insert to OrbitalElements_<CollectionName> DB table if requested
    print(f'\nSaving from DataFrame to {table_name_elt}...')
    df2db(df=df_elt, schema=schema, table=table_name_elt, truncate=truncate, chunksize=chunksize, verbose=verbose)

# ********************************************************************************************************************* 
def load_csv_batch(calc_type_cd: str):
    """
    Loads CSV files of the requested data type into database.
    INPUTS:
        calc_type_cd:   Code describing calculation type.
                        Must be one of 'vec', 'elt' for state vectors, orbital elements.
    OUTPUTS:
        None.  Loads tables into the database.
    """
    # Get calculation type
    table_prefix: str = calc_type_tbl[calc_type_cd]

    # Name of DB table
    schema: str = 'KS'
    table: str = f'{table_prefix}'

    # Delegate to csvs2db
    csvs2db(schema=schema, table=table)

# ********************************************************************************************************************* 
def load_csv_batches():
    """Program flow when in load_csv mode"""
    # List of the calculation types to save
    calc_type_cds = ['vec', 'elt']

    # Status
    print()
    print_stars()
    print(f'Loading CSVs for calculation types {calc_type_cds}.')

    # Iterate over desired calculation type
    for calc_type_cd in calc_type_cds:
        load_csv_batch(calc_type_cd=calc_type_cd)

# ********************************************************************************************************************* 
def main():
    """Integrate the orbits of the planets and major moons"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description='Integrate planets and selected asteroids in Rebound '
            'with initial conditions from Horizons (planets) and asteroid orbital elements.')
    parser.add_argument('n0', nargs='?', metavar='n0', type=int, default=0,
                        help='the first asteroid number to process')
    parser.add_argument('n_ast', nargs='?', metavar='B', type=int, default=1000,
                        help='the number of asteroids to process in this batch'),
    parser.add_argument('--epoch', nargs='?', metavar='EP', type=int, default=59000,
                        help='epoch of the base simulation that is integrated forward and backwards, as an MJD')
    parser.add_argument('--mjd0', nargs='?', metavar='t0', type=int, default=48000, # originally 40400
                        help='epoch of the first date in the integration, as an MJD.')
    parser.add_argument('--mjd1', nargs='?', metavar='t1', type=int, default=63000, # originally 77600
                        help='epoch of the last date in the integration, as an MJD.')
    parser.add_argument('--steps_per_day', nargs='?', metavar='SPD', type=int, default=1,
                        help='the (max) number of steps per day taken by the integrator')
    # parser.add_argument('--load_csv', dest='load_csv', action='store_const', const=True, default=False,
    #                     help="Don\'t do another solar system integration, just load cached CSV into DB.")
    parser.add_argument('--truncate', dest='truncate', action='store_const', const=True, default=False,
                        help='Whether to truncate tables before inserting.')
    parser.add_argument('--single_thread', dest='single_thread', action='store_const', const=True, default=False,
                        help='run in single threaded mode; better if running many parallel instances.')
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

    # Flags
    truncate: bool = args.truncate
    # load_csv: bool = args.load_csv
    single_thread: bool = args.single_thread
    progbar: bool = not args.quiet
    dry_run: bool = args.dry_run

    # Date range for integration
    mjd0: int = args.mjd0
    mjd1: int = args.mjd1
    steps_per_day: int = args.steps_per_day
    # Epoch as a date for reporting
    epoch_dt = mjd_to_date(epoch)
    mjd0_dt = mjd_to_date(mjd0)
    mjd1_dt = mjd_to_date(mjd1)
    width_yrs: float = (mjd1 - mjd0) / 365.25
    times_saved: int = (mjd1-mjd0) * steps_per_day

    # Integrator settings
    integrator: str = 'ias15'
    epsilon: float = 2.0**-32

    # Report arguments and integrator settings
    print_stars()
    print(f'*n0             : {n0:06d}')
    print(f'*n1             : {n1:06d}')
    print(f'*epoch          : {epoch} ({epoch_dt})')
    print(f' date range mjd : {mjd0} to {mjd1}')
    print(f'*steps_per_day  : {steps_per_day}')
    print(f' times to save  : {times_saved}')
    # print(f'*load_csv       : {load_csv}')
    print(f'*truncate       : {truncate}')
    print(f'*dry_run        : {dry_run}')

    # Quit early if it was a dry run
    if dry_run:
        print('\n This was a dry run.  Bye!')
        sys.exit()

    # # If we are just loading CSVs, don't do the integration, just try to reload them into DB and quit
    # if load_csv:
    #     load_csv_batches()
    #     # Now quit early after the CSVs are loaded
    #     sys.exit()

    # Set chunk_size for writing out DataFrame to database; DE435 export with ~700m rows crashed
    chunk_size: int = 2**19
    
    # Simulation with initial configuration for planets
    sim = make_sim_asteroids(epoch=epoch, n0=n0, n1=n1)
  
    # Delegate to process_sim
    process_sim(sim=sim, n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1, steps_per_day=steps_per_day, 
                truncate=truncate, single_thread=single_thread, progbar=progbar)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
