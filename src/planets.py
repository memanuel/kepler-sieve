"""
Trajectories for Planets
Compute in rebound, save to database.
Example calls:
$ python planets.py --collection a --epoch 59000 --mjd0 57000 --mjd1 61000
$ python planets.py --collection a --epoch 59000 --half_width 2000
$ python planets.py --collection a --epoch 59000
$ python planets.py --epoch 59000
$ python planets.py

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
from rebound_sim import make_sim_planets, make_sim_de435
from rebound_integrate import integrate_df
from db_utils import df2db, csvs2db, sp2df, truncate_table

# Typing
from typing import List

# ********************************************************************************************************************* 
# Collections of bodies to integrate: Planets or DE435
collection_name_tbl = {
    'P' : 'Planets',
    'D' : 'DE435'
}

# Type of calculation: state vectors or orbital elements
calc_type_tbl = {
    'vec' : 'StateVectors',
    'elt' : 'OrbitalElements'
}

# ********************************************************************************************************************* 
def process_sim(sim, collection_cd: str, mjd0: int, mjd1: int, steps_per_day: int, 
                save_elements: bool, truncate: bool):
    """
    Integrate a simulation and save it to database
    INPUTS:
        sim:            Simulation for the desired collection of bodies as of the start epoch
        collection_cd:  Code of the collection of bodies, e.g. 'P' for Planets or 'D' for DE435
        mjd0:           First date to process
        mjd1:           Last date to process
        save_elements:  Whether to save orbital elements (True) or just state vectors (False)
        truncate:       Flag indicating whether to truncate DB tables
    """
    # Columns for state vectors and orbital element frames
    cols_vec = ['TimeID', 'BodyID', 'mjd', 'qx', 'qy', 'qz', 'vx', 'vy', 'vz']
    cols_elt = ['TimeID', 'BodyID', 'mjd', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'M']

    # Mapping to rename columns to match DB; MariaDB has case insensitive column names :(
    elt_col_map = {
        'Omega':'Omega_node', 
        'omega': 'omega_peri',
    }

    # Generate collection name from collection code
    collection_name = collection_name_tbl[collection_cd]
    
    # Generate table names; add suffix with collection name for state vectors,
    schema = 'KS'
    table_name_vec = f'StateVectors_{collection_name}'
    table_name_elt = f'OrbitalElements_{collection_name}'

    # Flags for building simulation archive    
    progbar: bool = True
    # Set chunksize for writing out DataFrame to database; DE435 export with ~700m rows crashed
    chunksize: int = 2**19
    # Verbosity for DF to DB upload
    verbose: bool = True

    # Status
    print()
    print_stars()
    print(f'Integrating {collection_name} from {mjd0} to {mjd1}...')

    # Run the simulation and save as a DataFrame
    df = integrate_df(sim_epoch=sim, mjd0=mjd0, mjd1=mjd1, steps_per_day=steps_per_day, 
                      save_elements=save_elements, progbar=progbar)

    # DataFrame with the state vectors
    df_vec = df[cols_vec]

    # DataFrame with the orbital elements; excludes the Sun (other elements are relative to Sun)
    if save_elements:
        mask = (df.BodyID != 10)
        df_elt = df[mask][cols_elt]
        df_elt.rename(columns=elt_col_map, inplace=True)

    # Status
    print()
    print_stars()
    print(f'Saving from DataFrame to {table_name_vec}...')

    # Insert to StateVectors_<CollectionName> DB table
    try:
        df2db(df=df_vec, schema=schema, table=table_name_vec, truncate=truncate, chunksize=chunksize, verbose=verbose)
    except:
        if save_elements:
            print("Problem with DB insertion... Continuing to save orbital elements.")

    # Insert to OrbitalElements_<CollectionName> DB table if requested
    if save_elements:
        print(f'\nSaving from DataFrame to {table_name_elt}...')
        df2db(df=df_elt, schema=schema, table=table_name_elt, truncate=truncate, chunksize=chunksize, verbose=verbose)

# ********************************************************************************************************************* 
def load_csv_batch(collection_cd: str, calc_type_cd: str):
    """
    Loads CSV files of the requested data type into database.
    INPUTS:
        collection_cd:  Code describing source of integration; 'P' for Planets or 'D' for DE435.
        calc_type_cd:   Code describing calculation type.
                        Must be one of 'vec', 'elt' for state vectors, orbital elements.
    OUTPUTS:
        None.  Loads tables into the database.
    """
    # Get calculation type and collection
    collection_name: str = collection_name_tbl[collection_cd]
    table_prefix: str = calc_type_tbl[calc_type_cd]

    # Name of DB table
    schema: str = 'KS'
    table: str = f'{table_prefix}_{collection_name}'

    # Delegate to csvs2db
    csvs2db(schema=schema, table=table)

# ********************************************************************************************************************* 
def load_csv_batches(collection_cd: str, save_elements: bool):
    """Program flow when in load_csv mode"""
    # List of the calculation types to save
    calc_type_cds = ['vec']
    if save_elements:
        calc_type_cds.append('elt')

    # Status
    print()
    print_stars()
    print(f'Loading CSVs for calculation types {calc_type_cds} and data collections {collection_cds}.')

    # Iterate over desired calculation type
    for calc_type_cd in calc_type_cds:
        load_csv_batch(collection_cd=collection_cd, calc_type_cd=calc_type_cd)

# ********************************************************************************************************************* 
def main():
    """Integrate the orbits of the planets and major moons"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description='Integrate solar system objects in Rebound '
            'with initial conditions from Horizons in the DE-435 integration.')
    parser.add_argument('--collection', nargs='?', metavar='COLL', type=str, default='p',
                        help='collection of bodies to integrate: p- planets or d- DE435')
    parser.add_argument('--epoch', nargs='?', metavar='EP', type=int, default=59000,
                        help='epoch of the base simulation that is integrated forward and backwards, as an MJD')
    parser.add_argument('--half_width', nargs='?', metavar='HW', type=int, default=0,
                        help='half of the width of the interval in days, which is symmetric about the epoch.')
    parser.add_argument('--mjd0', nargs='?', metavar='t0', type=int, default=48000, # originally 40400
                        help='epoch of the first date in the integration, as an MJD.')
    parser.add_argument('--mjd1', nargs='?', metavar='t1', type=int, default=63000, # originally 77600
                        help='epoch of the last date in the integration, as an MJD.')
    parser.add_argument('--steps_per_day', nargs='?', metavar='SPD', type=int, default=0,
                        help='the (max) number of steps per day taken by the integrator')
    parser.add_argument('--skip_elements', dest='skip_elements', action='store_const', const=True, default=False,
                        help='Skip orbital elements and only save state vectors.')
    parser.add_argument('--load_csv', dest='load_csv', action='store_const', const=True, default=False,
                        help="Don\'t do another solar system integration, just load cached CSV into DB.")
    parser.add_argument('--truncate', dest='truncate', action='store_const', const=True, default=False,
                        help='Whether to truncate tables before inserting.')
    # parser.add_argument('--regen_diff', dest='regen_diff', action='store_const', const=True, default=False,
    #                    help='Whether to regenerate DB table IntegrationDiff.')
    parser.add_argument('--dry_run', dest='dry_run', action='store_const', const=True, default=False,
                        help='Dry run: report parsed inputs then quit.')
    args = parser.parse_args()
    
    # Unpack command line arguments
    collection_cd: str = args.collection.upper()
    collection_name: str = collection_name_tbl[collection_cd]
    if collection_cd not in ('P', 'D'):
        raise ValueError("Collection code must be one of 'P' (Planets) or 'D' (DE435)")
    epoch: int = args.epoch
    
    # Default steps_per_day is 288 for Planets, 48 for DE435
    spd_dflt_tbl = {
        'P' : 288,
        'D' : 48,
    }
    steps_per_day: int = args.steps_per_day if args.steps_per_day else spd_dflt_tbl[collection_cd]

    # Flags
    save_elements: bool = not args.skip_elements
    truncate: bool = args.truncate
    # regen_diff: bool = args.regen_diff
    load_csv: bool = args.load_csv
    dry_run: bool = args.dry_run

    # Date range for testing
    mjd0: int = args.mjd0
    mjd1: int = args.mjd1
    # Override date inputs if half_width was specified
    if args.half_width > 0:
        mjd0 = epoch - args.half_width
        mjd1 = epoch + args.half_width
    
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
    print(f'*collection     : {collection_cd} ({collection_name})')
    print(f'*epoch          : {epoch} ({epoch_dt})')
    print(f' date range mjd : {mjd0} to {mjd1}')
    print(f' date range     : {mjd0_dt} to {mjd1_dt}')
    print(f' full width     : {width_yrs:3.1f} years')
    print(f'*steps_per_day  : {steps_per_day}')
    print(f' times to save  : {times_saved}')
    print(f'*skip_elements  : {not save_elements}')
    print(f'*load_csv       : {load_csv}')
    print(f'*truncate       : {truncate}')
    # print(f'regen int diff : {regen_diff}')
    print(f'*dry_run        : {dry_run}')

    # Quit early if it was a dry run
    if dry_run:
        print('\n This was a dry run.  Bye!')
        sys.exit()

    # If we are just loading CSVs, don't to the integration, just try to reload them into DB and quit
    if load_csv:
        load_csv_batches(save_elements=save_elements, run_planets=run_planets, run_de435=run_de435)
        # Now quit early after the CSVs are loaded
        sys.exit()

    # Set chunk_size for writing out DataFrame to database; DE435 export with ~700m rows crashed
    chunk_size: int = 2**19
    
    # Pick function to make the initial simulation based on the collection code
    sim_func_tbl = {
        'P' : make_sim_planets,
        'D' : make_sim_de435,
    }
    make_sim_func = sim_func_tbl[collection_cd]

    # Simulation with initial configuration for planets
    sim = make_sim_func(epoch=epoch, integrator=integrator, epsilon=epsilon, 
                        steps_per_day=steps_per_day, load_file=False)
    # Delegate to process_sim
    process_sim(sim=sim, collection_cd=collection_cd, mjd0=mjd0, mjd1=mjd1, steps_per_day=steps_per_day, 
                save_elements=save_elements, truncate=truncate)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
