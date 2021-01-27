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

# Plotting
import matplotlib.pyplot as plt

# MSE imports
from utils import plot_style, print_stars
from astro_utils import mjd_to_date
from rebound_utils import make_sim_planets, make_sim_de435, integrate_df
from db_utils import df2db, df2db_chunked, sp2df

# Typing
from typing import List

# ********************************************************************************************************************* 
# Set plot style
plot_style()

# ********************************************************************************************************************* 
def process_sim(sim, collection_cd: str, mjd0: int, mjd1: int, steps_per_day: int, truncate: bool):
    """
    Integrate a simulation and save it to database
    INPUTS:
        sim:            Simulation for the desired collection of bodies as of the start epoch
        collection_cd:  Code of the collection of bodies, e.g. 'P' for Planets or 'D' for DE435
        mjd0:           First date to process
        mjd1:           Last date to process
        truncate:       Flag indicating whether to truncate DB tables
    """
    # Columns for state vectors and orbital element frames
    cols_vec = ['TimeID', 'BodyID', 'MJD', 'qx', 'qy', 'qz', 'vx', 'vy', 'vz']
    cols_elt = ['TimeID', 'BodyID', 'MJD', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'M']

    # Mapping to rename columns to match DB; MariaDB has case insensitive column names :(
    elt_col_map = {
        'Omega':'Omega_node', 
        'omega': 'omega_peri',
    }

    # Generate collection name from collection code
    collection_tbl = {
        'P' : 'Planets',
        'D' : 'DE435'
    }
    collection_name = collection_tbl[collection_cd]
    # Are we building the default collection with the Planets?
    is_planets: bool = (collection_cd == 'P')
    
    # Generate table names; add suffix with collection name for state vectors,
    # but only save orbital elements for the faster Planets integration
    table_name_vec = f'Integration_{collection_name}'
    table_name_elt = 'OrbitalElements'

    # Compute time_step from steps_per_day
    time_step: np.float64 = np.float64(1.0 / steps_per_day)

    # Flags for building simulation archive    
    save_elements: bool = is_planets  # save elements when we are running Planets integration, but not DE435
    progbar: bool = True

    # Set chunk_size for writing out DataFrame to database; DE435 export with ~700m rows crashed
    chunk_size: int = 2**19

    # Status
    print()
    print_stars()

    # Run the simulation and save as a DataFrame
    df = integrate_df(sim_epoch=sim, mjd0=mjd0, mjd1=mjd1, time_step=time_step, 
                        save_elements=save_elements, progbar=progbar)

    # DataFrame with the state vectors
    df_vec = df[cols_vec]

    # DataFrame with the orbital elements; excludes the Sun (other elements are relative to Sun)
    if save_elements:
        mask = (df.BodyName != 'Sun')
        df_elt = df[mask][cols_elt]
        df_elt.rename(columns=elt_col_map, inplace=True)

    # Save to Integration_<CollectionName> DB table
    print()
    print_stars()
    df2db_chunked(df=df_vec, schema='KS', table=table_name_vec, 
                    chunk_size=chunk_size, truncate=truncate, progbar=True)

    # Save to StateVectors table if we are running the Planets integration
    if is_planets:
        df2db_chunked(df=df_vec, schema='KS', table='StateVectors', 
                        chunk_size=chunk_size, truncate=truncate, progbar=True)
        
    # Save to OrbitalElements table if requested (if / only running planets)
    if save_elements:
        df2db_chunked(df=df_elt, schema='KS', table=table_name_elt, 
                        chunk_size=chunk_size, truncate=truncate, progbar=True)

# ********************************************************************************************************************* 
def main():
    """Integrate the orbits of the planets and major moons"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description='The solar system integrations against Horizons results.')
    parser.add_argument('--collection', nargs='?', metavar='COLL', type=str, default='p',
                        help='collection of bodies to integrate: p- planets; d- DE435; a-all (both planets and DE435)')
    parser.add_argument('--epoch', nargs='?', metavar='EP', type=int, default=59000,
                        help='epoch of the base simulation that is integrated forward and backwards, as an MJD')
    parser.add_argument('--half_width', nargs='?', metavar='HW', type=int, default=0,
                        help='half of the width of the interval in days, which is symmetric about the epoch.')
    parser.add_argument('--mjd0', nargs='?', metavar='t0', type=int, default=40400,
                        help='epoch of the first date in the integration, as an MJD.')
    parser.add_argument('--mjd1', nargs='?', metavar='t1', type=int, default=77600,
                        help='epoch of the last date in the integration, as an MJD.')
    parser.add_argument('--steps_per_day', nargs='?', metavar='SPD', type=int, default=24,
                        help='the (max) number of steps per day taken by the integrator')
    parser.add_argument('--truncate', dest='truncate', action='store_const', const=True, default=False,
                        help='Whether to truncate tables before inserting.')
    parser.add_argument('--regen_diff', dest='regen_diff', action='store_const', const=False, default=True,
                        help='Whether to regenerate DB table IntegrationDiff.')
    args = parser.parse_args()
    
    # Unpack command line arguments
    epoch: int = args.epoch                 # MJD 59000 = 2020-05-31
    steps_per_day: int = args.steps_per_day
    truncate: bool = args.truncate
    regen_diff: bool = args.regen_diff

    # Flags for planets and de435 integrations
    run_planets: bool = args.collection in ('p', 'a')
    run_de435: bool = args.collection in ('d', 'a')

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
    print(f'collection     : {args.collection} (run_planets = {run_planets}, run_de435 = {run_de435})')
    print(f'epoch          : {epoch} ({epoch_dt})')
    print(f'date range mjd : {mjd0} to {mjd1}')
    print(f'date range     : {mjd0_dt} to {mjd1_dt}')
    print(f'full width     : {width_yrs:3.1f} years')
    print(f'steps_per_day  : {steps_per_day}')
    print(f'times to save  : {times_saved}')
    print(f'truncate       : {truncate}')
    # print(f'regen int diff : {regen_diff}')

    # Set chunk_size for writing out DataFrame to database; DE435 export with ~700m rows crashed
    chunk_size: int = 2**19

    # If planets were requested, run the simulation and test the results
    if run_planets:
        # Simulation with initial configuration for planets
        sim = make_sim_planets(epoch=epoch, integrator=integrator, epsilon=epsilon, steps_per_day=steps_per_day)
        # Delegate to process_sim
        process_sim(sim=sim, collection_cd='P', mjd0=mjd0, mjd1=mjd1, steps_per_day=steps_per_day, truncate=truncate)

    # If DE435 was requested, run the simulation and test the results
    if run_de435:
        # Simulation with initial configuration for DE435
        sim = make_sim_de435(epoch=epoch, integrator=integrator, epsilon=epsilon, steps_per_day=steps_per_day)
        # Delegate to process_sim
        process_sim(sim=sim, collection_cd='D', mjd0=mjd0, mjd1=mjd1, steps_per_day=steps_per_day, truncate=truncate)

    # Regenerate IntegrationDiff table if requested
    if regen_diff:
        print(f'\nRegenerating table KS.Integration_Diff...')
        sp2df(sp_name='KS.MakeTable_IntegrationDiff', params=dict())

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
