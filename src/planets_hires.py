"""
Hi resolution trajectory for Earth and the Sun.
Compute in rebound, save to database.
Example calls:
$ python planets_hires.py
$ python planets_hires.py --epoch 59000 --mjd0 57000 --mjd1 61000

Michael S. Emanuel
2021-05-26
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
from rebound_sim import make_sim_planets
from rebound_integrate import integrate_df
from db_utils import df2db

# Typing
from typing import List

# ********************************************************************************************************************* 
def process_sim(sim, mjd0: int, mjd1: int, steps_per_day: int):
    """
    Integrate a simulation and save it to database
    INPUTS:
        sim:            Simulation for the planets as of the start epoch
        mjd0:           First date to process
        mjd1:           Last date to process
    """
    # Convert steps per day to interval_p and interval_q
    interval_p: int = 1
    interval_q: int = steps_per_day
    # Flags for building simulation archive    
    save_elements: bool = False
    progbar: bool = True
    # Set chunksize for writing out DataFrame to database
    chunksize: int = 2**19
    # Verbosity for DF to DB upload
    verbose: bool = True

    # Status
    print()
    print_stars()
    print(f'Integrating planets from {mjd0} to {mjd1}...')

    # Run the simulation and save as a DataFrame
    df = integrate_df(sim_epoch=sim, mjd0=mjd0, mjd1=mjd1, interval_p=interval_p, interval_q=interval_q, 
                      save_elements=save_elements, progbar=progbar)

    # Mask for earth and sun; columns for state vectors
    mask_earth = (df.BodyID == 399)
    mask_sun = (df.BodyID == 10)
    cols = ['TimeID', 'mjd', 'qx', 'qy', 'qz', 'vx', 'vy', 'vz']

    # The masked DataFrame for earth and sun
    df_earth = df[mask_earth][cols]
    df_sun = df[mask_sun][cols]

    # Table names to write output
    schema = 'KS'
    table_name_earth = 'StateVectors_Earth'
    table_name_sun = 'StateVectors_Sun'

    # Status
    print()
    print_stars()
    print(f'Saving from DataFrame to {table_name_earth} and {table_name_sun}...')

    # Insert to StateVectors_Earth DB table
    df2db(df=df_earth, schema=schema, table=table_name_earth, chunksize=chunksize, verbose=verbose)
    # Insert to StateVectors_Sun DB table
    df2db(df=df_sun, schema=schema, table=table_name_sun, chunksize=chunksize, verbose=verbose)

    # Save data frames as HDF5 files
    fname_h5 = '../data/planets/StateVectors_HiRes.h5'
    key_earth = 'df_earth'
    key_sun = 'df_sun'
    print(f'Saving data frames to {fname_h5} with key names {key_earth} and {key_sun}...')
    df_earth.to_hdf(fname_h5, key=key_earth, mode='w')
    df_sun.to_hdf(fname_h5, key=key_sun, mode='a')

# ********************************************************************************************************************* 
def main():
    """Integrate the orbits of the planets and major moons"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description='Integrate solar system objects in Rebound '
            'with initial conditions from Horizons in the DE-435 integration.')
    parser.add_argument('--epoch', nargs='?', metavar='EP', type=int, default=59000,
                        help='epoch of the base simulation that is integrated forward and backwards, as an MJD')
    parser.add_argument('--half_width', nargs='?', metavar='HW', type=int, default=0,
                        help='half of the width of the interval in days, which is symmetric about the epoch.')
    parser.add_argument('--mjd0', nargs='?', metavar='t0', type=int, default=48000, # originally 40400
                        help='epoch of the first date in the integration, as an MJD.')
    parser.add_argument('--mjd1', nargs='?', metavar='t1', type=int, default=63000, # originally 77600
                        help='epoch of the last date in the integration, as an MJD.')
    parser.add_argument('--dry_run', dest='dry_run', action='store_const', const=True, default=False,
                        help='Dry run: report parsed inputs then quit.')
    args = parser.parse_args()
    
    # Unpack command line arguments
    epoch: int = args.epoch
    dry_run: bool = args.dry_run

    # The steps per day is always consistent with one minute intervals    
    steps_per_day: int = 24*60

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
    print(f'*epoch          : {epoch} ({epoch_dt})')
    print(f' date range mjd : {mjd0} to {mjd1}')
    print(f' date range     : {mjd0_dt} to {mjd1_dt}')
    print(f' full width     : {width_yrs:3.1f} years')
    print(f'*steps_per_day  : {steps_per_day}')
    print(f' times to save  : {times_saved}')
    print(f'*dry_run        : {dry_run}')

    # Quit early if it was a dry run
    if dry_run:
        print('\n This was a dry run.  Bye!')
        sys.exit()

    # Set chunk_size for writing out DataFrame to database; DE435 export with ~700m rows crashed
    chunk_size: int = 2**19
    
    # Simulation with initial configuration for planets
    sim = make_sim_planets(epoch=epoch, integrator=integrator, epsilon=epsilon, 
                        steps_per_day=steps_per_day, load_file=False)
    # Delegate to process_sim
    process_sim(sim=sim, mjd0=mjd0, mjd1=mjd1, steps_per_day=steps_per_day)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
