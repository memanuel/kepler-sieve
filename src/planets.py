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

# Utility
import argparse

# Plotting
import matplotlib.pyplot as plt

# MSE imports
from utils import plot_style, print_stars
from astro_utils import mjd_to_date
from rebound_utils import make_sim_planets, make_sim_de435
from rebound_archive import make_archive
from rebound_test import test_integration

# Typing
from typing import List

# ********************************************************************************************************************* 
# Set plot style
plot_style()

# ********************************************************************************************************************* 
def main():
    """Integrate the orbits of the planets and major moons"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description='The solar system integrations against Horizons results.')
    parser.add_argument('--collection', nargs='?', metavar='COLL', type=str, default='p',
                        help='collection of bodies to integrate: p- planets; d- DE435; a-all (both planets and DE435)')
    parser.add_argument('--epoch', nargs='?', metavar='EP', type=int, default=59000,
                        help='epoch of the base simulation that is integrated forward and backwards, as an MJD')
    parser.add_argument('--half_width', nargs='?', metavar='HW', type=int, default=2000,
                        help='half of the width of the interval in days, which is symmetric about the epoch.')
    parser.add_argument('--mjd0', nargs='?', metavar='t0', type=int, default=None,
                        help='epoch of the first date in the integration, as an MJD.')
    parser.add_argument('--mjd1', nargs='?', metavar='t1', type=int, default=None,
                        help='epoch of the last date in the integration, as an MJD.')
    parser.add_argument('--steps_per_day', nargs='?', metavar='SPD', type=int, default=-1,
                        help='the (max) number of steps per day taken by the integrator')
    args = parser.parse_args()
    
    # Unpack command line arguments
    epoch: int = args.epoch                 # MJD 59000 = 2020-05-31
    steps_per_day: int = args.steps_per_day if args.steps_per_day >= 0 else 16

    # Flags for planets and de435 integrations
    run_planets: bool = args.collection in ('p', 'a')
    run_de435: bool = args.collection in ('d', 'a')

    # Date range for testing
    mjd0: int = args.mjd0 or epoch - args.half_width
    mjd1: int = args.mjd1 or epoch + args.half_width
    
    # Epoch as a date for reporting
    epoch_dt = mjd_to_date(epoch)
    mjd0_dt = mjd_to_date(mjd0)
    mjd1_dt = mjd_to_date(mjd1)
    width_yrs: float = (mjd1 - mjd0) / 365.25
    dates_saved: int = (mjd1-mjd0) * steps_per_day

    # Integrator settings
    integrator: str = 'ias15'
    epsilon: float = 2.0**-32

    # Report arguments and integrator settings
    print_stars()
    print(f'collection     : {args.collection} (run_planets = {run_planets}, run_de435 = {run_de435})')
    print(f'epoch          : {epoch} ({epoch_dt})')
    # print(f'half_width     : {args.half_width} days')
    print(f'date range mjd : {mjd0} to {mjd1}')
    print(f'date range     : {mjd0_dt} to {mjd1_dt}')
    print(f'full width     : {width_yrs:3.1f} years')
    print(f'steps_per_day  : {steps_per_day}')
    print(f'dates to save  : {dates_saved}')

    # Shared time_step and save_step
    time_step: int = 1
    save_step: int = 1

    # Flags for building simulation archive
    save_elements: bool = False
    progbar: bool = True

    # If planets were requested, run the simulation and test the results
    if run_planets:
        print()
        print_stars()
        # Simulation with initial configuration for planets
        sim_planets = make_sim_planets(epoch=epoch, integrator=integrator, epsilon=epsilon, steps_per_day=steps_per_day)
        # Run the planets simulation
        sa_planets = make_archive(sim_epoch=sim_planets, 
                 mjd0=mjd0, mjd1=mjd1, time_step=time_step, save_step=save_step,
                 save_elements=save_elements, progbar=progbar)

    # If DE435 was requested, run the simulation and test the results
    if run_de435:
        print()
        print_stars()
        # Simulation with initial configuration for planets and DE435
        sim_de435 = make_sim_de435(epoch=epoch, integrator=integrator, epsilon=epsilon, steps_per_day=steps_per_day)
        # Run the DE435 simulation
        sa_de435 = make_archive(sim_epoch=sim_de435, 
                 mjd0=mjd0, mjd1=mjd1, time_step=time_step, save_step=save_step,
                 save_elements=save_elements, progbar=progbar)
        # Test the DE435 integration if requested
        if args.test: 
            pos_err_de435, ang_err_de435 = \
                test_integration(sa=sa_de435, test_bodies=test_bodies, 
                sim_name='DE435', test_name=test_name, test_step=test_step,
                verbose=verbose, make_plot=make_plot)        

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
