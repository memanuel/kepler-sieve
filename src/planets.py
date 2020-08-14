"""
Harvard IACS Masters Thesis
Trajectories for Known Asteroids

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
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
from rebound_utils import make_sim_planets, make_sim_de435, extend_sim_ast
from rebound_utils import make_archive, get_asteroids
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
    parser.add_argument('--collection', nargs='?', metavar='COLL', type=str, default='a',
                        help='collection of bodies to integrate: p- planets; d- DE435; a-all (both planets and DE435)')
    parser.add_argument('--epoch', nargs='?', metavar='EP', type=int, default=59000,
                        help='epoch of the base simulation that is integrated forward and backwards, as an MJD')
    parser.add_argument('--half_width', nargs='?', metavar='HW', type=int, default=2000,
                        help='half of the width of the interval in days, which is symmetric about the epoch 59000')
    parser.add_argument('--steps_per_day', nargs='?', metavar='SPD', type=int, default=-1,
                        help='the (max) number of steps per day taken by the integrator')
    parser.add_argument('--test', default=False, action='store_true',
                        help='run in test mode')
    args = parser.parse_args()
    
    # Unpack command line arguments
    epoch: int = args.epoch                 # MJD 59000 = 2020-05-31
    steps_per_day: int = args.steps_per_day if args.steps_per_day >= 0 else 16
    is_test_mode: bool = args.test

    # Flags for planets and de435 integrations
    run_planets: bool = args.collection in ('p', 'a')
    run_de435: bool = args.collection in ('d', 'a')

    # Collection of test asteroids - first n_ast numbered asteroids
    n_ast: int = 100
    test_name: str = f'{n_ast} Asteroids'
    ast: pd.DataFrame = get_asteroids()
    body_names_ast: List[str] = list(ast.BodyName[0:n_ast].values)
    test_bodies = ['Earth'] + body_names_ast

    # Date range for testing
    mjd0: int = epoch - args.half_width
    mjd1: int = epoch + args.half_width
    # Interval in days for testing snapshots
    test_step: int = 50
    
    # Epoch as a date for reporting
    epoch_dt = mjd_to_date(epoch)
    mjd0_dt = mjd_to_date(mjd0)
    mjd1_dt = mjd_to_date(mjd1)
    width_yrs: float = (mjd1 - mjd0) / 365.25

    # Integrator settings
    integrator: str = 'ias15'
    epsilon: float = 2.0**-32

    # Report arguments and integrator settings
    print(f'collection     : {args.collection} (run_planets = {run_planets}, run_de435 = {run_de435})')
    print(f'epoch          : {epoch} ({epoch_dt})')
    print(f'half_width     : {args.half_width} days')
    # print(f'date range mjd : {mjd0} to {mjd1}')
    print(f'date range     : {mjd0_dt} to {mjd1_dt}')
    print(f'full width     : {width_yrs:3.1f} years')
    print(f'steps_per_day  : {steps_per_day}')
    print(f'run tests      : {is_test_mode}')
    if is_test_mode:
        print(f'test_step      : {test_step}')
    print('')

    # Simulation with initial configuration for planets and DE435
    sim_planets = make_sim_planets(epoch=epoch, integrator=integrator, epsilon=epsilon, steps_per_day=steps_per_day)
    sim_de435 = make_sim_de435(epoch=epoch, integrator=integrator, epsilon=epsilon, steps_per_day=steps_per_day)

    # If running in test mode, need to add the first n_ast asteroids to the base simulation
    if is_test_mode:
        # We want to add the additional test asteroids as test particles
        add_as_test: bool = True
        extend_sim_ast(sim=sim_planets, n_ast=n_ast, add_as_test=add_as_test)
        extend_sim_ast(sim=sim_de435, n_ast=n_ast, add_as_test=add_as_test)

    # Shared time_step and save_step
    time_step: int = 1
    save_step: int = 1

    # Flags for building simulation archive
    save_elements: bool = False
    progbar: bool = True

    # Flags for testing integration
    verbose: bool = False
    make_plot: bool = False

    # If planets were requested, run the simulation and test the results
    if run_planets:
        print()
        print_stars()
        # Run the planets simulation
        sa_planets = make_archive(sim_epoch=sim_planets, 
                 mjd0=mjd0, mjd1=mjd1, time_step=time_step, save_step=save_step,
                 save_elements=save_elements, progbar=progbar)
        # Test the planets integration if requested
        if args.test:
            pos_err_planets, ang_err_planets = \
                test_integration(sa=sa_planets, test_bodies=test_bodies, 
                sim_name='Planets', test_name=test_name, test_step=test_step,
                verbose=verbose, make_plot=make_plot)

    # If DE435 was requested, run the simulation and test the results
    if run_de435:
        print()
        print_stars()
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

    # Plot comparison of results for different integrations
    plot_comparison: bool = (args.test and run_planets and run_de435)
    if plot_comparison:
        # Times to be tested - every step days in simulation range
        t_test = np.arange(mjd0, mjd1+1, test_step).astype(np.int32)

        # Epochs to be tested - generate from test times
        test_epochs = t_test + mjd0

        # Plot position error
        fig, ax = plt.subplots(figsize=[16,10])
        ax.set_title(f'Position Error on {n_ast} Test Asteroids')
        ax.set_ylabel('RMS Position Error in AU')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0,))
        if run_planets:
            ax.plot(test_epochs, pos_err_planets, label='Planets', marker='o', color='blue')
        if run_de435:
            ax.plot(test_epochs, pos_err_de435, label='DE435', marker='+', color='red')
        ax.grid()
        ax.legend()
        fig.savefig(fname=f'../figs/integration_test/planets/sim_pos_error_comp.png', bbox_inches='tight')

        # Plot angle error
        fig, ax = plt.subplots(figsize=[16,10])
        ax.set_title(f'Angle Error on {n_ast} Test Asteroids')
        ax.set_ylabel('RMS Angle Error vs. Earth in Arcseconds')
        if run_planets:
            ax.plot(test_epochs, ang_err_planets, label='Planets', marker='+', color='blue')
        if run_de435:
            ax.plot(test_epochs, ang_err_de435, label='DE435', marker='o', color='red')
        ax.grid()
        ax.legend()
        fig.savefig(fname=f'../figs/integration_test/planets/sim_ang_error_comp.png', bbox_inches='tight')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
