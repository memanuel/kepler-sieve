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

# MSE imports
from utils import plot_style, print_stars
from rebound_utils import make_sim_planets, make_sim_de435
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
    parser.add_argument('type', nargs='?', metavar='type', type=str, default='a',
                        help='type of integration: p- planets; d- DE435; a-all')
    parser.add_argument('steps_per_day', nargs='?', metavar='SPD', type=int, default=-1,
                        help='the (max) number of steps per day taken by the integrator')
    parser.add_argument('--test', default=False, action='store_true',
                        help='run in test mode')
    args = parser.parse_args()
    
    # Unpack command line arguments
    integration_type = args.type
    steps_per_day: int = args.steps_per_day
    
    # Flags for planets and de435 integrations
    run_planets: bool = integration_type in ('p', 'a')
    run_de435: bool = integration_type in ('d', 'a')

    # Collection of test asteroids - first 100 numbered asteroids
    n_ast: int = 100
    test_name: str = f'{n_ast} Asteroids'
    ast = get_asteroids()
    body_names_add = list(ast.BodyName[0:n_ast].values)
    test_bodies = ['Earth'] + body_names_add

    # Date range for testing
    epoch: int = 59000      # 2020-05-31
    # half_width: int = 500
    half_width: int = 7300
    mjd0: int = epoch - half_width
    mjd1: int = epoch + half_width
    # Interval in days for testing snapshots
    test_step: int = 50

    # Integrator settings
    integrator: str = 'ias15'
    epsilon: float = 2.0**-32

    # Integrator time step
    steps_per_day: int = steps_per_day if steps_per_day >= 0 else 16

    # Simulation with initial configuration for planets and DE435
    sim_planets = make_sim_planets(epoch=epoch, body_names_add=body_names_add, 
                                   integrator=integrator, epsilon=epsilon, steps_per_day=steps_per_day)
    sim_de435 = make_sim_de435(epoch=epoch, body_names_add=body_names_add,
                               integrator=integrator, epsilon=epsilon, steps_per_day=steps_per_day)    

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
        ax.set_title(f'Position Error on 100 Test Asteroids')
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
        ax.set_title(f'Angle Error on 100 Test Asteroids')
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
