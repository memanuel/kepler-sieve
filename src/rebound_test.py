"""
Harvard IACS Masters Thesis
Utilities for working with Rebound simulations and archives

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

# Core
import numpy as np
import pandas as pd

# Astronomy
import rebound

# Utilities
from datetime import datetime
import argparse

# Plot
import matplotlib.pyplot as plt

# MSE imports
from utils import rms, plot_style, print_stars
from astro_utils import mjd_to_datetime, datetime_to_year, cart_to_sph
from horizons import make_sim_horizons
from rebound_utils import make_sim_planets, make_sim_de435, extend_sim
from rebound_utils import make_archive, get_asteroids

# Typing
from typing import List, Tuple, Dict, Set, Optional

# ********************************************************************************************************************* 
# Set plot style
plot_style()

# ********************************************************************************************************************* 
def sim_cfg_array(sim: rebound.Simulation, body_names: Optional[List[str]]=None) -> np.array:
    """Extract the Cartesian configuration of each body in the simulation"""
    # Allocate array of configurations
    num_objects: int = sim.N if body_names is None else len(body_names)
    cfgs: np.array = np.zeros(shape=(num_objects, 6))

    # Iterate through particles
    ps = sim.particles

    # Create list of (i, p) pairs
    if body_names is None:
        ips = enumerate(ps)
    else:
        ips = [(i, ps[body_name]) for i, body_name in enumerate(body_names)]

    # Extract the configuration of each particle
    for i, p in ips:
        cfgs[i] = np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz])

    return cfgs

# ********************************************************************************************************************* 
def sim_elt_array(sim: rebound.Simulation, body_names=None) -> np.array:
    """Extract the orbital elements of each body in the simulation"""
    # Allocate array of elements
    num_objects: int = sim.N-1 if body_names is None else len(body_names)
    elts: np.array = np.zeros(shape=(num_objects, 9))
    
    # Iterate through particles AFTER the primary
    ps: List[rebound.Particles] = sim.particles
    primary: rebound.Particle = ps[0]
    if body_names is None:
        for i, p in enumerate(ps[1:]):
            orb = p.calculate_orbit(primary=primary)
            elts[i] = np.array([orb.a, orb.e, orb.inc, orb.Omega, orb.omega, 
                                orb.f, orb.M, orb.pomega, orb.l])
    else:
        for i, body_name in enumerate(body_names):
            # look up the particle with this name
            p = ps[body_name]
            orb = p.calculate_orbit(primary=primary)
            elts[i] = np.array([orb.a, orb.e, orb.inc, orb.Omega, orb.omega, 
                                orb.f, orb.M, orb.pomega, orb.l])
    return elts

# ********************************************************************************************************************* 
def report_sim_difference(sim0: rebound.Simulation, sim1: rebound.Simulation, 
                          body_names: Optional[List[str]] = None, verbose: bool=False) -> \
                          Tuple[np.array, np.array]:
    """Report the difference between two simulations on a summary basis"""
    # Get the body names
    if body_names is None:
        body_names = sim0.body_names

    # Extract configuration arrays for the two simulations
    cfg0: np.array = sim_cfg_array(sim0, body_names)
    cfg1: np.array = sim_cfg_array(sim1, body_names)
    
    # Displacement of each body to earth
    earth_idx: int = np.argmax(body_names == 'Earth')
    q0: np.array = cfg0[:, 0:3] - cfg0[earth_idx, 0:3]
    q1: np.array = cfg1[:, 0:3] - cfg1[earth_idx, 0:3]
    
    # Right Ascension and Declination
    r0, asc0, dec0 = cart_to_sph(q0)
    r1, asc1, dec1 = cart_to_sph(q1)

    # Error in asc and dec; convert from radians to arcseconds
    asc_err: np.array = np.degrees(np.abs(asc1-asc0)) * 3600
    dec_err: np.array = np.degrees(np.abs(dec1-dec0)) * 3600

    # Take differences
    cfg_diff: np.array = (cfg1 - cfg0)
    pos_diff: np.array = cfg_diff[:, 0:3]
    vel_diff: np.array = cfg_diff[:, 3:6]

    # Error in position and velocity in barycentric coordinates; skip the sun
    pos_err: np.array = np.linalg.norm(pos_diff, axis=1)
    pos_err_den: np.array = np.linalg.norm(cfg0[:, 0:3], axis=1)
    pos_err_den[0] = 1.0
    pos_err_rel: np.array = pos_err / pos_err_den
    vel_err: np.array = np.linalg.norm(vel_diff, axis=1)
    vel_err_den: np.array = np.linalg.norm(cfg0[:, 3:6], axis=1)
    vel_err_den[0] = 1.0
    vel_err_rel: np.array = vel_err / vel_err_den

    if verbose:
        print(f'\nPosition difference - absolute & relative')
        print(f'(Angle errors in arcseconds, position in AU)')
        print(f'Body       : Phi     : Theta   : Pos AU  : Pos Rel : Vel Rel')
        body_names_short: List[str] = [nm.replace(' Barycenter', '') for nm in body_names]
        for i, nm in enumerate(body_names_short):
            print(f'{nm:10} : {asc_err[i]:5.2e}: {dec_err[i]:5.2e}: {pos_err[i]:5.2e}: '
                  f'{pos_err_rel[i]:5.2e}: {vel_err_rel[i]:5.2e}')
        print(f'Overall    : {rms(asc_err):5.2e}: {rms(dec_err):5.2e}: {rms(pos_err):5.2e}: '
              f'{rms(pos_err_rel):5.2e}: {rms(vel_err_rel):5.2e}')

    # Extract orbital element arrays from the two simulations
    elt0: np.array = sim_elt_array(sim0, body_names[1:])
    elt1: np.array = sim_elt_array(sim1, body_names[1:])

    # Take differences
    elt_diff: np.array = (elt1 - elt0)
    # Angle differences are mod two pi
    two_pi: float = 2.0 * np.pi
    elt_diff[:,2:] = (elt_diff[:,2:] +np.pi ) % two_pi - np.pi
    
    # Compute RMS difference by orbital element
    elt_rms: np.array = rms(elt_diff, axis=0)
    elt_err: np.array = np.abs(elt_diff)

    # Names of selected elements
    elt_names: List[str] = ['a', 'e', 'inc', 'Omega', 'omega', 'f', 'M', 'pomega', 'long']

    # Report RMS orbital element differences
    if verbose:
        print(f'\nOrbital element errors:')
        print(f'elt    : RMS      : worst      : max_err  : HRZN        : REB')
        for j, elt in enumerate(elt_names):
            idx = np.argmax(elt_err[:, j])
            worse = body_names_short[idx+1]
            print(f'{elt:6} : {elt_rms[j]:5.2e} : {worse:10} : {elt_err[idx, j]:5.2e} : '
                  f'{elt0[idx, j]:11.8f} : {elt1[idx, j]:11.8f}')
        print(f'RMS (a, e, inc) =          {rms(elt_diff[:,0:3]):5.2e}')
        print(f'RMS (f, M, pomega, long) = {rms(elt_diff[:,5:9]):5.2e}')      

    # One summary error statistic
    ang_err: np.array = rms(np.array([asc_err, dec_err]))
    
    # Return the RMS position error and angle errors
    return pos_err, ang_err
    
# ********************************************************************************************************************* 
def test_integration(sa: rebound.SimulationArchive,
                     test_bodies: Optional[List[str]], 
                     sim_name: str, 
                     test_name: str,
                     test_step: int = 20, 
                     verbose: bool = False,
                     make_plot: bool = False) -> \
                     Tuple[np.array, np.array]:
    """Test the integration of the planets against Horizons data"""
    # Extract body_collection from the simulation archive
    body_collection: str = sa.body_collection

    # Test bodies is input manually in case it's a strict subset of the bodies in the collection
    if test_bodies is None:
        test_bodies = sa.sim_epoch.body_names

    # Extract mjd0, mjd1 from simulation archive
    mjd0: int = sa.mjd0
    mjd1: int = sa.mjd1

    # Times to be tested - every step days in simulation range
    t_test = np.arange(np.int32(sa.tmin), np.int32(sa.tmax)+1, np.int32(test_step)).astype(np.int32)

    # Epochs to be tested - generate from test times
    test_epochs = t_test + mjd0
    sz: int = test_epochs.shape[0]
    # Generate verbose output on the last test epoch if requested
    verbose_epochs = test_epochs[sz-1:sz]

    # Convert epochs to dates and years for plotting
    test_dates = [mjd_to_datetime(mjd) for mjd in test_epochs]
    test_years = [datetime_to_year(test_date) for test_date in test_dates]

    # Errors on these dates
    ang_errs: List[np.array] = []
    pos_errs: List[np.array] = []
    
    # Test the specified times
    print(f'Testing accuracy of {sim_name} integration vs. Horizons on {sz} dates from {mjd0} to {mjd1}...')
    for t in t_test:
        # The epoch of this test time
        epoch: int = mjd0 + t
        # The reference simulation from Horizons, including test bodies
        sim0: rebound.Simulation = make_sim_horizons(body_collection=body_collection, epoch=epoch)
        extend_sim(sim=sim0, body_names=test_bodies, add_as_test=True)
        # The test simulation from the simulation archive
        sim1: rebound.Simulation = sa.getSimulation(t=t, mode='exact')
        # Verbosity flag and screen print if applicable
        report_this_date: bool = (epoch in verbose_epochs) and verbose
        if report_this_date:
            dt_t: datetime = mjd_to_datetime(epoch)
            print(f'\nDifference on {dt_t}:')
        # Run the test
        pos_err, ang_err = report_sim_difference(sim0=sim0, sim1=sim1, 
                                                 body_names=test_bodies, verbose=report_this_date)
        # Save position and angle errors
        pos_errs.append(pos_err)
        ang_errs.append(ang_err)
    
    # Plot error summary
    pos_err_rms: np.array = np.array([rms(x) for x in pos_errs])
    ang_err_rms: np.array = np.array([rms(x) for x in ang_errs])
        
    if make_plot:
        # Chart titles
        sim_name_chart = sim_name
        test_name_chart = test_name

        # Error in the position
        fig, ax = plt.subplots(figsize=[16,10])
        ax.set_title(f'Position Error of {test_name_chart} in {sim_name_chart} Integration')
        ax.set_ylabel(f'RMS Position Error in AU')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0,))
        # ax.plot(test_years, pos_err_rms, marker='o', color='red')
        ax.plot(test_epochs, pos_err_rms, marker='o', color='red')
        ax.grid()
        fname: str = f'../figs/integration_test/sim_error_{sim_name}_{test_name}_pos.png'
        fig.savefig(fname=fname, bbox_inches='tight')

        # Error in the angle
        fig, ax = plt.subplots(figsize=[16,10])
        ax.set_title(f'Angle Error of {test_name_chart} in {sim_name_chart} Integration')
        ax.set_ylabel(f'RMS Angle Error vs. Earth in Arcseconds')
        # ax.plot(test_years, ang_err_rms, marker='o', color='blue')
        ax.plot(test_epochs, ang_err_rms, marker='o', color='blue')
        ax.grid()
        fname: str = f'../figs/integration_test/sim_error_{sim_name}_{test_name}_angle.png'
        fig.savefig(fname=fname, bbox_inches='tight')
    
    if verbose:
        print(f'\nError by Date:')
        print('DATE       : ANG   : AU  ')
        for i, dt_t in enumerate(test_dates):
            print(f'{dt_t.date()} : {ang_err_rms[i]:5.3f} : {pos_err_rms[i]:5.3e}')
    
    # Compute average error
    mean_ang_err = np.mean(ang_err_rms)
    mean_pos_err = np.mean(pos_err_rms)
    print(f'\nMean RMS error in {sim_name} integration of {test_name} test objects:')
    print(f'AU   : {mean_pos_err:5.3e}')
    print(f'angle: {mean_ang_err:5.3f}')
    
    # Return summary of errors in position and angles
    return pos_err_rms, ang_err_rms
