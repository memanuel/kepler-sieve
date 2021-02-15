"""
Synchronize the orbital elements for asteroids quoted by JPL.
Data is loaded in DB table JPL.AsteroidElement.
Saved calculations are synchronized to the desired date (MJD 59000) 
and saved in DB table KS.AsteroidElement_Ref

Example call:
$ python asteroid_elements_synch.py
$ python asteroid_elements_synch.py --epoch 59000

Michael S. Emanuel
2021-02-10
"""

# Core
import numpy as np
import pandas as pd

# Astronomy
import rebound

# Utility
import argparse
import sys
from tqdm.auto import tqdm as tqdm_auto

# Local imports
from utils import print_stars
from db_utils import sp2df, df2db
from rebound_sim import make_sim_planets
from asteroid_element import get_ast_ref_elts_jpl, add_asteroid_elts, update_asteroid_elements

# Radians in a circle
tau = 2.0 * np.pi

# ********************************************************************************************************************* 
def process_dates(epoch, elt_dates, max_dates: int = None, progbar: bool=False):
    """Process all the orbital elements to the reference epoch"""
    # Save the original input argument to epoch_out
    epoch_out: int = epoch
    # The first date to process
    epoch = np.int32(elt_dates.epoch[0])
    # Get structure of the elements DataFrame, but no rows
    elts = get_ast_ref_elts_jpl(epoch=0)

    # Dates to process
    epochs: np.array = elt_dates.epoch.values
    N_date: int = epochs.shape[0]

    # Set up tqdm iterator for the dates
    i_max: int = N_date-1 if max_dates is None else min(max_dates, N_date-1)
    ii = list(range(i_max))
    if progbar:
        ii = tqdm_auto(ii)

    # Iterate through the dates
    for i in ii:
        # Two epochs: current one (with quoted elements) and next one
        epoch: int = epochs[i]
        epoch_next: int = epochs[i+1]

        # Fetch the new orbital elements at this epoch
        elts_new = get_ast_ref_elts_jpl(epoch=epoch)
        # Append the new elements to the main elts frame
        elts = pd.concat([elts, elts_new]).reset_index(drop=True)

        # Initialize a planets simulation as of epoch
        sim = make_sim_planets(epoch=epoch, load_file=True)
        sim.t = np.float64(epoch)
        sim.N_active = sim.N
        # Add asteroids previously processed with their current orbital elements
        asteroid_ids = add_asteroid_elts(sim=sim, elts=elts)

        # Integrate to the next epoch
        sim.integrate(epoch_next, exact_finish_time=1)
        # Update the asteroid elements as of this epoch
        update_asteroid_elements(sim=sim, elts=elts, epoch=epoch)

    # Integrate one last time if necessary
    if sim.epoch != epoch_out:
        sim.integrate(epoch_out, exact_finish_time=1)
        update_asteroid_elements(sim=sim, elts=elts, epoch=epoch_out)

    return sim, elts

# ********************************************************************************************************************* 
def calc_ref_elt(sim: rebound.Simulation, progbar: bool=False):
    """
    Caclulate reference orbital elements for insertion into KS.AsteroidElement_Ref
    INPUTS:
        sim:      The simulation the asteroid will be added to.
    OUTPUTS:
        elts:   New DataFrame with all required elements for inserting into DB Table KS.AsteroidElement_Ref
    """
    # Number of asteroids and active particles
    N: int = sim.N - sim.N_active
    N_active: int = sim.N_active

    # The epoch and TimeID as arrays of length N
    epoch = np.repeat(sim.t, N)
    TimeID = np.int32(np.rint(epoch * 24 * 60))

    # Grab the AsteroidID array from the simulation object
    AsteroidID = sim.asteroid_ids

    # Empty vectors for remaining columns
    # AsteroidID = np.zeros(N)
    a = np.zeros(N)
    e = np.zeros(N)
    inc = np.zeros(N)
    Omega_node = np.zeros(N)
    omega_peri = np.zeros(N)
    f = np.zeros(N)
    M = np.zeros(N)
    d = np.zeros(N)
    v = np.zeros(N)
    h = np.zeros(N)
    period = np.zeros(N)
    mean_motion = np.zeros(N)
    T_peri = np.zeros(N)
    pomega = np.zeros(N)

    # Calculate all orbits with the Sun as primary
    sim.calculate_orbits(primary=sim.particles[0])

    # Set up asteroid iterator
    ii = list(range(N))
    if progbar:
        ii = tqdm_auto(ii)

    # Loop through asteroids in the simulation
    for i in ii:
        # The index j in the simulation is offset by number of active particles
        j: int = N_active + i
        # The active particle
        p = sim.particles[j]
        # Save all the columns out from the particle to the arrays
        # Save angles in the standard interval [0, 2 pi)
        # AsteroidID[i] = np.int32(p.hash)
        a[i] = p.a
        e[i] = p.e
        inc[i] = p.inc
        Omega_node[i] = p.Omega % tau
        omega_peri[i] = p.omega % tau
        f[i] = p.f % tau
        M[i] = p.M % tau
        d[i] = p.d 
        v[i] = p.v
        h[i] = p.h
        period[i] = p.P
        mean_motion[i] = p.n
        T_peri[i] = p.T
        pomega[i] = p.pomega % tau

    # Wrap the arrays into a Python dictionary
    data_tbl = {
        'AsteroidID': AsteroidID,
        'TimeID': TimeID,
        'epoch': epoch,
        'a': a,
        'e': e,
        'inc': inc,
        'Omega_node': Omega_node,
        'omega_peri': omega_peri,
        'f': f,
        'M': M,
        'd': d,
        'v': v,
        'h': h,
        'period': period,
        'mean_motion': mean_motion,
        'T_peri': T_peri,
        'pomega': pomega,
    }

    # Convert the Python dict into a DataFrame and return it
    ref_elts = pd.DataFrame(data_tbl)
    return ref_elts

# ********************************************************************************************************************* 
def main():
    """Synchronize all orbital elements"""
    # Process command line arguments
    parser = argparse.ArgumentParser(description='Synchronize orbital elements in JPL.AsteroidElement '
            'to desired epoch and save output into KS.AsteroidElement_Ref.')
    parser.add_argument('--epoch', nargs='?', metavar='EP', type=int, default=59000,
                        help='epoch to which quoted orbital elements are synchronized.')
    parser.add_argument('--batch_size', nargs='?', metavar='BS', type=int, default=256,
                        help='the number of dates to process in each batch')
    args = parser.parse_args()

    # Unpack command line arguments
    epoch: int = args.epoch
    batch_size: int = args.batch_size

    # Count the number of dates
    elt_dates = sp2df(sp_name='JPL.GetAsteroidRefElementDates', params={'epoch':epoch})
    date_count: int = elt_dates.shape[0]
    ast_count: int = elt_dates.AsteroidCount.sum()
    print(f'Found {date_count} dates to process with {ast_count} asteroids.')

    # Set up batches    
    batch_count: int = int(np.ceil(date_count / batch_size))
    print(f'Running {batch_count} batches of size {batch_size} dates...')

    for i in range(batch_count):
        # Status
        print_stars()
        print(f'Starting batch {i}.')
        # Synchronize the elements
        print(f'Synchronizing orbital elements to epoch {epoch}...')
        sim, elts = process_dates(epoch=epoch, elt_dates=elt_dates, max_dates=batch_size, progbar=True)
        # Extract reference elements DataFrame
        print(f'Extracting reference orbital elements for DB insertion...')
        ref_elts = calc_ref_elt(sim, progbar=True)
        # Insert into KS.AsteroidElement_Ref
        df2db(df=ref_elts, schema='KS', table='AsteroidElement_Ref', single_thread=True, verbose=True)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
