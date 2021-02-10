"""
Synchronize the orbital elements for asteroids quoted by JPL.
Data is loaded in DB table JPL.AsteroidElement.
Saved calculations are synchronized to the desired date (MJD 59000) 
and saved in DB table KS.AsteroidElement_Ref

Example call:
$ python planets.py --epoch 59000

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

# Local imports
from db_utils import sp2df
from rebound_sim import make_sim_planets
from asteroid_element import get_ast_ref_elts_jpl, add_asteroid_elts, update_asteroid_elements

# ********************************************************************************************************************* 
def process_dates(epoch, max_dates: int = None):
    """Process all the orbital elements to the reference epoch"""
    # Save the original input argument to epoch_out
    epoch_out: int = epoch
    # All the reference dates with orbital elements that must be brought forward
    elt_dates = sp2df(sp_name='KS.GetAsteroidRefElementDates', params={'epoch':epoch_out})
    # The first date to process
    epoch = np.int32(elt_dates.epoch[0])
    # Get structure of the elements DataFrame, but no rows
    elts = get_ast_ref_elts_jpl(epoch=0)

    # Number of dates to process
    epochs: np.array = elt_dates.epoch.values
    N_date: int = epochs.shape[0]

    # Iterate through the remaining dates
    i_max: int = N_date if max_dates is None else max_dates
    for i in range(i_max):
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
def main():
    """Synchronize all orbital elements"""
    # Process command line arguments
    parser = argparse.ArgumentParser(description='Synchronize orbital elements in JPL.AsteroidElement '
            'to desired epoch and save output into KS.AsteroidElement_Ref.')
    parser.add_argument('--epoch', nargs='?', metavar='EP', type=int, default=59000,
                        help='epoch to which quoted orbital elements are synchronized.')
    args = parser.parse_args()

    # Unpack command line arguments
    epoch: int = args.epoch

    print(f'Synchronizing orbital elements to epoch {epoch}.')
    sim, elts = process_dates(epoch=epoch, max_dates=3)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
