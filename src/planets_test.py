"""
Test integrated path of planets vs. Horizons results in DB

Example calls:
$ python planets_test.py --mjd0 55350 --mjd1 61650
$ python planets_test.py --epoch 59000 --half_width 3650

Michael S. Emanuel
2021-01-26
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
from db_utils import sp2df

# ********************************************************************************************************************* 
# Set plot style
plot_style()

# ********************************************************************************************************************* 
def get_integration_diff(BodyCollectionCD: str, mjd0: int, mjd1: int):
    """
    Get integration error for a collection of bodies over a date range.
    INPUTS:
        BodyCollectionCD: Collection of bodies integrated in Rebound, e.g. 'P' for Planets, 'D' for DE435
        mjd0: First date in range
        mjd1: Last date in range
    OUTPUTS:
        df:   Pandas DataFrame with TimeID, MJD, dq_rel, dv_rel
    """

    # Wrap arguments to GetIntegrationDiffByDate stored procedure
    sp_name = 'KS.GetIntegrationDiffByDate'
    params = {
        'BodyCollectionCD': BodyCollectionCD,
        'mjd0': mjd0,
        'mjd1': mjd1,
    }

    # Run SQL and return as a DataFrame
    df = sp2df(sp_name=sp_name, params=params)
    return df

# ********************************************************************************************************************* 
def main():
    """Integrate the orbits of the planets and major moons"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description='The solar system integrations against Horizons results.')
    parser.add_argument('--epoch', nargs='?', metavar='EP', type=int, default=59000,
                        help='epoch of the base simulation that is integrated forward and backwards, as an MJD')
    parser.add_argument('--half_width', nargs='?', metavar='HW', type=int, default=0,
                        help='half of the width of the interval in days, which is symmetric about the epoch.')
    parser.add_argument('--mjd0', nargs='?', metavar='t0', type=int, default=55350,
                        help='epoch of the first date in the integration, as an MJD.')
    parser.add_argument('--mjd1', nargs='?', metavar='t1', type=int, default=62650,
                        help='epoch of the last date in the integration, as an MJD.')
    args = parser.parse_args()
    
    # Unpack command line arguments
    epoch: int = args.epoch                 # MJD 59000 = 2020-05-31

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

    # Report arguments and integrator settings
    print_stars()
    print(f'epoch          : {epoch} ({epoch_dt})')
    print(f'date range mjd : {mjd0} to {mjd1}')
    print(f'date range     : {mjd0_dt} to {mjd1_dt}')
    print(f'full width     : {width_yrs:3.1f} years')

    # Run error on planets
    print()
    print_stars()
    df = get_integration_diff(BodyCollectionCD='P', mjd0=mjd0, mjd1=mjd1)
    dq_rel = np.mean(df['dq_rel'])
    print('Mean Relative Error - Integration with Planets:')
    print(f'{dq_rel:5.3e}')

    # Run error on DE435
    print()
    print_stars()
    df = get_integration_diff(BodyCollectionCD='D', mjd0=mjd0, mjd1=mjd1)
    dq_rel = np.mean(df['dq_rel'])
    print('Mean Relative Error - Integration with DE435:')
    print(f'{dq_rel:5.3e}')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
