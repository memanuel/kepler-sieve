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
def get_integration_diff(BodyCollectionCD: str, mjd0: int, mjd1: int, by_date: bool):
    """
    Get integration error for a collection of bodies over a date range.
    INPUTS:
        BodyCollectionCD: Collection of bodies integrated in Rebound, e.g. 'P' for Planets, 'D' for DE435
        mjd0: First date in range
        mjd1: Last date in range
        by_date: Flag; true to get summary by date, false for detail by planet
    OUTPUTS:
        df:   Pandas DataFrame with TimeID, MJD, dq_rel, dv_rel
    """

    # Wrap arguments to GetIntegrationDiffByDate stored procedure
    suffix: str = 'ByDate' if by_date else ''
    sp_name: str = f'KS.GetIntegrationDiff{suffix}'
    params = {
        'BodyCollectionCD': BodyCollectionCD,
        'mjd0': mjd0,
        'mjd1': mjd1,
    }

    # Run SQL and return as a DataFrame
    df = sp2df(sp_name=sp_name, params=params)
    return df

# ********************************************************************************************************************* 
def plot_errors(df_p: pd.DataFrame, df_d: pd.DataFrame, window: int):
    """
    Generate plot of relative errors.
    INPUTS:
        df_p: DataFrame of integration errors for Planets.
        df_d: DataFrame of integration errors for DE435.
    OUTPUTS:
        Saves a figure to 
    """

    # Array of dates from MJDs
    dt = np.array([mjd_to_date(mjd) for mjd in df_p.MJD])

    # Set window for rolling average
    half_window: int = window // 2
    window = 2*half_window
    n = df_p.shape[0]

    # Generate moving average error
    plot_x = dt[half_window:n-half_window]
    plot_y_p = df_p.dq_rel.rolling(window).mean()[window:]
    plot_y_d = df_d.dq_rel.rolling(window).mean()[window:]

    # Plot position error
    fig, ax = plt.subplots(figsize=[16,10])
    ax.set_title('Relative Position Error on Sun & Planets')
    ax.set_ylabel('Relative Position Error $\\frac{|\Delta q|}{\sigma_{q}}$')
    ax.set_xlabel(f'End Date ({window} day moving average)')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0,))
    ax.plot(plot_x, plot_y_p, label='Planets', color='blue')
    # ax.plot(plot_x, plot_y_d, label='DE435', color='red')
    ax.grid()
    ax.legend()
    fig.savefig(fname='../figs/integration_test/planets/pos_error.png', bbox_inches='tight')    

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
    df_p = get_integration_diff(BodyCollectionCD='P', mjd0=mjd0, mjd1=mjd1, by_date=True)
    mean_err_p = np.mean(df_p['dq_rel'])
    print('Mean Relative Error - Integration with Planets:')
    print(f'{mean_err_p:5.3e}')

    # Run error on DE435
    print()
    print_stars()
    df_d = get_integration_diff(BodyCollectionCD='D', mjd0=mjd0, mjd1=mjd1, by_date=True)
    mean_err_d = np.mean(df_d['dq_rel'])
    print('Mean Relative Error - Integration with DE435:')
    print(f'{mean_err_d:5.3e}')

    # Generate plots
    plot_errors(df_p=df_p, df_d=df_d, window=180)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
