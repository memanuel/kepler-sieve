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
import sys

# Plotting
import matplotlib.pyplot as plt

# MSE imports
from utils import plot_style, print_stars
from astro_utils import mjd_to_date
from rebound_sim import make_sim_planets
from rebound_integrate import integrate_df
from db_utils import sp2df

# ********************************************************************************************************************* 
# Set plot style
plot_style()

# ********************************************************************************************************************* 
def get_integration_diff(body_collection: str, mjd0: int, mjd1: int, by_date: bool):
    """
    Get integration error for a collection of bodies over a date range.
    INPUTS:
        body_collection: Collection of bodies integrated in Rebound, e.g. 'Planets', or 'DE435'
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
        'BodyCollectionName': body_collection,
        'mjd0': mjd0,
        'mjd1': mjd1,
    }

    # Run SQL and return as a DataFrame
    df = sp2df(sp_name=sp_name, params=params)
    return df

# ********************************************************************************************************************* 
def live_integration_diff(df: pd.DataFrame):
    """Test an integration data frame against Horizons"""
    # Test this integration on the last epoch in the DataFrame
    epoch = np.max(df.MJD)

    # Mask the data frame to match the desired epoch
    mask = (df.TimeID == epoch*24*60)
    # Extract this date and the fields to test
    df_k = df[mask][['BodyID', 'qx', 'qy', 'qz', 'vx', 'vy', 'vz']]
    # The position according to the live integration
    q1 = df_k[['qx', 'qy', 'qz']].values

    
    # Look up the configuration on Horizons
    params={'BodyCollectionName': 'Planets', 'epoch': epoch}
    df_h = sp2df('JPL.GetHorizonsStateCollection', params)
    # Filter down to just the relevant columns
    df_h[['BodyID', 'qx', 'qy', 'qz', 'vx', 'vy', 'vz']]
    # The position according to Horizons
    q2 = df_h[['qx', 'qy', 'qz']].values

    # Return the mean relative position error
    dq = np.linalg.norm(q1-q2, axis=1)
    dq_rel = dq / np.linalg.norm(q1, axis=1)
    mean_dq_rel = np.mean(dq_rel)
    return mean_dq_rel

# ********************************************************************************************************************* 
def report_error(body_collection: str, mjd0: int, mjd1: int):
    """Report errors for one integration"""
    # Get DataFrames of errors by date and for each body
    df = get_integration_diff(body_collection=body_collection, mjd0=mjd0, mjd1=mjd1, by_date=False)
    dfd = get_integration_diff(body_collection=body_collection, mjd0=mjd0, mjd1=mjd1, by_date=True)
    
    # Error on all 11 bodies (Sun, 9 Planets, Moon)
    mean_dq = np.mean(dfd['dq'])
    mean_dq_rel = np.mean(dfd['dq_rel'])

    # Error on 11 bodies excluding Mercury
    mask = (df.BodyID != 1)
    dfd_xm = df[mask].groupby(dfd.TimeID).mean()
    mean_dq_xm = np.mean(dfd_xm['dq'])
    mean_dq_rel_xm = np.mean(dfd_xm['dq_rel'])

    # Report results
    print()
    print_stars()
    print(f'Mean Error - Integration with {body_collection}:')
    print('                All Bodies : Ex Mercury')
    print(f'Absolute Error: {mean_dq:5.3e}  : {mean_dq_xm:5.3e}')
    print(f'Relative Error: {mean_dq_rel:5.3e}  : {mean_dq_rel_xm:5.3e}')

    return df, dfd, dfd_xm

# ********************************************************************************************************************* 
def plot_errors(dfd_p: pd.DataFrame, dfd_d: pd.DataFrame, window: int):
    """
    Generate plot of relative errors.
    INPUTS:
        dfd_p: DataFrame of integration errors for Planets grouped by date.
        dfd_d: DataFrame of integration errors for DE435 grouped by date.
    OUTPUTS:
        Saves a figure to 
    """

    # Array of dates from MJDs
    dt = np.array([mjd_to_date(mjd) for mjd in dfd_p['MJD']])

    # Set window for rolling average
    half_window: int = window // 2
    window = 2*half_window
    n = dfd_p.shape[0]

    # Generate moving average error
    plot_x = dt[half_window:n-half_window]
    plot_y_p = dfd_p.dq_rel.rolling(window).mean()[window:]
    # plot_y_d = dfd_d.dq_rel.rolling(window).mean()[window:]

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
    parser.add_argument('--test', dest='test', action='store_const', const=True, default=False,
                        help='Quick test of a live integration, then quit.')
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

    # First check if we are doing a live test
    if args.test:
        print('\nPerforming live integration test on planets:')
        sim_epoch = make_sim_planets(epoch=epoch)
        df: pd.DataFrame = \
            integrate_df(sim_epoch=sim_epoch, mjd0=mjd0, mjd1=mjd1, interval_p=5, interval_q=1, 
                         save_elements=False, progbar=True)
        mean_dq_rel: float = live_integration_diff(df)
        isOK: bool = mean_dq_rel < 1.0E-5
        print(f'\nMean relative position error on epoch {mjd1}: {mean_dq_rel:5.3e}')
        msg = 'PASS' if isOK else 'FAIL'
        print(f'***** {msg} *****')
        sys.exit()

    # Run error on Planets
    df_p, dfd_p, dfd_xm_p = report_error(body_collection='Planets', mjd0=mjd0, mjd1=mjd1)

    # Run error on DE435
    # df_d, dfd_d, dfd_xm_d = report_error(body_collection='DE435', mjd0=mjd0, mjd1=mjd1)

    # Generate plots
    plot_errors(dfd_p=dfd_p, dfd_d=dfd_p, window=180)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
