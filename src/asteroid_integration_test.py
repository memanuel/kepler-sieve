"""
Test integrated path of asteroids vs. Horizons results in DB

Example call:
$ python planets_test.py

Michael S. Emanuel
2021-02-11
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
def get_diff_by_date():
    """
    Get difference in integration of asteroids between JPL and MSE from database. 
    Returns df: Pandas DataFrame with MJD, dq, dq_rel, match_count
    """
    # Wrap arguments to stored procedure
    sp_name: str = 'KS.GetAsteroidIntegrationDiffByDate'
    # Run SQL and return as a DataFrame
    df = sp2df(sp_name=sp_name)
    return df

# ********************************************************************************************************************* 
def get_diff_by_ast(mjd0: int=48000, mjd1: int = 63000):
    """
    Get difference in integration of asteroids between JPL and MSE from database. 
    Returns df: Pandas DataFrame with MJD, dq, dq_rel, match_count
    """
    # Wrap arguments to GetAsteroidIntegrationDiffstored procedure
    sp_name: str = 'KS.GetAsteroidIntegrationDiffByAst'
    params = {
        'mjd0': mjd0,
        'mjd1': mjd1
    }
    # Run SQL and return as a DataFrame
    df = sp2df(sp_name=sp_name, params=params)
    return df

# ********************************************************************************************************************* 
def report_error(df: pd.DataFrame):
    """Report errors for one integration"""
    # Mean error on all asteroids that matched
    mean_dq = np.mean(df['dq'])
    mean_dq_rel = np.mean(df['dq_rel'])
    # Number of matching asteroids
    match_count = np.rint(np.mean(df['match_count']))

    # Report results
    print()
    print_stars()
    print(f'Mean Error on {match_count} Matched Asteroids:')
    print(f'Absolute Error: {mean_dq:5.3e}')
    print(f'Relative Error: {mean_dq_rel:5.3e}')

    return df

# ********************************************************************************************************************* 
def plot_errors_by_date(df: pd.DataFrame, window: int):
    """
    Generate plot of relative errors.
    INPUTS:
        df:     DataFrame of integration errors for Asteroids grouped by date.
    OUTPUTS:
        Saves 2 figures to /figs/integration_test/asteroids
    """

    # Array of dates from MJDs
    dt = np.array([mjd_to_date(mjd) for mjd in df.MJD])

    # Set window for rolling average
    half_window: int = window // 2
    window = 2*half_window
    n = dfd_p.shape[0]

    # Generate moving average error
    plot_x = dt[half_window:n-half_window]
    plot_dq = df.dq.rolling(window).mean()[window:]
    plot_dq_rel = df.dq_rel.rolling(window).mean()[window:]

    # Plot absolute position error
    fig, ax = plt.subplots(figsize=[16,10])
    ax.set_title('Absolute Position Error on Asteroids')
    ax.set_ylabel('Absolute Position Error $|\Delta q|}$')
    ax.set_xlabel(f'End Date ({window} day moving average)')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0,))
    ax.plot(plot_x, plot_dq_rel, label='Asteroids', color='red')
    ax.grid()
    ax.legend()
    fig.savefig(fname='../figs/integration_test/asteroids/pos_error.png', bbox_inches='tight')    

    # Plot relative position error
    fig, ax = plt.subplots(figsize=[16,10])
    ax.set_title('Relative Position Error on Asteroids')
    ax.set_ylabel('Relative Position Error $\\frac{|\Delta q|}{a}$')
    ax.set_xlabel(f'End Date ({window} day moving average)')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0,))
    ax.plot(plot_x, plot_dq_rel, label='Asteroids', color='blue')
    ax.grid()
    ax.legend()
    fig.savefig(fname='../figs/integration_test/asteroids/pos_error_rel.png', bbox_inches='tight')    

# ********************************************************************************************************************* 
def main():
    """Calculate, report and plot error on asteroids integration"""

    # Get error
    df_dt = get_diff_by_date()
    df_ast = get_diff_by_ast()

    # Report error
    report_error(df=df_dt)

    # Generate plots
    plot_errors_by_date(df=df_dt, window=15)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
