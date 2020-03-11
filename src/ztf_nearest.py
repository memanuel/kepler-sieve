"""
Harvard IACS Masters Thesis
ztf_nearest.py
Search ZTF DataFrame against calculated asteroid orbits and implied direction from Palomar.
Find the nearest asteroid to each ZTF observation.

Michael S. Emanuel
11-Mar-2020
"""

# Libary imports
import numpy as np
import pandas as pd
import argparse

# MSE imports
from ztf_data import load_ztf_det_all, ztf_nearest_ast


# ********************************************************************************************************************* 
def main():
    """Main routine for integrating the orbits of known asteroids"""
    # Process command line arguments
    parser = argparse.ArgumentParser(description='Find the nearest asteroid to ZTF data.')
    parser.add_argument('n0', nargs='?', metavar='n0', type=int, default=0,
                        help='the first asteroid number to process')
    parser.add_argument('n_ast', nargs='?', metavar='B', type=int, default=1000,
                        help='the number of asteroids to process in this batch')
    parser.add_argument('--progress', default=False, action='store_true',
                        help='display progress bar')
    # parser.add_argument('--test', default=False, action='store_true',
    #                    help='run in test mode')
    args = parser.parse_args()
    
    # # If run in test mode, run tests without processing any asteroid trajectories
    # if args.test:
    #     # Test  that initial orbital elements recovered from the JPL file
    #     print_header(f'Testing recovery of initial orbital elements with JPL text file vs. Horizons')
    #     test_element_recovery(verbose=True)

    #     # Test the integration vs. Horizons
    #     print_header(f'Testing asteroid integration vs. Horizons')
    #     test_asteroid_sim(verbose=True, make_plot=True)
        
    #     # Test numpy arrays
    #     print_header(f'Testing Numpy array vs. simulation archive:')
    #     test_numpy(verbose=True)
        
    #     # Quit early in test mode: don't want to do any integrations
    #     print()
    #     exit()

    # Unpack command line arguments
    n0: int = args.n0
    n1: int = n0 + args.n_ast
    progbar: bool = args.progress

    # File name for 

    # Load ZTF detections as a DataFrame (no data about nearest asteroids yet)
    ztf, mjd_unq = load_ztf_det_all()

    # Date range in ZTF data
    # mjd_min = np.min(mjd_unq)
    # mjd_max = np.max(mjd_unq)
    # dt_min = mjd_to_date(mjd_min)
    # dt_max = mjd_to_date(mjd_max)

    # Observatory site for ZTF data
    site_name = 'palomar'

    # Build splined positions and observations against unique observation times
    ast_pos, earth_pos, ast_dir = spline_ast_vec_dir(n0=n0, n1=n1, mjd=mjd_unq, site_name=site_name)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
