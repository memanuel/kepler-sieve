"""
Direction from known asteroids to Earth center implied by integrated trajectories.
Example calls:
$ python asteroid_direction.py 0 1000

Functions in this module:
main()

Michael S. Emanuel
2021-06-01
"""

# Core
import numpy as np
import pandas as pd

# Commandline arguments
import argparse

# ********************************************************************************************************************* 
def main():
    """Integrate the orbits of the planets and selected batch of asteroids"""

    # Process command line arguments
    parser = argparse.ArgumentParser(description='Calculated direction from known asteroids to Earth center '
    'implied by rebound integration.  Populates DB table KS.AsteroidDirections.')
    parser.add_argument('n0', nargs='?', metavar='n0', type=int, default=0,
                        help='the first asteroid number to process')
    parser.add_argument('n_ast', nargs='?', metavar='B', type=int, default=1000,
                        help='the number of asteroids to process in this batch'),
    
    # Unpack command line arguments
    args = parser.parse_args()
    
    # Block of asteroids to integrate and epoch
    n0: int = args.n0
    n1: int = n0 + args.n_ast

    # Report arguments
    print(f'Processing asteroid directions for asteroid number {n0} <= AsteroidID < {n1}...')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
