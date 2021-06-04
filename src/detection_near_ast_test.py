"""
Test calculation of asteroid detections near to known asteroid trajectories.
Example call:
$ python detection_near_ast_test.py

Functions in this module:
direction_diff(u0, u1)
report_direction_diff(u_pre, u_vec, u_elt)
main()

Michael S. Emanuel
2021-06-02
"""

# Core
import numpy as np
import pandas as pd

# Local imports
from detection_near_ast import get_data_detections, asteroid_batch_prelim, \
    asteroid_batch_vec, asteroid_batch_elt, light_time_adj, cols_u_ast, arcmin_margin
from astro_utils import dist2sec
from utils import print_stars

# ********************************************************************************************************************* 
def direction_diff(u0: np.ndarray, u1: np.ndarray):
    """Compute difference in directions in arc seconds"""
    du = u1 - u0
    delta = np.sqrt(np.sum(np.square(du), axis=-1))
    # Change in arc seconds
    delta_sec = dist2sec(delta)
    return delta_sec
    
# ********************************************************************************************************************* 
def report_direction_diff(u_pre, u_vec, u_elt: pd.DataFrame):
    """
    Calculate and report difference between three methods of computing directions 
    to a known asteroid at a detection time
    """
    delta_pre_vec = np.max(direction_diff(u0=u_pre, u1=u_vec))
    delta_pre_elt = np.max(direction_diff(u0=u_pre, u1=u_elt))
    delta_vec_elt = np.max(direction_diff(u0=u_vec, u1=u_elt))

    # Report results
    print(f'Max change in direction between different methods:')
    print(f'pre vs. vec: {delta_pre_vec:5.2f} arc seconds')
    print(f'pre vs. elt: {delta_pre_elt:5.2f} arc seconds')
    print(f'vec vs. elt: {delta_vec_elt:5.2e} arc seconds')    

# ********************************************************************************************************************* 
def test_direction():
    """Test direction calculation on a small batch of asteroids"""
    # Inputs
    # d0 = 1
    # d1 = 1000
    d0 = 37501
    d1 = 38501
    n0 = 1
    n1 = 1000
    deg_max = 15.0
    arcmin_max = deg_max * 60.0

    # Status
    print_stars()
    print('Test of direction calculation with three methods:')
    print(f'(1) spline direction from KS.AsteroidDirections')
    print(f'(2) spline asteroid position from KS.AsteroidVectors; iteratively solve for light time')
    print(f'(3) spline asteroid position from KS.AsteroidElements; iteratively solve for light time')
    print(f'Inputs used:')
    print(f'd0: {d0}')
    print(f'd1: {d1}')
    print(f'n0: {n0}')
    print(f'n1: {n1}')
    print(f'deg_max: {deg_max}')

    # Get detections
    print(f'\nLoading detections with DetectionID between {d0} and {d1}...')
    det = get_data_detections(d0=d0, d1=d1)

    # First pass at detections near asteroids
    print(f'Calculating preliminary asteroid directions with AsteroidID between {n0} and {n1}...')
    dna_pre = asteroid_batch_prelim(det=det, n0=n0, n1=n1, arcmin_max=arcmin_max)
    row_count = dna_pre.shape[0]
    print(f'Found {row_count} rows in preliminary search with angular distance < {deg_max} degrees.')

    # Sharpen estimate using splined asteroid vectors
    print(f'\nRefining search using splined asteroid vectors...')
    dna_vec = dna_pre.copy()
    asteroid_batch_vec(dna=dna_vec, n0=n0, n1=n1, arcmin_max=arcmin_max)    
    # Report light time error with vectors
    _  = light_time_adj(df=dna_vec, verbose=True)

    # Sharpen estimate using splined asteroid elements
    print(f'\nRefining search using splined orbital elements...')
    dna_elt = dna_pre.copy()
    asteroid_batch_elt(dna=dna_elt, n0=n0, n1=n1, arcmin_max=arcmin_max)
    # Report light time error with vectors
    _  = light_time_adj(df=dna_vec, verbose=True)

    # Direction from three methods
    u_pre = dna_pre[cols_u_ast].values
    u_vec = dna_vec[cols_u_ast].values
    u_elt = dna_elt[cols_u_ast].values
    print()
    report_direction_diff(u_pre=u_pre, u_vec=u_vec, u_elt=u_elt)

    # Condition for passing
    thresh_pre_vec: float = arcmin_margin * 60.0
    thresh_vec_elt: float = 0.1
    delta_pre_vec = np.max(direction_diff(u0=u_pre, u1=u_vec))
    delta_vec_elt = np.max(direction_diff(u0=u_vec, u1=u_elt))
    is_ok: bool = (delta_pre_vec < thresh_pre_vec) and (delta_vec_elt < thresh_vec_elt)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')

# ********************************************************************************************************************* 
def main():
    """Run all tests in this module"""
    print('Testing detection_near_ast.')
    test_direction()

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
