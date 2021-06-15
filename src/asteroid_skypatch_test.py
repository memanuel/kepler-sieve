"""
Calculate sky patch occupied by known asteroids over time and populate DB table KS.AsteroidSkyPatch.
Example calls:
$ python asteroid_skypatch.py 0 1000

Functions in this module:

main()

Michael S. Emanuel
2021-06-15
"""

# Core
import numpy as np

# Local
from asteroid_skypatch import calc_ast_skypatch
from asteroid_spline import make_spline_ast_dir
from sky_patch import SkyPatchID2dir
from astro_utils import dist2deg
from utils import print_stars

# ********************************************************************************************************************* 
# SkyPatch grid size
N_sp = 1024

# Minutes in one day
mpd: int = 1440

# Threshold
thresh_du: float = 120.0

# ********************************************************************************************************************* 
def test_skypatch_dir(name1: str, name2: str, u1: np.ndarray, u2: np.ndarray, verbose: bool=False) -> float:
    """
    Report 
    INPUTS:
        name1: Descriptive name of the first source, e.g. 'spline'
        name2: Descriptive name of the second source, e.g. 'skypatch'
        u1:    Array of directions from source 1; shape Nx3
        u2:    Array of directions from source 2; shape Nx3
        verbose: Whether to report results to console
    """
    # Difference in unit directions
    du = u2 - u1
    du_norm = np.sqrt(np.sum(np.square(du), axis=-1))
    du_deg = dist2deg(du_norm)
    du_sec = du_deg * 3600.0

    # Calculate mean, median and max difference in degrees
    du_mean = np.mean(du_sec)
    du_median = np.median(du_sec)
    du_max = np.max(du_sec)

    if verbose:
        print(f'Angle Difference: {name2} vs. {name1} in arc seconds')
        print(f'*Mean  : {du_mean:9.1f}*')
        print(f' Median: {du_median:9.1f}')
        print(f' Max   : {du_max:9.1f}')

    # Return the difference of direction vectors in seconds of arc
    return du_mean
    
# ********************************************************************************************************************* 
def test_skypatch_spline(n0: int, n1: int, mjd0: int, mjd1: int, interval_min: int):
    """
    Compare the direction splined at the midpoint of each segment to the center of the SkyPatch in that segment.
    INPUTS:
        n0:     First asteroid to test.
        n1:     Last asteroid to test.
        mjd0:   First date to test.
        mjd1:   Last date to test.
    """
    # Build DataFrame with skypatches
    df = calc_ast_skypatch(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1, interval_min=interval_min)

    # Time at the midpoint of each interval
    t_obs_0 = df.TimeID_0.values / mpd
    t_obs_1 = df.TimeID_1.values / mpd
    t_obs = (t_obs_0 + t_obs_1)/2.0

    # Extract arrays for asteroid_id and sky_patch_id from DataFrame
    asteroid_id = df.AsteroidID.values
    sky_patch_id = df.SkyPatchID.values

    # Spline for asteroid directions vs. time
    spline_u = make_spline_ast_dir(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Evaluate splined asteroid directions
    u_ast, _ = spline_u(ts=t_obs, asteroid_id=asteroid_id)

    # Direction at the center of each sky patch
    u_sky = SkyPatchID2dir(SkyPatchID=sky_patch_id, N=N_sp)

    # Compare splined asteroid direction to center of skypatch
    print()
    print_stars()
    print('Test skypatch segments by comparing two directions:')
    print('(1) spline  : Take the midpoint time of each segment, and spline the direction of the asteroid '
          'using output of make_spline_ast_dir().')
    print('(2) skypatch: Take the center of the reported skypatch in that segment.\n')
    du = test_skypatch_dir(name1='spline', name2='skypatch', u1=u_ast, u2=u_sky, verbose=True)

    # Test result
    is_ok: bool = (du < thresh_du)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# ********************************************************************************************************************* 
def main():

    # Block of asteroids to test
    n0: int = 1
    n1: int = 100

    # Time range for asteroid directions
    mjd0: int = 58000
    mjd1: int = 60000
    interval_min: int = 15

    # Keep track of test results
    is_ok: bool = True

    # Test center of sky patch vs. splined direction
    is_ok &= test_skypatch_spline(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1, interval_min=interval_min)

    # Overall test result
    msg: str = 'PASS' if is_ok else 'FAIL'
    print()
    print_stars()
    print('Overall test of asteroid skypatch calculations:')
    print(f'**** {msg} ****')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
