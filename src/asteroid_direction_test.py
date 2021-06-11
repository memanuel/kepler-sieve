"""
Test direction from known asteroids to Earth center implied by integrated trajectories.
Example calls:
$ python asteroid_direction.py 0 1000

Functions in this module:
main()

Michael S. Emanuel
2021-06-04
"""

# Core
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

# Astronomy
from astropy.units import deg

# Local imports
from asteroid_spline import make_spline_df
from asteroid_direction import calc_distance, calc_dir_linear, calc_dir_linear_topos, calc_dir_spline, \
    calc_dir_ast2obs_df
from ra_dec import radec2dir
from db_utils import sp2df, df2db
from astro_utils import dist2sec
from ra_dec_test import direction_diff
from utils import print_stars

# ********************************************************************************************************************* 
# Column names for astrometric and apparent direction
cols_u_ast = ['ux_ast', 'uy_ast', 'uz_ast']
cols_u_app = ['ux_app', 'uy_app', 'uz_app']

# Column names for testing 
cols_q_obs = ['qObs_x', 'qObs_y', 'qObs_z']
cols_q_ast = ['qAst_x', 'qAst_y', 'qAst_z']
cols_v_ast = ['vAst_x', 'vAst_y', 'vAst_z']
cols_u = ['ux', 'uy', 'uz']

# ********************************************************************************************************************* 
def jpl_ast_dir_populate():
    """
    Populate DB table JPL.AsteroidDirection from JPL.AsteroidDirectionImport
    """
    # Get the imported RA/DEC from JPL
    sp_name = 'JPL.GetAsteroidDirection_Import'
    params = {
        'n0': 0,
        'n1': 1000000,
    }
    df = sp2df(sp_name=sp_name, params=params)

    # The observation time is just a number array
    obstime_mjd = df.mjd.values
    # Convert numpy arrays to angles in degrees - astrometric direction
    ra_ast = df.RA_ast.values * deg
    dec_ast = df.DEC_ast.values * deg
    # Convert numpy arrays to angles in degrees - apparentdirection
    ra_app = df.RA_app.values * deg
    dec_app = df.DEC_app.values * deg

    # Convert RA/DEC to a direction for both astrometric and apparent
    u_ast = radec2dir(ra=ra_ast, dec=dec_ast, obstime_mjd=obstime_mjd)
    u_app = radec2dir(ra=ra_app, dec=dec_app, obstime_mjd=obstime_mjd)

    # Save directions to DataFrame
    df[cols_u_ast] = u_ast
    df[cols_u_app] = u_app

    # Write asteroid directions to new DB table JPL.AsteroidDirections
    columns = ['AsteroidID', 'ObservatoryID', 'TimeID', 'mjd', 
               'RA_ast', 'DEC_ast', 'ux_ast', 'uy_ast', 'uz_ast', 
               'RA_app', 'DEC_app', 'ux_app', 'uy_app', 'uz_app', 
               'LightTime', 'r', 'rDot', 'delta', 'deltaDot']
    df2db(df=df, schema='JPL', table='AsteroidDirection', columns=columns)

# ********************************************************************************************************************* 
def test_dir_vectors(mjd0: int, mjd1: int):
    """
    Test for consistency between JPL and MSE integrations in the implied direction from earth geocenter to asteroid.
    INPUTS:
        mjd0: The first date to include in the test
        mjd1: The last date to include in the test
    """
    # Directions from two different sources
    df_jpl = sp2df(sp_name='JPL.AstrometricDirectionTest')
    df_mse = sp2df(sp_name='JPL.AstrometricDirectionAndVectorTest')

    # Mask for selected date range
    mask_range_jpl = (mjd0 <= df_jpl.tObs) & (df_jpl.tObs <= mjd1)
    mask_range_mse = (mjd0 <= df_mse.tObs) & (df_mse.tObs <= mjd1)

    # Mask for overlapping dates (multiples of 20)
    interval = 20*1440
    mask_overlap_jpl = (df_jpl.TimeID % interval) == 0
    mask_overlap_mse = (df_mse.TimeID % interval) == 0

    # Combined mask satisfies both requirements
    mask_jpl = mask_range_jpl & mask_overlap_jpl
    mask_mse = mask_range_mse & mask_overlap_mse

    # Mask both frames down to selected rows
    df_jpl = df_jpl[mask_jpl].reindex()
    df_mse = df_mse[mask_mse].reindex()    

    # State vectors according to JPL integration
    q_ast_jpl = df_jpl[cols_q_ast].values
    v_ast_jpl = df_jpl[cols_v_ast].values

    # State vectors according to MSE integration
    q_ast_mse = df_mse[cols_q_ast].values
    v_ast_mse = df_mse[cols_v_ast].values

    # The observatory position according to JPL and MSE
    q_obs_jpl = df_jpl[cols_q_obs].values
    q_obs_mse = df_mse[cols_q_obs].values

    # Directions calculated in linear model with both state vectors
    u_jpl, delta_jpl = calc_dir_linear(q_tgt=q_ast_jpl, v_tgt=v_ast_jpl, q_obs=q_obs_jpl)
    u_mse, delta_mse = calc_dir_linear(q_tgt=q_ast_mse, v_tgt=v_ast_mse, q_obs=q_obs_mse)

    # Compare the two directions
    print('Compare directions from Geocenter the first 16 asteroids between two sources of state vectors:')
    print(f'Data compared between dates {mjd0} and {mjd1}.')
    print('JPL: State vectors sourced from JPL using Horizons.')
    print('MSE: Integration done in rebound.')
    print('Both sources have directions calculated using calc_dir_linear().')
    err = direction_diff(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, verbose=True)
    # Test result
    thresh: float = 0.5
    is_ok: bool = (err < thresh)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# ********************************************************************************************************************* 
def test_dir_linear(state_vec_src: str, mjd0: int, mjd1: int):
    """
    Test MSE astrometric direction calculation vs. JPL using the linear model.
    INPUTS:
        state_vec_src: One of 'JPL' or 'MSE' - the source of the state vectors used
        mjd0: The first date to include in the test
        mjd1: The last date to include in the test
    """
    # Check state_vec_src
    if state_vec_src not in ('JPL', 'MSE'):
        raise ValueError("Bad state_vec_src! Must be one of 'JPL' or 'MSE'.")

    # Table of SP name keyed by state_vec_src
    sp_name_tbl = {
        'JPL': 'JPL.AstrometricDirectionTest',
        'MSE': 'JPL.AstrometricDirectionAndVectorTest' 
    }

    # Import the JPL asteroid directions
    sp_name = sp_name_tbl[state_vec_src]
    df = sp2df(sp_name=sp_name)

    # Mask data down to the selected date range
    mask =  (mjd0 <= df.tObs) & (df.tObs <= mjd1)
    df = df[mask].reindex()

    # Position and velocity according to selected source of state vectors
    q_obs = df[cols_q_obs].values
    q_ast = df[cols_q_ast].values
    v_ast = df[cols_v_ast].values

    # Direction according to JPL
    u_jpl = df[cols_u].values

    # MSE calculation of the astrometric direction in linear model
    u_mse, delta_mse = calc_dir_linear(q_tgt=q_ast, v_tgt=v_ast, q_obs=q_obs)

    # Compare the two directions
    print(f'Compare directions from Geocenter to the first 16 asteroids using {state_vec_src} state vectors in linear model:')
    print(f'Data compared between dates {mjd0} and {mjd1}.')
    print('JPL: Directions quoted by JPL.  Astrometric RA/DEC converted to vectors using radec2dir.')
    print('MSE: MSE calculations using linear model with function calc_dir_linear in asteroid_direction.py.\n')
    err = direction_diff(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, verbose=True)

    # Test result
    thresh: float = 1.0
    is_ok: bool = (err < thresh)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# ********************************************************************************************************************* 
def test_dir_linear_topos(site_name: str, mjd0: int, mjd1: int):
    """
    Test for consistency between JPL and MSE integrations in the implied direction from earth geocenter to asteroid.
    INPUTS:
        site_name: String ID for observatory in Skyfield.  One of 'geocenter' or 'palomar'
        mjd0: The first date to include in the test
        mjd1: The last date to include in the test
    """
    # Check site_name
    if site_name not in ('geocenter', 'palomar'):
        raise ValueError("Bad site_name! Must be one of 'geocenter' or 'palomar'.")

    # Table of ObservatoryShortName keyed by site_name
    obs_name_tbl = {
        'geocenter': 'Geocenter',
        'palomar': 'ZTF'
    }
    obs_name: str = obs_name_tbl[site_name]

    # Directions from two different sources
    params_jpl = {
        'ObservatoryShortName': obs_name,
        'n0': 0,
        'n1': 100,
    }
    df_jpl = sp2df(sp_name='JPL.GetAsteroidDirection', params=params_jpl)
    df_mse = sp2df(sp_name='JPL.AstrometricDirectionAndVectorTest')

    # Mask for selected date range
    mask_range_jpl = (mjd0 <= df_jpl.tObs) & (df_jpl.tObs <= mjd1)
    mask_range_mse = (mjd0 <= df_mse.tObs) & (df_mse.tObs <= mjd1)

    # Mask for overlapping dates (multiples of 20)
    interval = 20*1440
    mask_overlap_jpl = (df_jpl.TimeID % interval) == 0
    mask_overlap_mse = (df_mse.TimeID % interval) == 0

    # Combined mask satisfies both requirements
    mask_jpl = mask_range_jpl & mask_overlap_jpl
    mask_mse = mask_range_mse & mask_overlap_mse

    # Mask both frames down to selected rows
    df_jpl = df_jpl[mask_jpl].reindex()
    df_mse = df_mse[mask_mse].reindex()    

    # Directions according to JPL; choose the astrometric direction, not the apparent direction
    u_jpl = df_jpl[cols_u_ast].values

    # State vectors according to MSE integration
    q_earth_mse = df_mse[cols_q_obs].values
    q_ast_mse = df_mse[cols_q_ast].values
    v_ast_mse = df_mse[cols_v_ast].values
    # The observation time array
    obstime_mjd_mse = df_mse.tObs.values

    # Directions calculated in linear model with topos correction with both state vectors
    u_mse, delta_mse = calc_dir_linear_topos(q_tgt=q_ast_mse, v_tgt=v_ast_mse, q_earth=q_earth_mse, 
                                                obstime_mjd=obstime_mjd_mse, site_name=site_name)

    # Compare the two directions
    print('Compare directions from Palomar to the first 16 asteroids between MSE and JPL:')
    print(f'Data compared between dates {mjd0} and {mjd1}.')
    print('JPL: Quoted direction from observatory site to asteroid.')
    print('MSE: State vectors from rebound integration. Direction with topos from calc_dir_linear_topos().')
    err = direction_diff(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, verbose=True)

    # Test result
    thresh: float = 1.0
    is_ok: bool = (err < thresh)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# ********************************************************************************************************************* 
def test_dir_spline(state_vec_src: str, mjd0: int, mjd1: int):
    """
    Test MSE astrometric direction calculation vs. JPL using the iterated spline model.
    INPUTS:
        state_vec_src: One of 'JPL' or 'MSE' - the source of the state vectors used
        mjd0: The first date to include in the test
        mjd1: The last date to include in the test
    """
    # Check state_vec_src
    if state_vec_src not in ('JPL', 'MSE'):
        raise ValueError("Bad state_vec_src! Must be one of 'JPL' or 'MSE'.")

    # Table of SP name keyed by state_vec_src
    sp_name_tbl = {
        'JPL': 'JPL.AstrometricDirectionTest',
        'MSE': 'JPL.AstrometricDirectionAndVectorTest' 
    }

    # Import the JPL asteroid directions
    sp_name = sp_name_tbl[state_vec_src]
    df = sp2df(sp_name=sp_name)

    # Mask data down to the selected date range
    mask =  (mjd0 <= df.tObs) & (df.tObs <= mjd1)
    df = df[mask].reindex()

    # Position according to selected source of state vectors
    q_obs = df[cols_q_obs].values
    # Direction according to JPL
    u_jpl = df[cols_u].values
    # Other arrays used by calc_dir_spline
    asteroid_id = df.AsteroidID.values
    t_obs = df.tObs.values
    light_time = df.LightTime.values

    # Build asteroid position spline from asteroid state vectors
    id_col = 'AsteroidID'
    time_col = 'tAst'
    spline_q_ast = make_spline_df(df=df, cols_spline=cols_q_ast, id_col=id_col, time_col=time_col)

    # MSE calculation of the astrometric direction in linear model
    u_mse, delta_mse = calc_dir_spline(spline_q_ast=spline_q_ast, q_obs=q_obs, t_obs=t_obs, 
                        asteroid_id=asteroid_id, light_time=light_time, iters=2)

    # Compare the two directions
    print(f'Compare directions from Geocenter to the first 16 asteroids using {state_vec_src} state vectors in spline model:')
    print(f'Data compared between dates {mjd0} and {mjd1}.')
    print('JPL: Directions quoted by JPL.  Astrometric RA/DEC converted to vectors using radec2dir.')
    print('MSE: MSE calculations using linear model with function calc_dir_spline in asteroid_direction.py.\n')
    err = direction_diff(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, verbose=True)

    # Test result
    thresh: float = 1.0
    is_ok: bool = (err < thresh)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# ********************************************************************************************************************* 
def main():
    # Build table JPL.AsteroidDirections from JPl.AsteroidDirections_Import
    # jpl_ast_dir_populate()
    # print('Regenerated DB table JPL.AsteroidDirection')

    # Set date range for testing
    epoch: int = 59000
    half_width: int = 2000
    mjd0: int = epoch - half_width
    mjd1: int = epoch + half_width

    # Initialize test result
    is_ok: bool = True

    # Test the MSE integration vs. JPL by comparing state vectors
    print()
    print_stars()
    is_ok &= test_dir_vectors(mjd0=mjd0, mjd1=mjd1)
    
    # Test the asteroid directions in linear model with JPL state vectors
    print()
    print_stars()
    is_ok &= test_dir_linear(state_vec_src='JPL', mjd0=mjd0, mjd1=mjd1)

    # Test the asteroid directions in linear model with MSE state vectors
    print()
    print_stars()
    is_ok &= test_dir_linear(state_vec_src='MSE', mjd0=mjd0, mjd1=mjd1)

    # Test the asteroid directions in linear model with MSE state vectors with topos
    print()
    print_stars()
    is_ok &= test_dir_linear_topos(site_name='palomar', mjd0=mjd0, mjd1=mjd1)
    
    # Test the asteroid directions in splined model with JPL state vectors
    print()
    print_stars()
    is_ok &= test_dir_spline(state_vec_src='JPL', mjd0=mjd0, mjd1=mjd1)

    # Overall test result
    msg: str = 'PASS' if is_ok else 'FAIL'
    print()
    print_stars()
    print('Overall test of asteroid directions:')
    print(f'**** {msg} ****')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
