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
# import numpy as np
# import pandas as pd

# Astronomy
from astropy.units import deg

# Local imports
from asteroid_spline import make_spline_df
from asteroid_direction import calc_dir_linear, calc_dir_linear_topos, calc_dir_spline, calc_dir_ast2obs
from ra_dec import radec2dir
from db_utils import sp2df, df2db
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

# DataFrames of JPL directions with two sources of state vectors
df_jpl = sp2df(sp_name='JPL.AsteroidDirectionTest')
df_mse = sp2df(sp_name='JPL.AsteroidDirectionVectorTest')

# Table of DataFrame keyed by state_vec_src
df_tbl_src = {
    'JPL': df_jpl,
    'MSE': df_mse,
}

# Mask JPL DataFrame to separate Geocenter and ZTF observatories with JPL state vectors
mask_jpl_geo = (df_jpl.ObservatoryID==0)
mask_jpl_ztf = (df_jpl.ObservatoryID==1)
df_jpl_geo = df_jpl[mask_jpl_geo].copy()
df_jpl_ztf = df_jpl[mask_jpl_ztf].copy()

# Mask MSE DataFrame to separate Geocenter and ZTF observatories with MSE state vectors
mask_mse_geo = (df_mse.ObservatoryID==0)
mask_mse_ztf = (df_mse.ObservatoryID==1)
df_mse_geo = df_mse[mask_mse_geo].copy()
df_mse_ztf = df_mse[mask_mse_ztf].copy()

# Table of DataFrame keyed by (state_vec_src, obs_name)
df_tbl = {
    ('JPL', 'geo'): df_jpl_geo,
    ('JPL', 'ztf'): df_jpl_ztf,
    ('MSE', 'geo'): df_mse_geo,
    ('MSE', 'ztf'): df_mse_ztf,
}

# Table of observatory long names
obs_name_tbl = {
    'geo': 'Geocenter', 
    'ztf': 'Palomar',
}

# Range of asteroids to test
n0: int = 1
n1: int = 17

# ********************************************************************************************************************* 
def jpl_ast_dir_populate():
    """
    Populate DB table JPL.AsteroidDirection from JPL.AsteroidDirectionImport
    """
    # Get the imported RA/DEC from JPL
    sp_name = 'JPL.GetAsteroidDirection_Import'
    params = {
        'n0': n1,
        'n1': n1,
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
def test_dir_vectors(obs_name: str, mjd0: int, mjd1: int):
    """
    Test for consistency between JPL and MSE integrations in the implied direction from earth geocenter to asteroid.
    INPUTS:
        obs_name:       One of 'geo' or 'ztf' - the location of the observatory
        mjd0:           The first date to include in the test
        mjd1:           The last date to include in the test
    """
    # Check obs_name
    if obs_name not in ('geo', 'ztf'):
        raise ValueError("Bad obs_name! Must be one of 'geo' or 'ztf'.")
    obs_long_name = obs_name_tbl[obs_name]

    # Directions from two different sources
    df_jpl = df_tbl[('JPL', obs_name)].copy()
    df_mse = df_tbl[('MSE', obs_name)].copy()

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
    print()
    print_stars()
    print(f'Compare directions from {obs_long_name} to the first 16 asteroids between two sources of state vectors:')
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
def test_topos(mjd0: int, mjd1: int):
    """
    Compare offset in directions from when observatory source switches from Geocenter to Palomar.
    INPUTS:
        mjd0:           The first date to include in the test
        mjd1:           The last date to include in the test
    """
    # Get directions from both observatories according to JPL along with JPL state vectors
    df_geo = df_tbl[('JPL', 'geo')].copy()
    df_ztf = df_tbl[('JPL', 'ztf')].copy()

    # Mask data down to the selected date range
    mask_geo =  (mjd0 <= df_geo.tObs) & (df_geo.tObs <= mjd1)
    mask_ztf =  (mjd0 <= df_ztf.tObs) & (df_ztf.tObs <= mjd1)
    df_geo = df_geo[mask_geo].reindex()
    df_ztf = df_ztf[mask_ztf].reindex()

    # State vectors according to JPL
    q_ast_geo = df_geo[cols_q_ast].values
    v_ast_geo = df_geo[cols_v_ast].values
    q_ast_ztf = df_ztf[cols_q_ast].values
    v_ast_ztf = df_ztf[cols_v_ast].values

    # The observatory position according to JPL
    q_obs_geo = df_geo[cols_q_obs].values
    q_obs_ztf = df_ztf[cols_q_obs].values

    # Directions quoted by JPL
    u_geo_jpl = df_geo[cols_u].values
    u_ztf_jpl = df_ztf[cols_u].values

    # Directions calculated in MSE linear model
    u_geo_mse, _ = calc_dir_linear(q_tgt=q_ast_geo, v_tgt=v_ast_geo, q_obs=q_obs_geo)
    u_ztf_mse, _ = calc_dir_linear(q_tgt=q_ast_ztf, v_tgt=v_ast_ztf, q_obs=q_obs_ztf)

    # Change in direction according to JPL and MSE
    du_jpl = u_ztf_jpl - u_geo_jpl
    du_mse = u_ztf_mse - u_geo_mse

    # Compare the two directions
    print()
    print_stars()
    print(f'Compare directions change between Geocenter and Palomar on the first 16 asteroids:')
    print('JPL: Quoted direction ZTF->Asteroid minus GEO->Asteroid.')
    print('MSE: Calculated direction ZTF->Asteroid minus GEO->Asteroid using JPL state vectors.')
    print('Both sources have directions calculated using calc_dir_linear().')
    err = direction_diff(name1='JPL', name2='MSE', u1=du_jpl, u2=du_mse, verbose=True)
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
        state_vec_src:  One of 'JPL' or 'MSE' - the source of the state vectors used
        mjd0:           The first date to include in the test
        mjd1:           The last date to include in the test
    """
    # Check state_vec_src
    if state_vec_src not in ('JPL', 'MSE'):
        raise ValueError("Bad state_vec_src! Must be one of 'JPL' or 'MSE'.")    

    # Import the JPL asteroid directions and state vectors from the selected source
    df = df_tbl_src[state_vec_src].copy()

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
    print()
    print_stars()
    print(f'Compare directions to the first 16 asteroids using {state_vec_src} state vectors in linear model:')
    print(f'Data compared between dates {mjd0} and {mjd1}.')
    print('JPL: Directions quoted by JPL.  Astrometric RA/DEC converted to vectors using radec2dir.')
    print('MSE: MSE calculations using linear model with function calc_dir_linear in asteroid_direction.py.')
    print('Both models tested combining JPL data for Geocenter and Palomar observatory locations.')
    err = direction_diff(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, verbose=True)

    # Test result
    thresh: float = 1.0
    is_ok: bool = (err < thresh)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# # ********************************************************************************************************************* 
# def test_dir_spline(state_vec_src: str, mjd0: int, mjd1: int):
#     """
#     Test MSE astrometric direction calculation vs. JPL using the iterated spline model and JPL state vectors.
#     INPUTS:
#         state_vec_src:  One of 'JPL' or 'MSE' - the source of the state vectors used
#         mjd0: The first date to include in the test
#         mjd1: The last date to include in the test
#     """
#     # Check state_vec_src
#     if state_vec_src not in ('JPL', 'MSE'):
#         raise ValueError("Bad state_vec_src! Must be one of 'JPL' or 'MSE'.")    

#     # Import the JPL asteroid directions and state vectors from the selected source
#     df = df_tbl_src[state_vec_src].copy()

#     # Mask data down to the selected date range
#     mask =  (mjd0 <= df.tObs) & (df.tObs <= mjd1)
#     df = df[mask].reindex()

#     # Position and velocity according to selected source of state vectors
#     q_obs = df[cols_q_obs].values
#     q_ast = df[cols_q_ast].values
#     v_ast = df[cols_v_ast].values

#     # Direction according to JPL
#     u_jpl = df[cols_u].values

#     # Build asteroid position spline from asteroid state vectors
#     id_col = 'AsteroidID'
#     time_col = 'tObs'
#     spline_q_ast = make_spline_df(df=df, cols_spline=cols_q_ast, id_col=id_col, time_col=time_col)

#     # MSE calculation of the astrometric direction in linear model
#     u_mse, delta_mse = calc_dir_spline(spline_q_ast=spline_q_ast, q_obs=q_obs, t_obs=t_obs, 
#                         asteroid_id=asteroid_id, light_time=light_time, iters=2)

#     # Compare the two directions
#     print()
#     print_stars()
#     print(f'Compare directions from Geocenter to the first 16 asteroids using JPL state vectors in spline model:')
#     print(f'Data compared between dates {mjd0} and {mjd1}.')
#     print('JPL: Directions quoted by JPL. Astrometric RA/DEC converted to vectors using radec2dir.')
#     print('MSE: MSE calculations using linear model with function calc_dir_spline.\n')
#     err = direction_diff(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, verbose=True)

#     # Test result
#     thresh: float = 1.0
#     is_ok: bool = (err < thresh)
#     msg: str = 'PASS' if is_ok else 'FAIL'
#     print(f'**** {msg} ****')
#     return is_ok

# ********************************************************************************************************************* 
def test_calc_dir_ast2obs(obs_name: str, mjd0: int, mjd1: int):
    """
    Test MSE end to end calculation of asteroid direction vs. JPL.
    INPUTS:
        site_name: String ID for observatory in Skyfield.  One of 'geocenter' or 'palomar'
        mjd0: The first date to include in the test
        mjd1: The last date to include in the test
    """
    # Check obs_name
    if obs_name not in ('geo', 'ztf'):
        raise ValueError("Bad obs_name! Must be one of 'geo' or 'ztf'.")
    obs_long_name = obs_name_tbl[obs_name]

    # Import the JPL asteroid directions
    df_jpl = df_tbl[('JPL', obs_name)].copy()

    # Mask data down to the selected date range
    mask =  (mjd0 <= df_jpl.tObs) & (df_jpl.tObs <= mjd1)
    df_jpl = df_jpl[mask].reset_index(drop=True)

    # Selected times and asteroids
    ts = df_jpl.tObs.values
    asteroid_id = df_jpl.AsteroidID.values
    # Direction according to JPL
    u_jpl = df_jpl[cols_u].values

    # MSE calculation of the direction in end to end model
    df_mse = calc_dir_ast2obs(ts=ts, asteroid_id=asteroid_id, site_name=site_name)
    u_mse = df_mse[cols_u].values

    # Compare the two directions
    print()
    print_stars()
    print(f'Compare directions from {site_name} to the first 16 asteroids using MSE end to end model:')
    print(f'Data compared between dates {mjd0} and {mjd1}.')
    print('JPL: Directions quoted by JPL.  Astrometric RA/DEC converted to vectors using radec2dir.')
    print('MSE: MSE calculations using linear model with function calc_dir_ast2obs.\n')
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

    # Test the MSE integration vs. JPL by comparing state vectors; use both Geocenter and Palomar
    is_ok &= test_dir_vectors(obs_name='geo', mjd0=mjd0, mjd1=mjd1)
    is_ok &= test_dir_vectors(obs_name='ztf', mjd0=mjd0, mjd1=mjd1)
    
    # Test topos adjustment
    is_ok &= test_topos(mjd0=mjd0, mjd1=mjd1)

    # Test the asteroid directions in linear model with JPL state vectors
    is_ok &= test_dir_linear(state_vec_src='JPL', mjd0=mjd0, mjd1=mjd1)

    # Test the asteroid directions in linear model with MSE state vectors
    is_ok &= test_dir_linear(state_vec_src='MSE', mjd0=mjd0, mjd1=mjd1)

    # # Test the asteroid directions in splined model with JPL state vectors
    # is_ok &= test_dir_spline(mjd0=mjd0, mjd1=mjd1)

    # # Test asteroid directions in end to end model
    # is_ok &= test_calc_dir_ast2obs(site_name='geocenter', mjd0=mjd0, mjd1=mjd1)

    # Overall test result
    msg: str = 'PASS' if is_ok else 'FAIL'
    print()
    print_stars()
    print('Overall test of asteroid directions:')
    print(f'**** {msg} ****')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
