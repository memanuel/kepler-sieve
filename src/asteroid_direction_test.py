"""
Test direction from known asteroids to Earth center implied by integrated trajectories.
Example call:
$ python asteroid_direction_test.py

Functions in this module:
jpl_ast_dir_populate(n0, n1)
test_ast_dir(name1, name2, u1, u2, lt1, lt2)
test_dir_vectors(obs_name, mjd0, mjd1)
test_topos(mjd0, mjd1)
test_dir_linear(state_vec_src, mjd0, mjd1)
test_calc_dir_ast2obs(mjd0, mjd1)
test_asteroid_dir_db(mjd0, mjd1)
test_asteroid_dir_linear_vs_itersp(mjd0, mjd1, n0, n1)
test_asteroid_dir_spline(mjd0, mjd1, n0, n1)
main()

Michael S. Emanuel
2021-06-04
"""

# Core
import numpy as np
# import pandas as pd

# Astronomy
from astropy.units import deg

# Local imports
from asteroid_spline import make_spline_ast_vec, make_spline_ast_dir
from asteroid_direction import c, calc_dir_linear, calc_dir_ast2obs, calc_dir_ast2obs_spline, prep_ast_block
from planets_interp import get_earth_pos
from ra_dec import radec2dir
from db_utils import sp2df, df2db
from astro_utils import dist2deg
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

# Test thresholds
thresh_u_vec: float = 0.1
thresh_u_topos: float = 0.05
thresh_u_linear: float = 1.0
thresh_u_tot: float = 1.0
thresh_lt: float = 1.0E-4

# ********************************************************************************************************************* 
def jpl_ast_dir_populate(n0: int, n1:int):
    """
    Populate DB table JPL.AsteroidDirection from JPL.AsteroidDirectionImport
    """
    # Get the imported RA/DEC from JPL
    sp_name = 'JPL.GetAsteroidDirection_Import'
    params = {
        'n0': n0,
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
               'LightTime', 'r', 'rDot', 'delta', 'deltaDot', 'Mag']
    df2db(df=df, schema='JPL', table='AsteroidDirection', columns=columns)

# *************************************************************************************************
def test_ast_dir(name1: str, name2: str, u1: np.ndarray, u2: np.ndarray, 
                 lt1: np.ndarray, lt2: np.ndarray, verbose: bool=False) -> float:
    """
    Report 
    INPUTS:
        name1: Descriptive name of the first source, e.g. 'JPL'
        name2: Descriptive name of the second source, e.g. 'MSE'
        u1:    Array of directions from source 1; shape Nx3
        u2:    Array of directions from source 2; shape Nx3
        lt1:   Array of light times from source 1; shape N
        lt2:   Array of light times from source 2; shape N
        verbose: Whether to report results to console
    """
    # Difference in unit directions
    du = u2 - u1
    du_norm = np.sqrt(np.sum(np.square(du), axis=-1))
    du_deg =  dist2deg(du_norm)
    du_sec =  du_deg * 3600.0

    # Calculate mean, median and max difference in degrees
    du_mean = np.mean(du_sec)
    du_median = np.median(du_sec)
    du_max = np.max(du_sec)

    # Calculate light time error
    dlt_mean = np.mean(np.abs(lt2-lt1))

    if verbose:
        print(f'Angle Difference: {name2} vs. {name1} in arc seconds')
        print(f'*Mean  : {du_mean:9.3f}*')
        print(f' Median: {du_median:9.3f}')
        print(f' Max   : {du_max:9.3f}')
        print(f'Light time difference in minutes:')
        print(f' Mean  : {dlt_mean:9.2e}')

    # Return the difference of direction vectors in seconds of arc
    return du_mean, dlt_mean
    
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
    # Compute light time
    lt_jpl = delta_jpl / c
    lt_mse = delta_mse / c

    # Compare the two directions
    print()
    print_stars()
    print(f'Compare directions from {obs_long_name} to the first 16 asteroids between two sources of state vectors:')
    print(f'Data compared between dates {mjd0} and {mjd1}.')
    print('JPL: State vectors sourced from JPL using Horizons.')
    print('MSE: Integration done in rebound.')
    print('Both sources have directions calculated using calc_dir_linear().')
    du, dlt = test_ast_dir(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, lt1=lt_jpl, lt2=lt_mse, verbose=True)

    # Test result
    is_ok: bool = (du < thresh_u_vec) and (dlt < thresh_lt)
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
    # Light time according to JPL
    lt_geo_jpl = df_geo.LightTime.values
    lt_ztf_jpl = df_ztf.LightTime.values

    # Directions calculated in MSE linear model
    u_geo_mse, delta_geo_mse = calc_dir_linear(q_tgt=q_ast_geo, v_tgt=v_ast_geo, q_obs=q_obs_geo)
    u_ztf_mse, delta_ztf_mse = calc_dir_linear(q_tgt=q_ast_ztf, v_tgt=v_ast_ztf, q_obs=q_obs_ztf)
    # Compute light time
    lt_geo_mse = delta_geo_mse / c
    lt_ztf_mse = delta_ztf_mse / c

    # Change in direction according to JPL and MSE
    du_jpl = u_ztf_jpl - u_geo_jpl
    du_mse = u_ztf_mse - u_geo_mse

    # Change in light time according to JPL and MSE
    dlt_jpl = lt_ztf_jpl - lt_geo_jpl
    dlt_mse = lt_ztf_mse - lt_geo_mse

    # Compare the two changes due to topos
    print()
    print_stars()
    print(f'Compare directions change between Geocenter and Palomar on the first 16 asteroids:')
    print('JPL: Quoted direction ZTF->Asteroid minus GEO->Asteroid.')
    print('MSE: Calculated direction ZTF->Asteroid minus GEO->Asteroid using JPL state vectors.')
    print('Both sources have directions calculated using calc_dir_linear().')
    du, dlt = test_ast_dir(name1='JPL', name2='MSE', u1=du_jpl, u2=du_mse, lt1=dlt_jpl, lt2=dlt_mse, verbose=True)
    # Test result
    is_ok: bool = (du < thresh_u_topos) and (dlt < thresh_lt)
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

    # Direction and light time according to JPL
    u_jpl = df[cols_u].values
    lt_jpl = df.LightTime.values

    # MSE calculation of the astrometric direction in linear model
    u_mse, delta_mse = calc_dir_linear(q_tgt=q_ast, v_tgt=v_ast, q_obs=q_obs)
    lt_mse = delta_mse / c

    # Compare the two directions
    print()
    print_stars()
    print(f'Compare directions to the first 16 asteroids using {state_vec_src} state vectors in linear model:')
    print(f'Data compared between dates {mjd0} and {mjd1}.')
    print('JPL: Directions quoted by JPL.  Astrometric RA/DEC converted to vectors using radec2dir.')
    print('MSE: MSE calculations using linear model with function calc_dir_linear in asteroid_direction.py.')
    print('Both models tested combining JPL data for Geocenter and Palomar observatory locations.')
    du, dlt = test_ast_dir(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, lt1=lt_jpl, lt2=lt_mse, verbose=True)

    # Test result
    is_ok: bool = (du < thresh_u_tot) and (dlt < thresh_lt)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# ********************************************************************************************************************* 
def test_calc_dir_ast2obs(mjd0: int, mjd1: int):
    """
    Test MSE end to end calculation of asteroid direction vs. JPL.
    INPUTS:
        mjd0:       The first date to include in the test
        mjd1:       The last date to include in the test
    """
    # Import the JPL asteroid directions and MSE state vectors
    df_jpl = df_tbl_src['MSE'].copy()

    # Mask data down to the selected date range
    mask =  (mjd0 <= df_jpl.tObs) & (df_jpl.tObs <= mjd1)
    df_jpl = df_jpl[mask].reset_index(drop=True)

    # Selected times and asteroids
    t_obs = df_jpl.tObs.values
    asteroid_id = df_jpl.AsteroidID.values
    # Direction and light time according to JPL
    u_jpl = df_jpl[cols_u].values
    lt_jpl = df_jpl.LightTime.values
    # Observatory location according to JPL (plus MSE topos correction)
    q_obs = df_jpl[cols_q_obs].values

    # MSE calculation of the direction in end to end model
    df_mse = calc_dir_ast2obs(t_obs=t_obs, asteroid_id=asteroid_id, q_obs=q_obs)
    u_mse = df_mse[cols_u].values
    lt_mse = df_mse.LightTime.values

    # Compare the two directions
    print()
    print_stars()
    print(f'Compare directions to the first 16 asteroids using MSE end to end model:')
    print(f'Data compared between dates {mjd0} and {mjd1}.')
    print('JPL: Directions quoted by JPL.  Astrometric RA/DEC converted to vectors using radec2dir.')
    print('MSE: MSE calculations using end to end model with function calc_dir_ast2obs.\n')
    du, dlt = test_ast_dir(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, lt1=lt_jpl, lt2=lt_mse, verbose=True)

    # Test result
    is_ok: bool = (du < thresh_u_tot) and (dlt < thresh_lt)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# ********************************************************************************************************************* 
def test_asteroid_dir_db(mjd0: int, mjd1: int):
    """Test asteroid direction calculations written to database"""
    # Run SP comparing JPL and MSE directions
    sp_name = 'KS.TestAsteroidDirection'
    params = {
        'mjd0': mjd0,
        'mjd1': mjd1,
    }
    df = sp2df(sp_name=sp_name, params=params)

    # Extract the directions according to JPL and MSE on overlapping dates (every 20th day)
    cols_jpl = [col_nm + '_jpl' for col_nm in cols_u]
    cols_mse = [col_nm + '_mse' for col_nm in cols_u]
    u_jpl = df[cols_jpl].values
    u_mse = df[cols_mse].values
    # Extract the light times
    lt_jpl = df.LightTime_jpl.values
    lt_mse = df.LightTime_mse.values

    # Compare the two directions
    print()
    print_stars()
    print(f'Compare directions between JPL and MSE data written to DB table KS.AsteroidDirection:')
    print('JPL: Directions quoted by JPL.  Astrometric RA/DEC converted to vectors using radec2dir.')
    print('MSE: MSE calculations using calc_dir_ast2obs and written to DB in program asteroid_direction.py.\n')
    du, dlt = test_ast_dir(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, lt1=lt_jpl, lt2=lt_mse, verbose=True)

    # Test result
    is_ok: bool = (du < thresh_u_tot) and (dlt < thresh_lt)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# ********************************************************************************************************************* 
def test_asteroid_dir_linear_vs_itersp(mjd0: int, mjd1: int, n0: int, n1: int):
    """Test asteroid directions in linear model vs. iterated spline model."""
    # Prepare asteroid block
    interval_min = 60   # test spline on directions sampled every hour
    t_obs, asteroid_id = prep_ast_block(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1, interval_min=interval_min)
    # Earth position at these times
    q_obs = get_earth_pos(ts=t_obs)

    # Asteroid position in the iterated spline model
    df_itersp = calc_dir_ast2obs_spline(t_obs=t_obs, asteroid_id=asteroid_id, q_obs=q_obs)
    # Asteroid position in the linear model
    df_linear = calc_dir_ast2obs(t_obs=t_obs, asteroid_id=asteroid_id, q_obs=q_obs)

    # Extract direction arrays
    u_ast_itersp = df_itersp[cols_u].values
    u_ast_linear = df_linear[cols_u].values
    # Extract light time arrays
    lt_itersp = df_itersp.LightTime.values
    lt_linear = df_linear.LightTime.values

    # Compare two different methods
    print()
    print_stars()
    print('Compare directions between iterated spline model and linear model.')
    print('Test compares every 1 hour to exercise both splines.')
    print('itersp: Asteroid position is evaluated via splined elements; iterated to solve for light time.')
    print('linear: Asteroid position and velocity evalute via splined elements, then passed to linear model.\n')
    du, dlt = test_ast_dir(name1='itersp', name2='linear', u1=u_ast_itersp, u2=u_ast_linear, 
                            lt1=lt_itersp, lt2=lt_linear, verbose=True)

    # Test result
    is_ok: bool = (du < thresh_u_tot) and (dlt < thresh_lt)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# ********************************************************************************************************************* 
def test_asteroid_dir_spline(mjd0: int, mjd1: int, n0: int, n1: int):
    """Test asteroid directions splined directly from saved directions vs. those built from splined vectors."""
    # Prepare asteroid block
    interval_min = 60   # test spline on directions sampled every hour
    t_obs, asteroid_id = prep_ast_block(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1, interval_min=interval_min)
    # Build spline of asteroid direction
    spline_u = make_spline_ast_dir(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)
    # Build spline of asteroid vectors
    spline_ast_vec = make_spline_ast_vec(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Evaluate spline of asteroid directions
    u_ast_spline, lt_spline = spline_u(t_obs, asteroid_id)
    # Evaluate asteroid position and velocity
    q_ast, v_ast = spline_ast_vec(ts=t_obs, asteroid_id=asteroid_id)
    # Earth position at these times
    q_obs = get_earth_pos(ts=t_obs)
    # Calculate asteroid directions from position and velocity
    u_ast_linear, delta_linear = calc_dir_linear(q_tgt=q_ast, v_tgt=v_ast, q_obs=q_obs)
    lt_linear = delta_linear / c

    # Compare two different methods
    print()
    print_stars()
    print('Compare directions between direct spline of saved directions and full calculation with splined elements.')
    print('Directions are saved every 4 days in DB; test compares every 1 hour to exercise the spline.')
    print('linear: Spline orbital elements into vectors; then compute direciton using linear model on q_ast, v_ast.')
    print('spline: Directly spline u_ast and light_time saved in DB table KS.AsteroidDirection.\n')
    du, dlt = test_ast_dir(name1='linear', name2='spline', u1=u_ast_linear, u2=u_ast_spline, 
                            lt1=lt_linear, lt2=lt_spline, verbose=True)

    # Test result
    is_ok: bool = (du < thresh_u_tot) and (dlt < thresh_lt)
    msg: str = 'PASS' if is_ok else 'FAIL'
    print(f'**** {msg} ****')
    return is_ok

# *************************************************************************************************
# Console program runs all tests.
# *************************************************************************************************

# ********************************************************************************************************************* 
def main():
    # Build table JPL.AsteroidDirections from JPl.AsteroidDirections_Import
    # jpl_ast_dir_populate(n0=0, n1=100)
    # print('Regenerated DB table JPL.AsteroidDirection')

    # Set date range for testing
    epoch: int = 59000
    half_width: int = 2000
    mjd0: int = epoch - half_width
    mjd1: int = epoch + half_width

    # Set asteroid range for testing splines
    n0: int = 1
    n1: int = 17

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

    # Test asteroid directions in end to end model
    is_ok &= test_calc_dir_ast2obs(mjd0=mjd0, mjd1=mjd1)

    # Test linear model vs. iterated spline model
    is_ok &= test_asteroid_dir_linear_vs_itersp(mjd0=mjd0, mjd1=mjd1, n0=n0, n1=n1)

    # Test asteroid directions written to DB
    is_ok &= test_asteroid_dir_db(mjd0=mjd0, mjd1=mjd1)

    # Test direct spline of saved asteroid directions vs. calculation via splined vectors
    is_ok &= test_asteroid_dir_spline(mjd0=mjd0, mjd1=mjd1, n0=n0, n1=n1)

    # Overall test result
    msg: str = 'PASS' if is_ok else 'FAIL'
    print()
    print_stars()
    print('Overall test of asteroid directions:')
    print(f'**** {msg} ****')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
