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
from asteroid_direction import calc_direction_linear, calc_dir_ast2obs, calc_distance
from ra_dec import radec2dir
from db_utils import sp2df, df2db
from astro_utils import dist2sec
from ra_dec_test import direction_diff
from utils import print_stars

# ********************************************************************************************************************* 
# Column names for astrometric and apparent direction
cols_dir_ast = ['ux_ast', 'uy_ast', 'uz_ast']
cols_dir_app = ['ux_app', 'uy_app', 'uz_app']

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
    df[cols_dir_ast] = u_ast
    df[cols_dir_app] = u_app

    # Write asteroid directions to new DB table JPL.AsteroidDirections
    columns = ['AsteroidID', 'ObservatoryID', 'TimeID', 'mjd', 
               'RA_ast', 'DEC_ast', 'ux_ast', 'uy_ast', 'uz_ast', 
               'RA_app', 'DEC_app', 'ux_app', 'uy_app', 'uz_app', 
               'LightTime']
    df2db(df=df, schema='JPL', table='AsteroidDirection', columns=columns)

# ********************************************************************************************************************* 
def test_dir_linear(state_vec_src: str):
    """
    Test MSE astrometric direction calculation vs. JPL using the linear model.
    INPUTS:
        state_vec_src: One of 'JPL' or 'MSE' - the source of the state vectors used
    """
    # Check state_vec_src
    if state_vec_src not in ('JPL', 'MSE'):
        raise ValueError("Bad state_vec_src! Must be one of 'JPL' or 'MSE'.")
    # Table of SP name
    sp_name_tbl = {
        'JPL': 'JPL.AstrometricDirectionTest',
        'MSE': 'JPL.AstrometricDirectionAndVectorTest' 
    }

    # Import the JPL asteroid directions
    sp_name = sp_name_tbl[state_vec_src]
    df = sp2df(sp_name=sp_name)

    # Position and velocity according to JPL
    q_obs = df[cols_q_obs].values
    q_ast = df[cols_q_ast].values
    v_ast = df[cols_v_ast].values

    # Direction according to JPL
    u_jpl = df[cols_u].values

    # MSE calculation of the astrometric direction in linear model
    u_mse = calc_direction_linear(q_tgt=q_ast, v_tgt=v_ast, q_obs=q_obs)

    # Compare the two directions
    print(f'Compare directions on the first 16 asteroids using {state_vec_src} state vectors:')
    print('JPL: Directions quoted by JPL.  Astrometric RA/DEC converted to vectors using radec2dir.')
    print('MSE: MSE calculations using linear model with function calc_direction_linear in asteroid_direction.py.\n')
    direction_diff(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, verbose=True)

# ********************************************************************************************************************* 
def test_ast_dir_vectors():
    """
    Test for consistency between JPL and MSE integrations in the implied direction from earth geocenter to asteroid.
    """
    # Directions from two different sources
    df_jpl = sp2df(sp_name='JPL.AstrometricDirectionTest')
    df_mse = sp2df(sp_name='JPL.AstrometricDirectionAndVectorTest')

    # Mask both frames down to overlapping dates (multiples of 20)
    mask_interval = 20*1440
    mask_jpl = (df_jpl.TimeID % mask_interval) == 0
    df_jpl = df_jpl[mask_jpl].reindex()
    mask_mse = (df_mse.TimeID % mask_interval) == 0
    df_mse = df_mse[mask_mse].reindex()    

    # State vectors according to JPL integration
    q_obs_jpl = df_jpl[cols_q_obs].values
    q_ast_jpl = df_jpl[cols_q_ast].values
    v_ast_jpl = df_jpl[cols_v_ast].values

    # State vectors according to MSE integration
    q_obs_mse = df_mse[cols_q_obs].values
    q_ast_mse = df_mse[cols_q_ast].values
    v_ast_mse = df_mse[cols_v_ast].values

    # Directions calculated in linear model with both state vectors
    u_jpl = calc_direction_linear(q_tgt=q_ast_jpl, v_tgt=v_ast_jpl, q_obs=q_obs_jpl)
    u_mse = calc_direction_linear(q_tgt=q_ast_mse, v_tgt=v_ast_mse, q_obs=q_obs_mse)

    # Compare the two directions
    print('Compare directions on the first 16 asteroids between two sources of state vectors:')
    print('JPL: State vectors sources from JPL using Horizons.')
    print('MSE: Integration done in rebound.')
    print('Both sources have directions calculated using calc_direction_linear().')
    direction_diff(name1='JPL', name2='MSE', u1=u_jpl, u2=u_mse, verbose=True)

# ********************************************************************************************************************* 
def test_ast_dir(dfh: pd.DataFrame, asteroid_id: int, iters: int):
    """
    Test MSE direction calculation vs. JPL for one asteroid
    INTPUTS:
        dfh:            The DataFrame of JPL test data.
        asteroid_id:    The asteroid whose direction is being tested
    """
    # Create a copy of dfh
    dfh = dfh.copy()
    # Mask down to selected asteroid
    mask = (dfh.AsteroidID == asteroid_id)
    dfh = dfh[mask]

    # Inputs to calc_dir_ast2obs
    n0: int = asteroid_id
    n1: int = n0+2

    # Calculate the direction and light time
    dfm = calc_dir_ast2obs(n0=n0, n1=n1, iters=iters)

    # Mask down to selected AsteroidID
    mask = (dfm.AsteroidID == asteroid_id)
    dfm = dfm[mask]

    # Two time vectors for MSE data
    t_ast_mse = dfm.tAst.values
    t_obs_mse = dfm.tObs.values

    # Build spline for MSE asteroid direction and light time keyed by tObs
    x = t_obs_mse
    spline_u = CubicSpline(x=x, y=dfm[cols_dir].values)
    spline_LT = CubicSpline(x=x, y=dfm.LightTime.values)

    # Horizons arrays: both times
    t_ast_jpl = dfh.tAst.values
    t_obs_jpl = dfh.tObs.values
    # Horizons arrays: astrometric and apparent directions
    u_ast_jpl = dfh[cols_dir_ast].values
    u_app_jpl = dfh[cols_dir_app].values
    # Horizons array: light time
    light_time_jpl = dfh.LightTime.values

    # Spline MSE directions and light time to hopefully match JPL arrays
    x_jpl = t_obs_jpl
    u_mse = spline_u(x=x_jpl)
    light_time_mse = spline_LT(x=x_jpl)

    # Calculate the distance between the two directions
    du_ast = calc_distance(q0=u_ast_jpl, q1=u_mse)
    du_app = calc_distance(q0=u_app_jpl, q1=u_mse)
    # Convert the distances to arc seconds
    du_ast_sec = dist2sec(du_ast)
    du_app_sec = dist2sec(du_app)

    # Difference in light time
    dLT = np.abs(light_time_mse - light_time_jpl)*60.0

    # Threshold
    thresh_sec: float = 1.0
    du_sec_max = np.max(du_ast_sec)
    is_ok: bool = (du_sec_max < thresh_sec)
    msg: str = 'PASS' if is_ok else 'FAIL'

    # Report results
    print(f'Test asteroid direction for AsteroidID={asteroid_id} with {iters} light-time iterations.')
    print('\nDistance between JPL astrometric direction and MSE directions in arc seconds:')
    print(f'Mean: {np.mean(du_ast_sec):8.2f}')
    print(f'Max : {np.max(du_ast_sec):8.2f}')
    print('\nDistance between JPL apparent direction and MSE directions in arc seconds:')
    print(f'Mean: {np.mean(du_app_sec):8.2f}')
    print(f'Max : {np.max(du_app_sec):8.2f}')
    print(f'\nMean light time difference in seconds: {np.mean(dLT):5.3e}')
    print(f'**** {msg} ****')

# ********************************************************************************************************************* 
def main():
    # Build table JPL.AsteroidDirections from JPl.AsteroidDirections_Import
    # jpl_ast_dir_populate()
    # print('Regenerated DB table JPL.AsteroidDirection')

    # Test the asteroid directions with JPL state vectors
    print()
    print_stars()
    test_dir_linear(state_vec_src='JPL')

    # Test the asteroid directions with MSE state vectors
    print()
    print_stars()
    test_dir_linear(state_vec_src='MSE')

    # Test the MSE integration vs. JPL by comparing state vectors
    print()
    print_stars()
    test_ast_dir_vectors()
    
# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
