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
from asteroid_direction import calc_dir_ast2obs, calc_distance
from ra_dec import radec2dir
from db_utils import sp2df, df2db
from astro_utils import dist2sec
from utils import print_stars

# ********************************************************************************************************************* 
# Column names for astrometric and apparent direction
cols_dir = ['ux', 'uy', 'uz']
cols_dir_ast = ['ux_ast', 'uy_ast', 'uz_ast']
cols_dir_app = ['ux_app', 'uy_app', 'uz_app']

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
    columns = ['AsteroidID', 'TimeID', 'mjd', 
               'RA_ast', 'DEC_ast', 'ux_ast', 'uy_ast', 'uz_ast', 
               'RA_app', 'DEC_app', 'ux_app', 'uy_app', 'uz_app', 
               'LightTime']
    df2db(df=df, schema='JPL', table='AsteroidDirection', columns=columns)

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

    # Take every 4th row to get data every 20 days; this matches exactly on MSE integration points
    dfh = dfh.iloc[::20]
    dfh.reset_index(inplace=True)

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
    jpl_ast_dir_populate()
    print('Regenerated DB table JPL.AsteroidDirection')

    # Get test DataFrame from JPL
    sp_name = 'JPL.AsteroidDirectionTest'
    dfh = sp2df(sp_name=sp_name)

    # Test the asteroid directions
    print_stars()
    test_ast_dir(dfh=dfh, asteroid_id=1, iters=2)
    pass

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
