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

    # Convert numpy arrays to angles in degrees
    ra = df.RA.values * deg
    dec = df.DEC.values * deg
    # The observation time is just a number array
    obstime_mjd = df.mjd.values

    # Convert RA/DEC to a direction
    u = radec2dir(ra=ra, dec=dec, obstime_mjd=obstime_mjd)

    # Save directions to DataFrame
    cols_dir = ['ux', 'uy', 'uz']
    df[cols_dir] = u

    # Write asteroid directions to new DB table JPL.AsteroidDirections
    columns = ['AsteroidID', 'TimeID', 'mjd', 'RA', 'DEC', 'ux', 'uy', 'uz', 
               'r', 'rDot', 'delta', 'deltaDot', 'LightTime', 'Mag']
    df2db(df=df, schema='JPL', table='AsteroidDirection', columns=columns)


# ********************************************************************************************************************* 
def test_ast_dir(dfh: pd.DataFrame, asteroid_id: int):
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
    dfm = calc_dir_ast2obs(n0=n0, n1=n1)

    # Mask down to selected AsteroidID
    mask = (dfm.AsteroidID == asteroid_id)
    dfm = dfm[mask]

    # Two time vectors for MSE data
    t_ast_mse = dfm.tAst.values
    t_obs_mse = dfm.tObs.values

    # Build spline for MSE asteroid direction and light time keyed by either tAst or tObs
    cols_dir =['ux', 'uy', 'uz']
    spline_u_tAst = CubicSpline(x=t_ast_mse, y=dfm[cols_dir].values)
    # spline_u_tObs = CubicSpline(x=t_obs_mse, y=dfm[cols_dir].values)

    # Build splines for MSE light time keyed by either tAst or tObs
    spline_LT_tAst = CubicSpline(x=t_ast_mse, y=dfm.LightTime.values)
    # spline_LT_tObs = CubicSpline(x=t_obs_mse, y=dfm.LightTime.values)

    # Extract vectors from Horizons test frame for both times, direction and light time
    t_ast_jpl = dfh.tAst.values
    t_obs_jpl = dfh.tObs.values
    u_jpl = dfh[cols_dir].values
    light_time_jpl = dfh.LightTime.values

    # Spline MSE directions and light time
    u_mse = spline_u_tAst(x=t_ast_jpl)
    light_time_mse = spline_LT_tAst(x=t_ast_jpl)

    # Calculate the distance between the two directions
    du = calc_distance(q0=u_jpl, q1=u_mse)
    # Convert the distance to arc seconds
    du_sec = dist2sec(du)
    # Mean and max
    du_sec_mean = np.mean(du_sec)
    du_sec_max = np.max(du_sec)

    # Difference in light time
    dLT = np.abs(light_time_mse - light_time_jpl)
    dLT_mean = np.mean(dLT)*60

    # Threshold
    thresh_sec: float = 1.0
    is_ok: bool = (du_sec_max < thresh_sec)
    msg: str = 'PASS' if is_ok else 'FAIL'

    # Report results
    print(f'Test asteroid direction for AsteroidID={asteroid_id}.')
    print('Distance between JPL and MSE directions in arc seconds:')
    print(f'Mean: {du_sec_mean:8.2f}')
    print(f'Max : {du_sec_max:8.2f}')
    print(f'Mean light time difference in seconds: {dLT_mean:5.3e}')
    print(f'**** {msg} ****')

# ********************************************************************************************************************* 
def main():
    # Build table JPL.AsteroidDirections from JPl.AsteroidDirections_Import
    # jpl_ast_dir_populate()

    # Get test DataFrame from JPL
    sp_name = 'JPL.AsteroidDirectionTest'
    dfh = sp2df(sp_name=sp_name)

    # Test the asteroid directions    
    test_ast_dir(AsteroidID=1)
    pass

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
