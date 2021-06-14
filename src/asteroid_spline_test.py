"""
Load integrated asteroid trajectories as Pandas DataFrames.
Add earth or sun position / state vectors to asteroid DataFrame.

Functions in this module:
get_df_shape(df, id_col)
make_spline_df(df, cols_spline, id_col, time_col)
spline_data(df, ts, cols_spline, id_col, time_col)
spline_ast_data(df_ast, ts, cols_spline)
spline_ast_vec(vec, ts)
spline_ast_elt(elt, ts)
test_ast_spline_elt()

Michael S. Emanuel
Sat Sep 21 10:38:38 2019
"""

# Core
import numpy as np

# Local imports
from planets_interp import get_sun_pos
from orbital_element import unpack_elt_df, unpack_elt_np, elt2pos, anomaly_M2f
from asteroid_data import load_ast_data, load_ast_pos, load_ast_vectors
from orbital_element_test import angle_distance, report_test
from asteroid_spline import make_spline_df, make_spline_ast_elt, make_spline_ast_pos, make_spline_ast_vec
from utils import print_stars

# ********************************************************************************************************************* 
# Column collections
cols_q = ['qx', 'qy', 'qz']
cols_v = ['vx', 'vy', 'vz']
cols_elt = ['a', 'e', 'inc', 'Omega', 'omega', 'f', 'M']

# Test thresholds
thresh_q_vec: float = 1.0E-6
thresh_q: float = 1.0E-8
thresh_v: float = 1.0E-10
thresh_angle: float = 1.0E-9
thresh_e: float = 1.0E-9

# ********************************************************************************************************************* 
def test_ast_spline_df(n0: int, n1: int):
    """Test splining of asteroid position directly and via orbital elements"""
    # Get test orbital elements from the first 10 asteroids
    df_ast = load_ast_data(n0=n0, n1=n1)

    # Array of times for output
    # Choose a set of interpolation times that exactly match the ORIGINAL sample times for testing
    ts = df_ast.mjd.values
    # Array of AsteroidIDs for output
    asteroid_id = df_ast.AsteroidID.values

    # Down-sample from every 4 days to every 8 days to get a real test
    interval_days: int = 8
    interval = 24*60*interval_days
    mask = (df_ast.TimeID % interval == 0)
    df_in = df_ast[mask].copy()

    # Spline asteroid position directly
    time_col = 'mjd'
    id_col = 'AsteroidID'
    spline_pos_direct = make_spline_df(df=df_in, cols_spline=cols_q, time_col=time_col, id_col=id_col)    
  
    # Spline vectors via orbital elements
    time_col = 'mjd'
    id_col = 'AsteroidID'
    spline_elt = make_spline_df(df=df_in, cols_spline=cols_elt, time_col=time_col, id_col=id_col)    

    # Wrap up a function that returns the position vector using the splined elements
    def spline_pos_from_elt(ts: np.ndarray, asteroid_id: np.ndarray):
        """Splineed position as a function of time and asteroid_id"""
        # Caclulate the splined orbital elements for the input times and asteroid_ids
        elt = spline_elt(ts, asteroid_id)
        # Unpack the elements
        a, e, inc, Omega, omega, f, M = unpack_elt_np(elt)        
        # Recover f from M
        f = anomaly_M2f(M=M, e=e)
        # Compute q in the heliocentric frame
        q_hel: np.ndarray = elt2pos(a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)
        # Position  of the sun
        q_sun: np.ndarray = get_sun_pos(ts=ts)
        # The position of the asteroid in the barycentric frame
        q: np.ndarray = q_hel + q_sun
        return q

    # Unpack / calculated the position from (1) original data (2) interpolated via vectors (3) interpolated via elements    
    q0 = df_ast[cols_q].values
    q1 = spline_pos_direct(ts, asteroid_id)
    q2 = spline_pos_from_elt(ts, asteroid_id)

    # Position error
    dq1: np.ndarray = q1 - q0
    dq2: np.ndarray = q2 - q0
    err_q1: np.ndarray = np.sqrt(np.sum(np.square(dq1), axis=-1))
    err_q2: np.ndarray = np.sqrt(np.sum(np.square(dq2), axis=-1))   

    # Report the results
    print()
    print_stars()
    print(f'Splined position using orbital elements; asteroids {n0} to {n1}.')
    print(f'Original data interval is 4 days.  Spline built using interval of {interval_days} days.')
    is_ok: bool = True
    is_ok &= report_test(err=err_q1, test_name='ast_spline_elt (position q, spline on vectors)', thresh=thresh_q_vec)
    is_ok &= report_test(err=err_q2, test_name='ast_spline_elt (position q, spline on elements)', thresh=thresh_q)

    return is_ok

# ********************************************************************************************************************* 
def test_ast_spline_elt(n0: int, n1: int, mjd0: int, mjd1: int):
    """Test end to end function splining asteroid orbital elements"""
    # Get test orbital elements from the selected asteroids and mask to selected date range
    df = load_ast_data(n0=n0, n1=n1)
    mask = (mjd0 <= df.mjd) & (df.mjd <= mjd1)
    df = df[mask]
    # Unpack the quoted elements
    a0, e0, inc0, Omega0, omega0, f0, M0 = unpack_elt_df(df)
    # Build the element spline
    spline_elt = make_spline_ast_elt(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # The times and asteroid IDs for the input; this is also where we want to evaluate the output spline
    ts = df.mjd.values
    asteroid_id = df.AsteroidID.values
    # Calculate the splined orbital elements for the input times and asteroid_ids
    elt = spline_elt(ts, asteroid_id)
    # Unpack the elements
    a1, e1, inc1, Omega1, omega1, f1, M1 = unpack_elt_np(elt)        
    # Recover f from M
    f1 = anomaly_M2f(M=M1, e=e1)

    # Evaluate errors on seven orbital elements
    err_a = np.abs(a1-a0)
    err_e = np.abs(e1-e0)
    err_inc = angle_distance(inc0, inc1)
    err_Omega = angle_distance(Omega0, Omega1)
    err_omega = angle_distance(omega0, omega1)
    err_f = angle_distance(f0, f1)
    err_M = angle_distance(M0, M1)

    # Report the results
    print()
    print_stars()
    print(f'Splined orbital elements vs. MSE calculated elements in DB; asteroids {n0} to {n1}.')
    is_ok: bool = True
    is_ok &= report_test(err=err_a, test_name='a', thresh=thresh_q)
    is_ok &= report_test(err=err_e, test_name='e', thresh=thresh_e)
    is_ok &= report_test(err=err_inc, test_name='inc', thresh=thresh_angle)
    is_ok &= report_test(err=err_Omega, test_name='Omega', thresh=thresh_angle)
    is_ok &= report_test(err=err_omega, test_name='omega', thresh=thresh_angle)
    is_ok &= report_test(err=err_f, test_name='f', thresh=thresh_angle)
    is_ok &= report_test(err=err_M, test_name='M', thresh=thresh_angle)

    return is_ok

# ********************************************************************************************************************* 
def test_ast_spline_pos(n0: int, n1: int, mjd0: int, mjd1: int):
    """Test function splining asteroid positions via orbital elements"""
    # Get test orbital positions from the selected asteroids and mask to selected date range
    df = load_ast_pos(n0=n0, n1=n1)
    mask = (mjd0 <= df.mjd) & (df.mjd <= mjd1)
    df = df[mask]
    # Unpack the quoted position
    q0 = df[cols_q].values
    # Build the position spline
    spline_pos = make_spline_ast_pos(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # The times and asteroid IDs for the input; this is also where we want to evaluate the output spline
    ts = df.mjd.values
    asteroid_id = df.AsteroidID.values
    # Calculate the splined asteroid position for the input times and asteroid_ids
    q1 = spline_pos(ts, asteroid_id)

    # Position error
    dq: np.ndarray = q1 - q0
    err_q: np.ndarray = np.sqrt(np.sum(np.square(dq), axis=-1))

    # Report the results
    print()
    print_stars()
    print(f'Splined position using orbital elements; asteroids {n0} to {n1}.')
    is_ok: bool = True
    is_ok &= report_test(err=err_q, test_name='ast_spline_pos', thresh=thresh_q)

    return is_ok

# ********************************************************************************************************************* 
def test_ast_spline_vec(n0: int, n1: int, mjd0: int, mjd1: int):
    """Test function splining asteroid state vectors via orbital elements"""
    # Get test state vectors from the selected asteroids and mask to selected date range
    df = load_ast_vectors(n0=n0, n1=n1)
    mask = (mjd0 <= df.mjd) & (df.mjd <= mjd1)
    df = df[mask]
    # Unpack the quoted state vectors (q and v)
    q0 = df[cols_q].values
    v0 = df[cols_v].values
    # Build the spline of state vectors
    spline_vec = make_spline_ast_vec(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # The times and asteroid IDs for the input; this is also where we want to evaluate the output spline
    ts = df.mjd.values
    asteroid_id = df.AsteroidID.values
    # Calculate the splined asteroid position for the input times and asteroid_ids
    q1, v1 = spline_vec(ts, asteroid_id)

    # Error in state vectors
    dq: np.ndarray = q1 - q0
    dv: np.ndarray = v1 - v0
    err_q: np.ndarray = np.sqrt(np.sum(np.square(dq), axis=-1))
    err_v: np.ndarray = np.sqrt(np.sum(np.square(dv), axis=-1))

    # Report the results
    print()
    print_stars()
    print(f'Splined state vectors using orbital elements; asteroids {n0} to {n1}.')
    is_ok: bool = True
    is_ok &= report_test(err=err_q, test_name='ast_spline_vec (q)', thresh=thresh_q)
    is_ok &= report_test(err=err_v, test_name='ast_spline_vec (v)', thresh=thresh_v)

    return is_ok

# ********************************************************************************************************************* 
if __name__ == '__main__':
    # Set range of asteroids to test
    n0: int = 1
    n1: int = 16
    # Set range of dates to test
    mjd0: int = 48000
    mjd1: int = 63000

    # Status update; track whether all tests pass
    print(f'Running test suite on asteroid_spline.py...')    
    is_ok: bool = True

    # Test generic asteroid spline back end
    is_ok &= test_ast_spline_df(n0=n0, n1=n1)

    # Test splining of orbital elements
    is_ok &= test_ast_spline_elt(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Test splining of asteroid position
    is_ok &= test_ast_spline_pos(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Test splining of asteroid state vectors
    is_ok &= test_ast_spline_vec(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1)

    # Report overall test results
    msg: str = 'PASS' if is_ok else 'FAIL'
    print()
    print_stars()
    print(f'Overall test result:')
    print(f'**** {msg} ****')
