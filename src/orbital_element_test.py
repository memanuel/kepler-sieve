"""
Utilities for working with orbital elements.
Reference: Solar System Dynamics, C.D. Murray and S.F. Dermott, 1999

Michael S. Emanuel
2021-01-29
"""

# Core
import numpy as np
import pandas as pd

# Local imports
from asteroid_element import get_ast_data
from planets_interp import get_sun_vectors
from orbital_element import unpack_elt_df, elt2pos, elt2vec
from orbital_element import anomaly_f2E, anomaly_E2f, anomaly_E2M, anomaly_M2E_danby, anomaly_M2E, anomaly_f2M
from orbital_element_table import anomaly_M2E_table

# ********************************************************************************************************************* 
def make_test_elements(n1: int=10000, epoch: int = 59000):
    """
    Return a set of test orbital elements
    INPUTS:
        n1:         Last asteroid number to include
        epoch:      Epoch as of which vectors and elements computed
                    Should be a multiple of 4 in (48000, 63000)
    OUTPUTS:
        elts:       DataFrame with state vectors and orbital elements from the saved MSE integration
    """
    # Build test elements; save them to disk; and return them
    elts = get_ast_data(n0=0, n1=n1, epoch=epoch)
    fname = '../data/test/asteroid_elements.h5'
    elts.to_hdf(fname, key='elts')
    return elts

# ********************************************************************************************************************* 
def get_test_elements():
    """
    Return a set of test orbital elements
    INPUTS:
        None
    OUTPUTS:
        elts:       DataFrame with orbital elements from the saved MSE integration
    """
    fname = '../data/test/asteroid_elements.h5'
    elts = pd.read_hdf(fname, key='elts')
    return elts

# ********************************************************************************************************************* 
def angle_distance(x, y):
    """
    Calculate a distance between two angles on the unit sphere
    This distance 
    INPUTS:
        x:      The first angle to be compared
        y:      The second angle to be compared
    OUTPUTS:
        r:      The distance 2 sin( (x-y)/2 )
    """
    # The difference in the angles
    tau = 2.0 * np.pi 
    d1 = (x - y) % tau
    d2 = (y - x) % tau
    d = np.minimum(d1, d2)
    # Use sine of the half angle
    return 2.0 * np.sin(0.5 * d)

# ********************************************************************************************************************* 
def report_test(err: np.array, test_name: str, thresh: float) -> bool:
    """
    Report results of one test
    INPUTS:
        err:            Array of error values
        test_name:      Name of this test
    """
    # Calculate RMS and max error
    err_rms = np.sqrt(np.mean(np.square(err)))
    err_max = np.max(err)

    # Is the result a pass?
    is_ok: bool = err_max < thresh
    result: str = 'PASS' if is_ok else 'FAIL'

    # Report results
    print(f'\nTest: {test_name}')
    print(f'RMS error: {err_rms:5.3e}')
    print(f'Max error: {err_max:5.3e}')
    print(f'*** {result} ***')
    return is_ok

# ********************************************************************************************************************* 
def test_E2f() -> bool:
    """Test round trip conversion between E and f"""
    # Get test elements and unpack them
    elts = get_test_elements()
    e = elts.e.values
    f = elts.f.values
    
    # Calculate E from f and e
    E = anomaly_f2E(f=f, e=e)
    # Recover f from E
    f2 = anomaly_E2f(E=E, e=e)

    # Calculate the distance between these two angles
    err = angle_distance(f, f2)

    # Report the results
    return report_test(err=err, test_name='E2f', thresh=1.0E-9)

# ********************************************************************************************************************* 
def test_E2M():
    """Test conversion from E to M"""
    # Get test elements and unpack them
    elts = get_test_elements()
    e = elts.e.values
    f = elts.f.values
    M = elts.M.values
    
    # Calculate E from f and e
    E = anomaly_f2E(f=f, e=e)
    # Recover M from E
    M2 = anomaly_E2M(E=E, e=e)

    # Calculate the distance between these two angles
    err = angle_distance(M, M2)

    # Report the results
    return report_test(err=err, test_name='E2M', thresh=1.0E-9)

# ********************************************************************************************************************* 
def test_M2E_danby(n: int = 3):
    """Test conversion from M to E using Danby iteration method"""    
    # Get test elements and unpack them
    elts = get_test_elements()
    e = elts.e.values
    f = elts.f.values
    M = elts.M.values
    
    # Calculate E from f and e
    E = anomaly_f2E(f=f, e=e)
    # Recover E from M
    E2 = anomaly_M2E_danby(M=M, e=e, n=n)

    # Calculate the distance between these two angles
    err = angle_distance(E, E2)

    # Report the results
    return report_test(err=err, test_name=f'M2E (Danby, {n} iterations)', thresh=1.0E-9)

# ********************************************************************************************************************* 
def test_M2E_table():
    """Test conversion from M to E using interpolation table"""
    # Get test elements and unpack them
    elts = get_test_elements()
    e = elts.e.values
    f = elts.f.values
    M = elts.M.values
    
    # Calculate E from f and e
    E = anomaly_f2E(f=f, e=e)
    # Recover E from M
    E2 = anomaly_M2E_table(M=M, e=e)

    # Calculate the distance between these two angles
    err = angle_distance(E, E2)

    # Report the results
    return report_test(err=err, test_name='M2E (interpolation table; then 2 iterations of Danby)', thresh=1.0E-9)

# ********************************************************************************************************************* 
def test_M2E():
    """Test conversion from M to E"""
    # Get test elements and unpack them
    elts = get_test_elements()
    e = elts.e.values
    f = elts.f.values
    M = elts.M.values
    
    # Calculate E from f and e
    E = anomaly_f2E(f=f, e=e)
    # Recover E from M
    E2 = anomaly_M2E(M=M, e=e)

    # Calculate the distance between these two angles
    err = angle_distance(E, E2)

    # Report the results
    return report_test(err=err, test_name='M2E', thresh=1.0E-9)

# ********************************************************************************************************************* 
def test_M2E_internal(n: int = 3):
    """Test conversion from M to E for internal consistency on a round trip"""
    # Get test elements and unpack them
    elts = get_test_elements()
    e = elts.e.values
    f = elts.f.values

    # Calculate M from f and e using Kepler's equation
    M = anomaly_f2M(f=f, e=e)

    # Calculate E from f and e
    E1 = anomaly_f2E(f=f, e=e)
    # Recover E from M
    # E2 = anomaly_M2E(M=M, e=e, n=n)
    E2 = anomaly_M2E_danby(M=M, e=e, n=n)

    # Calculate the distance between these two angles
    err = angle_distance(E1, E2)

    # Report the results
    return report_test(err=err, test_name=f'M2E round trip ({n} iterations)', thresh=1.0E-9)

# ********************************************************************************************************************* 
def test_elt2pos():
    """Test conversion from orbital elements to position"""
    # Get test elements and unpack them
    elts = get_test_elements()

    # The position according to the integration
    cols_q = ['qx', 'qy', 'qz']
    q: np.array = elts[cols_q].values

    # Unpack orbital elements
    a, e, inc, Omega, omega, f = unpack_elt_df(elts)
    
    # Compute q in the heliocentric frame
    q_hel: np.array = elt2pos(a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)

    # Position of the sun
    ts = elts.mjd.values
    q_sun, v_sun = get_sun_vectors(ts)

    # The recovered position of the asteroid
    q2: np.array = q_hel + q_sun

    # Position reconstruction error
    dq: np.array = q2 - q
    err: np.array = np.sqrt(np.sum(np.square(dq), axis=1))
    # Report the results
    return report_test(err=err, test_name='elt2pos', thresh=1.0E-9)

# ********************************************************************************************************************* 
def test_elt2vec():
    """Test conversion from orbital elements to state vectors"""
    # Get test elements and unpack them
    elts = get_test_elements()

    # The state vectors according to the integration
    cols_q = ['qx', 'qy', 'qz']
    cols_v = ['vx', 'vy', 'vz']
    q: np.array = elts[cols_q].values
    v: np.array = elts[cols_v].values

    # Unpack orbital elements
    a, e, inc, Omega, omega, f = unpack_elt_df(elts)

    # Compute q and v in the heliocentric frame
    q_hel, v_hel = elt2vec(a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)

    # Position and velocity of the sun
    ts = elts.mjd.values
    q_sun, v_sun = get_sun_vectors(ts)

    # The recovered position and velocity of the asteroid
    q2: np.array = q_hel + q_sun
    v2: np.array = v_hel + v_sun

    # Position reconstruction error
    dq: np.array = q2 - q
    dv: np.array = v2 - v
    err_q: np.array = np.sqrt(np.sum(np.square(dq), axis=-1))
    err_v: np.array = np.sqrt(np.sum(np.square(dv), axis=-1))

    # Report the results
    is_ok_q = report_test(err=err_q, test_name='elt2vec (position q)', thresh=1.0E-9)
    is_ok_v = report_test(err=err_v, test_name='elt2vec (velocity v)', thresh=1.0E-9)
    return (is_ok_q and is_ok_v)

# ********************************************************************************************************************* 
def test_all():
    """Running test suite on orbital elements"""
    # Maintain overall flag with result
    is_ok: bool = True

    # Run all tests 
    is_ok &= test_E2f()
    is_ok &= test_E2M()
    # is_ok &= test_M2E_danby()
    # is_ok &= test_M2E_table()
    is_ok &= test_M2E()
    is_ok &= test_M2E_internal(n=3)
    is_ok &= test_elt2pos()
    is_ok &= test_elt2vec()

    # Report results
    result: str = 'PASS' if is_ok else 'FAIL'
    print(f'\nOverall results for orbital elements:')
    print(f'*** {result} ***')
    return is_ok

# ********************************************************************************************************************* 
if __name__ == '__main__':
    print(f'Running test suite on orbital_element.py...')
    test_all()
