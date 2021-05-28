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
from asteroid_data import ast_add_sun_vectors
from orbital_element import anomaly_f2E, anomaly_E2f, anomaly_E2M, anomaly_M2E, elt2pos

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
    d = np.abs(x - y)
    # Use sine of the half angle
    return 2.0 * np.sin(0.5 * d)

# ********************************************************************************************************************* 
def report_test(err: np.array, test_name: str, thresh: float):
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

# ********************************************************************************************************************* 
def test_E2f():
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
    report_test(err=err, test_name='E2f', thresh=1.0E-9)

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
    report_test(err=err, test_name='E2M', thresh=1.0E-9)

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
    report_test(err=err, test_name='M2E', thresh=1.0E-9)

# ********************************************************************************************************************* 
def test_elt2pos():
    """Test conversion from orbital elements to position"""
    # Get test elements and unpack them
    elts = get_test_elements()

    # The position according to the integration
    cols_q = ['qx', 'qy', 'qz']
    q = elts[cols_q].values

    # Unpack orbital elements
    a = elts.a.values
    e = elts.e.values
    inc = elts.inc.values
    Omega = elts.Omega.values
    omega = elts.omega.values
    f = elts.f.values
    
    # Add sun vectors; need this later to convert q from heliocentric to barycentric
    ast_add_sun_vectors(elts)
    # Compute q in the heliocentric frame
    q_hel = elt2pos(a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)
    # Position of the sun
    cols_q_sun = ['sun_qx', 'sun_qy', 'sun_qz']
    q_sun = elts[cols_q_sun].values

    # The recovered position of the asteroid
    q2 = q_hel + q_sun

    # Position reconstruction error
    dq = q2 - q
    err = np.sqrt(np.sum(np.square(dq), axis=1))
    # Report the results
    report_test(err=err, test_name='elt2pos', thresh=1.0E-9)

# ********************************************************************************************************************* 
def test_all():
    """Running test suite on orbital elements"""
    test_E2f()
    test_E2M()
    test_M2E()
    test_elt2pos()

# ********************************************************************************************************************* 
if __name__ == '__main__':
    print(f'Running test suite on orbital_element.py...')
    test_all()
