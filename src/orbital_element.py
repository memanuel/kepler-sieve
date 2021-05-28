"""
Utilities for working with orbital elements.
Reference: Solar System Dynamics, C.D. Murray and S.F. Dermott, 1999

Michael S. Emanuel
2021-01-29
"""

# Core
import numpy as np
import pandas as pd

# Utility
from collections import namedtuple

# ********************************************************************************************************************* 
# Named tuple data type for orbital elements with just the six 
# a, e, inc, Omega, omega, f
OrbitalElement_aeiOof = namedtuple('OrbitalElement', 'a e inc Omega omega f')

# Named tuple data type for orbital elements that includes the mean anomaly M
# a, e, inc, Omega, omega, f, M
OrbitalElement_aeiOofM = namedtuple('OrbitalElement', 'a e inc Omega omega f M')

# ********************************************************************************************************************* 
# Functions for converting between anomalies: mean, true, and eccentric
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def anomaly_E2f(E: np.array, e: np.array) -> np.array:
    """
    Convert the eccentric anomaly E to the true anomaly f
    INPUTS:
        E: The eccentric anomaly
        e: The eccentricity
    OUTPUTS:
        f: The true anomaly
    """
    # SSD equation 2.46
    tan_half_f: np.array = np.sqrt((1.0 + e) / (1.0 - e)) * np.tan(E/2.0)
    f: np.array = np.arctan(tan_half_f) * 2.0
    return f

# ********************************************************************************************************************* 
def anomaly_f2E(f: np.array, e: np.array) -> np.array:
    """
    Convert the true anomaly f to the eccentric anomaly E
    INPUTS:
        f: The true anomaly
        e: The eccentricity
    OUTPUTS:
        E: The eccentric anomaly
    """
    # SSD equation 2.46
    tan_half_f: np.array = np.tan(f * 0.5)
    tan_half_E: np.array = tan_half_f / np.sqrt((1.0 + e) / (1.0 - e))
    E: np.array = np.arctan(tan_half_E) * 2.0
    return E

# ********************************************************************************************************************* 
def anomaly_E2M(E: np.array, e: np.array) -> np.array:
    """
    Convert the eccentric anomaly E to the mean anomaly M
    INPUTS:
        E: The eccentric anomaly
        e: The eccentricity
    OUTPUTS:
        M: The mean anomaly
    """
    # SSD equation 2.52
    M: np.array = E - e * np.sin(E)
    return M

# ********************************************************************************************************************* 
def danby_iteration(M: np.array, e: np.array, E: np.array) -> np.array:
    """
    Perform one iteration of the Danby algorithm for computing E from M
    See SSD equation 2.62 on page 36
    INPUTS:
        M:      The mean anomaly
        e:      The eccentricity
        E:      The current estimate of the eccentric anomaly E
    OUTPUTS:
        E_next: The improved estimate of the eccentric anomaly E

    """
    # The objective function that is set to zero using Newton-Raphson is
    # f(E) = E - e Sin E - M
    
    # Save two intermediate arrays that are reused
    eSinE = e * np.sin(E)
    eCosE = e * np.cos(E)
    
    # Save the value of the function and its first three derivatives
    f0 = E - eSinE - M
    f1 = 1.0 - eCosE
    f2 = eSinE
    f3 = eCosE

    # The three delta adjustments; see SSD Equation 2.62
    d1 = -f0 / f1
    d2 = -f0 / (f1 + 0.5*d1 * f2)
    d3 = -f0 / (f1 + 0.5*d2 * f2 + (1.0/6.0)*np.square(d2)*f3)
        
    E_next = E + d3
    return E_next

# ********************************************************************************************************************* 
def danby_guess(M: np.array, e: np.array) -> np.array:
    """
    Initial guess E0 for iterative calculation of E from M
    See SSD equation 2.64 on page 36
    """
    k = 0.85
    E0 = M + np.sign(np.sin(M))*k*e
    return E0

# ********************************************************************************************************************* 
def anomaly_M2E(M: np.array, e: np.array) -> np.array:
    """
    Convert the mean anomaly M to the eccentric anomaly E using Danby iterations
    INPUTS:
        M: The mean anomaly
        e: The eccentricity
    OUTPUTS:
        E: The eccentric anomaly
    """
    # The initial guess
    E0 = danby_guess(M=M, e=e)
    # Three iterations of Danby algorithm
    E1 = danby_iteration(M=M, e=e, E=E0)
    E2 = danby_iteration(M=M, e=e, E=E1)
    E3 = danby_iteration(M=M, e=e, E=E2)
    return E3

# ********************************************************************************************************************* 
# Convert from orbital elements to state vectors
# ********************************************************************************************************************* 

def elt2pos(a: np.array, e: np.array, inc: np.array, Omega: np.array, omega: np.array, f: np.array):
    """
    Convert from orbital elements to state vectors.
    See SSD page 51, equation 2.122.
    INPUTS:
        a:          semimajor axis
        e:          eccentricity
        inc:        inclination
        Omega:      longitude of ascending node
        omega:      argument of percenter
        f:          true anomaly
    OUTPUTS:
        q:          position; shaped same as elements with one additional axis
    """

    # Calculate the distance from the center, r; SSD equation 2.20
    r = a * (1.0 - np.square(e)) / (1 + e * np.cos(f))

    # Calculate intermediate results used for angular rotations

    # The angle in the elliptic plane, measured from the reference direction
    theta = omega + f 
    # Trigonometric functions of the angles
    cos_inc = np.cos(inc)
    sin_inc = np.sin(inc)
    cos_Omega = np.cos(Omega)
    sin_Omega = np.sin(Omega)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    # The cartesian position coordinates; see SSD equation 2.122
    qx = r * (cos_Omega*cos_theta - sin_Omega*sin_theta*cos_inc)
    qy = r * (sin_Omega*cos_theta + cos_Omega*sin_theta*cos_inc)
    qz = r * (sin_theta*sin_inc)
    q = np.stack([qx, qy, qz], axis=1)

    return q

# ********************************************************************************************************************* 
def elt2vec(a: np.array, e: np.array, inc: np.array, Omega: np.array, omega: np.array, f: np.array):
    """
    Convert from orbital elements to state vectors.
    See SSD page 51, equation 2.122.
    INPUTS:
        a:          semimajor axis
        e:          eccentricity
        inc:        inclination
        Omega:      longitude of ascending node
        omega:      argument of percenter
        f:          true anomaly
    OUTPUTS:
        qx:         position; x
        qy:         position; y
        qz:         position; z
        vx:         velocity; x
        vy:         velocity; y
        vz:         velocity; z
    """

    # Calculate the distance from the center, r; SSD equation 2.20
    r = a * (1.0 - np.square(e)) / (1 + e * np.cos(f))

    # Calculate intermediate results used for angular rotations
    cos_inc = np.cos(inc)
    sin_inc = np.sin(inc)
    cos_Omega = np.cos(Omega)
    sin_Omega = np.sin(Omega)
    cos_omega = np.cos(omega)
    sin_omega = np.sin(omega)


# ********************************************************************************************************************* 
# Convert from state vectors to orbital elements
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
# Testing
# ********************************************************************************************************************* 

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
def test_all():
    """Running test suite on orbital elements"""
    test_E2f()
    test_E2M()
    test_M2E()

# ********************************************************************************************************************* 
if __name__ == '__main__':
    print(f'Running test suite on orbital_element.py...')
    test_all()