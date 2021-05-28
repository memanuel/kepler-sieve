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
def anomaly_M2E(M: np.array, e: np.array) -> np.array:
    """
    Convert the mean anomaly M to the eccentric anomaly E
    INPUTS:
        M: The mean anomaly
        e: The eccentricity
    OUTPUTS:
        E: The eccentric anomaly
    """
    # SSD equation 2.63
    # TODO: implement this
    raise NotImplementedError

# ********************************************************************************************************************* 
# Convert from orbital elements to configuration vectors
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
# Convert from configuration vectors to orbital elements
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
    Calculate a distance between two angles using one minus the cosine as a metric
    This distance is 0.0 for angles on top of each other (up to 2 pi) 
    and a maximum of 2 for angles that are opposite (pi apart)
    INPUTS:
        x:      The first angle to be compared
        y:      The second angle to be compared
    OUTPUTS:
        r:      The distance 1 - cos(x-y)
    """
    # The difference in the angles
    d = x - y
    # Use one minus cosine
    return 1.0 - np.cos(d)

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
    print(f'Test: {test_name}')
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
