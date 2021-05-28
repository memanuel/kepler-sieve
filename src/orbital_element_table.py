"""
Utilities for working with orbital elements.
Reference: Solar System Dynamics, C.D. Murray and S.F. Dermott, 1999

Michael S. Emanuel
2021-01-29
"""

# Core
import numpy as np
from scipy import interpolate

# Local
from orbital_element import danby_iteration, danby_guess, anomaly_M2E_danby

# Name the constant for radians in a circle (two pi)
tau = 2.0 * np.pi

# ********************************************************************************************************************* 
# Functions for converting between anomalies: mean, true, and eccentric
# ********************************************************************************************************************* 

# *************************************************************************************************
def anomaly_M2E_newton(M: np.ndarray, e: np.ndarray) -> np.ndarray:
    """
    Compute the eccentric anomaly E from the mean anomaly M and eccentricity e using and Newton's Method.
    This is the implementation that does not depend on a table of initial guesses
    See: https://en.wikipedia.org/wiki/Eccentric_anomaly
    INPUTS:
        M: The mean anomaly
        e: The eccentricity
    OUTPUTS:
        E: The Eccentric anomaly

    Kepler's Equation gives us
        M = E - e sin(E)
    which implies that E - e sin(E) - M = 0
    Think of this as a function f(E) = E - e sin(E) - M
    with derivative f'(E) = 1 - e cos(E)
    we solve for f(E) = 0 using Newton's Method
    """

    # Put M in the interval [0, 2*pi)
    M %= tau

    # Maximum number of iterations for Newton-Raphson
    max_iter: int = 40

    # Tolerance for maximum error
    err_tol: np.float64 = 2.0**-49

    # Initial guess is just a copy of M
    k = 0.85
    E = danby_guess(M=M, e=e)

    # Perform at most max_iter iterations of Newton's method; quit early if tolerance achieved    
    for i in range(max_iter):
        # The current function value f(E)
        f = E - e * np.sin(E) - M
        # Is the max error below the tolerance? If so, quit early
        max_err = np.max(np.abs(f))
        if max_err < err_tol:
            # print(f'Converged with error {max_err:5.2e} after {i+1} iterations.')
            break
        # The derivative f'(E)
        fp = 1.0 - e * np.cos(E)
        # Update E using Newton's method
        E -= f / fp
        # Restore E to range [0, 2 pi) to handle nasty corner cases where very high eccentricity
        # and high M cause a naive iteration to diverge!
        E %= tau

    # Return the (hopefully) converged eccentric anomaly E
    return E

# ********************************************************************************************************************* 
def anomaly_M2E_make_interp(N_M: int = 256, N_e: int = 256) -> interpolate.RectBivariateSpline:
    """Create an interpolation function for inital guesses to compute E from M and e."""
    # Rows of M: N_M+1 evenly spaced angles in [0, 2 pi], e.g. [0, 45, 90, ... 315, 360] degrees for N_M = 8
    # Include both endpoints for stability in the interpolation
    M = tau * np.linspace(0.0, 1.0, N_M+1)
    # Repeat the rows of M N_e times
    M_grid = np.tile(M, (N_e, 1))

    # Columns of e: N_e log spaced eccentricities in [0, 1 - 2^-48]
    # Spacing gets closer together as eccentricity approaches 1.0; log(1.0 - e) is uniformly spaced
    # log_1me = np.linspace(0.0, -48.0*np.log(2.0), N_e)
    log_1me = np.linspace(0.0, -6.0*np.log(2.0), N_e) 
    e = 1.0 - np.exp(log_1me)
    # Repeat the columns of e N_m+1 times
    e_grid = np.tile(e, (N_M+1, 1)).T

    # Compute the eccentric anomaly E using the implementation version of anomaly_M2E
    E = anomaly_M2E_danby(M_grid, e_grid)
    # E = anomaly_M2E_newton(M_grid, e_grid)
    # Overwrite the last column, for M= 360 degrees, with E= 360 degrees.  
    # This improves interpolation accuracy for angles very close to 360 degrees
    E[:, -1] = tau
    # Take transpose so shape is N_M x N_e to match API of interpolate.RectBivariateSpline
    E = E.T

    # Create a 2D interpolator; degree 5 for M, degree 3 for e
    interp = interpolate.RectBivariateSpline(x=M, y=e, z=E, kx=5, ky=3)

    return interp

# *************************************************************************************************
# Create one copy of the interpolation function for the module
interp_M2E: interpolate.RectBivariateSpline = anomaly_M2E_make_interp(N_M=2**12, N_e=2**8)

# *************************************************************************************************
def anomaly_M2E_table(M: np.ndarray, e: np.ndarray):
    """
    Compute the eccentric anomaly E from the mean anomaly M and eccentricity e
    Uses interpolation table from anomaly_M2E_make_interp for initial guess and Danby iteration
    INPUTS:
        M: The mean anomaly
        e: The eccentricity
    OUTPUTS:
        E: The Eccentric anomaly
    """   
    # Compute initial guess using interpolation function
    E0 = interp_M2E.ev(M, e)

    # Delegate to danby_iteration for Newton's Method
    E1 = danby_iteration(M, e, E0)
    E2 = danby_iteration(M, e, E1)
    
    return E2
