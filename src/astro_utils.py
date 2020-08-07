"""
Harvard IACS Masters Thesis
Astronomy Utilities

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

# Core
import numpy as np
from scipy import interpolate
from datetime import date, datetime, timedelta

# Astronomy
import astropy
from astropy.units import deg, au, km, meter, day, minute, second
from astropy.coordinates import SkyCoord, ICRS, GCRS, BarycentricMeanEcliptic

# Typing
from typing import Tuple, Optional

# *************************************************************************************************
# Constant with the base date for julian day conversions
julian_base_date: date = date(1899,12,31)
julian_base_datetime: datetime = datetime(1899, 12, 31, 12, 0, 0)
# The Julian date of 1899-12-31 12:00  is 2415020; this is the "Dublin JD"
# https://en.wikipedia.org/wiki/Julian_day
julian_base_number: int = 2415020

# The modified Julian Date is 2400000.5 less than the Julian Base Number
modified_julian_offset: float = 2400000.5

# Equivalently, it is the number of days from the epoch beginning 1858-11-178 0:00 (midnight)
# http://scienceworld.wolfram.com/astronomy/ModifiedJulianDate.html
modified_julian_base_date: date = date(1858, 11, 17)
modified_julian_base_datetime: datetime = datetime(1858, 11, 17, 0, 0, 0)

# Number of seconds in one day
day2sec: float = 24.0 * 3600.0
sec2day: float = 1.0 / day2sec

# Radians in a circle
tau = 2.0 * np.pi

# *************************************************************************************************
def date_to_jd(t: date) -> int:
   """Convert a Python date to a Julian day"""
   # Compute the number of days from Julian Base Date to date t
   dt = t - julian_base_date
   # Add the julian base number to the number of days from the julian base date to date t
   return julian_base_number + dt.days

# *************************************************************************************************
def date_to_mjd(t: date) -> int:
   """Convert a Python datetime to a Modified Julian day"""
   # Compute the number of days from Julian Base Date to date t
   dt = t - modified_julian_base_date
   return dt.days

# *************************************************************************************************
def jd_to_date(jd: int) -> date:
   """Convert an integer julian date to a Python date"""
   dt = timedelta(days=jd - julian_base_number)
   return julian_base_date + dt

# *************************************************************************************************
def mjd_to_date(mjd: int) -> date:
   """Convert an integer modified julian date to a Python date"""
   dt = timedelta(days=mjd)
   return modified_julian_base_date + dt

# *************************************************************************************************
def datetime_to_jd(t: datetime) -> float:
    """Convert a Python datetime to a Julian day"""
    # Compute the number of days from Julian Base Date to date t
    dt = t - julian_base_datetime
    # Add the julian base number to the number of days from the julian base date to date t
    return julian_base_number + dt.days + sec2day * dt.seconds 

# *************************************************************************************************
def datetime_to_mjd(t: datetime, epoch = modified_julian_base_datetime) -> float:
    """Convert a Python datetime to a Modified Julian day"""
    # Compute the number of days from Julian Base Date to date t
    dt = t - epoch
    return dt.days + sec2day * dt.seconds

# *************************************************************************************************
def jd_to_datetime(jd: float) -> date:
    """Convert a floating point julian date to a Python datetime"""
    interval = jd - julian_base_number
    dt = timedelta(seconds=day2sec*interval)
    return julian_base_datetime + dt

# *************************************************************************************************
def mjd_to_datetime(mjd: float) -> date:
    """Convert an integer modified julian date to a Python date"""
    dt = timedelta(seconds=day2sec*mjd)
    return modified_julian_base_datetime + dt

# *************************************************************************************************
def jd_to_mjd(jd: float) -> date:
    """Convert a floating point julian date to a modified julian date"""
    return jd - modified_julian_offset

# *************************************************************************************************
def mjd_to_jd(mjd: float) -> date:
    """Convert a floating point modified julian date to a julian date"""
    return mjd + modified_julian_offset

# ********************************************************************************************************************* 
def dist2rad(dist):
    """Convert a cartesian distance on unit sphere in [0, 2] to radians in [0, pi]"""
    x_rad = np.arcsin(0.5 * dist) * 2.0
    return x_rad

# ********************************************************************************************************************* 
def rad2dist(x_rad):
    """Convert a distance on unit sphere from radians in [0, pi] to cartesian distance in [0, 2]"""
    return np.sin(0.5 * x_rad) * 2.0

# ********************************************************************************************************************* 
def dist2deg(dist):
    """Convert a cartesian distance on unit sphere in [0, 2] to degrees in [0, 180]"""
    x_rad = dist2rad(dist)
    return np.rad2deg(x_rad)

# ********************************************************************************************************************* 
def deg2dist(x_deg):
    """Convert a distance on unit sphere from degrees in [0, 180] to cartesian distance in [0, 2]"""
    x_rad = np.deg2rad(x_deg)
    return rad2dist(x_rad)

# ********************************************************************************************************************* 
def dist2sec(dist):
    """Convert a cartesian distance on unit sphere in [0, 2] to arc seconds in [0, 180*3600]"""
    x_rad = dist2rad(dist)
    return np.rad2deg(x_rad) * 3600.0

# *************************************************************************************************
def xyz_to_sph(x: np.array, y: np.array, z: np.array):
    """
    Convert a Cartesian coordinates x, y, z of a displacement vector to 
    spherical coordinates r, alt, az
    Used only for error checking, not RA/DEC calculations.
    See ra_dec.py for conversions between Cartesian and RA/DEC coordinates.
    """
    # The distance R
    r = np.sqrt(x*x + y*y + z*z)

    # The azimuth
    az = np.arctan2(y, x)

    # The altitude; use mask to avoid divide by zero when r=0
    alt = np.zeros_like(z)
    mask = r>0
    alt[mask] = np.arcsin(z[mask] / r[mask])

    return r, alt, az

# *************************************************************************************************
def cart_to_sph(q: np.array):
    """
    Convert a Cartesian coordinates q with shape (N,3) to spherical coordinates r, alt, az"""
    # Unpack x, y, z
    x = q[:, 0]
    y = q[:, 1]
    z = q[:, 2]
    # Delegate to xyz_to_sph
    return xyz_to_sph(x, y, z)

# ********************************************************************************************************************* 
def reverse_velocity(sim):
    """Reverse the velocities in a simulation for backwards time integration; modifies sim in place"""
    for p in sim.particles:
        vx, vy, vz = p.vx, p.vy, p.vz
        p.vx = -vx
        p.vy = -vy
        p.vz = -vz

# *************************************************************************************************
def anomaly_M2E_impl(M: np.ndarray, e: np.ndarray, E0: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute the eccentric anomaly E from the mean anomaly M and eccentricity e using and Newton's Method.
    This is the implementation that does not depend on a table of initial guesses
    See: https://en.wikipedia.org/wiki/Eccentric_anomaly
    INPUTS:
        M: The mean anomaly
        e: The eccentricity
        E0: Initial guess
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
    
    # Use the initial guess E0 if provided; otherwise use M
    E = E0 if E0 is not None else M.copy()

    # Maximum number of iterations for Newton-Raphson
    max_iter: int = 40

    # Tolerance for maximum error
    err_tol: np.float64 = 2.0**-49

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

# *************************************************************************************************
def anomaly_M2E_make_interp(N_M: int = 256, N_e: int = 256) -> interpolate.RectBivariateSpline:
    """Create an interpolation function for inital guesses to compute E from M and e."""
    # Rows of M: N_M+1 evenly spaced angles in [0, 2 pi], e.g. [0, 45, 90, ... 315, 360] degrees for N_M = 8
    # Include both endpoints for stability in the interpolation
    M = tau * np.linspace(0.0, 1.0, N_M+1)
    # Repeat the rows of M N_e times
    M_grid = np.tile(M, (N_e, 1))

    # Columns of e: N_e log spaced eccentricities in [0, 1 - 2^-48]
    # Spacing gets closer together as eccentricity approaches 1.0; log(1.0 - e) is uniformly spaced
    log_1me = np.linspace(0.0, -48.0*np.log(2.0), N_e)
    e = 1.0 - np.exp(log_1me)
    # Repeat the columns of e N_m+1 times
    e_grid = np.tile(e, (N_M+1, 1)).T

    # Compute the eccentric anomaly E using the implementation version of anomaly_M2E
    E = anomaly_M2E_impl(M_grid, e_grid)
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
def anomaly_M2E(M: np.ndarray, e: np.ndarray):
    """
    Compute the eccentric anomaly E from the mean anomaly M and eccentricity e
    Uses interpolation table from anomaly_M2E_make_interp for initial guess and
    Newton's Method as implemented by anomaly_M2E_impl to refine guess to convergence.
    See: https://en.wikipedia.org/wiki/Eccentric_anomaly
    INPUTS:
        M: The mean anomaly
        e: The eccentricity
    OUTPUTS:
        E: The Eccentric anomaly
    """   
    # Compute initial guess using interpolation function
    E0 = interp_M2E.ev(M, e) % tau
    
    # Delegate to anomaly_M2E_impl for Newton's Method
    E = anomaly_M2E_impl(M, e, E0)
    
    return E

# *************************************************************************************************
def anomaly_E2f(E: np.ndarray, e: np.ndarray) -> np.ndarray:
    """
    Convert the eccentric anomaly E to the true anomaly f
    INPUTS:
        E: The Eccenctric anomaly
        e: The eccentricity
    OUTPUTS:
        f: The true anomaly    
    """
    # See https://en.wikipedia.org/wiki/Eccentric_anomaly
    # Use formula f = 2 arg(sqrt(1-e) cos(E/2), sqrt(1+e) sin(E/2))
    
    # Apply formula
    ecc_ratio = np.sqrt( (1.0 + e) / (1.0 - e))
    half_E = E * 0.5
    f = 2.0 * np.arctan(ecc_ratio * np.tan(half_E))
    # Shift angle returned to the interval [0, 2 pi)
    return f % tau
    
# *************************************************************************************************
# Testing
# *************************************************************************************************
def test_julian_day():
    """Test Julian Day conversions"""
    # Known conversion from Small-Body browser
    t = datetime(2019, 4, 27, 0)
    jd = 2458600.5
    mjd = 58600

    # Compute recovered time and julian date    
    t_rec = jd_to_datetime(jd)
    jd_rec = datetime_to_jd(t)
    mjd_rec = datetime_to_mjd(t)
    
    # Errors vs. known dates
    err_t = t_rec - t
    err_jd = jd_rec - jd
    err_mjd = mjd_rec - mjd
    
    # Did the test pass?
    isOK: bool = (err_t == timedelta(seconds=0)) and (err_jd == 0.0) and (err_mjd == 0.0)
    msg: str = 'PASS' if isOK else 'FAIL'
    
    # Test results to screen
    print('Known Synchronous Times from NASA Small Body Browser:')
    print(f't = {t}')
    print(f'jd = {jd}')
    print(f'mjd = {mjd}')
    # print(f't_rec = {t_rec}')
    # print(f'jd_rec = {jd_rec}')
    # print(f'jd_rec = {jd_rec}')
    print(f'Error in t: {err_t}')
    print(f'Error in jd: {err_jd}')
    print(f'Error in mjd: {err_mjd}')
    print(f'*** {msg} ***')
    
# *************************************************************************************************
def test_anomaly_M2E(N_test = 10**6):
    """Test conversion from mean anomaly M to eccentric anomaly E"""
    # Randomly generate test cases
    np.random.seed(42)
    M = tau * np.random.rand(N_test)
    e = np.random.rand(N_test)
    
    # Compute eccentric anomaly
    E = anomaly_M2E(M, e)
    
    # Recover M using Kepler's equation
    M_rec = E - e * np.sin(E)
    
    # Get max error between input and recovered mean anomaly M
    err = np.abs(M_rec - M)
    idx = np.argmax(err)
    max_err = err[idx]    
    
    # Check against tolerance
    tol = 1.0E-8
    isOK: bool = (max_err < tol)
    msg: str = 'PASS' if isOK else 'FAIL'
    
    # Test results to screen
    print('Conversion of mean anomaly M to eccentric anomaly E with anomaly_M2E:')
    print(f'Test {N_test} randomly generated input pairs (M, e).')
    print('Compare input and recovered M from Kepler''s Equation M = E - e sin E')
    print(f'Max error = {max_err:5.3e} = 2^ {np.log2(max_err):6.2f}')
    print(f'M = {M[idx]}, e = {e[idx]}, E = {E[idx]}')
    print(f'*** {msg} ***')

# *************************************************************************************************
def main():
    test_julian_day()
    test_anomaly_M2E()

if __name__ == '__main__':
    main()
