"""
Utilities for working with orbital elements.
Reference: Solar System Dynamics, C.D. Murray and S.F. Dermott, 1999

Michael S. Emanuel
2021-01-29
"""

# Core
import numpy as np

# Utility
from collections import namedtuple

# Typing
from typing import Tuple

# ********************************************************************************************************************* 
# Named tuple data type for orbital elements with just the six 
# a, e, inc, Omega, omega, f
OrbitalElement_aeiOof = namedtuple('OrbitalElement', 'a e inc Omega omega f')

# Named tuple data type for orbital elements that includes the mean anomaly M
# a, e, inc, Omega, omega, f, M
OrbitalElement_aeiOofM = namedtuple('OrbitalElement', 'a e inc Omega omega f M')

# Name the constant for radians in a circle (two pi)
tau = 2.0 * np.pi

# The gravitational constant in unit system (days, AU, Msun)
# see rebound documentation for exact value        
# sim = make_sim_planets(epoch=59000)
# G_ = sim.G_
G_ = 2.959122082855910945e-04

# mu is the gravitational field strength: mu = G * (m0 + m1)
# here m0 = 1.0 (Sun) and we assume m1 is light, i.e. m0 = 0.0
mu = G_ * 1.0

# ********************************************************************************************************************* 
# Pack and unpack arrays
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def unpack_elt_df(df: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Unpack a DataFrame of orbital elements
    INPUTS:
        df:     A DataFrame of orbital elements
                Must include columns named a, e, inc, Omega, omega, f
    OUTPUTS:
        Tuple of 6 numpy arrays for these elements.
    """
    a = df.a.values
    e = df.e.values
    inc = df.inc.values
    Omega = df.Omega.values
    omega = df.omega.values
    f = df.f.values
    return a, e, inc, Omega, omega, f

# ********************************************************************************************************************* 
def unpack_vector(v: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unpack three spatial indices of a vector"""
    vx = v[:, 0]
    vy = v[:, 1]
    vz = v[:, 2]
    return vx, vy, vz

# ********************************************************************************************************************* 
# Functions for converting between anomalies: mean, true, and eccentric
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def anomaly_E2f(E: np.ndarray, e: np.ndarray) -> np.ndarray:
    """
    Convert the eccentric anomaly E to the true anomaly f
    INPUTS:
        E: The eccentric anomaly
        e: The eccentricity
    OUTPUTS:
        f: The true anomaly
    """
    # SSD equation 2.46; solve for f in terms of E
    tan_half_f: np.ndarray = np.sqrt((1.0 + e) / (1.0 - e)) * np.tan(E/2.0)
    f: np.ndarray = np.arctan(tan_half_f) * 2.0
    return f

# ********************************************************************************************************************* 
def anomaly_f2E(f: np.ndarray, e: np.ndarray) -> np.ndarray:
    """
    Convert the true anomaly f to the eccentric anomaly E
    INPUTS:
        f: The true anomaly
        e: The eccentricity
    OUTPUTS:
        E: The eccentric anomaly
    """
    # SSD equation 2.46; solve for E in terms of f
    tan_half_f: np.ndarray = np.tan(f * 0.5)
    tan_half_E: np.ndarray = tan_half_f / np.sqrt((1.0 + e) / (1.0 - e))
    E: np.ndarray = np.arctan(tan_half_E) * 2.0
    return E

# ********************************************************************************************************************* 
def anomaly_E2M(E: np.ndarray, e: np.ndarray) -> np.ndarray:
    """
    Convert the eccentric anomaly E to the mean anomaly M
    INPUTS:
        E: The eccentric anomaly
        e: The eccentricity
    OUTPUTS:
        M: The mean anomaly
    """
    # SSD equation 2.52; this is Kepler's Equation
    M: np.ndarray = E - e * np.sin(E)
    return M

# ********************************************************************************************************************* 
def anomaly_f2M(f: np.ndarray, e: np.ndarray) -> np.ndarray:
    """
    Convert the true anomaly f to the mean anomaly M
    INPUTS:
        f: The true anomaly
        e: The eccentricity
    OUTPUTS:
        E: The eccentric anomaly
    """
    # Delegate to anomaly_f2E
    E = anomaly_f2E(f=f, e=e)
    # Delegate to anomaly_E2M
    M = anomaly_E2M(E=E, e=e)
    return M

# ********************************************************************************************************************* 
def danby_iteration(M: np.ndarray, e: np.ndarray, E: np.ndarray) -> np.ndarray:
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
    # This is just the error term if we take Kepler's Equation and subtract M from both sides
    # This function will be zero at the correct value of E
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
    d2 = -f0 / (f1 + 0.5*d1*f2)
    d3 = -f0 / (f1 + 0.5*d2*f2 + (1.0/6.0)*np.square(d2)*f3)

    E_next = E + d3
    return E_next

# ********************************************************************************************************************* 
def danby_guess(M: np.ndarray, e: np.ndarray) -> np.ndarray:
    """
    Initial guess E0 for iterative calculation of E from M
    See SSD equation 2.64 on page 36
    """
    k = 0.85
    E0 = M + np.sign(np.sin(M))*k*e
    return E0

# ********************************************************************************************************************* 
def anomaly_M2E_danby(M: np.ndarray, e: np.ndarray, n: int=3) -> np.ndarray:
    """
    Convert the mean anomaly M to the eccentric anomaly E using Danby iterations
    INPUTS:
        M: The mean anomaly
        e: The eccentricity
    OUTPUTS:
        E: The eccentric anomaly
    """
    # The initial guess
    E = danby_guess(M=M, e=e)
    # n iterations of Danby algorithm
    for k in range(n):
        E = danby_iteration(M=M, e=e, E=E)
    return E

# ********************************************************************************************************************* 
def anomaly_M2E(M: np.ndarray, e: np.ndarray) -> np.ndarray:
    """
    Convert the mean anomaly M to the eccentric anomaly E.
    Implemention can be changed behind the scenes later.
    INPUTS:
        M: The mean anomaly
        e: The eccentricity
    OUTPUTS:
        E: The eccentric anomaly
    """
    # Delegate to anomaly_M2E_danby with n=3
    E = anomaly_M2E_danby(M=M, e=e, n=3)
    return E

# ********************************************************************************************************************* 
def anomaly_M2f(M: np.ndarray, e: np.ndarray) -> np.ndarray:
    """
    Convert the mean anomaly M to the true anomaly f.
    Done by converting first to E, then to f.
    INPUTS:
        M: The mean anomaly
        e: The eccentricity
    OUTPUTS:
        E: The eccentric anomaly
    """
    # Delegate to anomaly_M2E
    E = anomaly_M2E(M=M, e=e)
    # Delegate to anomaly_M2f
    f = anomaly_E2f(E=E, e=e)
    return f

# ********************************************************************************************************************* 
# Convert from orbital elements to state vectors
# ********************************************************************************************************************* 

def elt2pos(a: np.ndarray, e: np.ndarray, inc: np.ndarray, Omega: np.ndarray, omega: np.ndarray, f: np.ndarray) -> np.ndarray:
    """
    Convert from orbital elements to position vector.
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
    The position q is relative to the primary, i.e. the sun for heliocentric elements
    """

    # Calculate the distance from the center, r; SSD equation 2.20
    r: np.ndarray = a * (1.0 - np.square(e)) / (1.0 + e * np.cos(f))

    # Calculate intermediate results used for angular rotations
    # The angle in the elliptic plane, measured from the reference direction
    theta: np.ndarray = omega + f 
    # Trigonometric functions of the angles
    cos_inc: np.ndarray = np.cos(inc)
    sin_inc: np.ndarray = np.sin(inc)
    cos_Omega: np.ndarray = np.cos(Omega)
    sin_Omega: np.ndarray = np.sin(Omega)
    cos_theta: np.ndarray = np.cos(theta)
    sin_theta: np.ndarray = np.sin(theta)

    # The cartesian position coordinates; see SSD equation 2.122
    qx: np.ndarray = r * (cos_Omega*cos_theta - sin_Omega*sin_theta*cos_inc)
    qy: np.ndarray = r * (sin_Omega*cos_theta + cos_Omega*sin_theta*cos_inc)
    qz: np.ndarray = r * (sin_theta*sin_inc)
    q: np.ndarray = np.stack([qx, qy, qz], axis=-1)

    return q

# ********************************************************************************************************************* 
def elt2vec(a: np.ndarray, e: np.ndarray, inc: np.ndarray, Omega: np.ndarray, omega: np.ndarray, f: np.ndarray) \
            -> Tuple[np.ndarray, np.ndarray]:
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
        v:          velocity; shaped same as elements with one additional axis
    Both q and v are relative to the primary;
    So for heliocentric elements, they give RELATIVE position and velocity to the sun.
    """

    # The equations used here were taken from the rebound library
    # The position calculations are equivalent to the ones above from SSD.
    # The velocity calculations are a bit more involved, and I did not see them with explicit equations in SSD.

    # sine and cosine of the angles inc, Omega, omega, and f
    ci = np.cos(inc)
    si = np.sin(inc)
    cO = np.cos(Omega)
    sO = np.sin(Omega)
    co = np.cos(omega)
    so = np.sin(omega)
    cf = np.cos(f)
    sf = np.sin(f)

    # Distance from center
    one_minus_e2 = 1.0 - np.square(e)
    one_plus_e_cos_f = 1.0 + e*np.cos(f)
    r = a * one_minus_e2 / one_plus_e_cos_f

    # Current speed
    v0 = np.sqrt(mu / a / one_minus_e2)

    # Position
    # qx = r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
    # qy = r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
    # qz = r*(so*cf+co*sf)*si
    # the term cos_omega*cos_f - sin_omega*sin_f appears 2 times
    # the term sin_omega*cos_f + cos_omega*sin_f appears 3 times
    cocf_sosf = co*cf - so*sf
    socf_cosf = so*cf + co*sf
    qx = r*(cO*cocf_sosf - sO*socf_cosf*ci)
    qy = r*(sO*cocf_sosf + cO*socf_cosf*ci)
    qz = r*(socf_cosf*si)
    # Wrap the position vector
    q: np.ndarray = np.stack([qx, qy, qz], axis=-1)
    
    # Velocity
    # vx = v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
    # vy = v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
    # vz = v0*((e+cf)*co*si - sf*si*so)
    # The term e+cf appears three times
    epcf = e + cf
    # The term cocO appears twice
    cocO = co*cO
    # The term cosO appears twice
    cosO = co*sO
    # The term so*sO appears twice
    sosO = so*sO
    # The terms socO appears twice
    socO = so*cO

    # Simplified expression for velocity with substitutions
    vx = v0*(epcf*(-ci*cosO - socO) - sf*(cocO - ci*sosO))
    vy = v0*(epcf*(ci*cocO - sosO)  - sf*(cosO + ci*socO))
    vz = v0*(epcf*co*si - sf*si*so)
    # Wrap the velocity vector
    v: np.ndarray = np.stack([vx, vy, vz], axis=-1)

    # Return tuple of position and velocity
    return q, v

# ********************************************************************************************************************* 
# Convert from state vectors to orbital elements
# ********************************************************************************************************************* 
