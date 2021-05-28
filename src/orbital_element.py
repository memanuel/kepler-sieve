"""
Utilities for working with orbital elements.
Reference: Solar System Dynamics, C.D. Murray and S.F. Dermott, 1999

Michael S. Emanuel
2021-01-29
"""

# Core
import numpy as np
from scipy import interpolate

# Utility
from collections import namedtuple

# Typing
from typing import Optional

# ********************************************************************************************************************* 
# Named tuple data type for orbital elements with just the six 
# a, e, inc, Omega, omega, f
OrbitalElement_aeiOof = namedtuple('OrbitalElement', 'a e inc Omega omega f')

# Named tuple data type for orbital elements that includes the mean anomaly M
# a, e, inc, Omega, omega, f, M
OrbitalElement_aeiOofM = namedtuple('OrbitalElement', 'a e inc Omega omega f M')

# Name the constant for radians in a circle (two pi)
tau = 2.0 * np.pi

# The gravitational constant in unit system (years, AU, Msun)
# numerical value close to 4 pi^2; see rebound documentation for exact value        
G_ = 39.476926421373
mu = G_ * 1.0

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
def anomaly_M2E_danby(M: np.array, e: np.array, n: int=3) -> np.array:
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
def anomaly_M2E(M: np.array, e: np.array) -> np.array:
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
# Convert from orbital elements to state vectors
# ********************************************************************************************************************* 

def elt2pos(a: np.array, e: np.array, inc: np.array, Omega: np.array, omega: np.array, f: np.array):
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
                    This position is relative to the primary, i.e. the sun for heliocentric elements
    """

    # Calculate the distance from the center, r; SSD equation 2.20
    r: np.array = a * (1.0 - np.square(e)) / (1.0 + e * np.cos(f))

    # Calculate intermediate results used for angular rotations
    # The angle in the elliptic plane, measured from the reference direction
    theta: np.array = omega + f 
    # Trigonometric functions of the angles
    cos_inc: np.array = np.cos(inc)
    sin_inc: np.array = np.sin(inc)
    cos_Omega: np.array = np.cos(Omega)
    sin_Omega: np.array = np.sin(Omega)
    cos_theta: np.array = np.cos(theta)
    sin_theta: np.array = np.sin(theta)

    # The cartesian position coordinates; see SSD equation 2.122
    qx: np.array = r * (cos_Omega*cos_theta - sin_Omega*sin_theta*cos_inc)
    qy: np.array = r * (sin_Omega*cos_theta + cos_Omega*sin_theta*cos_inc)
    qz: np.array = r * (sin_theta*sin_inc)
    q: np.array = np.stack([qx, qy, qz], axis=1)

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
    one_plus_e_cos_f = 1.0 + e * np.cos(f)
    r = a * one_minus_e2 / one_plus_e_cos_f
    
    # Current speed
    v0 = np.sqrt(mu / a / one_minus_e2)

    # Position
    # qx = r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
    # qy = r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
    # qz = r*(so*cf+co*sf)*si
    # the term cos_omega*cos_f - sin_omega*sin_f appears 2 times
    # the term sin_omega*cos_f + cos_omega*sin_f appears 3 times
    cocf_sosf = co*cf-so*sf
    socf_cosf = so*cf+co*sf
    qx = r*(cO*cocf_sosf - sO*socf_cosf*ci)
    qy = r*(sO*cocf_sosf + cO*socf_cosf*ci)
    qz = r*socf_cosf*si
    
    # Velocity
    # vx = v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
    # vy = v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
    # vz = v0*((e+cf)*co*si - sf*si*so)
    # The term e+cf appears three times
    epcf = e + cf
    # vx = v0*(epcf*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
    # vy = v0*(epcf*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
    # vz = v0*(epcf*co*si - sf*si*so)
    # The term cocO appears twice
    cocO = co * cO
    # The term cosO appears twice
    cosO = co * sO
    # The term so*sO appears twice
    sosO = so * sO
    # The terms socO appears twice
    socO = so*cO

    # Simplified expression for velocity with substitutions
    vx = v0*(epcf*(-ci*cosO - socO) - sf*(cocO - ci*sosO))
    vy = v0*(epcf*(ci*cocO - sosO)  - sf*(cosO + ci*socO))
    vz = v0*(epcf*co*si - sf*si*so)

# ********************************************************************************************************************* 
# Convert from state vectors to orbital elements
# ********************************************************************************************************************* 