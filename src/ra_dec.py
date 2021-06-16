"""
Astronomy: Transformation of Right Ascension & Declination (RA/DEC) to Ecliptic Directions

Michael S. Emanuel
Fri Aug 23 16:13:28 2019

Functions in this module:
radec2dir(ra, dec, t_obs)
dir2radec(u, t_obs)
"""

# Core
import numpy as np

# Astronomy
import astropy
from astropy.units import au, deg
from astropy.coordinates import SkyCoord, ICRS, BarycentricMeanEcliptic

# MSE imports
from astro_utils import infer_shape

# Typing
from typing import Tuple

# *************************************************************************************************
# Convert between a direction as a RA/DEC and a unit vector (ux, uy, uz) in the B.M.E. frame
# *************************************************************************************************

# ********************************************************************************************************************* 
def radec2dir(ra: np.ndarray, dec: np.ndarray, mjd: np.ndarray) -> np.array:
    """
    Convert a RA and DEC as of observation time to a unit displacement vector u = (ux, uy, uz) in ecliptic plane.
    INPUTS:
        ra:     An astrometric Right Ascension in the ICRF; array in degrees
        dec:    An astromentric Declination in the ICRF;  array in degrees
        mjd:    The observation time as a modified julian day
    RETURNS:
        u:      An array [ux, uy, uz] on the unit sphere in the the ecliptic frame
    EXAMPLE:
        u = radec2dir(ra=76.107414227, dec=23.884882701, t_obs=58600.0)
        (this is Mars viewed from Earth at mjd 58600 / 2019-04-27 with JPL RA and DEC)
    """
    # Build the observation as a SkyCoord in the ICRS (oriented with earth, origin at barycenter)
    obstime = astropy.time.Time(mjd, format='mjd')
    # Add units of degrees to ra and dec
    ra_deg = ra * deg
    dec_deg = dec * deg
    # The observation as an astropy object in the ICRS
    obs_icrs = SkyCoord(ra=ra_deg, dec=dec_deg, obstime=obstime, frame=ICRS)
    # Convert to the barycentric ecliptic frame (oriented with ecliptic, origin at barycenter)
    obs_ecl = obs_icrs.transform_to(BarycentricMeanEcliptic)
    # Extract the direction vector
    u = obs_ecl.cartesian.xyz
    # Return as a numpy array of shape Nx3 
    # (astropy default shapes the array as 3xN)
    return np.moveaxis(u.value, source=0, destination=-1)

# ********************************************************************************************************************* 
def dir2radec(u: np.array, mjd: np.array) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space.
    INPUTS:
        u:      An array [ux, uy, uz] on the unit sphere in the the ecliptic frame
        mjd:    Observation time on Earth as a modified Julian day
    RETURNS:
        ra:     Right Ascension; with units (astropy degrees)
        dec:    Declination; with units (astropy degrees)
    """
    # The observation time
    obstime = astropy.time.Time(mjd, format='mjd')

    # infer data shape
    data_axis, space_axis, shape = infer_shape(u)
    # Unpack position and convert to 1 au distance
    if space_axis==0:
        x, y, z = u * au
    else:
        x, y, z = u.T * au
    
    # The observation in the BarycentricMeanEcliptic coordinate frame
    frame = BarycentricMeanEcliptic
    obs_frame = SkyCoord(x=x, y=y, z=z, obstime=obstime,  representation_type='cartesian', frame=frame)
    # Transform the observation to the ICRS frame
    obs_icrs = obs_frame.transform_to(ICRS)
    # Extract the RA and DEC in degrees
    ra = obs_icrs.ra.deg
    dec = obs_icrs.dec.deg
    # Return (ra, dec) as a tuple
    return (ra, dec)
