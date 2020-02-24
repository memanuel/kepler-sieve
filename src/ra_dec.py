"""
Harvard IACS Masters Thesis
Astronomy: Transformation of Right Ascension & Declination (RA/DEC) to Ecliptic Directions

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

import numpy as np
import astropy
from astropy.units import deg, au, km, meter, day, minute, second, arcsec
from astropy.coordinates import SkyCoord, ICRS, GCRS, BarycentricMeanEcliptic
from datetime import date, datetime, timedelta
from typing import Tuple

# Constants for distance of zero AU and velocity of zero AU / day
zero_au = 0.0 * au
zero_au_day = 0.0 * au / day
zero_km_sec = 0.0 * km / second
# Speed of light; express this in AU / minute
light_speed = astropy.constants.c.to(au / minute)

# ********************************************************************************************************************* 
def radec2dir(ra: float, dec: float, obstime_mjd: float) -> np.array:
    """
    Convert a RA and DEC as of observation time to a unit displacement vector u = (ux, uy, uz) in ecliptic plane.
    INPUTS:
        ra: An astrometric Right Ascension in the ICRF; passed with units (default degrees)
        dec: An astromentric Declination in the ICRF;  passed with units (default degrees)
        obtime_mjd: The observation time as a modified julian day
    RETURNS:
        u: An array [ux, uy, uz] on the unit sphere in the the ecliptic frame
    EXAMPLE:
        u = radec2dir(ra=76.107414227*deg, dec=23.884882701*deg, obstime_mjd=58600.0)
        (this is Mars viewed from Earth at mjd 58600 / 2019-04-27 with JPL RA and DEC)
    """
    # Build the observation as a SkyCoord in the ICRS (oriented with earth, origin at barycenter)
    obstime = astropy.time.Time(obstime_mjd, format='mjd')
    obs_icrs = SkyCoord(ra=ra, dec=dec, obstime=obstime, frame=ICRS)
    # Convert to the barycentric ecliptic frame (oriented with ecliptic, origin at barycenter)
    obs_ecl = obs_icrs.transform_to(BarycentricMeanEcliptic)
    u = obs_ecl.cartesian.xyz
    return u.value

# ********************************************************************************************************************* 
def radec_app2dir(ra: float, dec: float, obstime_mjd: float) -> np.array:
    """
    Convert an apparent RA and DEC as of observation time to a unit displacement 
    vector u = (ux, uy, uz) in ecliptic plane.
    INPUTS:
        ra: An apparent Right Ascension in the ICRF; passed with units (default degrees)
        dec: An apparent Declination in the ICRF;  passed with units (default degrees)
        obtime_mjd: The observation time as a modified julian day
    RETURNS:
        u: An array [ux, uy, uz] on the unit sphere in the the ecliptic frame
    EXAMPLE:
        u = radec_app2dir(ra=76.391533*deg, dec=23.90881*deg, obstime_mjd=58600.0)
        (this is Mars viewed from Earth at mjd 58600 / 2019-04-27 with JPL RA and DEC)
    """
    # this is a placeholder function declaration.
    # probably don't need this
    pass

# ********************************************************************************************************************* 
def qv2dir(q_body: np.ndarray, v_body: np.ndarray, q_earth: np.ndarray) -> np.ndarray:
    """
    Compute the direction of displacement from earth to a space body as a unit displacement vector
    u = (ux, uy, uz) in the ecliptic plane.
    INPUTS:
        q_body: position of this body in ecliptic coordinate frame; passed with units (default AU)
        v_body: velocity of this body in ecliptic coordinate frame; passed with units (default AU / day)
        q_earth: position of earth in ecliptic coordinate frame; passed with units (default AU)
    RETURNS:
        u: An array [ux, uy, uz] on the unit sphere in the ecliptic frame
    EXAMPLE:
        u = qv2dir(q_body=np.array([-0.328365, 1.570624, 0.040733])*au, 
                   v_body=np.array([-0.013177, -0.001673, 0.000288])*au/day,
                   q_earth=np.array([-0.813785, -0.586761, -0.000003])*au)
    """
    # displacement from earth to body; in AU
    q_rel = (q_body - q_earth)

    # distance; in AU
    r = np.linalg.norm(q_rel, axis=0) * au
    
    # light time in minutes
    light_time = (r / light_speed)
    
    # adjusted relative position, accounting for light time
    dq_lt = light_time * v_body.to(au/minute)     # convert velocity to au / minute b/c time is in minutes
    q_rel_lt = q_rel - dq_lt
    
    # adjusted direction
    r_lt = np.linalg.norm(q_rel_lt, axis=0) * au
    u = q_rel_lt / r_lt
    return u.value

# ********************************************************************************************************************* 
def dir2radec(u: np.array, obstime_mjd: np.array) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space.
    INPUTS:
        u: An array [ux, uy, uz] on the unit sphere in the the ecliptic frame
        obstime_mjd: Observation time on Earth as a modified Julian day
    RETURNS:
        ra: Right Ascension; with units (astropy degrees)
        dec: Declination; with units (astropy degrees)
    """
    # The observation time
    obstime = astropy.time.Time(obstime_mjd, format='mjd')

    # Unpack position and convert to 1 au distance
    x, y, z = u * au
    
    # The observation in the given coordinate frame (usually BarycentricMeanEcliptic)
    frame = BarycentricMeanEcliptic
    obs_frame = SkyCoord(x=x, y=y, z=z, obstime=obstime,  representation_type='cartesian', frame=frame)
    # Transform the observation to the ICRS frame
    obs_icrs = obs_frame.transform_to(ICRS)
    # Extract the RA and DEC in degrees
    ra = obs_icrs.ra.deg * deg
    dec = obs_icrs.dec.deg * deg
    # Return (ra, dec) as a tuple
    return (ra, dec)

# *************************************************************************************************
def dir2radec(u: np.array, obstime_mjd: np.array) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space.
    INPUTS:
        u: An array [ux, uy, uz] on the unit sphere in the the ecliptic frame
        obstime_mjd: Observation time on Earth as a modified Julian day
    RETURNS:
        ra: Right Ascension; with units (astropy degrees)
        dec: Declination; with units (astropy degrees)
    """
    # The observation time
    obstime = astropy.time.Time(obstime_mjd, format='mjd')

    # Unpack position and convert to 1 au distance
    x, y, z = u * au
    
    # The observation in the given coordinate frame (usually BarycentricMeanEcliptic)
    frame = BarycentricMeanEcliptic
    obs_frame = SkyCoord(x=x, y=y, z=z, obstime=obstime,  representation_type='cartesian', frame=frame)
    # Transform the observation to the ICRS frame
    obs_icrs = obs_frame.transform_to(ICRS)
    # Extract the RA and DEC in degrees
    ra = obs_icrs.ra.deg * deg
    dec = obs_icrs.dec.deg * deg
    # Return (ra, dec) as a tuple
    return (ra, dec)

# *************************************************************************************************
def radec_diff(name1, name2, ra1, dec1, ra2, dec2, obstime_mjd, verbose: bool=False):
    """Report difference in RA/DEC calculation according to two methods"""
    # difference in RA/DEC between skyfield and JPL
    ra_diff = np.abs(ra2 - ra1)
    dec_diff = np.abs(dec2 - dec1)
    # convert to arcseconds
    ra_diff_sec = ra_diff.to(arcsec).value
    dec_diff_sec = dec_diff.to(arcsec).value

    # Convert both to unit directions
    u1 = radec2dir(ra=ra1, dec=dec1, obstime_mjd=obstime_mjd)
    u2 = radec2dir(ra=ra2, dec=dec2, obstime_mjd=obstime_mjd)

    # Difference in unit directions
    u_diff = u2 - u1
    u_diff_norm = np.linalg.norm(u_diff, axis=0) 
    u_diff_deg =  np.rad2deg(u_diff_norm)
    u_diff_sec =  u_diff_deg * 3600

    # The name for the differences
    name_diff = 'Diff'

    # Is input data arrays?
    is_array: bool = isinstance(ra1, np.ndarray)
    is_scalar: bool = not is_array
        
    if is_array:
        ra_diff_mean = np.mean(ra_diff.value)
        dec_diff_mean = np.mean(dec_diff.value)
        u_diff_mean = np.mean(u_diff_deg)
        u_diff_median = np.median(u_diff_deg)
        u_diff_max = np.max(u_diff_deg)
        ra_diff_mean_sec = ra_diff_mean * 3600
        dec_diff_mean_sec = dec_diff_mean * 3600
        u_diff_mean_sec = u_diff_mean * 3600
        u_diff_median_sec = u_diff_median * 3600
        u_diff_max_sec = u_diff_max * 3600

    if verbose:
        # report RA
        if is_scalar:
            print(f'Difference in RA: {name2} - {name1}')
            print(f'{name2:12}: {ra2:8.6f} deg')
            print(f'{name1:12}: {ra1:8.6f} deg')
            print(f'{name_diff:12}: {ra_diff:8.6f} deg ({ra_diff_sec:5.2f} seconds)')
        else:
            # print(f'{name_diff:12}: {ra_diff_mean:8.6f} deg ({ra_diff_mean_sec:5.2f} seconds)')
            pass

        # report DEC
        if is_scalar:
            print(f'\nDifference in DEC: {name2} - {name1}')
            print(f'{name2:12}: {dec2:8.6f} deg')
            print(f'{name1:12}: {dec1:8.6f} deg')
            print(f'{name_diff:12}: {dec_diff:8.6f} deg ({dec_diff_sec:5.2f} seconds)')
        else:
            # print(f'{name_diff:12}: {dec_diff_mean:8.6f} deg ({dec_diff_mean_sec:5.2f} seconds)')
            pass

        # report direction vectors
        if is_scalar:
            print(f'\nUnit direction vectors:')
            print(f'{name2:12}:', u2)
            print(f'{name1:12}:', u1)
            print(f'{name_diff:12}:', u_diff)
            print(f'AngleDiff   : {u_diff_deg:8.6f} deg ({u_diff_sec:5.3f} seconds)')
        else:
            print(f'Angle Difference: {name2} vs. {name1}')
            print(f'Mean  : {u_diff_mean:10.6f} deg ({u_diff_mean_sec:8.3f} seconds)')
            print(f'Median: {u_diff_median:10.6f} deg ({u_diff_median_sec:8.3f} seconds)')
            print(f'Max   : {u_diff_max:10.6f} deg ({u_diff_max_sec:8.3f} seconds)')

    # Return the difference of direction vectors in seconds of arc
    return u_diff_sec

# ********************************************************************************************************************* 
def qv2obs(q: np.array, v: np.array, mjd: np.array, frame=BarycentricMeanEcliptic) -> SkyCoord:
    """
    Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space.
    INPUTS:
        q: Position of this body in given coordinate frame; passed with units (default AU)
        v: Velocity of this body in given coordinate frame; passed with units (default AU / day)
        mjd: Observation time on Earth as a modified Julian day
        frame: Astropy coordinate frame; defaults to BarycentricMeanEcliptic
    RETURNS:
        obs_gcrs: SkyCoordinate instance for this observation in the Geocentric frame (GCRS)
    """
    # The observation time
    obstime = astropy.time.Time(mjd, format='mjd')

    # Unpack position and velocity
    x, y, z = q
    v_x, v_y, v_z = v
    
    # The observation in the given coordinate frame (usually BarycentricMeanEcliptic)
    obs_frame = SkyCoord(x=x, y=y, z=z, v_x=v_x, v_y=v_y, v_z=v_z, obstime=obstime, 
                         representation_type='cartesian', differential_type='cartesian',
                         frame=frame)
    # Transform the observation to the Geocentric frame
    obs_gcrs = obs_frame.transform_to(GCRS)
    return obs_gcrs

# # ********************************************************************************************************************* 
# def qv2radec(q: np.array, v: np.array, mjd: np.array, frame=BarycentricMeanEcliptic) -> Tuple[np.ndarray, np.ndarray]:
#     """
#     Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space.
#     INPUTS:
#         q: Position of this body in given coordinate frame; passed with units (default AU)
#         v: Velocity of this body in given coordinate frame; passed with units (default AU / day)
#         mjd: Observation time on Earth as a modified Julian day
#         frame: Astropy coordinate frame; defaults to BarycentricMeanEcliptic
#     RETURNS:
#         ra: Right Ascension in degrees
#         dec: Declination in degrees
#         r: Distance in AU
#     """
#     # Delegate to qv2obs to get the observation in the GCRS frame
#     obs_gcrs = qv2obs(q=q, v=v, mjd=mjd, frame=frame)
#     # Extract the RA and DEC in degrees
#     ra = obs_gcrs.ra.deg
#     dec = obs_gcrs.dec.deg
#     r = obs_gcrs.distance.au
#     # Return (ra, dec, r) as a tuple
#     return (ra, dec, r)

# # ********************************************************************************************************************* 
# def qvrel2obs(q_body: np.array, v_body: np.array, 
#               q_earth: np.array, v_earth: np.array, 
#               mjd: np.array, frame=BarycentricMeanEcliptic) -> SkyCoord:
#     """
#     Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space
#     vs. integrated position and velocity of earth.
#     INPUTS:
#         q_body: Position of this body in given coordinate frame; passed with units (default AU)
#         v_body: Velocity of this body in given coordinate frame; passed with units (default AU / day)
#         q_earth: Position of earth in given coordinate frame; passed with units (default AU)
#         v_earth: Velocity of earth in given coordinate frame; passed with units (default AU / day)
#         mjd: Observation time on Earth as a modified Julian day
#         frame: Astropy coordinate frame; defaults to BarycentricMeanEcliptic
#     RETURNS:
#         obs: SkyCoordinate instance for this observation in the Geocentric frame (GCRS)
#     """
#     # The observation time
#     obstime = astropy.time.Time(mjd, format='mjd')
    
#     # Correct shape of zero_au and zero_km_sec because constants appear not to broadcast in SkyCoord constructor
#     if isinstance(mjd, np.ndarray):
#         shape = obstime.shape
#         zero_au = np.zeros(shape) * au
#         zero_km_sec = np.zeros(shape) * km / second

#     # Represent the earth in the GCRS frame.  Easy because it has position & velocity zero!
#     earth_gcrs = SkyCoord(x=zero_au, y=zero_au, z=zero_au, v_x=zero_km_sec, v_y=zero_km_sec, v_z=zero_km_sec,
#                           representation_type='cartesian', differential_type='cartesian',
#                           obstime=obstime, frame=GCRS)
#     # Represent earth in BarycentricMeanEcliptic frame, using Cartesian representation
#     earth_bme = earth_gcrs.transform_to(BarycentricMeanEcliptic)
#     earth_bme.representation_type = 'cartesian'
#     earth_bme.differential_type = 'cartesian'
    
#     # Position and speed of earth according to input in the given frame
#     x, y, z = q_earth
#     v_x, v_y, v_z = v_earth
#     earth_input_frame = SkyCoord(x=x, y=y, z=z, v_x=v_x, v_y=v_y, v_z=v_z,
#                                  representation_type='cartesian', differential_type='cartesian',
#                                  obstime=obstime, frame=frame)
#     # Compute input position of earth in BME; do transformation if input wasn't in BME
#     if frame != BarycentricMeanEcliptic:
#         # Convert earth to spherical representation so it can transform into BarycentricMeanEcliptic
#         earth_input_frame.representation_type='spherical'
#         earth_input_frame.differential_type='spherical'
#         # Position and speed of earth in the BarycentricMeanEcliptic frame
#         earth_input_bme = earth_input_frame.transform_to(BarycentricMeanEcliptic)
#     else:
#         earth_input_bme = earth_input_frame 
    
#     # Correction factor to add to position and velocity of earth in input frame so it matches BME
#     dq = earth_bme.cartesian.xyz - earth_input_bme.cartesian.xyz
#     dv = earth_bme.velocity.d_xyz - earth_input_bme.velocity.d_xyz
    
#     # Unpack position and velocity of space body in the given frame
#     x, y, z = q_body
#     v_x, v_y, v_z = v_body
#     body_frame = SkyCoord(x=x, y=y, z=z, v_x=v_x, v_y=v_y, v_z=v_z,
#                           representation_type='cartesian', differential_type = 'cartesian',
#                           obstime=obstime, frame=frame)
#     # Compute input position of body in BME; do transformation if input wasn't in BME
#     if frame != BarycentricMeanEcliptic:
#         # Convert body to spherical representation so it can transform into BarycentricMeanEcliptic
#         body_frame.representation_type='spherical'
#         body_frame.differential_type='spherical'
#         # Position and speed of space body in the BarycentricMeanEcliptic frame
#         body_bme = body_frame.transform_to(BarycentricMeanEcliptic)
#         # Represent the space body in Cartesian coordinates so we can extract q and v
#         body_bme.representation_type = 'cartesian'
#         body_bme.differential_type = 'cartesian'
#     else:
#         body_bme = body_frame
    
#     # Generate the corrected position and velocity of space body in BME
#     x, y, z = body_bme.cartesian.xyz + dq
#     v_x, v_y, v_z = body_bme.velocity.d_xyz + dv
#     # Observe space body in the BarcycentricMeanEcliptic frame, with the corrections
#     obs_bme = SkyCoord(x=x, y=y, z=z, v_x=v_x, v_y=v_y, v_z=v_z, obstime=obstime,
#                        representation_type='cartesian', differential_type='cartesian',
#                        frame=BarycentricMeanEcliptic)
#     # Transform observation to GCRS
#     obs_gcrs = obs_bme.transform_to(GCRS)
#     return obs_gcrs

# # ********************************************************************************************************************* 
# def qvrel2radec(q_body: np.array, v_body: np.array, 
#                 q_earth: np.array, v_earth: np.array, 
#                 mjd: np.array, frame=BarycentricMeanEcliptic) -> Tuple[np.array, np.array]:
#     """
#     Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space
#     vs. integrated position and velocity of earth.
#     INPUTS:
#         q_body: Position of this body in given coordinate frame; passed with units (default AU)
#         v_body: Velocity of this body in given coordinate frame; passed with units (default AU / day)
#         q_earth: Position of earth in given coordinate frame; passed with units (default AU)
#         v_earth: Velocity of earth in given coordinate frame; passed with units (default AU / day)
#         mjd: Observation time on Earth as a modified Julian day
#         frame: Astropy coordinate frame; defaults to BarycentricMeanEcliptic
#     RETURNS:
#         ra: Right Ascension in degrees
#         dec: Declination in degrees
#         r: Distance in AU
#     """
#     # Delegate to qvrel2obs to get the observation in the GCRS frame
#     obs_gcrs = qvrel2obs(q_body=q_body, v_body=v_body, q_earth=q_earth, v_earth=v_earth, mjd=mjd, frame=frame)
#     # Extract the RA and DEC in degrees
#     ra = obs_gcrs.ra.deg
#     dec = obs_gcrs.dec.deg
#     r = obs_gcrs.distance.au
#     # Return [ra, dec, r] as a tuple
#     return (ra, dec, r)

