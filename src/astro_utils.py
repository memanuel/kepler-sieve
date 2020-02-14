"""
Harvard IACS Masters Thesis
Astronomy Utilities

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

import numpy as np
import astropy
from astropy.units import deg, au, km, meter, day, minute, second
from astropy.coordinates import SkyCoord, ICRS, GCRS, BarycentricMeanEcliptic
from datetime import date, datetime, timedelta
from typing import Tuple

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

# Constants for distance of zero AU and velocity of zero AU / day
zero_au = 0.0 * au
zero_au_day = 0.0 * au / day
zero_km_sec = 0.0 * km / second
# Speed of light
light_speed = astropy.constants.c

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
def radec2dir(ra: float, dec: float, mjd: float) -> np.array:
    """
    Convert a RA and DEC as of observation time to a unit displacement vector u = (ux, uy, uz) in ecliptic plane.
    INPUTS:
    ra: An astrometric Right Ascension in the ICRF
    dec: An Astromentric Declination in the ICRF
    mjd: The observation time as a modified julian day
    RETURNS:
    u: An array [ux, uy, uz] on the unit sphere in the the Ecliptic frame
    EXAMPLE:
    u = radec2dir(ra=76.107414227, dec=23.884882701, mjd=58600.0)
    (this is Mars viewed from Earth at mjd 58600 / 2019-04-27 with JPL RA and DEC)
    """
    # Build the observation as a SkyCoord in the ICRS (oriented with earth, origin at barycenter)
    obstime = astropy.time.Time(mjd, format='mjd')
    obs_icrs = SkyCoord(ra=ra*deg, dec=dec*deg, obstime=obstime, frame=ICRS)
    # Convert to the barycentric ecliptic frame (oriented with ecliptic, origin at barycenter)
    obs_ecl = obs_icrs.transform_to(BarycentricMeanEcliptic)
    u = obs_ecl.cartesian.xyz
    return u.value

# ********************************************************************************************************************* 
def qv2obs(q: np.array, v: np.array, mjd: np.array, 
           frame=BarycentricMeanEcliptic, light_lag: bool = False) -> SkyCoord:
    """
    Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space.
    INPUTS:
        q: Position of this body in given coordinate frame; passed with units (default AU)
        v: Velocity of this body in given coordinate frame; passed with units (default AU / day)
        mjd: Observation time on Earth as a modified Julian day
        frame: Astropy coordinate frame; defaults to BarycentricMeanEcliptic
        light_lag: Flag, indicating whether to account for time for light to arrive; defaults to False
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
    
    # If method was called with light_lag flag, adjust for the speed of light
    if light_lag:
        # The light time; units will be in seconds
        light_time = obs_gcrs.distance.to(meter) / light_speed
        # The position adjustment; product of velocity in AU / day and light time in days
        dq_lt = v.to(au/day) * light_time.to(day)
        # Delegate second call to qv2obs, with this adjustment
        return qv2obs(q=q-dq_lt, v=v, mjd=mjd, frame=frame, light_lag=False)
    else:
        return obs_gcrs

# ********************************************************************************************************************* 
def qv2radec(q: np.array, v: np.array, mjd: np.array, 
             frame=BarycentricMeanEcliptic, light_lag: bool = False) -> Tuple[np.array, np.array]:
    """
    Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space.
    INPUTS:
        q: Position of this body in given coordinate frame; passed with units (default AU)
        v: Velocity of this body in given coordinate frame; passed with units (default AU / day)
        mjd: Observation time on Earth as a modified Julian day
        frame: Astropy coordinate frame; defaults to BarycentricMeanEcliptic
        light_lag: Flag, indicating whether to account for time for light to arrive; defaults to False
    RETURNS:
        ra: Right Ascension in degrees
        dec: Declination in degrees
        r: Distance in AU
    """
    # Delegate to qv2obs to get the observation in the GCRS frame
    obs_gcrs = qv2obs(q=q, v=v, mjd=mjd, frame=frame, light_lag=light_lag)
    # Extract the RA and DEC in degrees
    ra = obs_gcrs.ra.deg
    dec = obs_gcrs.dec.deg
    r = obs_gcrs.distance.au
    # Return (ra, dec, r) as a tuple
    return (ra, dec, r)

# ********************************************************************************************************************* 
def qvrel2obs(q_body: np.array, v_body: np.array, 
              q_earth: np.array, v_earth: np.array, 
              mjd: np.array, 
              frame=BarycentricMeanEcliptic, light_lag: bool = False) -> SkyCoord:
    """
    Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space
    vs. integrated position and velocity of earth.
    INPUTS:
        q_body: Position of this body in given coordinate frame; passed with units (default AU)
        v_body: Velocity of this body in given coordinate frame; passed with units (default AU / day)
        q_earth: Position of earth in given coordinate frame; passed with units (default AU)
        v_earth: Velocity of earth in given coordinate frame; passed with units (default AU / day)
        mjd: Observation time on Earth as a modified Julian day
        frame: Astropy coordinate frame; defaults to BarycentricMeanEcliptic
        light_lag: Flag, indicating whether to account for time for light to arrive; defaults to False
    RETURNS:
        obs: SkyCoordinate instance for this observation in the Geocentric frame (GCRS)
    """
    # The observation time
    obstime = astropy.time.Time(mjd, format='mjd')

    # Represent the earth in the GCRS frame.  Easy because it has position & velocity zero!
    earth_gcrs = SkyCoord(x=zero_au, y=zero_au, z=zero_au, v_x=zero_km_sec, v_y=zero_km_sec, v_z=zero_km_sec,
                          representation_type='cartesian', obstime=obstime, frame=GCRS)
    # Represent earth in BarycentricMeanEcliptic frame, using Cartesian representation
    earth_bme = earth_gcrs.transform_to(BarycentricMeanEcliptic)
    earth_bme.representation_type = 'cartesian'
    earth_bme.differential_type = 'cartesian'
    
    # Position and speed of earth according to input in the given frame
    x, y, z = q_earth
    v_x, v_y, v_z = v_earth
    earth_input_frame = SkyCoord(x=x, y=y, z=z, v_x=v_x, v_y=v_y, v_z=v_z,
                                 representation_type='cartesian', differential_type='cartesian',
                                 obstime=obstime, frame=frame)
    # Compute input position of earth in BME; do transformation if input wasn't in BME
    if frame != BarycentricMeanEcliptic:
        # Convert earth to spherical representation so it can transform into BarycentricMeanEcliptic
        earth_input_frame.representation_type='spherical'
        earth_input_frame.differential_type='spherical'
        # Position and speed of earth in the BarycentricMeanEcliptic frame
        earth_input_bme = earth_input_frame.transform_to(BarycentricMeanEcliptic)
    else:
        earth_input_bme = earth_input_frame 
    
    # Correction factor to add to position and velocity of earth in input frame so it matches BME
    dq = earth_bme.cartesian.xyz - earth_input_bme.cartesian.xyz
    dv = earth_bme.velocity.d_xyz - earth_input_bme.velocity.d_xyz
    
    # Unpack position and velocity of space body in the given frame
    x, y, z = q_body
    v_x, v_y, v_z = v_body
    body_frame = SkyCoord(x=x, y=y, z=z, v_x=v_x, v_y=v_y, v_z=v_z,
                          representation_type='cartesian', differential_type = 'cartesian',
                          obstime=obstime, frame=frame)
    # Compute input position of body in BME; do transformation if input wasn't in BME
    if frame != BarycentricMeanEcliptic:
        # Convert body to spherical representation so it can transform into BarycentricMeanEcliptic
        body_frame.representation_type='spherical'
        body_frame.differential_type='spherical'
        # Position and speed of space body in the BarycentricMeanEcliptic frame
        body_bme = body_frame.transform_to(BarycentricMeanEcliptic)
        # Represent the space body in Cartesian coordinates so we can extract q and v
        body_bme.representation_type = 'cartesian'
        body_bme.differential_type = 'cartesian'
    else:
        body_bme = body_frame
    
    # Generate the corrected position and velocity of space body in BME
    x, y, z = body_bme.cartesian.xyz + dq
    v_x, v_y, v_z = body_bme.velocity.d_xyz + dv
    # Observe space body in the BarcycentricMeanEcliptic frame, with the corrections
    obs_bme = SkyCoord(x=x, y=y, z=z, v_x=v_x, v_y=v_y, v_z=v_z, obstime=obstime,
                       representation_type='cartesian', differential_type='cartesian',
                       frame=BarycentricMeanEcliptic)
    # Transform observation to GCRS
    obs_gcrs = obs_bme.transform_to(GCRS)

    # If method was called with light_lag flag, adjust for the speed of light
    if light_lag:
        # The light time; units will be in seconds
        light_time = obs_gcrs.distance.to(meter) / light_speed
        # The position adjustment; product of velocity in AU / day and light time in days
        dq_lt = v.to(au/day) * light_time.to(day)
        # Delegate second call to qv2obs, with this adjustment
        return qv2obs(q=q-dq_lt, v=v, mjd=mjd, frame=frame, light_lag=False)
    else:
        return obs_gcrs

# ********************************************************************************************************************* 
def qvrel2radec(q_body: np.array, v_body: np.array, 
                q_earth: np.array, v_earth: np.array, 
                mjd: np.array, frame=BarycentricMeanEcliptic) -> Tuple[np.array, np.array]:
    """
    Compute a RA and DEC from Earth Geocenter given position and velocity of a body in space
    vs. integrated position and velocity of earth.
    INPUTS:
        q_body: Position of this body in given coordinate frame; passed with units (default AU)
        v_body: Velocity of this body in given coordinate frame; passed with units (default AU / day)
        q_earth: Position of earth in given coordinate frame; passed with units (default AU)
        v_earth: Velocity of earth in given coordinate frame; passed with units (default AU / day)
        mjd: Observation time on Earth as a modified Julian day
        frame: Astropy coordinate frame; defaults to BarycentricMeanEcliptic
        light_lag: Flag, indicating whether to account for time for light to arrive; defaults to False
    RETURNS:
        ra: Right Ascension in degrees
        dec: Declination in degrees
        r: Distance in AU
    """
    # Delegate to qvrel2obs to get the observation in the GCRS frame
    obs_gcrs = qvrel2obs(q_body=q_body, v_body=v_body, q_earth=q_earth, v_earth=v_earth, mjd=mjd,
                         frame=frame, light_lag=light_lag)
    # Extract the RA and DEC in degrees
    ra = obs_gcrs.ra.deg
    dec = obs_gcrs.dec.deg
    r = obs_gcrs.distance.au
    # Return [ra, dec, r] as a tuple
    return (ra, dec, r)

# *************************************************************************************************
def report_radec_diff(name1, name2, ra1, dec1, ra2, dec2, obstime_mjd):
    """Report difference in RA/DEC calculation according to two methods"""
    # difference in RA/DEC between skyfield and JPL
    ra_diff = ra2 - ra1
    dec_diff = dec2 - dec1
    # convert to seconds
    ra_diff_sec = ra_diff * 3600
    dec_diff_sec = dec_diff * 3600

    # Convert both to unit directions
    u1 = radec2dir(ra1, dec1, obstime_mjd)
    u2 = radec2dir(ra2, dec2, obstime_mjd)

    # Difference in unit directions
    u_diff = u2 - u1
    u_diff_norm = np.linalg.norm(u_diff) 
    u_diff_sec =  np.rad2deg(u_diff_norm) * 3600
    
    # The name for the differences
    name_diff = 'Diff'

    # report RA
    print(f'Difference in RA: {name2} - {name1}')
    print(f'{name2:12}: {ra2:8.6f} deg')
    print(f'{name1:12}: {ra1:8.6f} deg')
    print(f'{name_diff:12}: {ra_diff:8.6f} deg ({ra_diff_sec:5.2f} seconds)')

    # report DEC
    print(f'\nDifference in DEC: {name2} - {name1}')
    print(f'{name2:12}: {dec2:8.6f} deg')
    print(f'{name1:12}: {dec1:8.6f} deg')
    print(f'{name_diff:12}: {dec_diff:8.6f} deg ({dec_diff_sec:5.2f} seconds)')

    # report direction vectors
    print(f'\nUnit direction vectors:')
    print(f'{name2:12}:', u2)
    print(f'{name1:12}:', u1)
    print(f'{name_diff:12}:', u_diff)
    print(f'AngleDiff   : {u_diff_sec:5.3f} seconds')
    
    return u_diff_sec

# *************************************************************************************************
def xyz_to_sph(x: np.array, y: np.array, z: np.array):
    """
    Convert a Cartesian coordinates x, y, z of a displacement vector to 
    spherical coordinates r, alt, az"""
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
    Convert a Cartesian coordinates q with shape (N,3)o f a displacement vector to 
    spherical coordinates r, alt, az"""
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
def main():
    test_julian_day()

if __name__ == '__main__':
    main()
