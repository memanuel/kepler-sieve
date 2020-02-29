"""
Harvard IACS Masters Thesis
Astronomy: Transformation of Right Ascension & Declination (RA/DEC) to Ecliptic Directions

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

# Core
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

# Astronomy
import astropy
from astropy.units import deg, au, km, meter, day, minute, second, arcsec
from astropy.coordinates import SkyCoord, ICRS, GCRS, BarycentricMeanEcliptic, EarthLocation
import skyfield
from skyfield.api import Loader as SkyfieldLoader
from skyfield.toposlib import Topos

# Utility
import os
from datetime import date, datetime

# MSE imports
from utils import range_inc
from astro_utils import jd_to_mjd, mjd_to_jd, date_to_mjd

# Typing
from typing import Tuple, Optional

# ********************************************************************************************************************* 
# Constants for distance of zero AU and velocity of zero AU / day
zero_au = 0.0 * au
zero_au_day = 0.0 * au / day
zero_km_sec = 0.0 * km / second
# Speed of light; express this in AU / minute
light_speed = astropy.constants.c.to(au / minute)

# Create Skyfield loader in preferred location
skyfield_load = SkyfieldLoader('../data/skyfield')

# Load Skyfield timescale
ts_sf = skyfield_load.timescale()

# Load planetary positions using de435
planets_sf = skyfield_load('de435.bsp')
earth_sf = planets_sf['earth']

# Table with saved splines for topos adjustment
fname_topos = '../data/skyfield/topos_tbl.npz'
try:
    with np.load(fname_topos, allow_pickle=True) as npz:
        topos_tbl = npz['topos_tbl'][0]
except:
    topos_tbl = dict()
    np.savez(fname_topos, topos_tbl=topos_tbl)

# Suppress fake pandas warnings
pd.options.mode.chained_assignment = None

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
    # Return as a numpy array of shape 3xN (easier b/c this is astropy default)
    return u.value

# # ********************************************************************************************************************* 
# def radec_app2dir(ra: float, dec: float, obstime_mjd: float) -> np.array:
#     """
#     Convert an apparent RA and DEC as of observation time to a unit displacement 
#     vector u = (ux, uy, uz) in ecliptic plane.
#     INPUTS:
#         ra: An apparent Right Ascension in the ICRF; passed with units (default degrees)
#         dec: An apparent Declination in the ICRF;  passed with units (default degrees)
#         obtime_mjd: The observation time as a modified julian day
#     RETURNS:
#         u: An array [ux, uy, uz] on the unit sphere in the the ecliptic frame
#     EXAMPLE:
#         u = radec_app2dir(ra=76.391533*deg, dec=23.90881*deg, obstime_mjd=58600.0)
#         (this is Mars viewed from Earth at mjd 58600 / 2019-04-27 with JPL RA and DEC)
#     """
#     # this is a placeholder function declaration.
#     # probably don't need this
#     pass

# ********************************************************************************************************************* 
def site2geoloc(site_name: str, verbose: bool = False):
    """
    Compute the geolocation of an object from its name
    INPUTS:
        site_name: Name of an observatory or geolocation site recognized by astropy, e.g. 'geocenter' or 'Palomar'
    OUTPUTS:
        geoloc: Geolocation of the observatory as an astropy EarthLocation object
    """
    # Construct the geoolocation
    zero_meter = 0.0 * meter
    if site_name == 'geocenter':
        geoloc = EarthLocation(x=zero_meter, y=zero_meter, z=zero_meter)
    else:
        geoloc = EarthLocation.of_site(site_name)
    
    # Report results if required
    if verbose:
        print(f'Geolocation of {site_name}:')
        print(f'cartesian = {geoloc}')
        print(f'geodetic  = {geoloc.geodetic}') 
    return geoloc

# ********************************************************************************************************************* 
def infer_shape(q):
    """
    Infer the axes with data and space dimensions
    INPUTS:
        q: a vector
    OUTPUTS:
        data_axis: the index of the axis for the N data entries, e.g. 0 when data is Nx3
        space_axis: the index of the axis for the 3 space dimensions, e.g. 1 when data is Nx3
        shape: the shape of the vectors, e.g. (-1, 3) when data is Nx3
    """
    # 
    if q.shape[0] == 3:
        data_axis, space_axis = 1, 0
        shape = (3, -1)
    elif q.shape[1] == 3:
        data_axis, space_axis = 0, 1
        shape = (-1, 3)
    else:
        raise ValueError(f'Bad data shape! q_earth has shape {q_earth.shape}, shoulde by Nx3 or 3xN.')

    return data_axis, space_axis, shape

# ********************************************************************************************************************* 
def calc_topos_impl(obstime_mjd: np.ndarray, obsgeoloc: EarthLocation) -> np.ndarray:
    """
    Compute the topos adjustment from the center of earth to an observer's geolocation.
    INPUTS:
        obstime_mjd: observation time as a modified julian date
        obsgeoloc: geolocation of the observatory as an astropy EarthLocation object
    RETURNS:
        dq_topos: An array of displacements in the ecliptic frame
    """
    # compute the correction due to the observatory of obstime_mjd and geoloc are passed
    # dq_obs is the displacement from geocenter to the observatory
    # the observation times as skyfield time objects
    obstime_jd = mjd_to_jd(obstime_mjd)
    obstime_sf = ts_sf.tt_jd(obstime_jd)
    # unpack the geodetic coordinates of this observer geolocation
    longitude, latitude, height = obsgeoloc.geodetic
    longitude_degrees = longitude.to(deg).value
    latitude_degrees = latitude.to(deg).value
    elevation_m = height.to(meter).value
    # construct a Skyfield topos instance
    topos = Topos(latitude_degrees=latitude_degrees, 
                  longitude_degrees=longitude_degrees, 
                  elevation_m=elevation_m)
    # query the topos object at the observation times
    dq_topos = topos.at(obstime_sf).ecliptic_position().au.T * au
    dv_topos = topos.at(obstime_sf).ecliptic_velocity().km_per_s.T * km / second
    dv_topos = dv_topos.to(au/day)

    return dq_topos, dv_topos

# TODO: build a spline of topos adjustments for palomar, for speed; rename this calc_topos_impl
# and have calc_topos just do a lookup from the spline

# ********************************************************************************************************************* 
def calc_topos(obstime_mjd, site_name: str):
    """
    Calculate topos adjustment using cached spline
    """
    # convert site_name to lower case
    site_name = site_name.lower()
    # only look up splined topos if obstime_mjd is passed and site name is not geocenter
    if obstime_mjd is not None and site_name != 'geocenter':
        try:
            # Try to find the spline for this site name
            topos_spline = topos_tbl[site_name]
        except KeyError:
            print(f'Generating Topos spline for {site_name}...')
            # sample on a very big range of dates to ensure spline lookup is good
            mjd0 = date_to_mjd(date(1980,1,1))
            mjd1 = date_to_mjd(date(2060,1,1))
            mjd_long = np.arange(start=mjd0, stop=mjd1, step=0.0625)
            # convert the site name to a geolocation
            obsgeoloc = site2geoloc(site_name)
            # do the actual topos calculation with calc_topos_impl; pass long array mjd_long
            dq_topos, dv_topos = calc_topos_impl(obstime_mjd=mjd_long, obsgeoloc=obsgeoloc)
            # build a cubic spline for this site
            x = mjd_long
            y = np.concatenate([dq_topos.value, dv_topos.value], axis=1)
            topos_spline = CubicSpline(x=x, y=y)
            # save this spline to the topos table
            topos_tbl[site_name] = topos_spline
            # save revised topos_tbl to disk
            np.savez(fname_topos, topos_tbl=topos_tbl)
        
        # We are now guaranteed that topos_spline is the spline for this site
        # Evaluate it at the desired observation times
        spline_out = topos_spline(obstime_mjd)
        # unpack into dq and dv, with units
        dq_topos = spline_out[0:3] * au
        dv_topos = spline_out[3:6] * au/day
    else:
        # default is to use geocenter if obstime and geoloc are not passed
        dq_topos = np.zeros(3) * au
        dv_topos = np.zeros(3) * au/day
    return dq_topos

# ********************************************************************************************************************* 
def astrometric_dir(q_body: np.ndarray, v_body: np.ndarray, q_obs: np.ndarray):
    """
    Compute the astrometric direction from earth to a space body as a unit displacement vector
    u = (ux, uy, uz) in the ecliptic plane.
    INPUTS:
        q_body: position of this body in ecliptic coordinate frame; passed with units (default AU)
        v_body: velocity of this body in ecliptic coordinate frame; passed with units (default AU / day)
        q_obs:  position of observer in ecliptic coordinate frame; passed with units (default AU)
                typically this will be computed as q_earth + dq_topos
    RETURNS:
        u: An array [ux, uy, uz] on the unit sphere in the ecliptic frame
    """
    # displacement from observer on earth to body; in AU
    q_rel = q_body - q_obs

    # distance; in AU
    r = np.linalg.norm(q_rel, axis=space_axis, keepdims=True)
    
    # light time in minutes
    light_time = (r / light_speed)
    
    # adjusted relative position, accounting for light time
    dq_lt = light_time * v_body.to(au/minute)     # convert velocity to au / minute b/c time is in minutes
    q_rel_lt = q_rel - dq_lt
    
    # adjusted direction
    r_lt = np.linalg.norm(q_rel_lt, axis=space_axis, keepdims=True) * au
    u = q_rel_lt / r_lt
    return u.value

# ********************************************************************************************************************* 
def qv2dir(q_body: np.ndarray, v_body: np.ndarray, q_earth: np.ndarray, 
           obstime_mjd: Optional[np.ndarray] = None, 
           site_name: str = 'geocenter') -> np.ndarray:
    """
    Compute the direction of displacement from earth to a space body as a unit displacement vector
    u = (ux, uy, uz) in the ecliptic plane.
    INPUTS:
        q_body: position of this body in ecliptic coordinate frame; passed with units (default AU)
        v_body: velocity of this body in ecliptic coordinate frame; passed with units (default AU / day)
        q_earth: position of earth in ecliptic coordinate frame; passed with units (default AU)
        obstime_mjd: observation time as a modified julian date; only required if passing obsgeoloc 
        obsgeoloc: geolocation of the observatory as an astropy EarthLocation object
    RETURNS:
        u: An array [ux, uy, uz] on the unit sphere in the ecliptic frame
    EXAMPLE:
        u = qv2dir(q_body=np.array([-0.328365, 1.570624, 0.040733])*au, 
                   v_body=np.array([-0.013177, -0.001673, 0.000288])*au/day,
                   q_earth=np.array([-0.813785, -0.586761, -0.000003])*au,
                   obsgeoloc=[-2410346.78217658, -4758666.82504051, 3487942.97502457] * meter)
    """
    # compute the correction due to the observatory of obstime_mjd and site_name are passed
    dq_topos = calc_topos(objstime_mjd=obstime_mjd, site_name=site_name)

    # reshape dq_topos to match q_earth
    data_axis, space_axis, shape = infer_shape(q_earth)
    dq_topos = dq_topos.reshape(shape)

    # position of the observer in space
    q_obs = q_earth + dq_topos

    # calculate astrometric direction
    u = astrometric_dir(q_body=q_body, v_body=v_body, q_obs=q_obs)

    return u

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
    ra = obs_icrs.ra.deg * deg
    dec = obs_icrs.dec.deg * deg
    # Return (ra, dec) as a tuple
    return (ra, dec)

# *************************************************************************************************
def direction_diff(name1: str, name2: str, u1: np.ndarray, u2: np.ndarray, verbose: bool=False) -> float:
    """
    Report the difference in directions in arc seconds
    INPUTS:
        name1: Descriptive name of the first source, e.g. 'JPL'
        name2: Descriptive name of the second source, e.g. 'MSE'
        u1:    Array of directions from source 1; shape Nx3 or 3xN
        u2:    Array of directions from source 2; shape Nx3 or 3xN
        verbose: Whether to report results to console
    """
    # infer data shape
    data_axis, space_axis, shape = infer_shape(u1)

    # Difference in unit directions
    u_diff = u2 - u1
    u_diff_norm = np.linalg.norm(u_diff, axis=space_axis) 
    u_diff_deg =  np.rad2deg(u_diff_norm)
    u_diff_sec =  u_diff_deg * 3600

    # The name for the differences
    name_diff = 'Diff'

    # Calculate mean, median and max difference in degrees
    u_diff_mean = np.mean(u_diff_deg)
    u_diff_median = np.median(u_diff_deg)
    u_diff_max = np.max(u_diff_deg)

    # Convert differences to arc seconds
    u_diff_mean_sec = u_diff_mean * 3600
    u_diff_median_sec = u_diff_median * 3600
    u_diff_max_sec = u_diff_max * 3600

    if verbose:
        print(f'Angle Difference: {name2} vs. {name1}')
        print(f'Mean  : {u_diff_mean:10.6f} deg ({u_diff_mean_sec:8.3f} seconds)')
        print(f'Median: {u_diff_median:10.6f} deg ({u_diff_median_sec:8.3f} seconds)')
        print(f'Max   : {u_diff_max:10.6f} deg ({u_diff_max_sec:8.3f} seconds)')

    # Return the difference of direction vectors in seconds of arc
    return u_diff_mean_sec

# *************************************************************************************************
def radec_diff(name1, name2, ra1, dec1, ra2, dec2, obstime_mjd, 
               verbose: bool=False, verbose_radec:bool =False):
    """
    Report difference in RA/DEC calculation according to two sources
    INPUTS:
        name1: Descriptive name of the first source, e.g. 'JPL'
        name2: Descriptive name of the second source, e.g. 'MSE'
        u1:    Array of directions from source 1; shape Nx3 or 3xN
        u2:    Array of directions from source 2; shape Nx3 or 3xN
        verbose: Whether to report results to console (angle between direction)
        verbose_radec: Whether to report differences as difference in RA and DEC separately
    """
    # Convert both to unit directions
    u1 = radec2dir(ra=ra1, dec=dec1, obstime_mjd=obstime_mjd)
    u2 = radec2dir(ra=ra2, dec=dec2, obstime_mjd=obstime_mjd)

    # mean difference in RA/DEC in degrees
    ra_diff = np.mean(np.abs(ra2 - ra1).value)
    dec_diff = np.mean(np.abs(dec2 - dec1).value)
    # convert to arcseconds
    ra_diff_sec = ra_diff * 3600
    dec_diff_sec = dec_diff * 3600

    # Delegate to direction_diff
    u_diff_mean_sec = direction_diff(name1=name1, name2=name2, u1=u1, u2=u2, verbose=verbose)
    # Report RA/DEC difference if requested
    if verbose_radec:
        print(f'\nRA/DEC:')
        print(f'RA Mean : {ra_diff:10.6f} deg ({ra_diff_sec:8.3f} seconds)')
        print(f'DEC Mean: {dec_diff:10.6f} deg ({dec_diff_sec:8.3f} seconds)')

    return u_diff_mean_sec

# ********************************************************************************************************************
def skyfield_observe(observer_sf, body_sf, obstime_sf):
    """
    Observe a body in Skyfield
    INPUTS:
        observer_sf: observer as a Skyfield object, e.g. observer_sf=planets['earth']
        body_sf:     body to be observed as a Skyfield object, e.g. body_sf=planets['mars']
        obstime_sf:  observation time as a Skyfield time array, e.g. ts_sf.tt_jd(obstime_mjd)
    OUTPUTS:
        (ra, dec, delta): tuple of numpy arrays of astropy angles and distances (degrees and AU)
    EXAMPLE:
        # Load ephemeris and timescale
        planets_sf = skyfield_load('../data/jpl/ephemeris/de435.bsp')
        earth_sf = planets_sf['earth']
        mars_sf = planets_sf['mars barycenter']
        ts_sf = skyfield_load.timescale()
        obstime_mars_jd = df_obs_mars_geo.mjd.values
        obstime_mars_sf = ts_sf.tt_jd(obstime_mars_jd)

        # Observe Mars from Earth geocenter
        ra_mars_geo_sf, dec_mars_geo_sf, delta_mars_geo_sf = \
            skyfield_observe(observer_sf=earth_sf, body_sf=mars_sf, obstime_sf=obstime_mars_sf)
        
        # Calculate topos of Palomar observatory
        geoloc_pal = site2geoloc('palomar', verbose=True)
        lon, lat, height = geoloc_pal.geodetic
        palomar_topos = skyfield.toposlib.Topos(latitude_degrees=lat.value, longitude_degrees=lon.value, elevation_m=height.value)

        # Observe Mars from Palomar observatory
        palomar_sf = earth_sf + palomar_topos
        ra_mars_pal_sf, dec_mars_pal_sf, delta_mars_pal_sf = \
            skyfield_observe(observer_sf=palomar_sf, body_sf=mars_sf, obstime_sf=obstime_mars_sf)
    """
    
    # Observe body from observer with Skyfield
    obs_sf = observer_sf.at(obstime_sf).observe(body_sf)

    # Build Skyfield angle arrays (RA, DEC) and distance array (delta)
    ra_aa, dec_aa, delta_da = obs_sf.radec()

    # Extract degrees and AU to get arrays of astropy angles and distances
    ra = ra_aa._degrees * deg
    dec = dec_aa._degrees * deg
    delta = delta_da.au * au
    
    # Return ra, dec, delta as astropy arrays
    return (ra, dec, delta)
