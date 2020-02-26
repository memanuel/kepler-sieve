"""
Harvard IACS Masters Thesis
Astronomy: Transformation of Right Ascension & Declination (RA/DEC) to Ecliptic Directions

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

import numpy as np
import pandas as pd
import astropy
from astropy.units import deg, au, km, meter, day, minute, second, arcsec
from astropy.coordinates import SkyCoord, ICRS, GCRS, BarycentricMeanEcliptic, EarthLocation
from scipy.interpolate import CubicSpline
from datetime import date, datetime, timedelta
import os
from typing import Tuple, Optional

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
def qv2dir(q_body: np.ndarray, v_body: np.ndarray, q_earth: np.ndarray, 
           obstime_mjd: Optional[np.ndarray] = None, 
           obsgeoloc: Optional[EarthLocation] = None) -> np.ndarray:
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
    # compute the correction due to the observatory of obstime_mjd and geoloc are passed
    # dq_obs is the displacement from geocenter to the observatory
    if (obstime_mjd is not None and obsgeoloc is not None):
        # the observation times as astropy time objects
        obstime = astropy.time.Time(obstime_mjd, format='mjd')
        # the displacement from the geocenter to the observatory in the ICRS frame
        x, y, z = obsgeoloc.geocentric
        obs_icrs = SkyCoord(x=x, y=y, z=z, obstime=obstime, frame=ICRS, representation_type='cartesian')
        # displacement from geocenter to observatory in the BME frame
        obs_bme = obs_icrs.transform_to(BarycentricMeanEcliptic)
        dq_topos = obs_bme.cartesian.xyz.to(au)   
    else:
        # default is to use geocenter if obstime and geoloc are not passed
        dq_topos = np.zeros((3)) * au

    # infer data shape
    data_axis, space_axis, shape = infer_shape(q_earth)

    # reshape dq_topos to match q_earth
    dq_topos = dq_topos.reshape(shape)

    # position of the observer in space
    q_obs = q_earth + dq_topos

    # displacement from observer on earth to body; in AU
    q_rel = q_body - q_obs

    # distance; in AU
    r = np.linalg.norm(q_rel, axis=space_axis, keepdims=True) * au
    
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
    return u_diff_sec
    
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
    u_diff_sec = direction_diff(name1=name1, name2=name2, u1=u1, u2=u2, verbose=verbose)
    # Report RA/DEC difference if requested
    if verbose_radec:
        print(f'\nRA/DEC:')
        print(f'RA Mean : {ra_diff:10.6f} deg ({ra_diff_sec:8.3f} seconds)')
        print(f'DEC Mean: {dec_diff:10.6f} deg ({dec_diff_sec:8.3f} seconds)')

# ********************************************************************************************************************
# TESTING
# ********************************************************************************************************************

# ********************************************************************************************************************
def load_planet_pos_jpl(planet_name: str, dir_name: str = '../data/jpl/testing/hourly'):
    """
    Construct a DataFrame with the position and velocity of a planet according to JPL
    INPUTS:
        planet_name: The name of the planet, used in the file names, e.g. 'earth'
        dir_name:    The directory where the data files are stored, e.g. '../data/jpl/testing/hourly'
    """
    # Load the JPL position as CSV
    filename: str = os.path.join(dir_name, f'vectors-{planet_name}.txt')
    df = pd.read_csv(filename, index_col=False)
    # Add column for mjd
    df['mjd'] = (df['JulianDate'] - 2400000.5)
    # integer column for the date/time; count number of hours
    df['time_key'] = np.int32(np.round(df['mjd']*24))    
    # Order columns
    columns = ['mjd', 'JulianDate', 'time_key', 'X', 'Y', 'Z', 'VX', 'VY', 'VZ', 'LT', 'RG', 'RR']
    df = df[columns]
    return df

# ********************************************************************************************************************
def add_cols_obs(df):
    """Add new columns to a DataFrame of observations with the MJD, time_key and direction"""
    # Add column for mjd
    df['mjd'] = (df['JulianDate'] - 2400000.5)
    # integer column for the date/time
    df['time_key'] = np.int32(np.round(df['mjd']*24))
    # Add columns for the direction in ecliptic
    u = radec2dir(ra=df.RA.values * deg, dec=df.DEC.values * deg, obstime_mjd = df['mjd'].values)
    df['ux_jpl'] = u[0]
    df['uy_jpl'] = u[1]
    df['uz_jpl'] = u[2]
    return df

# ********************************************************************************************************************
def load_planet_obs_jpl(planet_name: str, observer_name: str, dir_name: str = '../data/jpl/testing/hourly'):
    """
    Construct a DataFrame with the observation of a planet according to JPL
    INPUTS:
        planet_name:   Name of the planet being observed, used in the file names, e.g. 'mars'
        observer_name: Name of the observer location, e.g. 'geocenter' or 'palomar'
        dir_name:      Directory where the data files are stored, e.g. '../data/jpl/testing/hourly'
    """
    # Load the JPL observations as CSV
    filename: str = os.path.join(dir_name, f'observer-{planet_name}-{observer_name}.txt')
    df = pd.read_csv(filename, index_col=False)
    
    # Add columns for time and directions
    add_cols_obs(df)

    # Order columns
    columns = ['mjd', 'JulianDate', 'time_key', 'ux_jpl', 'uy_jpl', 'uz_jpl', 
               'RA', 'DEC', 'RA_apparent', 'DEC_apparent',
               'delta', 'delta_dot', 'light_time',]
    df = df[columns]
    return df

# ********************************************************************************************************************
def obs_add_interp_qv(df_obs: pd.DataFrame, 
                      df_body: pd.DataFrame,
                      df_earth: pd.DataFrame,
                      source_name: str) -> None:
    """
    Add interpolated position and velocity to a DataFrame of observations
    INPUTS:
        df_obs:   DataFrame of observations
        df_body:  DataFrame with position of the body in BME frame
        df_earth: DataFrame with position of earth in BME frame
        source_name: Name of the source of the position data, e.g. 'jpl' or 'mse'
    OUTPUTS:
        Modifies df_obs in place
    """
    # Alias source_name to src for legibility
    src: str = source_name

    # Add new columns for position and velocity
    body_cols = [f'body_x_{src}', f'body_y_{src}', f'body_z_{src}', 
                 f'body_vx_{src}', f'body_vy_{src}', f'body_vz_{src}']
    earth_cols = [f'earth_x_{src}', f'earth_y_{src}', f'earth_z_{src}']
    new_cols = body_cols + earth_cols
    for col in new_cols:
        df_obs[col] = 0.0
        
    # Interpolator for earth position
    interp_t = df_earth.mjd.values
    interp_q = df_earth[['X', 'Y', 'Z']].values
    interp_earth = CubicSpline(x=interp_t, y=interp_q)
    
    # Set interpolated position of earth
    earth_q = interp_earth(df_obs.mjd.values)
    df_obs[earth_cols] = earth_q

    # Add interpolated position and velocity to df_obs
    # Build an interpolator for the body position and velocity
    interp_t = df_body.mjd.values
    interp_qv = df_body[['X', 'Y', 'Z', 'VX', 'VY', 'VZ']].values
    interp_ast = CubicSpline(x=interp_t, y=interp_qv)

    # The times to be interpolated
    obs_t = df_obs.mjd.values

    # Evaluate the interpolated position / velocity at the observation times
    obs_qv = interp_ast(obs_t)

    # Assign interpolated qv to the df_obs dataframe on the mask
    df_obs[body_cols] = obs_qv

# ********************************************************************************************************************
def obs_add_calc_dir(df_obs: pd.DataFrame, site_name: str, source_name: str) -> None:
    """
    Add calculated direction and RA/DEC to an observation DataFrame
    INPUTS:
        df_obs:      DataFrame of observations
        site_name:   Name of the location of the observer, e.g. 'geocenter' or 'palomar'
        source_name: Name of the source of the position data, e.g. 'jpl' or 'mse'
    OUTPUTS:
        Modifies df_obs in place
    """
    # Alias source_name to src for legibility
    src: str = source_name

    # Columns for position and velocity from this source
    q_body_cols = [f'body_x_{src}', f'body_y_{src}', f'body_z_{src}']
    v_body_cols = [f'body_vx_{src}', f'body_vy_{src}', f'body_vz_{src}']
    q_earth_cols = [f'earth_x_{src}', f'earth_y_{src}', f'earth_z_{src}']

    # Extract position and velocity of space body; build as Nx3 array with astropy units
    q_body = df_obs[q_body_cols].values * au
    v_body = df_obs[v_body_cols].values * au / day

    # Extract position of earth; build as Nx3 array with astropy units
    q_earth = df_obs[q_earth_cols].values * au

    # Geolocation of this site
    obsgeoloc = site2geoloc(site_name=site_name, verbose=False)

    # Compute the direction of the body from earth using qv2dir
    u = qv2dir(q_body=q_body, v_body=v_body, q_earth=q_earth, obsgeoloc=obsgeoloc)    

    # Columns for computed direction from this source
    u_cols = [f'ux_calc_{src}', f'uy_calc_{src}', f'uz_calc_{src}']

    # Create new columns in DataFrame
    for col in u_cols:
        df_obs[col] = 0.0

    # Set computed direction
    df_obs[u_cols] = u

# ********************************************************************************************************************
def obs_direction_diff(df_obs: pd.DataFrame, src1: str, src2: str, verbose: str = False):
    """
    Compute and report the difference in directions between two sources on an observation DataFrame.
    INPUTS:
        df_obs: DataFrame of observations with directions from multiple sources
        src1:   First source in comparison, e.g. 'jpl' for direction calculated from quoted RA/DEC
        src2:   Second source in comparison, e.g. 'calc_jpl' for MSE qv2dir() calculation of JPL positions
        verbose: Wheth
    RETURNS:
        u_mean_diff_sec: Mean difference between sources measured in arc seconds
    """
    # Columns for direction according to these sources
    u1_cols = [f'ux_{src1}', f'uy_{src1}', f'uz_{src1}']
    u2_cols = [f'ux_{src2}', f'uy_{src2}', f'uz_{src2}']

    # Assemble directions into arrays u1 and u2; shape will be Nx3
    u1 = df_obs[u1_cols].values
    u2 = df_obs[u2_cols].values

    # Delegate to direction_diff
    u_mean_diff_sec = direction_diff(name1=src1, name2=src2, u1=u1, u2=u2, verbose=verbose)