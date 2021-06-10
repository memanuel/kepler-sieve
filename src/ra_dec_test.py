"""
Test direction from known asteroids to Earth center implied by integrated trajectories.
Example calls:
$ python asteroid_direction.py 0 1000

Functions in this module:
main()

Michael S. Emanuel
2021-06-04
"""

# Core
import numpy as np

# Astronomy
from astropy.units import au, deg

# Local imports
from astro_utils import infer_shape
from ra_dec import radec2dir

# *************************************************************************************************
# Functions for testing / comparing these calculations to JPL and SkyField
# *************************************************************************************************

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
        ra1:   Array of right ascension from source 1
        dec1:  Array of declination from source 1
        ra2:   Array of right ascension from source 2
        dec2:  Array of declination from source 2
        verbose: Whether to report results to console (angle between direction)
        verbose_radec: Whether to report differences as difference in RA and DEC separately
    """
    # Convert both RA/DEC inputs to unit directions
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

# ********************************************************************************************************************* 
def main():
    pass

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
