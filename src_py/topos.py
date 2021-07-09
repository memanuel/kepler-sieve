"""
Astronomy: Calculate topos (position on Earth's surface)

Michael S. Emanuel
Fri Aug 23 16:13:28 2019

Functions in this module:
site2geoloc(site_name, verbose)
calc_topos_impl(obstime_mjd, obsgeoloc)
calc_topos(obstime_mjd, site_name)
"""

# Core
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

# Astronomy
from astropy.units import au, km, meter, day, second, deg
from astropy.coordinates import EarthLocation
from skyfield.api import Loader as SkyfieldLoader
from skyfield.toposlib import Topos

# Utility
from datetime import date

# Local
from astro_utils import mjd_to_jd, date_to_mjd
from db_utils import df2db

# ********************************************************************************************************************* 
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
        topos_tbl = npz['topos_tbl'].item()
    # print(f'loaded topos table {fname_topos}')
except:
    topos_tbl = dict()
    np.savez(fname_topos, topos_tbl=topos_tbl)

# Table keying site_name to ObservatoryID
site_name_tbl = {
    'geocenter': 0,
    'palomar': 1,
}

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
def calc_topos_impl(t_obs: np.ndarray, obsgeoloc: EarthLocation) -> np.ndarray:
    """
    Compute the topos adjustment from the center of earth to an observer's geolocation.
    INPUTS:
        t_obs:      observation time as a modified julian date
        obsgeoloc:  geolocation of the observatory as an astropy EarthLocation object
    RETURNS:
        dq_topos: An array of displacements in the ecliptic frame
    """
    # compute the correction due to the observatory of obstime_mjd and geoloc are passed
    # dq_obs is the displacement from geocenter to the observatory
    # the observation times as skyfield time objects
    t_obs_jd = mjd_to_jd(t_obs)
    t_obs_sf = ts_sf.tt_jd(t_obs_jd)
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
    dq_topos = topos.at(t_obs_sf).ecliptic_position().au.T * au
    dv_topos = topos.at(t_obs_sf).ecliptic_velocity().km_per_s.T * km / second
    # Convert the topos to flat numpy arrays in units of AU and AU / day
    dq_topos = dq_topos.to(au).value
    dv_topos = dv_topos.to(au/day).value

    return dq_topos, dv_topos

# ********************************************************************************************************************* 
def calc_topos(t_obs: np.ndarray, site_name: str):
    """
    Calculate topos adjustment using cached spline
    INPUTS:
        t_obs:      Observation times as a numpy array of mjds
        site_name:  Name of an observatory or geolocation site recognized by astropy, e.g. 'geocenter' or 'Palomar'
    """
    # convert site_name to lower case
    site_name = site_name.lower()
    # only look up splined topos if t_obs is passed and site name is not geocenter
    if t_obs is not None and site_name != 'geocenter':
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
            dq_topos, dv_topos = calc_topos_impl(t_obs=mjd_long, obsgeoloc=obsgeoloc)
            # build a cubic spline for this site
            x = mjd_long
            y = np.concatenate([dq_topos, dv_topos], axis=1)
            topos_spline = CubicSpline(x=x, y=y)
            # save this spline to the topos table
            topos_tbl[site_name] = topos_spline
            # save revised topos_tbl to disk
            np.savez(fname_topos, topos_tbl=topos_tbl)
        
        # We are now guaranteed that topos_spline is the spline for this site
        # Evaluate it at the desired observation times
        spline_out = topos_spline(t_obs)
        # unpack into dq and dv; units are in au and au/day
        dq_topos = spline_out[:, 0:3]
        dv_topos = spline_out[:, 3:6]
    else:
        # default is to use geocenter if t_obs and site_name are not passed
        dq_topos = np.zeros(3)
        dv_topos = np.zeros(3)
    return dq_topos, dv_topos

# ********************************************************************************************************************* 
def make_topos_df(site_name: str):
    """
    Generate DataFrame of topos adjustments
    INPUTS:
        site_name:  Name of an observatory or geolocation site recognized by astropy, e.g. 'geocenter' or 'Palomar'
    """
    # Date range to process
    mjd0: int = 48000
    mjd1: int = 63000
    TimeID_0: int = mjd0 * 1440
    TimeID_1: int = mjd1 * 1440

    # Array of TimeID and ts (in form of MJDs)
    time_id = np.arange(TimeID_0, TimeID_1+1, dtype=np.int64)
    # Convert array of time_ids into obstime_mjd
    t_obs: np.ndarray = time_id / 1440.0

    # The ObservatoryID
    observatory_id: int = site_name_tbl[site_name]
    observatory_ids = np.repeat(observatory_id, time_id.shape)

    # Calculate the arrays using calc_topos
    dq_topos, dv_topos = calc_topos(t_obs=t_obs, site_name=site_name)

    # Wrap this up into DataFrame
    df_tbl = {
        'ObservatoryID': observatory_ids,
        'TimeID': time_id,
        'mjd': t_obs,
    }
    df = pd.DataFrame(df_tbl)

    # Add displacement
    cols_dq = ['qx', 'qy', 'qz',]
    cols_dv = ['vx', 'vy', 'vz',]
    df[cols_dq] = dq_topos
    df[cols_dv] = dv_topos

    return df

# ********************************************************************************************************************* 
def main():
    """
    Populate DB table StateVectors_Topos
    """

    # Sites to process
    site_names = ['geocenter', 'palomar']

    # DB Table
    schema = 'KS'
    table = 'StateVectors_Topos'
    columns = ['ObservatoryID', 'TimeID', 'mjd', 'qx', 'qy', 'qz', 'vx', 'vy', 'vz']

    # Process the requested site names
    for site_name in site_names:
        df = make_topos_df(site_name=site_name)
        df2db(df=df, schema=schema, table=table, columns=columns)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()