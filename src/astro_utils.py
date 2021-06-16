"""
Astronomy Utilities

Michael S. Emanuel
Fri Aug 23 16:13:28 2019

Functions defined in this module:
date_to_jd(t)
date_to_mjd(t)
jd_to_date(jd)
mjd_to_date(mjd)
datetime_to_jd(t)
datetime_to_mjd(t, epoch)
datetime_to_year(t)
jd_to_datetime(jd)
mjd_to_datetime(mjd)
jd_to_mjd(jd)
mjd_to_jd(mjd)
dist2rad(dist)
rad2dist(x_rad)
dist2deg(dist)
deg2dist(x_deg)
dist2sec(dist)
"""

# Core
import numpy as np
from datetime import date, datetime, timedelta
from collections import namedtuple

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
# Named tuple data type for orbital elements
OrbitalElement = namedtuple('OrbitalElement', 'a e inc Omega omega f')

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
# The MJD of January 1, 2000; used for converting datetime to year (float)
mjd_01Jan2000: float = date_to_mjd(date(2000,1,1))

# *************************************************************************************************
def datetime_to_year(t: datetime):
    """Convert a datetime object to a year (float)"""
    mjd = datetime_to_mjd(t)
    return 2000.0 + (mjd - mjd_01Jan2000) / 365.25

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

# ********************************************************************************************************************* 
def reverse_velocity(sim):
    """Reverse the velocities in a simulation for backwards time integration; modifies sim in place"""
    for p in sim.particles:
        vx, vy, vz = p.vx, p.vy, p.vz
        p.vx = -vx
        p.vy = -vy
        p.vz = -vz

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
    if q.shape[0] == 3:
        data_axis, space_axis = 1, 0
        shape = (3, -1)
    elif q.shape[1] == 3:
        data_axis, space_axis = 0, 1
        shape = (-1, 3)
    else:
        raise ValueError(f'Bad data shape! q_earth has shape {q.shape}, should by Nx3 or 3xN.')

    return data_axis, space_axis, shape

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
def main():
    test_julian_day()
    # test_anomaly_M2E()

# *************************************************************************************************
if __name__ == '__main__':
    main()
