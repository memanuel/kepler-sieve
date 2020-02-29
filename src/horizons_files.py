"""
Harvard IACS Masters Thesis
JPL Horizons Testing:
Utilities for handling downloads from Horizons obtained at
https://ssd.jpl.nasa.gov/horizons.cgi#top
Used for testing especially RA/DEC transformations in ra_dec.py

Michael S. Emanuel
27-Feb-2020
"""

# Library imports
import os
import numpy as np
import pandas as pd
import astropy
from astropy.units import deg, au, km, meter, day, minute, second, arcsec
import skyfield
from scipy.interpolate import CubicSpline

# MSE imports
from utils import range_inc
from ra_dec import radec2dir, dir2radec, qv2dir, site2geoloc, direction_diff

# ********************************************************************************************************************
# Create Skyfield loader in preferred location
skyfield_load = skyfield.api.Loader('../data/skyfield')

# Load Skyfield timescale
ts_sf = skyfield_load.timescale()

# Load planetary positions using de435
planets_sf = skyfield_load('de435.bsp')
earth_sf = planets_sf['earth']

# Suppress fake pandas warnings
pd.options.mode.chained_assignment = None

# ********************************************************************************************************************
def load_pos_jpl(body_name: str, dir_name: str = '../data/jpl/testing/hourly'):
    """
    Construct a DataFrame with the position and velocity of a body according to JPL
    INPUTS:
        body_name: The name of the body, used in the file names, e.g. 'earth' or 'asteroid-001'
        dir_name:  The directory where the data files are stored, e.g. '../data/jpl/testing/hourly'
    """
    # Load the JPL position as CSV
    filename: str = os.path.join(dir_name, f'vectors-{body_name}.txt')
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
def load_ast_jpl(ast_num0: int, ast_num1: int, dir_name: str = '../data/jpl/testing/hourly'):
    """
    Construct a DataFrame with the position and velocity of a batch of asteroids according to JPL.
    INPUTS:
        ast_num0: the first asteroid number to process, e.g. 1
        ast_num1: the last asteroid number to process, e.g. 16
        dir_name: the directory where the data files are stored, e.g. '../data/jpl/testing/hourly'
    """

    # List of dataframes; one per asteroid
    df_ast_list = []

    # Load the JPL position of asteroids one at a time
    for j in range_inc(ast_num0, ast_num1):
        body_name = f'asteroid-{j:03d}'
        df_ast_j = load_pos_jpl(body_name=body_name, dir_name=dir_name)
        # Add column for the asteroid_num
        df_ast_j.insert(loc=0, column='asteroid_num', value=j)
        # Add this to list of frames
        df_ast_list.append(df_ast_j)

    # Concatenate dataframes
    df_ast = pd.concat(df_ast_list)
    return df_ast

# ********************************************************************************************************************
def add_cols_obs(df):
    """Add new columns to a DataFrame of observations with the MJD, time_key and direction"""
    # Add column for mjd
    df['mjd'] = (df['JulianDate'] - 2400000.5)
    # integer column for the date/time
    df['time_key'] = np.int32(np.round(df['mjd']*24))
    # compute direction in the BME (BarycentricMeanEcliptic) frame
    ra = df.RA_jpl.values * deg
    dec = df.DEC_jpl.values * deg
    obstime_mjd = df['mjd'].values
    u = radec2dir(ra=ra, dec=dec, obstime_mjd = obstime_mjd)
    # Add columns for the direction in BME
    # the direction u has shape 3xN
    df['ux_jpl'] = u[0]
    df['uy_jpl'] = u[1]
    df['uz_jpl'] = u[2]
    return df

# ********************************************************************************************************************
def load_obs_jpl(body_name: str, observer_name: str, dir_name: str = '../data/jpl/testing/hourly'):
    """
    Construct a DataFrame with the observation of a body according to JPL
    INPUTS:
        body_name:   Name of the body being observed, used in the file names, e.g. 'mars'
        observer_name: Name of the observer location, e.g. 'geocenter' or 'palomar'
        dir_name:      Directory where the data files are stored, e.g. '../data/jpl/testing/hourly'
    """
    # Load the JPL observations as CSV
    filename: str = os.path.join(dir_name, f'observer-{body_name}-{observer_name}.txt')
    df = pd.read_csv(filename, index_col=False)

    # Rename RA and DEC to RA_jpl and DEC_jpl
    df.rename(columns={'RA': 'RA_jpl', 'DEC': 'DEC_jpl'}, inplace=True)
    
    # Add columns for time and directions
    add_cols_obs(df)

    # Order columns
    columns = ['mjd', 'JulianDate', 'time_key',  
               'RA_jpl', 'DEC_jpl', 'ux_jpl', 'uy_jpl', 'uz_jpl',
               'RA_apparent', 'DEC_apparent',
               'delta', 'delta_dot', 'light_time',]
    df = df[columns]
    return df

# ********************************************************************************************************************
def load_obs_ast_jpl(ast_num0: int, ast_num1: int, observer_name: str, 
                     dir_name: str = '../data/jpl/testing/hourly'):
    """
    Construct a DataFrame with the observation of a batch of asteroids according to JPL
    INPUTS:
        ast_num0: the first asteroid number to process, e.g. 1
        ast_num1: the last asteroid number to process, e.g. 16
        observer_name: Name of the observer location, e.g. 'geocenter' or 'palomar'
        dir_name:      Directory where the data files are stored, e.g. '../data/jpl/testing/hourly'
    """
    
    # List of dataframes; one per asteroid
    df_obs_list = []

    # Load the JPL observations of asteroids one at a time
    for j in range_inc(ast_num0, ast_num1):
        body_name = f'asteroid-{j:03d}'
        df_obs_j = load_obs_jpl(body_name=body_name, observer_name=observer_name, dir_name=dir_name)
        # Add column for the asteroid_num
        df_obs_j.insert(loc=0, column='asteroid_num', value=j)       
        # Add this to list of frames
        df_obs_list.append(df_obs_j)

    # Concatenate dataframes
    df_obs = pd.concat(df_obs_list)
    return df_obs
    
# ********************************************************************************************************************
def obs_add_radec(df_obs: pd.DataFrame, ra: np.ndarray, dec: np.ndarray, source: str) -> None:
    """
    Add columns to an observation DataFrame with a RA/DEC quoted by an alternate source (e.g. Skyfield)
    INPUTS:
        df_obs: Dataframe to be modified
        ra:     Right Ascension, as a numpy array of astropy angles
        dec:    Declination, as a numpy array of astropy angles
        source: Name of this source used as column suffix, e.g. 'jpl', 'sf' or 'mse'
    """
    # Alias source to src for legibility
    src = source

    # Add columns for RA and DEC
    df_obs[f'RA_{src}'] = ra
    df_obs[f'DEC_{src}'] = dec
    
    # Create empty columns for direction
    cols_dir = [f'ux_{src}', f'uy_{src}', f'uz_{src}']
    for col in cols_dir:
        df_obs[col] = 0.0
        
    # Compute direction using radec2dir()
    obstime_mjd = df_obs.mjd.values
    u = radec2dir(ra=ra, dec=dec, obstime_mjd=obstime_mjd)

    # Add directions to observation DataFrame
    # Need transpose u.T b/c shape of u is 3xN
    df_obs[cols_dir] = u.T
    
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
def obs_ast_add_interp_qv(df_obs: pd.DataFrame, df_ast: pd.DataFrame, df_earth: pd.DataFrame, 
                          source_name: str) -> None:
    """
    Add interpolated position and velocity to a DataFrame of asteroid observations.
    INPUTS:
        df_obs:     DataFrame of asteroid observations
        df_ast:     DataFrame of asteroid position and velocity
        df_earth:   DataFrame of earth position
        source_name: Name of the source for positions, e.g. 'jpl'
    OUTPUTS:
        Modifies df_obs in place        
    """
    # Alias source name for legibility
    src = source_name

    # List of distinct asteroid numbers
    ast_nums = sorted(set(df_obs.asteroid_num))

    # Add new columns
    body_cols = [f'body_x_{src}', f'body_y_{src}', f'body_z_{src}', 
                 f'body_vx_{src}', f'body_vy_{src}', f'body_vz_{src}']
    earth_cols = [f'earth_x_{src}', f'earth_y_{src}', f'earth_z_{src}']
    cols = body_cols + earth_cols
    for col in cols:
        df_obs[col] = 0.0    
    
    # Augment the asteroids one at a time
    for j in ast_nums:
        mask_obs = df_obs.asteroid_num == j
        mask_pos = df_ast.asteroid_num == j
        df_obs_j = df_obs.loc[mask_obs]
        df_ast_j = df_ast.loc[mask_pos]
        obs_add_interp_qv(df_obs=df_obs_j, df_body=df_ast_j, df_earth=df_earth, source_name='jpl')    
        df_obs.loc[mask_obs, cols] = df_obs_j[cols]
 
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

    # Observation times and geolocation of this site
    obstime_mjd = df_obs.mjd.values
    obsgeoloc = site2geoloc(site_name=site_name, verbose=False)

    # Compute the direction of the body from earth using qv2dir
    u, delta = qv2dir(q_body=q_body, v_body=v_body, q_earth=q_earth, obstime_mjd=obstime_mjd, site_name=site_name)

    # Compute the RA and DEC from the direction
    obstime_mjd = df_obs.mjd.values
    ra, dec = dir2radec(u=u, obstime_mjd=obstime_mjd)    

    # Save the RA and DEC
    df_obs[f'RA_calc_{src}'] = ra
    df_obs[f'DEC_calc_{src}'] = dec

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
    return u_mean_diff_sec
