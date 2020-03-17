"""
Harvard IACS Masters Thesis
ZTF Data
Utilities for acquiring data from the ZTF2 database using the Alerce data broker.

Michael S. Emanuel
27-Feb-2020
"""

# Libraries for getting Alerce data out of ZTF2 database
import json
import psycopg2
from alerce.api import AlerceAPI

# Standard libraries
import numpy as np
import pandas as pd

# Astronomy related
from astropy.units import deg
from scipy.interpolate import CubicSpline

# Utility
import os
import datetime
from datetime import date
from tqdm.auto import tqdm

# MSE imports
from utils import range_inc
from astro_utils import date_to_mjd, deg2dist
from ra_dec import radec2dir
from asteroid_dataframe import spline_ast_vec_dir

# Typing
from typing import Optional

# ********************************************************************************************************************* 
# Global variables

# Credentials for ZTF2 connection
credentials_file = "../alerce/alercereaduser.json"
with open(credentials_file) as fh:
    cred = json.load(fh)["params"]

# Connect to ZTF2 database; this is a shared resource, just one for the module
try:
    conn = psycopg2.connect(dbname=cred['dbname'], user=cred['user'], host=cred['host'], password=cred['password'])
except:
    print(f'Unable to get ZTF2 connection.')
    conn = None

# Directory for files with nearest asteroid to ztf observations
ztf_ast_dir_name = '../data/ztf_ast'

# ********************************************************************************************************************* 
def load_ztf_det(mjd0: float, mjd1: float):
    """
    Load all the ZTF detections classified as asteroids in a date range.
    INPUTS:
        mjd0: Modified Julian Date (mjd) to start processing, e.g. 55197.0 for 2010-01-01
        mjd1: Modified Julian Date (mjd) to end processing, e.g. 55197.0+7.0 for one week
    RETURNS:
        df: Pandas DataFrame of detections.  Includes columns ObjectID, CandidateID, mjd, ra, dec, asteroid_prob
    """
    
    # SQL with parameters
    sql = \
    """
    select
        obj.oid as ObjectID,
        det.candid as CandidateID,
        det.mjd,
        det.ra,
        det.dec,
        det.magpsf as mag_psf,
        det.magap as mag_app,
        det.magnr as nag_nr,
        det.sigmara as sigma_ra,
        det.sigmadec as sigma_dec,
        obj.pclassearly as asteroid_prob
    from 
        detections as det
        inner join objects as obj on obj.oid = det.oid
    where
        %(mjd0)s <= det.mjd and det.mjd < %(mjd1)s and
        obj.classearly = %(asteroid_class)s
    """
    
    # Set up parameters
    params = {'mjd0':mjd0, 'mjd1':mjd1, 'asteroid_class':21}
    
    # Run query and return DataFrame
    df = pd.read_sql_query(sql=sql, con=conn, params=params)

    # Create new ObjectID column, with dtype Pandas string ('|S') rather than object
    df.insert(loc=0, column='ObjectID', value=df.objectid.astype('|S'))    
    # Drop old objectid column, which is now redundant
    df.drop(columns='objectid', inplace=True)
    
    # Create new CandidateID column, with dtype Pandas np.int64 rather than object
    df.insert(loc=1, column='CandidateID', value=df.candidateid.astype(np.int64))    
    # Drop old candidateid column, which is now redundant
    df.drop(columns='candidateid', inplace=True)

    return df

# ********************************************************************************************************************* 
def load_ztf_det_year(year: int, save_dir: str = '../data/ztf'):
    """
    Load all ZTF detections in a year
    INPUTS:
        year: the year to be processed, e.g. 2019
        save_dir: directory to save the DataFrame
    OUTPUTS:
        df:   one combined DataFrame with all ZTF2 asteroid detections in this year
              same columns as load_ztf_det.
    """
    # Check if file already exists.  If so, load it from memory and return early
    file_name = os.path.join(save_dir, f'ztf-detections-{year}.h5')
    try:
        df = pd.read_hdf(file_name)
        print(f'Loaded {file_name} from disk.')
        return df
    except:
        print(f'Querying ZTF2 for data...')
    
    # Start and end date
    mjd0 = date_to_mjd(date(year,1,1))
    mjd1 = date_to_mjd(date(year+1,1,1))

    # Sample dates at weekly intervals, with one long "week" at the end
    mjd = np.arange(mjd0, mjd1, 7, dtype=np.float64)
    mjd[-1] = mjd1
    # Number of weeks (will be 52)
    n: int = len(mjd)-1

    # Array of DataFrames
    df_list = np.empty(shape=n, dtype=np.object)
    
    # Process one week at a time
    for i in tqdm(range(n)):
        # Start and end of week i
        mjd0_week, mjd1_week = mjd[i:i+2]
        # print(f'i={i}; mjd0={mjd0_week}, mjd1={mjd1_week}')
        # Process this week
        df_i = load_ztf_det(mjd0=mjd0_week, mjd1=mjd1_week)
        df_list[i] = df_i
        # Save this week into the HDF5; append mode adds one week at a time
        # df_i.to_hdf(file_name, key='df', mode='a', append=True)
        
    # Assemble the weeks into one DataFrame and return it
    df = pd.concat(df_list)
    # Save the combined DataFrame into an HDF5 file
    df.to_hdf(file_name, key='df', mode='w')
    return df

# ********************************************************************************************************************* 
def ztf_det_add_dir(df: pd.DataFrame, file_name: str, dir_name: str='../data/ztf'):
    """
    Add calculated directions to DataFrame of ZTF observations
    INPUTS:
        df: DataFrame of ZTF observations including ra, dec and mjd columns
        file_name: Name to save DataFrame on disk, e.g. 'ztf-detections.h5'
        dir_name:  Directory to save DataFrame on disk, e.g. '../data/ztf'
    OUTPUTS:
        Modifies df in place and saves it to disk.
    """
    # Extract mjd, ra, and dec as vectors of astropy angles
    mjd = df.mjd.values
    ra = df.ra.values * deg
    dec = df.dec.values * deg

    # Compute the directions using radec2dir()
    u = radec2dir(ra=ra, dec=dec, obstime_mjd=mjd)    

    # Add these directions to the DataFrame in three columns after dec
    col_num_dec = df.columns.get_loc('dec')
    df.insert(loc=col_num_dec+1, column='ux', value=u[0])
    df.insert(loc=col_num_dec+2, column='uy', value=u[1])
    df.insert(loc=col_num_dec+3, column='uz', value=u[2])

    # Save modified DataFrame to disk
    path: str = os.path.join(dir_name, file_name)
    df.to_hdf(path, key='df', mode='w')

# ********************************************************************************************************************* 
def load_ztf_det_all(verbose: bool = False):
    """
    Load all available ZTF detections
    INPUTS:
        verbose:  Whether to print status messages to console
    RETURNS:
        df: Pandas DataFrame of detections.  Includes columns ObjectID, CandidateID, mjd, ra, dec, asteroid_prob
    """
    # Check if file already exists.  If so, load it from memory and return early
    save_dir = '../data/ztf'
    file_path = os.path.join(save_dir, f'ztf-detections.h5')
    try:
        df = pd.read_hdf(file_path)
        if verbose:
            print(f'Loaded {file_path} from disk.')
        # Generate the unique times
        mjd_unq = np.unique(df.mjd.values)
        return df, mjd_unq
    except:
        if verbose:
            print(f'Assembling {file_name} from ZTF detections by year...')   
        # Load data for all years with available data
        dfs = []
        year0 = 2018  # first year of ZTF2 data
        year1 = datetime.datetime.today().year
        for year in range_inc(year0, year1):
            df = load_ztf_det_year(year=year)
            dfs.append(df)

        # Combine frames into one DataFrame
        df = pd.concat(dfs)

        # Add calculated directions; save file to disk; and return it
        ztf_det_add_dir(df=df, file_name='ztf-detections.h5', dir_name='../data/ztf')

    # Drop the superfluous magnitude columns to avoid confusion
    df.drop(columns=['mag_psf', 'mag_nr'], inplace=True)
    # Rename the column 'mag_app' to 'mag'
    df.rename(columns={'mag_app':'mag'})

    # Generate the TimeStampID column.  This is integer key counting the unique times of observations.
    mjd_unq, time_stamp_id = np.unique(df.mjd.values, return_inverse=True)
    # Convert time_stamp_id to 32 bit integer
    time_stamp_id = time_stamp_id.astype(np.int32)
    # Add TimeStampID to the DataFrame
    col_num_ts = df.columns.get_loc('CandidateID') + 1
    df.insert(loc=col_num_ts, column='TimeStampID', value=time_stamp_id)

    # Save assembled DataFrame to disk
    df.to_hdf(file_path, key='ztf', mode='w')

    # Return the DataFrame of observations and vector of unique observation times
    return df, mjd_unq

# ********************************************************************************************************************* 
def ztf_ast_file_name(n0: int, n1: Optional[int]):
    """File name for data file with ZTF observations and nearest asteroid."""
    # Handle special case that we want all asteroids (n0=0, n1=None)
    if (n0==0 and n1 is None):
        file_name = 'ztf-nearest-ast.h5'
    # Handle special case that we want all asteroids starting from n0 (n1=None)
    elif (n1 is None):
        file_name = f'ztf-nearest-ast-{n0:06d}.h5'
    # General case: n0 and n1 both specified
    else:
        file_name = f'ztf-nearest-ast-{n0:06d}-{n1:06d}.h5'
    return file_name

# ********************************************************************************************************************* 
def ztf_ast_file_path(n0: int, n1: int):
    """File path for data file with ZTF observations and nearest asteroid."""
    file_name = ztf_ast_file_name(n0=n0, n1=n1)   
    file_path = os.path.join(ztf_ast_dir_name, file_name)
    return file_path

# ********************************************************************************************************************* 
def ztf_calc_nearest_ast(ztf: pd.DataFrame, 
                         ast_dir: pd.DataFrame, 
                         thresh_deg: float = 180.0,
                         progbar: bool = False,
                         verbose: bool = False):
    """
    Calculate the nearest asteroid to each observation in the ZTF data.
    INPUTS:
        ztf:     DataFrame of ZTF observations.  
                 Columns used include TimeStampID, ux, uy, uz, nearest_ast_num, nearest_ast_dist
        ast_dir: DataFrame of theoretical asteroid directions from observatory with observations.
                 mjd must match up with thos of mjd_unq (unique, sorted time stamps in ztf)
        thresh:  Threshold in degrees to consider an observation close to an asteroid
        progbar: Whether to display a tqdm progress bar
        verbose:  Whether to print status messages to console
    OUTPUTS:
        DataFrame ztf with the columns for nearest asteroid
    """
    # Distinct asteroid numbers in ast_dir
    asteroid_nums = np.unique(ast_dir.asteroid_num)
    # Range of asteroid numbers
    # ast_min: int = np.min(asteroid_nums)
    # ast_max: int = np.max(asteroid_nums)

    # Convert threshold from degrees to magnitude of direction difference
    thresh_dist = np.sin(np.deg2rad(thresh_deg/2.0))*2.0

    # Add columns for nearest asteroid number and distance if not already present
    if 'nearest_ast_num' not in ztf.columns:
        ztf.insert(loc=ztf.columns.size, column='nearest_ast_num', value=np.int32(-1))
    if 'nearest_ast_dist' not in ztf.columns:
        ztf.insert(loc=ztf.columns.size, column='nearest_ast_dist', value=2.0)
    # Add columns for RA, DEC, and direction of nearest asteroid
    cols_nearest_ast_radec = ['ast_ra', 'ast_dec']
    cols_nearest_ast_dir = ['ast_ux', 'ast_uy', 'ast_uz']
    for col in (cols_nearest_ast_radec + cols_nearest_ast_dir):
        if col not in ztf.columns:
            ztf.insert(loc=ztf.columns.size, column=col, value=np.nan)
        
    # Extract directions of the ZTF observations as an Mx3 array
    cols_radec = ['ra', 'dec']
    cols_dir = ['ux', 'uy', 'uz']
    # radec_ztf = ztf[cols_radec].values
    u_ztf = ztf[cols_dir].values

    # Extract TimeStampID as M array
    row_num = ztf.TimeStampID.values

    # Check asteroid distance one at a time
    iterates = tqdm(asteroid_nums) if progbar else asteroid_nums
    for asteroid_num in iterates:
        # Mask for this asteroid in the ast_dir DataFrame
        mask_ast = (ast_dir.asteroid_num == asteroid_num)

        # Directions of this asteroid at the unique time stamps
        radec_ast_unq = ast_dir.loc[mask_ast, cols_radec].values
        u_ast_unq = ast_dir.loc[mask_ast, cols_dir].values

        # Difference bewteen ztf direction and this asteroid direction
        dist_i = np.linalg.norm(u_ztf - u_ast_unq[row_num], axis=1)

        # Is this the closest seen so far?
        mask_close = (dist_i < thresh_dist) & (dist_i < ztf.nearest_ast_dist)
        # Row numbers corresponding to close observations
        row_num_close = row_num[mask_close]
        # Update nearest asteroid number and distance columns
        ztf.nearest_ast_num[mask_close] = asteroid_num
        ztf.nearest_ast_dist[mask_close] = dist_i
        # Save the RA, DEC and direction columns
        ztf.loc[mask_close, cols_nearest_ast_radec] = radec_ast_unq[row_num_close]
        ztf.loc[mask_close, cols_nearest_ast_dir] = u_ast_unq[row_num_close]

    return ztf

# ********************************************************************************************************************* 
def load_ztf_nearest_ast(n0: int=0, 
                         n1: int=None,
                         dir_name: str = '../data/ztf_ast'):
    """
    Load the nearest asteroid to each observation in the ZTF data.
    INPUTS:
        n0:       First asteroid number to process, inclusive (e.g. 0)
        n1:       Last asteroid number to process, exclusive (e.g. 1000)
        dir_name: Directory with h5 file
    OUTPUTS:
        DataFrame ztf with the columns for nearest asteroid
    """
    # File path
    file_path = ztf_ast_file_path(n0=n0, n1=n1)

    # Load the data file if its available, and regeneration was not requested
    try:
        df = pd.read_hdf(file_path)
    except:
        raise ValueError(f'ztf_load_nearest_ast: unable to load file {file_path}.')
    return df

# ********************************************************************************************************************* 
def ztf_nearest_ast(ztf: pd.DataFrame, 
                    n0: int, 
                    n1: int,
                    thresh_deg: float = 180.0,
                    dir_name: str = '../data/ztf_ast',
                    regen: bool = False,
                    progbar: bool = False,
                    verbose: bool = False):
    """
    Load or calculate the nearest asteroid to each observation in the ZTF data.
    INPUTS:
        ztf:      DataFrame of ZTF observations.  
                  Columns used include TimeStampID, ux, uy, uz, nearest_ast_num, nearest_ast_dist
        n0:       First asteroid number to process, inclusive (e.g. 0)
        n1:       Last asteroid number to process, exclusive (e.g. 1000)
        thresh:   Threshold in degrees to consider an observation close to an asteroid
        dir_name: Directory with h5 file
        regen:    Flag; force regeneration of file whether or not on disk
        progbar:  Whether to display a tqdm progress bar
        verbose:  Whether to print status messages to console
    OUTPUTS:
        DataFrame ztf with the columns for nearest asteroid
    """
    # File name and path
    file_name = ztf_ast_file_name(n0=n0, n1=n1)
    file_path = os.path.join(dir_name, file_name)

    # Load the data file if its available, and regeneration was not requested
    if not regen:
        try:
            df = pd.read_hdf(file_path)
            if verbose:
                print(f'Loaded {file_path} from disk.')
            return df
        except:
            if verbose:
                print(f'Unable to load {file_path}, computing nearest asteroids from {n0} to {n1}...')
    else:
        if verbose:
            print(f'Regenerating {file_path}, computing nearest asteroids from {n0} to {n1}...')
    
    # If we get here, we need to build the ztf_ast DataFrame by calculation

    # Unique times in ztf data
    mjd_unq = np.unique(ztf.mjd)

    # Observatory for ZTF data is always Palomar Mountain
    site_name = 'palomar'

    # Build splined positions and observations at unique observation times
    ast_pos, earth_pos, ast_dir = \
        spline_ast_vec_dir(n0=n0, n1=n1, mjd=mjd_unq, site_name=site_name, progbar=progbar)

    # Calculate nearest asteroid in this block with  ztf_calc_nearest_ast
    ztf_ast = ztf_calc_nearest_ast(ztf=ztf, ast_dir=ast_dir, thresh_deg=thresh_deg, 
                                   progbar=progbar, verbose=verbose)

    # Save assembled DataFrame to disk and return it
    ztf_ast.to_hdf(file_path, key='ztf_ast', mode='w')
    return ztf_ast

# ********************************************************************************************************************* 
def calc_hit_freq(ztf, thresh_sec: float):
    """
    Calculate number of close hits by asteroid number
    INPUTS:
        ztf: DataFrame with distance to nearest asteroid
        thresh_sec: Threshold in arc seconds for close observations
    """
    # Threshold distance and flag indicating whether observations are within threshold
    thresh_deg = thresh_sec / 3600.0
    thresh_dist = deg2dist(thresh_deg)
    is_close = ztf.nearest_ast_dist < thresh_dist
    # View of ztf limited to close observations
    ztfc = ztf[is_close]

    # Group close observations by asteroid number
    close_by_ast = ztfc.groupby(ztfc.nearest_ast_num)
    close_by_ast_count = close_by_ast.size()
    ast_num = close_by_ast_count.index.values
    hit_count = close_by_ast_count.values

    # Return numbers and counts
    return ast_num, hit_count

# ********************************************************************************************************************* 
def make_ztf_easy_batch(batch_size: int = 64, thresh_sec: float = 10.0):
    """
    Generate an "easy batch" to prototype asteroid search algorithm.
    The easy batch consists of all ZTF observations whose nearest asteroid is
    one of the 64 asteroids with the most hits at a 2.0 arc second threshold.
    """
    # Load all ZTF observations including nearest asteroid
    ztf = load_ztf_nearest_ast()

    # Asteroid numbers and hit counts
    ast_num, hit_count = calc_hit_freq(ztf=ztf, thresh_sec=2.0)

    # Sort the hit counts in descending order and find the top batch_size
    idx = np.argsort(hit_count)[::-1][0:batch_size]

    # Extract the asteroid number and hit count for this batch
    ast_num_batch = ast_num[idx]
    hit_count_batch = hit_count[idx]

    # Mask for entries in this batch
    is_in_batch = (ztf.nearest_ast_num == 0)
    for ast_num in ast_num_batch:
        is_in_batch = is_in_batch | (ztf.nearest_ast_num == ast_num)
    ztf_batch = ztf[is_in_batch].copy()

    # Threshold to be close enough
    thresh_deg = thresh_sec / 3600.0
    thresh_dist = deg2dist(thresh_deg)
    is_close = ztf_batch.nearest_ast_dist < thresh_dist
    ztf_batch = ztf_batch[is_close]

    # Reset index so rows counted from zero
    ztf_batch.reset_index()

    return ztf_batch

