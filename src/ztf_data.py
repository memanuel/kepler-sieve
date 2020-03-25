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
from astro_utils import date_to_mjd, deg2dist, dist2deg
from ra_dec import radec2dir
from asteroid_dataframe import calc_ast_data, spline_ast_vec_df, calc_ast_dir, spline_ast_vec_dir
from candidate_element import orbital_element_batch

# Typing
from typing import Optional, Dict

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
def ztf_reindex(ztf: pd.DataFrame, offset: np.int32 = 0) -> pd.DataFrame:
    """Reindex a ZTF DataFrame so it is sorted by (mjd, CandidateID)"""
    # Extract the two columns used for the sort: mjd and candidate_id
    mjd = ztf.mjd.values
    candidate_id = ztf.CandidateID.values
    # Perform a lexical sort: first by mjd, then by candidate_id
    index = np.lexsort(keys=(mjd, candidate_id))
    # Need to reset the index because there may be duplicates in the input
    ztf.reset_index(drop=True, inplace=True)
    # Reindex the dataframe
    ztf = ztf.reindex(index=index)
    # Sort the dataframe in place by the new index
    ztf.sort_index(inplace=True)
    # Add the offset to the index
    if offset != 0:
        ztf.index += offset
    # Return the sorted DataFrame
    return ztf

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
        
    # Assemble the weeks into one DataFrame
    df = pd.concat(df_list)
    # Reindex the DataFrame
    df = ztf_reindex(df)
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
        # Reindex the combined DataFrame; otherwise will have duplicated indices
        df = ztf_reindex(df)
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
def make_ztf_near_elt(ztf: pd.DataFrame, df_dir: pd.DataFrame, thresh_deg: float, 
                      R_deg: float = 0.20, progbar: bool=False) \
                      -> Dict[np.int32, pd.DataFrame]:
    """
    Assemble Python dict of DataFrames with ZTF observations near a batch of orbital elements.
    INPUTS:
        ztf:        DataFrame of candidate ZTF observations
        df_dir:     DataFrame of splined directions for elements at the unique observation times
        thresh_deg: Threshold for a close observation in degrees
        progbar:    Whether to display a progress bar
    """
    # Get unique element IDs
    element_ids = df_dir.element_id.values
    element_ids_unq = np.unique(element_ids)

    # Column collections used on ztf
    cols_catalog = ['ObjectID', 'CandidateID', 'TimeStampID', 'mjd']
    cols_radec = ['ra', 'dec']
    cols_dir = ['ux', 'uy', 'uz']
    cols_elt_dir = ['elt_ux', 'elt_uy', 'elt_uz']
    # Add nearest asteroid data to columns if available
    cols_nearest_ast_all = ['nearest_ast_num', 'nearest_ast_dist', 
                            'ast_ra', 'ast_dec', 'ast_ux', 'ast_uy', 'ast_uz']
    cols_nearest_ast = [col for col in cols_nearest_ast_all if col in ztf.columns]
    # All the output columns
    cols_out = cols_catalog + cols_radec + cols_dir + cols_nearest_ast
    # Flag indicating whether we have nearest asteroid data
    is_near_ast: bool = 'nearest_ast_num' in ztf.columns

    # Extract TimeStampID as (M,) array; name it row_num to emphasize that we use it to index into u_elt
    row_num = ztf.TimeStampID.values

    # Extract directions of the ZTF observations as an Mx3 array
    u_ztf = ztf[cols_dir].values

    # Threshold as Cartesian distance
    thresh_dist = deg2dist(thresh_deg)    

    # Threshold value of z = 1- s^2 / 2
    thresh_z = 1.0 - thresh_dist**2 / 2.0

    # Theshold for hits is 2.0 arc seconds (completely separate from threshold above, for inclusion)
    thresh_hit_sec = 2.0

    # Dictionary of DataFrames that are close
    ztf_tbl = dict()

    # Iterate over distinct element IDs
    iterates = tqdm(element_ids_unq) if progbar else element_ids_unq
    for element_id in iterates:
        # projected directions with this element id
        mask_elt = (df_dir.element_id == element_id)
        # Directions of this candidate element at the unique time stamps
        u_elt = df_dir.loc[mask_elt, cols_dir].values
        r_elt = df_dir.loc[mask_elt, 'delta'].values
        # Difference bewteen ztf direction and this element's direction
        dist_i = np.linalg.norm(u_ztf - u_elt[row_num], axis=1)
        # Which rows are within the threshold distance?
        mask_close = (dist_i < thresh_dist)
        # Row numbers corresponding to close observations
        row_num_close = row_num[mask_close]
        # Copy this slice
        ztf_i = ztf.loc[mask_close, cols_out].copy()
        # Insert columns with the ztf_id and element_id
        ztf_i.insert(loc=0, column='ztf_id', value=ztf_i.index.values)
        ztf_i.insert(loc=1, column='element_id', value=element_id)
        # Insert columns with the predicted directions of the elements
        for col in cols_elt_dir:
            ztf_i.insert(loc=ztf_i.columns.size, column=col, value=0.0)
        ztf_i[cols_elt_dir] = u_elt[row_num_close]
        ztf_i['elt_r'] = r_elt[row_num_close]
        # Save the distance in a column named s
        ztf_i['s'] = dist_i[mask_close]
        # Convert distance between element and observation to degrees
        # dist_deg = dist2deg(dist_i[mask_close])
        dist_deg = dist2deg(ztf_i.s)
        # Save distance in arc seconds
        ztf_i['dist_sec'] = dist_deg * 3600.0
        # Save the quantity z = 1 - s^2 / 2
        ztf_i['z'] = 1.0 - ztf_i.s**2 / 2.0
        # Save the quantity v = (1 - z) / (1 - thresh_z); for random observations, v will uniform on [0, 1]
        ztf_i['v'] = (1.0 - ztf_i.z) / (1.0 - thresh_z)
        # Compute score function with given resolution
        score_arg = -0.5 * (dist_deg / R_deg)**2
        ztf_i['score'] = np.exp(score_arg)
        # Flag for hits; an entry is a hit if these elements are within threshold of the ZTF observation
        ztf_i['is_hit'] = (ztf_i.dist_sec < thresh_hit_sec)
        # Flag for matches; a match is when the element_id matches the nearest asteroid
        # ztf_i['is_match'] = False
        if is_near_ast:
            ztf_i.is_match = (ztf_i.nearest_ast_num == ztf_i.element_id)
        # Save ztf_i to the table (dict) of elements
        ztf_tbl[element_id] = ztf_i

    # Combine rows into one big dataframe
    ztf_batch = pd.concat(ztf_tbl.values())

    # Reset the index so it does not have duplicate indices
    ztf_batch.reset_index(inplace=True, drop=True)

    return ztf_batch

# ********************************************************************************************************************* 
def make_ztf_batch(elts: pd.DataFrame, thresh_deg: float = 1.0, near_ast: bool = False):
    """
    Generate a ZTF batch with all ZTF observations within a threshold of the given elements
    INPUTS:
        elts:       Dataframe including element_id; 6 orbital elements
        thresh_deg: Threshold in degrees; only observations this close to elements are returned
        near_ast:   Whether to include data on the neareast asteroid
    """
    # Load all ZTF observations; include nearest asteroid data if requested
    ztf: pd.DataFrame
    if near_ast:
        # ZTF detections w/ nearest ast; need to regenerate unique times
        ztf = load_ztf_nearest_ast() 
        mjd = np.unique(ztf.mjd)
    else:
        # just the ZTF detections w/o nearest asteroid data
        ztf, mjd = load_ztf_det_all()

    # element_id is the unique identifier for each orbital element considered
    element_id = elts.element_id.values

    # Compute mjd0 and mjd1 from mjd_unq
    mjd0 = np.floor(np.min(mjd))
    mjd1 = np.ceil(np.max(mjd))

    # Calculate positions in this date range, sampled daily
    df_ast_daily, df_earth_daily, df_sun_daily = calc_ast_data(elts=elts, mjd0=mjd0, mjd1=mjd1, element_id=element_id)

    # Spline positions at ztf times
    df_ast, df_earth, df_sun = spline_ast_vec_df(df_ast=df_ast_daily, df_earth=df_earth_daily, df_sun=df_sun_daily, 
                                                 mjd=mjd, include_elts=False)
    # Direction from palomar
    df_dir = calc_ast_dir(df_ast=df_ast, df_earth=df_earth, site_name='palomar')

    # Calculate subset of ZTF data within threshold of this batch
    progbar = True
    ztf_batch = make_ztf_near_elt(ztf=ztf, df_dir=df_dir, thresh_deg=thresh_deg, progbar=progbar)

    return ztf_batch

# ********************************************************************************************************************* 
def elt_hash(elts: pd.DataFrame, thresh_deg: float = 1.0, near_ast: bool = False):
    """
    Load or generate a ZTF batch with all ZTF observations within a threshold of the given elements
    INPUTS:
        elts:       Dataframe including element_id; 6 orbital elements
        thresh_deg: Threshold in degrees; only observations this close to elements are returned
        near_ast:   Whether to include data on the neareast asteroid
    OUTPUTS:
        hash_id:    Unique ID for these inputs
    """
    # Columns of the Dataframe to hash
    cols_to_hash = ['element_id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch']
    # Tuple of int64; one per orbital element candidate
    hash_df = tuple((pd.util.hash_pandas_object(elts[cols_to_hash])).values)
    # Combine the element hash tuple with the threshold and near_ast flags
    thresh_int = int(thresh_deg*2**48)
    hash_id = hash(hash_df + (thresh_int, near_ast,))

    return hash_id
# ********************************************************************************************************************* 
def load_ztf_batch(elts: pd.DataFrame, thresh_deg: float = 1.0, near_ast: bool = False):
    """
    Load or generate a ZTF batch with all ZTF observations within a threshold of the given elements
    INPUTS:
        elts:       Dataframe including element_id; 6 orbital elements
        thresh_deg: Threshold in degrees; only observations this close to elements are returned
        near_ast:   Whether to include data on the neareast asteroid
    OUTPUTS:
        ztf_elt:    Dataframe including ZTF observation (mjd, ra, dec, ...)
                    elt_ux, elt_uy, elt_uz, 
    """
    # Get hash of arguments
    hash_id = elt_hash(elts=elts, thresh_deg=thresh_deg, near_ast=near_ast)

    # Name of file
    file_path = f'../data/ztf_elt/ztf_elt_{hash_id}.h5'

    # Try to load file if available
    try:
        ztf_elt = pd.read_hdf(file_path)
    # Generate it on the fly if it's not available
    except FileNotFoundError:
        ztf_elt = make_ztf_batch(elts=elts, thresh_deg=thresh_deg, near_ast=near_ast)
        ztf_elt.to_hdf(file_path, key='ztf_elt', mode='w')
    
    return ztf_elt

# ********************************************************************************************************************* 
def make_ztf_easy_batch(batch_size: int = 64, thresh_deg: float = 1.0):
    """
    Generate an "easy batch" to prototype asteroid search algorithm.
    The easy batch consists of all ZTF observations whose nearest asteroid is
    one of the 64 asteroids with the most hits at a 2.0 arc second threshold.
    """
    # Load all ZTF observations including nearest asteroid
    ztf_ast = load_ztf_nearest_ast()

    # Asteroid numbers and hit counts
    ast_num, hit_count = calc_hit_freq(ztf=ztf_ast, thresh_sec=2.0)

    # Sort the hit counts in descending order and find the top batch_size
    idx = np.argsort(hit_count)[::-1][0:batch_size]

    # Extract the asteroid number and hit count for this batch
    ast_num_best = ast_num[idx]
    # hit_count_best = hit_count[idx]

    # Orbital elements for best asteroids (dict of numpy arrays)
    element_id = np.sort(ast_num_best)
    elts = orbital_element_batch(element_id)

    # Delegate to make_ztf_batch
    ztf_batch = make_ztf_batch(elts=elts, thresh_deg=thresh_deg, near_ast=True)
    return ztf_batch, elts

# ********************************************************************************************************************* 
def load_ztf_easy_batch(batch_size: int = 64, thresh_deg: float = 1.0):
    """
    Load the ztf easy batch if available. Otherwise generate it.
    """

    # Name of the file
    thresh_sec = np.round(thresh_deg*3600.0).astype(np.int32)
    file_path = f'../data/ztf/ztf-easy-batch-n={batch_size}-thresh={thresh_sec}sec.h5'

    # Try to load file if available, otherwise generate it on the fly and save it
    try:
        ztf_batch = pd.read_hdf(file_path, key='ztf_batch')
        elts = pd.read_hdf(file_path, key='elts')
    except:
        ztf_batch, elts = make_ztf_easy_batch(batch_size=batch_size, thresh_deg=thresh_deg)
        ztf_batch.to_hdf(file_path, key='ztf_batch', mode='w')
        elts.to_hdf(file_path, key='elts', mode='a')
    return ztf_batch, elts

# ********************************************************************************************************************* 
def report_ztf_score(ztf, display: bool = True):
    """Report number of hits and breakdown of score between hits and noise on a ZTF batch"""
    # Number of elements in this data set
    num_elts = np.unique(ztf.element_id).size
    num_rows = ztf.shape[0]
    mean_rows = num_rows / num_elts

    # Count hits
    num_hits = np.sum(ztf.is_hit)
    mean_hits = num_hits / num_elts
    hit_pct = num_hits / num_rows * 100.0

    # Count matches
    num_matches = np.sum(ztf.is_match)
    mean_matches = num_matches / num_elts
    match_pct = num_matches / num_rows * 100.0

    # Mean score by category
    mean_score = np.sum(ztf.score) / num_elts
    mean_score_hits = np.sum(ztf[ztf.is_hit].score) / num_elts
    mean_score_match = np.sum(ztf[ztf.is_match].score) / num_elts
    mean_score_noise = np.sum(ztf[~ztf.is_match].score) / num_elts

    # Score breakdown by type
    score_hit_pct = mean_score_hits / mean_score * 100.0
    score_match_pct = mean_score_match / mean_score * 100.0
    score_noise_pct = mean_score_noise / mean_score * 100.0

    # Aggregation by element_id
    ztf['score_match'] = ztf.score * ztf.is_match
    ztf['score_noise'] = ztf.score * (~ztf.is_match)
    cols_ztf = ['score', 'score_match', 'score_noise', 'is_hit', 'is_match']
    score_by_elt = ztf[cols_ztf].groupby(by=ztf.element_id).sum()
    score_agg= score_by_elt.agg(['mean', 'std'])

    # Report
    if display:
        print(score_agg)
        print(f'\nHits and Matches per element; Percentage of Rows:')
        print(f'Hits:    {mean_hits:6.2f} / {hit_pct:5.2f}%')
        print(f'Matches: {mean_matches:6.2f} / {match_pct:5.2f}%')
        print(f'\nScore Contribution: Match vs. Noise')
        print(f'Match: {mean_score_match:6.2f} / {score_match_pct:5.2f}%')
        print(f'Noise: {mean_score_noise:6.2f} / {score_noise_pct:5.2f}%')

    return score_by_elt, score_agg