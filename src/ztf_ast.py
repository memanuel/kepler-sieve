"""
Harvard IACS Masters Thesis
ZTF Data
Calculate distances between ZTF observations and known asteroids.
Load DataFrame of ZTF observations augmented with nearest asteroid information.

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
from astro_utils import date_to_mjd, deg2dist, dist2deg, dist2sec
from ra_dec import radec2dir
from asteroid_dataframe import calc_ast_data, spline_ast_vec_df, calc_ast_dir, spline_ast_vec_dir

# Typing
from typing import Optional, Dict

# ********************************************************************************************************************* 
# Global variables

# Directory for files with nearest asteroid to ztf observations
ztf_ast_dir_name = '../data/ztf_ast'

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
