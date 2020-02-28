"""
Harvard IACS Masters Thesis
ZTF2 (Alerce) Data
Utilities for acquiring data from the Alerce data broker ZTF2 database.

Michael S. Emanuel
27-Feb-2020
"""

# Standard libraries
import pandas as pd
import numpy as np
from astropy.units import deg
import os
from datetime import date
# import time
from tqdm.auto import tqdm

# Libraries for getting Alerce data out of ZTF2 database
import json
import psycopg2
from alerce.api import AlerceAPI

# MSE imports
from astro_utils import date_to_mjd
from ra_dec import radec2dir

# ********************************************************************************************************************* 
# Get credentials for ZTF2 connection
credentials_file = "../alerce/alercereaduser.json"
with open(credentials_file) as fh:
    cred = json.load(fh)["params"]

# Connect to ZTF2 database; this is a shared resource, just one for the module
conn = psycopg2.connect(dbname=cred['dbname'], user=cred['user'], host=cred['host'], password=cred['password'])

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
        # Create an empty h5 file to save into
        # print(f'Created empty {file_name} on disk; querying ZTF2 for data...')
        # pd.DataFrame().to_hdf(file_name, key='df', mode='w')
    
    # Start and end date
    mjd0 = date_to_mjd(date(year,1,1))
    mjd1 = date_to_mjd(date(year+1,1,1))

    # Sample dates at weekly intervals, with one long "week" at the end
    mjd = np.arange(mjd0, mjd1, 7, dtype=np.float64)
    if mjd[-1] < mjd1:
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
