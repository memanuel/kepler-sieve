"""
Harvard IACS Masters Thesis
ZTF Data
Utilities for acquiring data from the ZTF2 database using the Alerce data broker.

Michael S. Emanuel
27-Feb-2020
"""

# Libraries for getting Alerce data out of ZTF2 database
# import psycopg2
import sqlalchemy

# Standard libraries
import numpy as np
import pandas as pd

# # Astronomy related
# from astropy.units import deg
# from scipy.interpolate import CubicSpline

# # Utility
# import os
# import datetime
# from datetime import date
# from tqdm.auto import tqdm

# # MSE imports
# from utils import range_inc
# from astro_utils import date_to_mjd, deg2dist, dist2deg, dist2sec
# from ra_dec import radec2dir

# # Typing
# from typing import Optional, Dict

# ********************************************************************************************************************* 
# Global variables

# Create database engine for ZTF2
def make_db_engine():
    # Credentials for ZTF2 connection
    username = 'alerceread'
    password = 'alerce2019'
    hostname = 'db.alerce.online:5432'
    dbname = 'ztf_v2'
    # Build DB url and build the DB engine
    db_url: str = f'postgresql+psycopg2://{username}:{password}@{hostname}/{dbname}'
    db_engine: sqlalchemy.engine = sqlalchemy.create_engine(db_url, pool_size=32)
    return db_engine

# Connect to ZTF2 database; this is a shared resource, just one for the module
db_engine: sqlalchemy.engine = make_db_engine()

# ********************************************************************************************************************* 
def load_ztf_classes():
    """Load contents of the ZTF class table"""
    # SQL statement
    sql_stmt = \
    """
    SELECT
        cls.id AS ClassID,
        cls.name AS ClassName
    FROM 
        public.class AS cls
    ORDER BY cls.id;
    """

    # Run query and return DataFrame
    with db_engine.connect() as conn:
        df = pd.read_sql(sql_stmt, conn)

    # Map column names (ZTF database has case insensitive column names)
    cols = ['ClassID', 'ClassName']
    mapper = {col.lower() : col for col in cols}
    df.rename(mapper=mapper, axis='columns', inplace=True)

    # Convert ObjectCD column from object to string data type
    df['ClassName'] = df['ClassName'].astype('|S')

    return df

# ********************************************************************************************************************* 
def load_ztf_objects(n0: int, sz: int):
    """
    Load entries from ZTF objects table
    INPUTS:
        n0: Start of batch
        sz: Number of records to load
    
    RETURNS:
        Pandas DataFrame of the objects.
    """
    # SQL statement
    sql_stmt = \
    """
    SELECT
        obj.oid AS ObjectCD,
        obj.nobs AS ObservationCount,
        obj.firstmjd AS mjd0,
        obj.lastmjd AS mjd1,
        obj.meanra AS MeanRA,
        obj.meandec AS MeanDEC,
        obj.mean_magap_g AS MeanMag_g,	
        obj.mean_magap_r AS MeanMag_r,	
        obj.mean_magpsf_g AS MeanMagPSF_g,	
        obj.mean_magpsf_r AS MeanMagPSF_r,	
        obj.pclassearly AS AsteroidProb
    FROM 
        public.objects AS obj
    WHERE
        -- Objects classified as asteroids
        obj.classearly=21 AND
        -- with at least 2 observations
        obj.nobs > 2 AND
        -- With at least 0.95 probability of being asteroids
        obj.pclassearly > 0.95
    ORDER BY obj.oid
    LIMIT %(limit)s OFFSET %(offset)s;
    """

    # Set up parameters
    params = {'limit':sz, 'offset':n0}

    # Run query and return DataFrame
    with db_engine.connect() as conn:
        df = pd.read_sql(sql_stmt, conn, params=params)

    # Map column names (ZTF database has case insensitive column names)
    cols = ['ObjectCD', 'ObservationCount', 'MeanRA', 'MeanDEC', 
            'MeanMag_g', 'MeanMag_r', 'MeanMagPSF_g', 'MeanMagPSF_r', 'AsteroidProb']
    mapper = {col.lower() : col for col in cols}
    df.rename(mapper=mapper, axis='columns', inplace=True)

    # Convert ObjectCD column from object to string data type
    df['ObjectCD'] = df['ObjectCD'].astype('|S')

    return df

# ********************************************************************************************************************* 
def load_ztf_det(n0: int, sz: int):
    """
    Load all the ZTF detections classified as asteroids in a date range.
    INPUTS:
        n0: Start of batch
        sz: Number of records to load
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
        det.ra as RA,
        det.dec as DEC,
        det.magpsf as MagPSF,
        det.magap as MagApp,
        det.magnr as MagNR,
        det.sigmara as Sigma_RA,
        det.sigmadec as Sigma_DEC,
        obj.pclassearly as AsteroidProb
    from 
        detections as det
        inner join objects as obj on obj.oid = det.oid
    where
        %(mjd0)s <= det.mjd and det.mjd < %(mjd1)s and
        obj.classearly = %(asteroid_class)s
    """
    
    # Set up parameters
    params = {'limit':sz, 'offset':n0}
    
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
