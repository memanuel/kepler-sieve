"""
ZTF Data Load
Load data from the ZTF2 database using the Alerce data broker into MSE database.

Michael S. Emanuel
27-Feb-2020
"""

# Libraries for getting Alerce data out of ZTF2 database
# import psycopg2
import sqlalchemy

# Standard libraries
import numpy as np
import pandas as pd

# Astronomy related
from astropy.units import deg

# Utility
import os
from itertools import count
from tqdm.auto import tqdm

# MSE imports
from db_utils import csv2db

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

# Location to save CSV files
data_dir = '../data/ztf'

# Columns for tables extracted
columns_cls = ['ObjectClassID', 'ObjectClassName']

columns_obj = \
        ['ObjectCD', 'ObservationCount', 'mjd0', 'mjd1', 'mean_ra', 'mean_dec', 
         'mean_mag_g', 'mean_mag_r', 'mean_mag_psf_g', 'mean_mag_psf_r', 'ObjectClassID', 'ClassificationProb']

columns_det = \
    ['DetectionID', 'ObjectCD', 'mjd', 'ra', 'dec', 'mag_psf', 'mag_app', 'mag_nr', 'sigma_ra', 'sigma_dec']

# ********************************************************************************************************************* 
def load_ztf_classes():
    """Load contents of the ZTF class table"""
    # SQL statement
    sql_stmt = \
    """
    SELECT
        cls.id AS ObjectClassID,
        cls.name AS ObjectClassName
    FROM 
        public.class AS cls
    ORDER BY ObjectClassID;
    """

    # Run query and return DataFrame
    with db_engine.connect() as conn:
        df = pd.read_sql(sql_stmt, conn)

    # Map column names (ZTF database has case insensitive column names)
    mapper = {col.lower() : col for col in columns_cls}
    df.rename(mapper=mapper, axis='columns', inplace=True)

    return df

# ********************************************************************************************************************* 
def load_ztf_objects(n0: int, sz: int, min_object_cd: str):
    """
    Load entries from ZTF objects table
    INPUTS:
        n0: Start of batch
        sz: Number of records to load
        min_object_cd: Minimum object code 
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
        obj.meanra AS mean_ra,
        obj.meandec AS mean_dec,
        obj.mean_magap_g AS mean_mag_g,	
        obj.mean_magap_r AS mean_mag_r,	
        obj.mean_magpsf_g AS mean_mag_psf_g,	
        obj.mean_magpsf_r AS mean_mag_psf_r,
        obj.classearly AS ObjectClassID,
        obj.pclassearly AS ClassificationProb        
    FROM 
        public.objects AS obj
    WHERE
        -- Only those greater than the last object code before we started this batch
        obj.oid > %(min_object_cd)s
    ORDER BY ObjectCD
    LIMIT %(limit)s OFFSET %(offset)s;
    """

    # Set up parameters
    params = {'limit':sz, 'offset':n0, 'min_object_cd':min_object_cd}

    # Run query and return DataFrame
    with db_engine.connect() as conn:
        df = pd.read_sql(sql_stmt, conn, params=params).fillna(-1.0)

    # Map column names (ZTF database has case insensitive column names)
    mapper = {col.lower() : col for col in columns_obj}
    df.rename(mapper=mapper, axis='columns', inplace=True)

    return df

# ********************************************************************************************************************* 
def load_ztf_detections(n0: int, sz: int, min_detection_id: int):
    """
    Load all the ZTF detections classified as asteroids in a date range.
    INPUTS:
        n0: Start of batch
        sz: Number of records to load
        min_object_id: Minimum detection ID
    RETURNS:
        df: Pandas DataFrame of detections.  Includes columns ObjectID, CandidateID, mjd, ra, dec, asteroid_prob
    """

    # SQL with parameters
    sql_stmt = \
    """
    SELECT
        det.candid as DetectionID,
        det.oid as ObjectCD,
        det.mjd,
        det.ra,
        det.dec,
        det.magpsf as mag_psf,
        det.magap as mag_app,
        det.magnr as mag_nr,
        det.sigmara as sigma_ra,
        det.sigmadec as sigma_dec
    FROM
        public.detections as det
    WHERE
        -- Only detections after the last one we've already loaded
        det.candid > %(min_detection_id)s
    ORDER BY DetectionID        
    LIMIT %(limit)s OFFSET %(offset)s;
    """

    # Set up parameters
    params = {'limit':sz, 'offset':n0, 'min_detection_id':min_detection_id}

    # Run query and return DataFrame
    with db_engine.connect() as conn:
        df = pd.read_sql(sql_stmt, conn, params=params).fillna(-1.0)

    # Map column names (ZTF database has case insensitive column names)
    mapper = {col.lower() : col for col in columns_det}
    df.rename(mapper=mapper, axis='columns', inplace=True)

    # Convert CandidateID column from object to np.int64
    df['DetectionID'] = df['DetectionID'].astype(np.int64)

    return df

# ********************************************************************************************************************* 
def main():
    """
    Main routine for console program
    """
    # Schema for permanent DB
    schema = 'ZTF'
    
    # File names for output
    fname_cls = os.path.join(data_dir, 'ObjectClass.csv')
    fname_obj = os.path.join(data_dir, 'Object.csv')
    fname_det = os.path.join(data_dir, 'Detection.csv')

    # For each table, three tasks:
    # (1) Extract the data
    # (2) Save as CSV file(s)
    # (3) Insert into KeplerDB.ZTF database schema

    # Process ObjectClass table
    df_cls = load_ztf_classes()
    df_cls.to_csv(fname_cls, index=False)
    csv2db(schema=schema, table='ObjectClass', columns=columns_cls, fname_csv=fname_cls)

    # The last object code already in the DB
    min_object_cd = str(sp2df('ZTF.GetObjectLast').LastObjectCD[0])
    print(f'Last ObjectCD found in ZTF.Object = {min_object_cd}.  Processing new objects...')

    # Process Object table in batches
    sz_obj = 100000
    # There are about 53.3 million objects according to the DB statistics
    nMax_obj: int = 55000000
    for n0 in tqdm(range(0, nMax_obj, sz_obj)):
        # df_obj = load_ztf_objects(n0=n0, sz=sz_obj, min_object_cd=min_object_cd)
        min_object_cd = str(sp2df('ZTF.GetObjectLast').LastObjectCD[0])
        df_obj = load_ztf_objects(n0=0, sz=sz_obj, min_object_cd=min_object_cd)
        df_obj.to_csv(fname_obj, index=False)
        csv2db(schema=schema, table='Object', columns=columns_obj, fname_csv=fname_obj)
        # Quit early if no rows in resultset
        if df_obj.shape[0] == 0:
            break

    # The last object code already in the DB
    min_detection_id= int(sp2df('ZTF.GetDetectionLast').LastDetectionID[0])
    print(f'Last DetectionID found in ZTF.Detection = {min_detection_id}.  Processing new objects...')

    # Process Detection table in batches
    sz_det = 100000
    # There are about 146.3 million detections according to the DB statistics
    nMax_det: int = 147000000
    for n0 in tqdm(range(0, nMax_det, sz_det)):
        # df_det = load_ztf_detections(n0=n0, sz=sz_det, min_detection_id=min_detection_id)
        min_detection_id= int(sp2df('ZTF.GetDetectionLast').LastDetectionID[0])
        df_det = load_ztf_detections(n0=0, sz=sz_det, min_detection_id=min_detection_id)
        df_det.to_csv(fname_det, index=False)
        csv2db(schema=schema, table='Detection', columns=columns_det, fname_csv=fname_det)
        # Quit early if no rows in resultset
        if df_det.shape[0] == 0:
            break

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()