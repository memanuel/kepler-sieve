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

# Astronomy related
from astropy.units import deg
# from scipy.interpolate import CubicSpline

# Utility
import os
from tqdm.auto import tqdm

# MSE imports
# from utils import range_inc
# from astro_utils import date_to_mjd, deg2dist, dist2deg, dist2sec
from ra_dec import radec2dir
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
columns_class = ['ClassID', 'ClassName']

columns_object = \
        ['ObjectCD', 'ObservationCount', 'MeanRA', 'MeanDEC', 
          'MeanMag_g', 'MeanMag_r', 'MeanMagPSF_g', 'MeanMagPSF_r', 'AsteroidProb']

columns_detection = \
    ['DetectionID', 'ObjectCD', 'mjd', 'RA', 'DEC', 
    'MagPSF', 'MagApp', 'MagNR', 'Sigma_RA', 'Sigma_DEC', 'AsteroidProb']

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
    ORDER BY ClassID;
    """

    # Run query and return DataFrame
    with db_engine.connect() as conn:
        df = pd.read_sql(sql_stmt, conn)

    # Map column names (ZTF database has case insensitive column names)
    mapper = {col.lower() : col for col in columns_class}
    df.rename(mapper=mapper, axis='columns', inplace=True)

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
    ORDER BY ObjectCD
    LIMIT %(limit)s OFFSET %(offset)s;
    """

    # Set up parameters
    params = {'limit':sz, 'offset':n0}

    # Run query and return DataFrame
    with db_engine.connect() as conn:
        df = pd.read_sql(sql_stmt, conn, params=params)

    # Map column names (ZTF database has case insensitive column names)
    mapper = {col.lower() : col for col in columns_object}
    df.rename(mapper=mapper, axis='columns', inplace=True)

    return df

# ********************************************************************************************************************* 
def load_ztf_detections(n0: int, sz: int):
    """
    Load all the ZTF detections classified as asteroids in a date range.
    INPUTS:
        n0: Start of batch
        sz: Number of records to load
    RETURNS:
        df: Pandas DataFrame of detections.  Includes columns ObjectID, CandidateID, mjd, ra, dec, asteroid_prob
    """

    # SQL with parameters
    sql_stmt = \
    """
    SELECT
        det.candid as DetectionID,
        obj.oid as ObjectCD,
        det.mjd,
        det.ra as RA,
        det.dec as DEC,
        det.magpsf as MagPSF,
        det.magap as MagApp,
        det.magnr as MagNR,
        det.sigmara as Sigma_RA,
        det.sigmadec as Sigma_DEC,
        obj.pclassearly as AsteroidProb
    FROM
        public.detections as det
        INNER JOIN public.objects AS obj ON obj.oid = det.oid
    ORDER BY DetectionID
    LIMIT %(limit)s OFFSET %(offset)s;
    """

    # Set up parameters
    params = {'limit':sz, 'offset':n0}

    # Run query and return DataFrame
    with db_engine.connect() as conn:
        df = pd.read_sql(sql_stmt, conn, params=params)

    # Map column names (ZTF database has case insensitive column names)
    mapper = {col.lower() : col for col in columns_detection}
    df.rename(mapper=mapper, axis='columns', inplace=True)

    # Convert CandidateID column from object to np.int64
    df['CandidateID'] = df['CandidateID'].astype(np.int64)

    return df

# ********************************************************************************************************************* 
def ztf_detection_add_dir(df: pd.DataFrame):
    """
    Add calculated directions to DataFrame of ZTF observations
    INPUTS:
        df: DataFrame of ZTF observations including ra, dec and mjd columns
    OUTPUTS:
        Modifies df in place.
    """
    # Extract mjd, ra, and dec as vectors of astropy angles
    mjd = df.mjd.values
    ra = df.ra.values * deg
    dec = df.dec.values * deg

    # Compute the directions using radec2dir()
    u = radec2dir(ra=ra, dec=dec, obstime_mjd=mjd)    

    # Add these directions to the DataFrame in three columns after dec
    col_num_dec = df.columns.get_loc('DEC')
    df.insert(loc=col_num_dec+1, column='ux', value=u[0])
    df.insert(loc=col_num_dec+2, column='uy', value=u[1])
    df.insert(loc=col_num_dec+3, column='uz', value=u[2])

# ********************************************************************************************************************* 
def main():
    """
    Main routine for console program
    """
    # Schema for permanent DB
    schema = 'ZTF'
    
    # File names for output
    fname_class = os.path.join(data_dir, 'ObjectClass.csv')

    # Column names for each type
    columns_class = ['']

    # Extract the object classes from Alerce's ZTF database, save as CSV, and insert into MSE database
    df_cls = load_ztf_classes()
    df_cls.to_csv(fname_class, index=False)
    # csv2db(schema=schema, table='ObjectClass', columns=columns_class, fname_csv=fname_class)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
