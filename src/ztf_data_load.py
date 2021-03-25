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
# from scipy.interpolate import CubicSpline

# Utility
import os
from itertools import count
from tqdm.auto import tqdm

# MSE imports
from db_utils import csv2db, df2db, sp2df
from astro_utils import mjd_to_datetime
from ra_dec import radec2dir, calc_topos
from rebound_sim import make_sim_planets
from rebound_integrate import integrate_mjds

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
        ['ObjectCD', 'ObservationCount', 'mjd0', 'mjd1', 'MeanRA', 'MeanDEC', 
         'MeanMag_g', 'MeanMag_r', 'MeanMagPSF_g', 'MeanMagPSF_r', 'ObjectClassID', 'ClassificationProb']

columns_det = \
    ['DetectionID', 'ObjectCD', 'mjd', 'RA', 'DEC', 'MagPSF', 'MagApp', 'MagNR', 'Sigma_RA', 'Sigma_DEC']

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
        obj.meanra AS MeanRA,
        obj.meandec AS MeanDEC,
        obj.mean_magap_g AS MeanMag_g,	
        obj.mean_magap_r AS MeanMag_r,	
        obj.mean_magpsf_g AS MeanMagPSF_g,	
        obj.mean_magpsf_r AS MeanMagPSF_r,
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
        det.ra as RA,
        det.dec as DEC,
        det.magpsf as MagPSF,
        det.magap as MagApp,
        det.magnr as MagNR,
        det.sigmara as Sigma_RA,
        det.sigmadec as Sigma_DEC
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
def calc_ztf_detection_times():
    """Update entries in table KS.DetectionTime associated with ZTF observations"""
    # Distinct ZTF detection times
    df = sp2df('ZTF.GetDetectionTimes')
    # Array of observations times as MJDs
    mjds = df['MJD'].values
    # Populate the CalendarDateTime field
    df['CalendarDateTime'] = np.array([mjd_to_datetime(mjd) for mjd in mjds])
    
    # All available data sources as a DataFrame
    ds = sp2df('KS.GetDataSources')
    # Populate DataSourceID and ObservatoryID fields
    df['DataSourceID'] = ds.DataSourceID[ds.DataSourceCD=='ZTF'].values[0]
    df['ObservatoryID'] = ds.ObservatoryID[ds.DataSourceCD=='ZTF'].values[0]

    # Integrate the planets saving outputs at these observation times
    print(f'Integrating planets on {df.shape[0]} distinct observation times...')
    sim_epoch = make_sim_planets(epoch=mjds[0])
    body_ids, body_names, q, v, elts = integrate_mjds(sim_epoch=sim_epoch, mjds=mjds, save_elements=False, progbar=True)

    # Earth position at these observation times
    earth_idx = np.argmax(body_names=='Earth')
    q_earth = q[:,earth_idx,:]
    v_earth = v[:,earth_idx,:]

    # Calculate topos adjustment
    dq_topos, dv_topos = calc_topos(obstime_mjd=mjds, site_name='Palomar')

    # The position and velocity of the observatory
    q_obs = q_earth + dq_topos.value # in AU
    v_obs = v_earth + dv_topos.value # in AU / day

    # Position of the Sun
    sun_idx = np.argmax(body_names=='Sun')
    q_sun = q[:,sun_idx,:]

    # Save positions of observatory and sun to DataFrame
    df[['qObs_x', 'qObs_y', 'qObs_z']] = q_obs
    df[['vObs_x', 'vObs_y', 'vObs_z']] = v_obs
    df[['qSun_x', 'qSun_y', 'qSun_z']] = q_sun

    # Save these observation times and positions to DB
    df2db(df=df, schema='KS', table='DetectionTime', verbose=False, progbar=False)

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

    # Call SQL procedure to add new rows to ZTF.DetectionTime from ZTF.Detection
    sp_run('ZTF.MakeTable_DetectionTime')    

    # Rebuild the KS.DetectionTime entries coming from ZTF
    calc_ztf_detection_times()

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
