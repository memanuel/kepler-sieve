"""
ZTF Data
Calculations to enrich ZTF data after it is loaded.

Michael S. Emanuel
24-Mar-2021
"""

# Standard libraries
import numpy as np
import pandas as pd

# Astronomy related
from astropy.units import deg

# UI
from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from db_utils import df2db, sp2df, sp_run
from astro_utils import mjd_to_datetime
from ra_dec import radec2dir, calc_topos
from rebound_sim import make_sim_planets
from rebound_integrate import integrate_mjds

# ********************************************************************************************************************* 
def calc_ztf_detection_times():
    """Update entries in table KS.DetectionTime associated with ZTF observations"""
    # Distinct ZTF detection times
    df = sp2df('ZTF.GetDetectionTimes')
    # Array of observations times as MJDs
    mjds = df['mjd'].values
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
def ztf_raw_detections(sz: int = 1000000):
    """Populate missing rows of table KS.RawDetections from ZTF.Detections"""
    # Get the number of missing rows
    rMax = sp2df('KS.MakeTable_RawDetection_ZTF_RowCount').RowCount[0]
    # Set up a tqdm index counter to process the rows in chunks of sz at a time
    iMax: int = rMax // sz + 1
    idx = tqdm_auto(np.arange(iMax))
    # Process the rows in chunks of sz
    params={'sz':sz}    
    print(f'KS.RawDetections missing {rMax} rows from ZTF.Detections. Processing now...')
    for i in idx:
        sp_run('KS.MakeTable_RawDetection_ZTF', params=params)

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
    col_num_dec = df.columns.get_loc('dec')
    df.insert(loc=col_num_dec+1, column='ux', value=u[0])
    df.insert(loc=col_num_dec+2, column='uy', value=u[1])
    df.insert(loc=col_num_dec+3, column='uz', value=u[2])

# ********************************************************************************************************************* 
def main():
    """
    Main routine for console program
    """

    # Call SQL procedure to add new rows to ZTF.DetectionTime from ZTF.Detection
    sp_run('ZTF.MakeTable_DetectionTime')    

    # Rebuild the KS.DetectionTime entries coming from ZTF
    calc_ztf_detection_times()

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
