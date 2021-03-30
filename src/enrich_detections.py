"""
Calculations to enrich raw detections data after they are loaded into KS.RawDetections.
RawDetection has: DetectionID, DetectionTimeID, mjd, ra, dec, mag.
Detection has: DetectionTimeID, SkyPatchID, k, DetectionID, ux, uy, uz, mag

Michael S. Emanuel
29-Mar-2021
"""

# Standard libraries
import numpy as np
import pandas as pd

# Astronomy related
from astropy.units import deg

# UI
from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from db_utils import df2db, sp2df
from ra_dec import radec2dir
from sky_patch import dir2SkyPatchID

# ********************************************************************************************************************* 
def detection_add_dir(df: pd.DataFrame):
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
def calc_detections(mjd0: int, mjd1: int, N_sky_patch: int):
    """
    Assemble a batch of detections for insertion to Database from RawDetection data.
    INPUTS:
        mjd0: First date of detections (inclusive)
        mjd1: Last date of detections (exclusive)
        N_sky_patch: Grid size used for SkyPatch
    RETURNS:
        df:   DataFrame ready to insert.
              Columns before enrichment: DetectionTimeID, DetectionID, mjd, ra, dec, mag
              Columns added: ux, uy, uz, SkyPatchID, k
    """

    # Get DataFrame of raw detections
    params = {'mjd0':mjd0, 'mjd1': mjd1}
    df = sp2df('KS.GetRawDetections', params)

    # Calculate direction
    detection_add_dir(df)

    # Nx3 array of directions
    dir = df[['ux', 'uy', 'uz']].values

    # Calculate the SkyPatchID
    SkyPatchID = dir2SkyPatchID(dir=dir, N=N_sky_patch)

    # Add the SkyPatchID
    df['SkyPatchID'] = SkyPatchID

    # Add the field k for the detection number in each SkyPatch
    k = df.groupby(['DetectionTimeID', 'SkyPatchID']).cumcount() +1
    df['k'] = k

    return df

# ********************************************************************************************************************* 
def main():
    """
    Main routine for console program
    """

    # Set grid size for SkyPatch
    N_sky_patch: int = 1024

    # Columns to insert to DB
    columns = ['DetectionTimeID', 'SkyPatchID', 'k', 'DetectionID', 'mjd', 'ux', 'uy', 'uz', 'mag']

    # Get date range to roll up
    dts = sp2df('KS.GetRawDetectionDates')    
    mjd0_all = dts.mjd0[0]
    mjd1_all = dts.mjd1[0]
    mjd_width = mjd1_all - mjd0_all

    # Loop through dates
    sz: int = 10
    mjd_width= sz
    for i in tqdm_auto(range(0, mjd_width, sz)):
        # Date range for this loop
        mjd0 = mjd0_all + i * sz
        mjd1 = mjd0 + sz
        # Calculate the detections
        df = calc_detections(mjd0=mjd0, mjd1=mjd1, N_sky_patch=N_sky_patch)
        # Insert into DB table
        df2db(df=df, schema='KS', table='Detection', columns=columns)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()


