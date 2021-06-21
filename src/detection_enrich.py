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

# UI
from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from db_utils import df2db, sp2df, sp_run
from ra_dec import radec2dir
from orbital_element import unpack_vector
from sky_patch import dir2SkyPatchID, N_sp

# ********************************************************************************************************************* 
# Number of minutes in one day
mpd: int = 1440

# ********************************************************************************************************************* 
def detection_add_dir(df: pd.DataFrame):
    """
    Add calculated directions to DataFrame of ZTF observations
    INPUTS:
        df: DataFrame of ZTF observations including ra, dec and mjd columns
    OUTPUTS:
        Modifies df in place.
    """
    # Extract mjd, ra, and dec as numpy arrays; ra and dec flat arrays in degrees
    mjd = df.mjd.values
    ra = df.ra.values
    dec = df.dec.values

    # Compute the directions using radec2dir() and unpack it
    u = radec2dir(ra=ra, dec=dec, mjd=mjd)
    ux, uy, uz = unpack_vector(u)

    # Add these directions to the DataFrame in three columns after dec
    col_num_dec = df.columns.get_loc('dec')
    df.insert(loc=col_num_dec+1, column='ux', value=ux)
    df.insert(loc=col_num_dec+2, column='uy', value=uy)
    df.insert(loc=col_num_dec+3, column='uz', value=uz)

# ********************************************************************************************************************* 
def calc_detections(did0: int, did1: int, N_sp: int):
    """
    Assemble a batch of detections for insertion to Database from RawDetection data.
    INPUTS:
        did0: First DetectionID (inclusive)
        did1: Last DetectionID (exclusive)
        N_sp: Grid size used for SkyPatch
    RETURNS:
        df:   DataFrame ready to insert.
              Columns before enrichment: DetectionID, DetectionTimeID, mjd, ra, dec, mag
              Columns added: SkyPatchID, TimeID, k, ux, uy, uz
    """

    # Get DataFrame of raw detections
    params = {
        'did0': did0, 
        'did1': did1
    }
    df = sp2df('KS.GetRawDetections', params)

    # Handle edge case where there are no rows
    if df.shape[0]==0:
        return pd.DataFrame(columns = list(df.columns) + ['SkyPatchID', 'TimeID', 'ux', 'uy', 'uz',])

    # Calculate TimeID and add it to DataFrame
    time_id = np.floor(df.mjd.values*mpd).astype(np.int)
    df['TimeID'] = time_id

    # Calculate direction and add it to DataFrame
    detection_add_dir(df)
    # Nx3 array of directions
    dir = df[['ux', 'uy', 'uz']].values

    # Calculate and add the SkyPatchID
    SkyPatchID = dir2SkyPatchID(dir=dir, N=N_sp)
    df['SkyPatchID'] = SkyPatchID

    return df

# ********************************************************************************************************************* 
def write_detections():
    """
    Write to KS.Detection table from KS.RawDetection
    """
    # Columns to insert to DB
    columns = ['DetectionID', 'DetectionTimeID', 'TimeID', 'SkyPatchID', 'mjd', 'ux', 'uy', 'uz', 'mag']

    # All available (raw) Detection IDs
    did_raw = sp2df('KS.GetRawDetectionIDs')
    # DetectionIDs that have already been processed
    did_proc = sp2df('KS.GetDetectionIDs')
    # First and last DetectionID to process on this job
    # Start with the last processed DetectionID, plus 1; handle corner case where table is empty
    did0_job: int = (did_proc.DetectionID_1[0] or 0)+1
    # End with the last available RawDetectionID
    did1_job: int = (did_raw.DetectionID_1[0]+1 or 1)

    # Set batch size: number of detections
    b: int = 100000

    # Array of starting detection ID for each batch
    did0s = np.arange(did0_job, did1_job, b, dtype=np.int)
    # Status update
    print('Rolling up from KS.RawDetection to KS.Detection.')
    print('Enriching quoted RA/DEC with astrometric direction vector and SkyPatchID.')
    print(f'Processing detection IDs from {did0_job} to {did1_job}...')

    # Loop through the batches with a progress bar
    for did0 in tqdm_auto(did0s):
        # End Detection ID for this batch
        did1 = min(did0+b, did1_job)
        # Build DataFrame of enriched detections in this DetectionID range
        df = calc_detections(did0=did0, did1=did1, N_sp=N_sp)
        # Write enriched detections to DB table KS.Detection
        df2db(df=df, schema='KS', table='Detection', columns=columns, progbar=False)
        
# ********************************************************************************************************************* 
def write_tracklets(sz: int, start: int = 0):
    """Write all tracklets from the KS.Detection table.  Input sz sets batch size."""

    # Get range of DetectionTimePairIDs 
    dtpr = sp2df('KS.GetDetectionTimePairRange')
    dtp_min = dtpr.DetectionTimePairID_min[0]
    dtp_max = dtpr.DetectionTimePairID_max[0]
    # Apply start argument
    dtp_min = max(dtp_min, start)

    # Status
    print('Rolling up from KS.Detection to KS.Tracklet.')
    print(f'Processing DetectionTimePairIDs from {dtp_min} to {dtp_max}...')

    # Range of DetectionTimePairID
    dtp_ids = np.arange(dtp_min, dtp_max, sz, dtype=np.int32)

    # Loop through the batches
    for DetectionTimePairID_0 in tqdm_auto(dtp_ids):
        # The second DetectionID for this batch
        DetectionTimePairID_1 = DetectionTimePairID_0 + sz
        # Populate the tracklets by calling the DB stored procedure
        params = {
            'DetectionTimePairID_0': DetectionTimePairID_0, 
            'DetectionTimePairID_1': DetectionTimePairID_1}
        sp_run('KS.MakeTable_Tracklet', params)

# ********************************************************************************************************************* 
def main():
    """
    Main routine for console program
    """

    # Write from KS.RawDetection to KS.Detection
    write_detections()

    # # Get last DetectionPairID already loaded onto KS.Tracklet
    # ldtp = sp2df('KS.GetLastTrackletTimePair').LastDetectionTimePairID[0]
    # print(f'Last DetectionTimePairID loaded into KS.Tracklet is {ldtp}.')

    # # Write from KS.Detection to KS.Tracklet
    # write_tracklets(sz=500, start=ldtp)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
