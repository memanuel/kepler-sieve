

from tqdm.auto import tqdm as tqdm_auto
from db_utils import sql_run

# SQL query
sql_str = \
"""
INSERT INTO KS.Detection_temp
(DetectionID, DetectionTimeID, SkyPatchID, TimeID, mjd, ux, uy, uz, mag)
SELECT
	det.DetectionID,
	det.DetectionTimeID,
	det.SkyPatchID,
	CAST(FLOOR(det.mjd*1440) AS INT) AS TimeID,
	det.mjd,
	det.ux,
	det.uy,
	det.uz,
	det.mag
FROM
	KS.Detection_v1 AS det
WHERE det.DetectionID BETWEEN :did0 AND :did1-1;
"""

# Batch size
sz = 100000
# Range of detection_id to process
did0_job = 56000000
did1_job = 159000000
# Range of n
n0 = did0_job // sz
n1 = did1_job // sz

# Iterate over batches
for n in tqdm_auto(range(n0, n1)):
    # Range of DetectionID
    did0 = n*sz
    did1 = did0 + sz
    # Bind parameters and run SQL
    params = {
        'did0': did0,
        'did1': did1,
    }
    sql_run(sql_str=sql_str, params=params)
