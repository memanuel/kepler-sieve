

from tqdm.auto import tqdm as tqdm_auto
from db_utils import sql_run

# SQL query
def make_sql(n: int):
	sql_str = \
	f"""
INSERT INTO KS.AsteroidSkyPatch
(AsteroidID, Segment, SkyPatchID, TimeID_0, TimeID_1)
SELECT
	asp.AsteroidID, 
	asp.Segment, 
	asp.SkyPatchID, 
	asp.TimeID_0, 
	asp.TimeID_1
FROM
	KS.AsteroidSkyPatch_Stage_{n:02d} AS asp;
"""
	return sql_str

# No parameters - empty map
params: dict = dict()

# Iterate over batches
for n in tqdm_auto(range(0, 39)):
    # SQL string for this i
	sql_str: str = make_sql(n=n)
	sql_run(sql_str=sql_str, params=params)
