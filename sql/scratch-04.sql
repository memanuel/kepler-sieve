SELECT 
	asp.AsteroidID,
	asp.Segment
FROM 
	-- Start with AsteroidSkyPatch
	KS.AsteroidSkyPatch AS asp
WHERE
	asp.AsteroidID = 1
LIMIT 100;

SELECT
	spn.SkyPatchID_1,
	spn.SkyPatchID_2,
	(360*3600/pi())*ASIN(spn.dr_mid/2) AS dr_sec
FROM
	KS.SkyPatchNeighbor AS spn
LIMIT 100;	

SELECT ASIN(1.0) AS x;
