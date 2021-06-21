EXPLAIN
SELECT 
	asp.AsteroidID,
	asp.Segment,
	spn.SkyPatchID_2 AS SkyPatchID,
	asp.TimeID_0 - 15 AS TimeID_0,
	asp.TimeID_1
FROM 
	-- Start with AsteroidSkyPatch
	KS.AsteroidSkyPatch AS asp
	-- Neighboring sky patches
	INNER JOIN KS.SkyPatchNeighbor AS spn ON
		spn.SkyPatchID_1 = asp.SkyPatchID
WHERE
	asp.AsteroidID = 1
LIMIT 100;