SELECT 
	asp.AsteroidID,
	det.DetectionID,
-- 	asp.SkyPatchID,
	det.mjd,
	det.ux,
	det.uy,
	det.uz
FROM 
	-- Start with AsteroidSkyPatch
	KS.AsteroidSkyPatch AS asp
	-- Neighboring sky patches
	INNER JOIN KS.SkyPatchNeighbor AS spn ON
		spn.SkyPatchID_1 = asp.SkyPatchID
	-- Matching detections
	INNER JOIN KS.Detection AS det ON
		det.SkyPatchID = spn.SkyPatchID_2 AND
		det.TimeID BETWEEN asp.TimeID_0-15 AND asp.TimeID_1
WHERE
	-- Only selected range of asteroids
	(asp.AsteroidID BETWEEN 1 AND 10);