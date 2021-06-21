EXPLAIN
SELECT 
	asp.AsteroidID,
	det.DetectionID,
	det.mjd,
	det.ux,
	det.uy,
	det.uz
FROM 
	-- Start with AsteroidSkyPatch
	KS.AsteroidSkyPatch AS asp
	-- Matching detections
	INNER JOIN KS.Detection AS det ON
		det.SkyPatchID = asp.SkyPatchID AND
		det.TimeID BETWEEN asp.TimeID_0-15 AND asp.TimeID_1
WHERE
	-- Only selected range of asteroids
	(asp.AsteroidID BETWEEN 1 AND 10);