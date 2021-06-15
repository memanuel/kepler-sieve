SELECT
	asp.AsteroidID,
	asp.Segment,
	asp.SkyPatchID,
	asp.TimeID_0,
	asp.TimeID_1
FROM
	-- Start with AsteroidSkyPatch
	KS.AsteroidSkyPatch AS asp
	-- Detection times in this interval
	-- INNER JOIN KS.DetectionTime AS dt ON dt.
WHERE
	-- Only selected range of asteroids
	asp.AsteroidID BETWEEN 1 AND 1;