CREATE OR REPLACE TABLE temp.DNA_Candidate(
	DetectionTimeID INT NOT NULL,
	SkyPatchID INT NOT NULL,
	AsteroidID INT NOT NULL,
	PRIMARY KEY (DetectionTimeID, SkyPatchID, AsteroidID)
);

INSERT INTO temp.DNA_Candidate
(DetectionTimeID, SkyPatchID, AsteroidID)
SELECT 
	dt.DetectionTimeID,
	-- asp.SkyPatchID,
	spn.SkyPatchID_2,
	asp.AsteroidID
FROM 
	-- Start with AsteroidSkyPatch
	KS.AsteroidSkyPatch AS asp
	-- Get all neighboring skypatches
	INNER JOIN KS.SkyPatchNeighbor AS spn ON spn.SkyPatchID_1 = asp.SkyPatchID
	-- Matching detection times
	INNER JOIN KS.DetectionTime AS dt ON dt.HiResTimeID BETWEEN asp.TimeID_0 AND asp.TimeID_1
WHERE
	-- Only selected range of asteroids
	(asp.AsteroidID BETWEEN 1 AND 1)
	AND asp.TimeID_0 < 58010;