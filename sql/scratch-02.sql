CREATE OR REPLACE TABLE temp.ASP_Candidate(
	AsteroidID INT NOT NULL,
	Segment INT NOT NULL,
	SkyPatchID INT NOT NULL,
	TimeID_0 INT NOT NULL,
	TimeID_1 INT NOT NULL,
	PRIMARY KEY (AsteroidID, Segment),
	INDEX (SkyPatchID, TimeID_0, TimeID_1)
) engine memory;


INSERT INTO temp.ASP_Candidate
(AsteroidID, Segment, SkyPatchID, TimeID_0, TimeID_1)
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
	asp.AsteroidID = 1;
