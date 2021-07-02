DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionNearAstCand(
    IN AsteroidID_0 INT,
    IN AsteroidID_1 INT
)
COMMENT "Get candidate detections near asteroids in the given range of asteroids.  Brute force implementation"

BEGIN 

	SELECT
	dc.AsteroidID,
	dc.DetectionID
FROM
	KS.DetectionNearAsteroidCandidate AS dc
WHERE
	dc.AsteroidID BETWEEN AsteroidID_0 AND (AsteroidID_1-1)
ORDER BY dc.AsteroidID, dc.DetectionID;
	
END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionNearAstCand_bf(
    IN AsteroidID_0 INT,
    IN AsteroidID_1 INT,
    IN jn INT
)
COMMENT "Get candidate detections near asteroids in the given range of asteroids.  Brute force implementation"

BEGIN 

-- Set variables for this session
SET max_heap_table_size = 16*(1024*1024*1024);
SET tmp_table_size = 16*(1024*1024*1024);

-- The table name
SET @table_name = CONCAT('KS.AsteroidSkyPatch_Stage_', LPAD(jn, 2, '0'));

-- Staging table for SkyPatchID and time range by asteroid
CREATE OR REPLACE TEMPORARY TABLE KS.ASP_Candidate(
	SkyPatchID INT NOT NULL,
	AsteroidID INT NOT NULL,
	Segment INT NOT NULL,
	TimeID_1 INT NOT NULL,
	TimeID_2 INT NOT NULL,
	PRIMARY KEY (SkyPatchID, AsteroidID, Segment)
) ENGINE=Memory;

-- Staging table for (AsteroidID, DetectionID) pairs
CREATE OR REPLACE TEMPORARY TABLE KS.DNA_Candidate(
	AsteroidID INT NOT NULL,
	DetectionID INT NOT NULL,
	PRIMARY KEY (AsteroidID, DetectionID)
) ENGINE=Memory;

-- Batch of candidate skypatches including neighbors
INSERT INTO KS.ASP_Candidate
SELECT DISTINCT
	asp.AsteroidID,
	spn.SkyPatchID_2 AS SkyPatchID,
	asp.Segment,
	asp.TimeID_0-15 AS TimeID_0,
	asp.TimeID_1
FROM
	-- Start with AsteroidSkyPatch
	KS.AsteroidSkyPatch_Stage_00 AS asp
	-- Neighboring SkyPatches
	INNER JOIN KS.SkyPatchNeighbor AS spn ON
		spn.SkyPatchID_1 = asp.SkyPatchID
WHERE
	-- Only selected range of asteroids
	asp.AsteroidID BETWEEN AsteroidID_0 AND (AsteroidID_1-1);

-- Candidate detections matching these skypatches
INSERT INTO KS.DNA_Candidate
(AsteroidID, DetectionID)
SELECT DISTINCT
	asp.AsteroidID,
	det.DetectionID
FROM
	-- Start with AsteroidSkyPatch candidates
	KS.ASP_Candidate AS asp
	-- Matching detections
	INNER JOIN KS.Detection AS det ON
		det.SkyPatchID = asp.SkyPatchID AND
		det.TimeID BETWEEN asp.TimeID_1 AND asp.TimeID_2;

-- Query the selected detections along with supplementary information
SELECT
	-- Key fields: DetectionID and AsteroidID
	dc.DetectionID,
	dc.AsteroidID,
	-- Time of this detection
	det.mjd AS tObs,
	-- The position of the observatory
	dt.qObs_x,
	dt.qObs_y,
	dt.qObs_z,
	-- The direction of this detection
	det.ux AS uObs_x,
	det.uy AS uObs_y,
	det.uz AS uObs_z
FROM
	-- Start with the candidate AsteroidID, DetectionID pairs
	KS.DNA_Candidate AS dc
	-- Join Detection to get the position
	INNER JOIN KS.Detection AS det ON det.DetectionID = dc.DetectionID
	-- Join DetectionTime to get the observatory position
	INNER JOIN KS.DetectionTime AS dt ON dt.DetectionTimeID = det.DetectionTimeID
	ORDER BY dc.AsteroidID, dc.DetectionID;

END
$$

-- ********************************************************************************
DELIMITER ;
