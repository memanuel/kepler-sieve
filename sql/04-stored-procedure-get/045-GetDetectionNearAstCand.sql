DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionNearAstCand(
    IN AsteroidID_0 INT,
    IN AsteroidID_1 INT
)
COMMENT "Get candidate detections near asteroids in the given range of asteroids."

BEGIN 

-- Set variables for this session
SET max_heap_table_size = 16*(1024*1024*1024);
SET tmp_table_size = 16*(1024*1024*1024);

-- Range of asteroid
-- SET @AsteroidID_0 = 0;
-- SET @AsteroidID_1 = 10;

-- The table name
SET @k = CAST(AsteroidID_0 DIV 25000 AS INT);
SET @table_name = CONCAT('KS.AsteroidSkyPatch_Stage_', LPAD(@k, 2, '0'));

-- Staging table for SkyPatchID and time range by asteroid
CREATE OR REPLACE TEMPORARY TABLE KS.ASP_Candidate(
	SkyPatchID INT NOT NULL,
	AsteroidID INT NOT NULL,
	Segment INT NOT NULL,
	TimeID_1 INT NOT NULL,
	TimeID_2 INT NOT NULL,
	PRIMARY KEY (SkyPatchID, AsteroidID, Segment)
) ENGINE=Memory;

-- SQL statement template
SET @sql_str =
"
INSERT INTO KS.ASP_Candidate
SELECT DISTINCT
	spn.SkyPatchID_2 AS SkyPatchID,
	asp.AsteroidID,
	asp.Segment,
	asp.TimeID_0-15 AS TimeID_0,
	asp.TimeID_1
FROM
	-- Start with AsteroidSkyPatch
	@table_name AS asp
	-- Neighboring SkyPatches
	INNER JOIN KS.SkyPatchNeighbor AS spn ON
		spn.SkyPatchID_1 = asp.SkyPatchID
WHERE
	-- Only selected range of asteroids
	asp.AsteroidID BETWEEN @AsteroidID_0 AND (@AsteroidID_1-1);
";

-- Bind parameter values
SET @sql_str = REPLACE(@sql_str, '@AsteroidID_0', AsteroidID_0); 
SET @sql_str = REPLACE(@sql_str, '@AsteroidID_1', AsteroidID_1);
SET @sql_str = REPLACE(@sql_str, '@table_name', @table_name);

-- Batch of candidate skypatches including neighbors
PREPARE stmt FROM @sql_str;
EXECUTE stmt;
DEALLOCATE PREPARE stmt;

-- Batch of candidate skypatches including neighbors
-- INSERT INTO KS.ASP_Candidate
-- SELECT DISTINCT
-- 	spn.SkyPatchID_2 AS SkyPatchID,
-- 	asp.AsteroidID,
-- 	asp.Segment,
-- 	asp.TimeID_0-15 AS TimeID_0,
-- 	asp.TimeID_1
-- FROM
-- 	-- Start with AsteroidSkyPatch
-- 	KS.AsteroidSkyPatch_Stage_00 AS asp
-- 	-- Neighboring SkyPatches
-- 	INNER JOIN KS.SkyPatchNeighbor AS spn ON
-- 		spn.SkyPatchID_1 = asp.SkyPatchID
-- WHERE
-- 	-- Only selected range of asteroids
-- 	asp.AsteroidID BETWEEN AsteroidID_0 AND (AsteroidID_1-1);

-- Candidate detections matching these skypatches
SELECT DISTINCT
	det.DetectionID,
	asp.AsteroidID
FROM
	-- Start with AsteroidSkyPatch candidates
	temp.ASP_Candidate AS asp
	-- Matching detections
	INNER JOIN KS.Detection_v2 AS det ON
		det.SkyPatchID = asp.SkyPatchID AND
		det.TimeID BETWEEN asp.TimeID_1 AND asp.TimeID_2
ORDER BY det.DetectionID, asp.AsteroidID;

END
$$

-- ********************************************************************************
DELIMITER ;
