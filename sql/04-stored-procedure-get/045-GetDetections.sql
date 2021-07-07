DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetMaxDetectionID()
COMMENT "Get ID of last asteroid detections; add 1 for compatibility with exclusive uper index."

BEGIN 

SELECT
	COALESCE(max(d.DetectionID), -1)+1 AS MaxDetectionID
FROM
	KS.Detection AS d;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetections(
    IN d0 INT,
    IN d1 INT
)
COMMENT "Get asteroid detections in a range of detection IDs."

BEGIN 

SELECT
	-- Integer ID fields
	d.DetectionID,
	d.SkyPatchID,
	d.TimeID,
	d.DetectionTimeID,
	-- Data payload
	d.mjd,
	d.ux,
	d.uy,
	d.uz,
	d.mag
FROM
	KS.Detection AS d
WHERE
	d.DetectionID BETWEEN d0 AND (d1-1)
ORDER BY d.DetectionID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionsObs(
    IN d0 INT,
    IN d1 INT
)
COMMENT "Get asteroid detections in a range of detection IDs; also return the observatory position."

BEGIN 

SELECT
	-- Integer ID fields
	d.DetectionID,
	d.SkyPatchID,
	d.TimeID,
	d.DetectionTimeID,
	-- Data payload
	d.mjd,
	d.ux,
	d.uy,
	d.uz,
	d.mag,
	-- Location of observatory - from join to DetectionTime
	dt.qObs_x,
	dt.qObs_y,
	dt.qObs_z
FROM
	KS.Detection AS d
	INNER JOIN KS.DetectionTime AS dt ON dt.DetectionTimeID = d.DetectionTimeID
WHERE
	d.DetectionID BETWEEN d0 AND (d1-1)
ORDER BY d.DetectionID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionCandidates(
    IN d0 INT,
    IN d1 INT
)
COMMENT "Get asteroid detections in a range of detection IDs; lightweight version with just integer IDs."

BEGIN 

SELECT
	-- Integer ID fields
	d.DetectionID,
	d.SkyPatchID,
	d.TimeID
FROM
	KS.Detection AS d
WHERE
	d.DetectionID BETWEEN d0 AND (d1-1)
ORDER BY d.DetectionID;

END
$$

-- ********************************************************************************
DELIMITER ;
