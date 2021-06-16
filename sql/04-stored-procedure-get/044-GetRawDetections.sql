DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetRawDetectionIDs()
COMMENT "Get ID range of raw detection that have not yet been rolled up."

BEGIN 

SELECT
	MIN(rd.DetectionID) AS DetectionID_0,
	MAX(rd.DetectionID) AS DetectionID_1
FROM
	KS.RawDetection AS rd;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionIDs()
COMMENT "Get ID range of detection that have been rolled up."

BEGIN 

SELECT
	MIN(det.DetectionID) AS DetectionID_0,
	MAX(det.DetectionID) AS DetectionID_1
FROM
	KS.Detection AS det;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetRawDetectionDates()
COMMENT "Get date range of raw detection that have not yet been rolled up."

BEGIN 

SELECT
	CAST(FLOOR(MIN(dt.mjd)) AS INT) AS mjd0,
	CAST(CEILING(Max(dt.mjd)) AS INT) AS mjd1
FROM
	KS.DetectionTime AS dt
WHERE NOT EXISTS (
	SELECT d.DetectionID
	FROM KS.Detection AS d
	WHERE d.DetectionTimeID = dt.DetectionTimeID
    );

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetRawDetections(
    IN DetectionID_0 INT,
    IN DetectionID_1 INT
)
COMMENT "Get raw detections in the specified DetectionID range.  Only return detections that have not been rolled up into KS.Detections."

BEGIN 

SELECT
	rd.DetectionID,
	rd.DetectionTimeID,
	dt.mjd,
	rd.ra,
	rd.`dec`,
	rd.mag
FROM
	-- Start with RawDetection table
	KS.RawDetection AS rd
	-- Join DetectionTime table to get mjd of this detection
	INNER JOIN KS.DetectionTime AS dt ON dt.DetectionTimeID = rd.DetectionTimeID
WHERE
	-- Only selected range of DetectionIDs
	DetectionID_0 <= rd.DetectionID AND rd.DetectionID < DetectionID_1
    -- Only detections that have not already been rolled up into KS.Detection
    AND NOT EXISTS (
        SELECT d.DetectionID
        FROM KS.Detection AS d
        WHERE d.DetectionID = rd.DetectionID
    );
END
$$

-- ********************************************************************************
DELIMITER ;
