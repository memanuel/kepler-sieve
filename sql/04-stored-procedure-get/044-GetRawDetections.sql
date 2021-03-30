DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetRawDetections(
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get raw detections in the specified date range.  Only return detections that have not been rolled up into KS.Detections."

BEGIN 

SELECT
	rd.DetectionTimeID,
	rd.DetectionID,
	dt.mjd,
	rd.ra,
	rd.`dec`,
	rd.mag
FROM
	KS.RawDetection AS rd
	INNER JOIN KS.DetectionTime AS dt ON dt.DetectionTimeID = rd.DetectionTimeID
WHERE
	mjd0 <= dt.mjd AND dt.mjd < mjd1
    -- Only detections that have not already been rolled up into KS.Detection
    AND NOT EXISTS (
        SELECT d.DetectionID
        FROM KS.Detection AS d
        WHERE d.DetectionID = rd.DetectionID
    );

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetRawDetectionDates()
COMMENT "Get date range of raw detection that have not yet been rolled up."

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

-- ********************************************************************************
DELIMITER ;
