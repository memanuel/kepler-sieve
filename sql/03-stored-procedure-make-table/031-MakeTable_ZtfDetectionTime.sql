DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE ZTF.MakeTable_DetectionTime()
COMMENT "Populate the ZTF.DetectionTime table from ZTF.Detection (all distinct times with detections)"
BEGIN 

-- Populate ZTF.DetectionTime from ZTF.Detection
INSERT INTO ZTF.DetectionTime
(mjd)
SELECT
    det.mjd
FROM
    ZTF.Detection AS det
    -- Corresponding time already on ZTF.DetectionTime table
    LEFT JOIN ZTF.DetectionTime AS dt ON dt.mjd = det.mjd
-- Only detection times not already present
WHERE dt.DetectionTimeID IS NULL
-- We want distinct detection times only
GROUP BY det.mjd
ORDER BY det.mjd;

-- Populate KS.DetectionTimeSlice using a fixed size of 1 minute
INSERT INTO KS.DetectionTimeSlice
(DetectionTimeSliceID, IntegrationTimeID, mjd, mjd0, mjd1)
WITH t1 AS (
SELECT
	CAST(FLOOR(dt.mjd*24*60) AS INT) AS DetectionTimeSliceID,
	CAST(FLOOR(dt.mjd*24*12)*5 AS INT) AS IntegrationTimeID,
	-- Midpoint of the interval is 30 seconds MJD0 and 30 seconds before MJD1
	(FLOOR(dt.mjd*24*60)+0.5)/(24*60) AS mjd,
	-- Start of the interval is the floor of the number of minutes
	FLOOR(dt.mjd*24*60)/(24*60) AS mjd0,
	(FLOOR(dt.mjd*24*60)+1)/(24*60) AS mjd1
FROM
	ZTF.DetectionTime AS dt
)
SELECT
	t1.DetectionTimeSliceID,
	t1.IntegrationTimeID,
	t1.mjd,
	t1.mjd0,
	t1.mjd1
FROM
	t1
WHERE NOT EXISTS(
	SELECT dts.DetectionTimeSliceID
	FROM KS.DetectionTimeSlice AS dts
	WHERE dts.DetectionTimeSliceID = t1.DetectionTimeSliceID
)
GROUP BY t1.DetectionTimeSliceID;

END
$$

DELIMITER ;
