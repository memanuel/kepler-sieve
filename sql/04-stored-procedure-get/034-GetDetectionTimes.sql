DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionTimes()
COMMENT "Get the first and last DetectionTimeID from the DetectionTime."

BEGIN 

SELECT
	COALESCE(MIN(dt.DetectionTimeID),0) AS DetectionTimeID_min,
	COALESCE(MAX(dt.DetectionTimeID),0) AS DetectionTimeID_max
FROM
	KS.DetectionTime AS dt;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetTrackletDetectionTimes()
COMMENT "Get the first and last DetectionTimeID written to the Tracklet table."

BEGIN 

SELECT
	COALESCE(MIN(dt.DetectionTimeID),0) AS DetectionTimeID_min,
	COALESCE(MAX(dt.DetectionTimeID),0) AS DetectionTimeID_max
FROM
	KS.Tracklet AS t
INNER JOIN KS.Detection AS d ON d.DetectionID = t.DetectionID_1
INNER JOIN KS.DetectionTime AS dt ON dt.DetectionTimeID = d.DetectionTimeID;

END
$$

-- ********************************************************************************
DELIMITER ;
