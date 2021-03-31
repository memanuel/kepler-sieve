DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionTimePairs()
COMMENT "Get all available pairs of detection times."

BEGIN 

SELECT
	dtp.DetectionTimePairID,
	dtp.DetectionTimeID_1,
	dtp.DetectionTimeID_2,
	dtp.DataSourceID,
	dtp.mjd_bar,
	dtp.dt,
	dtp.mjd1,
	dtp.mjd2
FROM
	KS.DetectionTimePair AS dtp
ORDER BY dtp.DetectionTimePairID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionTimePairRange()
COMMENT "Get range detection time pair IDs."

BEGIN 

SELECT
	MIN(dtp.DetectionTimePairID) AS DetectionTimePairID_min,
	MAX(dtp.DetectionTimePairID) AS DetectionTimePairID_max
FROM
	KS.DetectionTimePair AS dtp;

END
$$

-- ********************************************************************************
DELIMITER ;
