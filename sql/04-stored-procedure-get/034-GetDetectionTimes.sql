DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionTimes()
COMMENT "Get all available detection times."

BEGIN 

SELECT
	dt.DetectionTimeID,
	dt.HiResTimeID AS TimeID,
	dt.mjd,
	dt.DataSourceID,
	dt.ObservatoryID,
	dt.qObs_x,
	dt.qObs_y,
	dt.qObs_z,
	dt.qSun_x,
	dt.qSun_y,
	dt.qSun_z	
FROM
	KS.DetectionTime AS dt
ORDER BY dt.DetectionTimeID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionTimeMaxID()
COMMENT "Get maximum DetectionTimeID."

BEGIN 

SELECT
	max(dt.DetectionTimeID) AS MaxDetectionTimeID	
FROM
	KS.DetectionTime AS dt;

END
$$

-- ********************************************************************************
DELIMITER ;
