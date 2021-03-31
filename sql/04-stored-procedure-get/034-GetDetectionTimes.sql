DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionTimes()
COMMENT "Get all available detection times."

BEGIN 

SELECT
	dt.DetectionID,
	dt.DetectionTimeSliceID,
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
ORDER BY dt.DetectionID;

END
$$

-- ********************************************************************************
DELIMITER ;
