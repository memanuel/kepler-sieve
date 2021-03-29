DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE ZTF.GetDetectionTimes()
COMMENT "Get all the distinct detection times ZTF.Detection table."

BEGIN 

SELECT
	dt.DetectionTimeID,
	dt.mjd,
	dts.DetectionTimeSliceID
FROM
	ZTF.DetectionTime AS dt
	INNER JOIN KS.DetectionTimeSlice AS dts ON
		dts.DetectionTimeSliceID = floor(dt.mjd*24*60)
ORDER BY dt.DetectionTimeID;

END
$$

DELIMITER ;
