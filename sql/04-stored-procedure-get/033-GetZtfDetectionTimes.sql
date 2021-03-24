DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE ZTF.GetDetectionTimes()
COMMENT "Get all the distinct detection times ZTF.Detection table."

BEGIN 

SELECT
	dt.DetectionTimeID,
	dt.MJD
FROM
	ZTF.DetectionTime AS dt
ORDER BY dt.DetectionTimeID;

END
$$

DELIMITER ;
