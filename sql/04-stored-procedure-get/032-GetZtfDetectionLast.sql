DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE ZTF.GetDetectionLast()
COMMENT "Get the last DetectionID in the ZTF.Detection table."

BEGIN 

SELECT
	COALESCE(max(det.DetectionID),0) AS LastDetectionID
FROM
	ZTF.Detection AS det;

END
$$

DELIMITER ;
