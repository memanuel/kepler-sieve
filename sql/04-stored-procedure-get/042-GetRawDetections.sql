DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetRawDetections(
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get raw detections in the specified date range."

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
	mjd0 <= dt.mjd AND dt.mjd < mjd1;

END
$$

DELIMITER ;
