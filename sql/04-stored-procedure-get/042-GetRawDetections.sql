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
	dt.MJD AS mjd,
	rd.RA AS ra,
	rd.`DEC` AS `dec`,
	rd.Mag AS mag
FROM
	KS.RawDetection AS rd
	INNER JOIN KS.DetectionTime AS dt ON dt.DetectionTimeID = rd.DetectionTimeID
WHERE
	mjd0 <= dt.MJD AND dt.MJD < mjd1;

END
$$

DELIMITER ;
