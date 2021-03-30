DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_DetectionTimePair()
COMMENT "Populate the KS.DetectionTimePair table from KS.DetectionTime"
BEGIN 

-- Populate KS.MakeTable_DetectionTimePair with pairs at most one minute apart
INSERT IGNORE INTO KS.DetectionTimePair
(DetectionTimeID_1, DetectionTimeID_2, DataSourceID, mjd1, mjd2)
SELECT
	dt1.DetectionTimeID,
	dt2.DetectionTimeID,
	dt1.DataSourceID,
	dt1.mjd AS mjd1,
	dt2.mjd AS mjd2
FROM
	KS.DetectionTime AS dt1
	INNER JOIN KS.DetectionTime AS dt2 ON 
		dt1.mjd < dt2.mjd AND dt2.mjd < dt1.mjd + 1.0 / (24*60) AND
        -- Only a pair of detection times from the same source
		dt2.DataSourceID = dt1.DataSourceID;

END
$$

DELIMITER ;
