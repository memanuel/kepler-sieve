DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_DetectionTimePair(
	IN width INT)
COMMENT "Populate the KS.DetectionTimePair table from KS.DetectionTime"
BEGIN 

-- Number of minutes in one day
SET @mpd = CAST(24*60 AS INT);
-- Convert width to a number of MJDs
SET @window_width = CAST(width / @mpd AS DOUBLE);
	
-- Populate KS.MakeTable_DetectionTimePair with pairs at most one minute apart
INSERT IGNORE INTO KS.DetectionTimePair
(DetectionTimeID_1, DetectionTimeID_2, DataSourceID, mjd1, mjd2, mjd, dt)
SELECT
	-- The two detection time IDs
	dt1.DetectionTimeID,
	dt2.DetectionTimeID,
	-- Shared data source
	dt1.DataSourceID,
	-- The two times
	dt1.mjd AS mjd1,
	dt2.mjd AS mjd2,
	-- Mean detection time
	(dt1.mjd + dt2.mjd)/2.0 AS mjd,
	-- Difference of detection times
	(dt2.mjd - dt1.mjd) AS dt
FROM
	KS.DetectionTime AS dt1
INNER JOIN KS.DetectionTime AS dt2 ON 
	dt1.mjd < dt2.mjd AND dt2.mjd < dt1.mjd + @window_width AND
    -- Only a pair of detection times from the same source
	dt2.DataSourceID = dt1.DataSourceID
ORDER BY dt1.DetectionTimeID, dt2.DetectionTimeID;

END $$

DELIMITER ;
