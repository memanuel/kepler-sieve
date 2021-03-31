DELIMITER $$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_DetectionTimeSlice(
	IN sz INT
)
COMMENT "Populate the KS.DetectionTimeSlice table from ZTF.DetectionTime"
BEGIN 

-- Set window size in minutes; MUST DIVIDE number of minutes in a day = 24*60
-- SET @sz = CAST(30 AS INT);
-- Number of minutes in one day
SET @mpd = CAST(24*60 AS INT);
-- Number of windows per day
SET @wpd = CAST( (@mpd DIV sz) AS INT);
-- Size of the window in days
SET @szd = CAST(sz / @mpd AS DOUBLE);

-- Populate KS.DetectionTimeSlice using the input TimeSlice
INSERT IGNORE INTO KS.DetectionTimeSlice
(DetectionTimeSliceID, IntegrationTimeID, mjd, mjd0, mjd1)
WITH t1 AS (
SELECT
	-- Integer ID is the floor number of minutes to multiples of the window size
	FLOOR(dt.mjd*@wpd)*sz AS DetectionTimeSliceID,
	-- Nearest IntegrationTime to the midpoint of the window; IntegrationTime has a width of 5 minutes
	ROUND(dt.mjd*24*12)*5 AS IntegrationTimeID,
	-- Midpoint of the interval as an MJD
	(FLOOR(dt.mjd*@wpd)+0.5) * @szd AS mjd,
	-- Start of the interval as an MJD
	(FLOOR(dt.mjd*@wpd)+0) * @szd AS mjd0,
	-- End of the interval as an MJD
	(FLOOR(dt.mjd*@wpd)+1) * @szd AS mjd1
FROM
	ZTF.DetectionTime AS dt
)
SELECT
	t1.DetectionTimeSliceID,
	t1.IntegrationTimeID,
	t1.mjd,
	t1.mjd0,
	t1.mjd1
FROM
	t1
GROUP BY t1.DetectionTimeSliceID;

END $$

