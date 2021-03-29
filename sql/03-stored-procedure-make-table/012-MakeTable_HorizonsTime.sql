DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_HorizonsTime()
COMMENT "Populate the HorizonsTime table from HorizonsImport"
BEGIN 

-- Empty table if necessary
TRUNCATE TABLE JPL.HorizonsTime;	
	
-- Distinct times; sufficient to query on just the sun because all major bodies were run on the same schedule.
INSERT INTO JPL.HorizonsTime
(TimeID, mjd, CalendarDate, CalendarDateTime, delta_T)
SELECT 
	FLOOR((hi.JD - 2400000.5) * 24 * 60) as TimeID,
	hi.JD - 2400000.5 as mjd,
	cast(hi.CalendarDateTime as Date) as Date,
	hi.CalendarDateTime as DateTime,
	hi.delta_T
FROM JPL.HorizonsImport as hi
WHERE hi.BodyNumber = 10 and hi.BodyTypeCD = 'S'
GROUP BY hi.JD
ORDER BY hi.JD;

END
$$

DELIMITER ;
