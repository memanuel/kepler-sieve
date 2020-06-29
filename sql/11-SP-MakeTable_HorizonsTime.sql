DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_HorizonsTime()
BEGIN 

-- Distinct times; sufficient to query on just one body because all bodies were run on the same schedule.
INSERT INTO JPL.HorizonsTime
(MinuteID, MJD, JD, CalendarDate, CalendarDateTime, delta_T)
SELECT 
	round((hi.JD - 2400000.5) * 24 * 60) as MinuteID,
	hi.JD - 2400000.5 as MJD,
	hi.JD,
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