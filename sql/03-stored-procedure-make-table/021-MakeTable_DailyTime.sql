DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_DailyTime()
COMMENT "Populate the DailyTime table from JPL.HorizonsTime"
BEGIN 

-- Populate IntegrationTime from HorizonsTime by joining pairs of dates to a minutes offset
INSERT INTO KS.DailyTime
(TimeID, MJD, CalendarDate, CalendarDateTime, delta_T)
SELECT
	ht.TimeID,
	ht.MJD,
	ht.CalendarDate,
	DATE_ADD(ht.CalendarDateTime, INTERVAL 0 MINUTE) AS CalendarDateTime,
	ht.delta_T 	
FROM
	JPL.HorizonsTime AS ht;

END
$$

DELIMITER ;
