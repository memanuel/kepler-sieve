DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_IntegrationTime()
COMMENT "Populate the IntegrationTime table from JPL.HorizonsTime"
BEGIN 

-- Populate IntegrationTime from HorizonsTime by joining pairs of dates to a minutes offset
	INSERT INTO KS.IntegrationTime
(TimeID, MJD, CalendarDate, CalendarDateTime, delta_T)
SELECT 
	ht.TimeID + minutes._ AS TimeID,
	ht.MJD + minutes.mjd_offset AS MJD,
	ht.CalendarDate,
	DATE_ADD(ht.CalendarDateTime, INTERVAL minutes._ MINUTE) AS CalendarDateTime,
	(minutes.wt0 * ht.delta_T) + (minutes.wt1 * ht1.delta_T) AS delta_T
FROM 
	-- The HorizonsTime record at the start of the day
	JPL.HorizonsTime as ht
	-- The HorizonsTime record one day in forward; used to interpolate delta_T
	INNER JOIN JPL.HorizonsTime AS ht1 ON
		ht1.TimeID = ht.TimeID + (24*60)
	-- All the minutes records
	CROSS JOIN KS.Minutes AS minutes;

-- Handle the last record in HorizonsTime
INSERT INTO KS.IntegrationTime
(TimeID, MJD, CalendarDate, CalendarDateTime, delta_T)
SELECT 
	ht.TimeID,
	ht.MJD,
	ht.CalendarDate,
	ht.CalendarDateTime,
	ht.delta_T
FROM 
	-- The HorizonsTime record for the very last day in HorizonsTime
	JPL.HorizonsTime as ht
WHERE
	ht.TimeID = 77600*24*60;

END
$$

DELIMITER ;
