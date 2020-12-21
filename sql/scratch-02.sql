INSERT INTO KS.TestTime
(IntegrationTimeID, MJD, JD, CalendarDate, CalendarDateTime, delta_T)
SELECT 
	ht.IntegrationTimeID + minutes._ AS IntegrationTimeID,
	ht.MJD + minutes.mjd_offset AS MJD,	
	ht.JD  + minutes.mjd_offset AS JD,
	ht.CalendarDate,
	DATE_ADD(ht.CalendarDateTime, INTERVAL minutes._ MINUTE) AS CalendarDateTime,
	(minutes.wt0 * ht.delta_T) + (minutes.wt1 * ht1.delta_T) AS delta_T
FROM 
	JPL.HorizonsTime as ht
	-- Time one day forward
	INNER JOIN JPL.HorizonsTime AS ht1 ON
		ht1.IntegrationTimeID = ht.IntegrationTimeID + 24*60
	INNER JOIN KS.Minutes AS minutes ON
		minutes._ = 0
LIMIT 10000;
