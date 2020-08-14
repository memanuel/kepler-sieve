DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_IntegrationTime()
COMMENT "Populate the IntegrationTime table from JPL.HorizonsTime"
BEGIN 

INSERT INTO KS.IntegrationTime
(MinuteID, MJD, JD, CalendarDate, CalendarDateTime, delta_T)
SELECT 
	FLOOR((ht.JD - 2400000.5) * 24 * 60) as MinuteID,
	ht.JD - 2400000.5 as MJD,
	ht.JD,
	cast(ht.CalendarDateTime as Date) as Date,
	ht.CalendarDateTime as DateTime,
	ht.delta_T
FROM JPL.HorizonsTime as ht;

# Populate the the IntegrationTimeID field on JPL.HorizonsTime
UPDATE 
	JPL.HorizonsTime AS ht
	INNER JOIN KS.IntegrationTime AS it ON it.MinuteID = ht.MinuteID
SET
	ht.IntegrationTimeID = it.IntegrationTimeID;

END
$$

DELIMITER ;
