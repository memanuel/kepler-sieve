DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_HiResTime()
COMMENT "Populate the IntegrationTime table from JPL.HorizonsTime"
BEGIN 

-- Populate IntegrationTime from HorizonsTime by joining pairs of dates to a minutes offset
INSERT INTO KS.HiResTime
(TimeID, mjd)
SELECT 
	it1.TimeID + m._ AS TimeID,
	it1.mjd + m.mjd_offset AS mjd
FROM 
	KS.IntegrationTime AS it1
	INNER JOIN KS.IntegrationTime AS it2 ON it2.TimeID = it1.TimeID+5
	CROSS JOIN KS.Minutes AS m ON m._ < 5;

-- Handle the last record in HorizonsTime
INSERT INTO KS.IntegrationTime
(TimeID, mjd)
SELECT 
	it.TimeID,
	it.mjd
FROM 
	KS.IntegrationTime as it
WHERE
	it.TimeID = 77600*24*60;

END
$$

DELIMITER ;
