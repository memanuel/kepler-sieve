DELIMITER $$

-- ********************************************************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetMassiveBodies()
COMMENT "Get a list of all massive bodies in DE-435 colleciton."

BEGIN 

SELECT
	mb.BodyID,
    mb.M,
    mb.GM
FROM
	KS.MassiveBody as mb
	INNER JOIN KS.Body AS b ON b.BodyID = mb.BodyID
ORDER BY mb.BodyID;

END
$$

-- ********************************************************************************************************************
DELIMITER ;
