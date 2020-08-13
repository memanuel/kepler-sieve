DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroids()
COMMENT "Get the state body IDs and names of asteroids"

BEGIN 

SELECT
	b.BodyID,
	b.BodyName
FROM
	KS.Body AS b
	INNER JOIN KS.BodyType AS bt ON bt.BodyTypeID = b.BodyTypeID
WHERE bt.BodyTypeCD = 'A'
ORDER BY b.BodyID;

END
$$

DELIMITER ;
