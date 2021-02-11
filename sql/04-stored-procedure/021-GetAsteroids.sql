DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroids()
COMMENT "Get the body IDs and names of asteroids"

BEGIN 

SELECT
	ast.AsteroidID,
	ast.BodyID,
	ast.AsteroidName,
	b.BodyName
FROM
	KS.Asteroid AS ast
	INNER JOIN KS.Body AS b ON b.BodyID = ast.BodyID
ORDER BY ast.AsteroidID;

END
$$

DELIMITER ;
