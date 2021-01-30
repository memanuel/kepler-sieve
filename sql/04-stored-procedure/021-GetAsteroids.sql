DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroids()
COMMENT "Get the body IDs and names of asteroids"

BEGIN 

SELECT
	ast.AsteroidID,
	ast.BodyID,
	ast.AsteroidName
FROM
	KS.Asteroid AS ast
ORDER BY ast.AsteroidID;

END
$$

DELIMITER ;
