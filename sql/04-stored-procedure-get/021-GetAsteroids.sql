DELIMITER $$

-- ********************************************************************************
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

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidIDs()
COMMENT "Get the AsteroidID field only for all numbered asteroids"

BEGIN 

SELECT
	ast.AsteroidID
FROM
	KS.Asteroid AS ast
ORDER BY ast.AsteroidID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidNumberRange()
COMMENT "Get ranges of asteroid IDs"

BEGIN 

SELECT
	ast.IsnumberedAsteroid,
	min(ast.AsteroidID) AS AsteroidID_0,
	max(ast.AsteroidID) AS AsteroidID_1
FROM
	KS.Asteroid AS ast
GROUP BY ast.IsnumberedAsteroid
ORDER BY ast.IsNumberedAsteroid DESC;

END
$$

DELIMITER ;
