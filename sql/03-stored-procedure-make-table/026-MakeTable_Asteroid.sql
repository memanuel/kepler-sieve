DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_Asteroid()
COMMENT "Populate the Asteroid table from JPL.AsteroidElement joined to KS.Body"
BEGIN 

INSERT INTO KS.Asteroid
(AsteroidID, AsteroidNumber, AsteroidName, BodyID, IsNumberedAsteroid, H, G)
SELECT
	ast.AsteroidNumber AS AsteroidID,
	ast.AsteroidNumber,
	ast.AsteroidName,
	b.BodyID,
	ast.IsNumberedAsteroid,
	ast.H,
	ast.G
FROM 
	JPL.AsteroidElement AS ast
	INNER JOIN KS.Body AS b ON b.BodyID = ast.AsteroidNumber+1000000
;

END
$$

DELIMITER ;
