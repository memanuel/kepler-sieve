DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidRefElements(
    IN epoch INT,
    IN n0 INT,
    IN n1 INT
)
COMMENT "Get all available reference orbital elements for asteroids on the given epoch."

BEGIN 

-- Compute TimeID from epoch
SET @TimeID = epoch * 24 * 60;

-- Select all available orbital elements at this time, sorted by AsteroidID.
SELECT
	elt.AsteroidID,
	elt.TimeID,
	ast.AsteroidName,
	elt.epoch,
	elt.a,
	elt.e,
	elt.inc,
	elt.Omega_node AS Omega,
	elt.omega_peri AS omega,
	elt.M,
	b.BodyID,
	b.BodyName
FROM
	KS.AsteroidElement_Ref AS elt
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidID = elt.AsteroidID
	INNER JOIN KS.Body AS b ON b.BodyID = ast.BodyID
WHERE
    (n0 <= ast.AsteroidID AND ast.AsteroidID < n1) 
    AND
    (elt.TimeID = @TimeID)
ORDER BY elt.AsteroidID;

/*
SELECT
	ast.AsteroidID,	
	dt.TimeID,
	ast.AsteroidName,
	elt.epoch,
	elt.a,
	elt.e,
	elt.inc,
	elt.Omega_node AS Omega,
	elt.omega_peri AS omega,
	elt.M AS M,
	b.BodyID,
	b.BodyName
FROM
	JPL.AsteroidElement AS elt
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt.AsteroidNumber
	INNER JOIN KS.Body AS b ON b.BodyID = ast.BodyID
	INNER JOIN KS.DailyTime AS dt ON dt.MJD = elt.epoch
WHERE
	dt.TimeID = @TimeID AND
	n0 <= ast.AsteroidID AND ast.AsteroidID < n1
ORDER BY ast.AsteroidID;
*/

END
$$

DELIMITER ;
