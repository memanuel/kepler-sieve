DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.GetAsteroidRefElements(
    IN epoch INT
)
COMMENT "Get reference orbital elements as directly quoted by JPL."

BEGIN 

-- Compute TimeID from epoch
SET @TimeID = epoch * 24 * 60;

-- Select all available orbital elements at this time, sorted by AsteroidID.
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
	LEFT JOIN KS.AsteroidElement_Ref AS relt ON
		relt.AsteroidID = ast.AsteroidID AND relt.TimeID = dt.TimeID
WHERE
	dt.TimeID = @TimeID AND relt.AsteroidID IS NULL
ORDER BY ast.AsteroidID;

END
$$

DELIMITER ;
