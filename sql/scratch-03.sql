DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.GetAsteroidRefElements(
    IN epoch INT
)
COMMENT "Get all available reference orbital elements as directly quoted by JPL."

BEGIN 

-- Compute TimeID from epoch
SET @TimeID = epoch * 24 * 60;

-- Select all available orbital elements at this time, sorted by AsteroidID.
SELECT
	ast.AsteroidID,
	it.TimeID,
	elt.epoch,
	elt.a,
	elt.e,
	elt.inc,
	elt.Omega_node AS Omega,
	elt.omega_peri AS omega,
	elt.M AS M
FROM
	JPL.AsteroidElement AS elt
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt.AsteroidNumber
	INNER JOIN KS.IntegrationTime AS it ON it.MJD = elt.epoch
WHERE
	it.TimeID = @TimeID
ORDER BY ast.AsteroidID;

END
$$

DELIMITER ;
