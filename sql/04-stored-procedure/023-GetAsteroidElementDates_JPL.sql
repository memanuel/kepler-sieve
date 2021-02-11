DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.GetAsteroidRefElementDates(
    IN epoch INT
)
COMMENT "Get all available reference orbital elements for asteroids on the given epoch."

BEGIN 

-- Select all distinct dates with orbital elements available
SELECT
	CAST(elt.epoch AS INT) AS epoch,
	COUNT(elt.AsteroidNumber) AS AsteroidCount
FROM
	JPL.AsteroidElement AS elt
	INNER JOIN KS.IntegrationTime AS it ON
		it.MJD = elt.epoch
-- Only take records that don't have a match
WHERE NOT EXISTS (
	SELECT relt.AsteroidID 
	FROM 
		KS.AsteroidElement_Ref AS relt
		INNER JOIN KS.Asteroid AS ast ON ast.AsteroidID = relt.AsteroidID
	WHERE ast.AsteroidNumber = elt.AsteroidNumber AND relt.epoch = epoch
	)
GROUP BY elt.epoch
ORDER BY elt.epoch;

END
$$

DELIMITER ;
