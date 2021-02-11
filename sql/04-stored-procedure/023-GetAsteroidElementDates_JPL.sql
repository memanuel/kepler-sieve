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
-- WHERE NOT EXISTS (
-- 	SELECT elt2.AsteroidNumber FROM JPL.AsteroidElement AS elt2 
-- 	WHERE elt2.AsteroidNumber = elt.AsteroidNumber AND elt2.epoch = epoch
-- 	)
GROUP BY elt.epoch
ORDER BY elt.epoch;

END
$$

DELIMITER ;
