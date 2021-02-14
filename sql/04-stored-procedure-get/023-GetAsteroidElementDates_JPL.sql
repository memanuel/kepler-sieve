DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.GetAsteroidRefElementDates(
    IN epoch INT
)
COMMENT "Get all available reference orbital elements for asteroids on the given epoch."

BEGIN 

-- Compute TimeID from epoch
SET @TimeID = epoch * 24 * 60;

	
-- Select all distinct dates with orbital elements available
SELECT
	dt0.TimeID,
	dt0.MJD AS epoch,
	COUNT(elt.AsteroidNumber) AS AsteroidCount
FROM
	JPL.AsteroidElement AS elt
	-- The matching asteroid and body to this quote from JPL
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt.AsteroidNumber
	INNER JOIN KS.Body AS b ON b.BodyID = ast.BodyID
	-- The epoch when the elements are quoted
	INNER JOIN KS.DailyTime AS dt0 ON dt0.MJD = elt.epoch
	-- The epoch as of which we want the results
	INNER JOIN KS.DailyTime AS dt1 ON dt1.TimeID = @TimeID
	-- Test if we already have a match on the AsteroidElement_Ref table
	LEFT JOIN KS.AsteroidElement_Ref AS relt ON
		relt.AsteroidID = ast.AsteroidID AND 
		relt.TimeID = dt1.TimeID
WHERE
	relt.AsteroidID IS NULL AND
	-- Only synchronize elements quoted PRIOR to the desired epoch
	dt0.TimeID < dt1.TimeID
GROUP BY elt.epoch
ORDER BY elt.epoch;

END
$$

DELIMITER ;
