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
	dt0.mjd AS epoch,
	COUNT(elt_jpl.AsteroidNumber) AS AsteroidCount
FROM
	-- Start with JPL reference elements
	JPL.AsteroidElement AS elt_jpl
	-- The matching asteroid and body to this quote from JPL
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt_jpl.AsteroidNumber
	INNER JOIN KS.Body AS b ON b.BodyID = ast.BodyID
	-- The epoch when the elements are quoted
	INNER JOIN KS.DailyTime AS dt0 ON dt0.mjd = elt_jpl.epoch
	-- The epoch as of which we want the results
	INNER JOIN KS.DailyTime AS dt1 ON dt1.TimeID = @TimeID
	-- Try to join AsteroidElement_Ref on the date we want, dt1
	LEFT JOIN KS.AsteroidElement_Ref AS elt_ks ON
		elt_ks.AsteroidID = ast.AsteroidID AND
		elt_ks.TimeID = dt1.TimeID
WHERE
	-- Only take rows that are not already in KS.AsteroidElement_Ref
	elt_ks.AsteroidID IS NULL AND
	-- Only synchronize elements quoted PRIOR to the desired epoch
	dt0.TimeID < dt1.TimeID
GROUP BY dt0.TimeID;

END
$$

DELIMITER ;
