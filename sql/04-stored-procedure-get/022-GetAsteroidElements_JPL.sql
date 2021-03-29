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
	dt0.TimeID,
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
	-- The matching asteroid and body to this quote from JPL
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt.AsteroidNumber
	INNER JOIN KS.Body AS b ON b.BodyID = ast.BodyID
	-- The epoch when the elements are quoted
	INNER JOIN KS.DailyTime AS dt0 ON dt0.mjd = elt.epoch
	-- The epoch as of which we want the results
	INNER JOIN KS.DailyTime AS dt1 ON dt1.TimeID = @TimeID
	-- Test if we already have a match on the AsteroidElement_Ref table
	LEFT JOIN KS.AsteroidElement_Ref AS relt ON
		relt.AsteroidID = ast.AsteroidID AND 
		relt.TimeID = dt1.TimeID
WHERE
	-- Only match the selected date
	elt.epoch = epoch
	-- Only asteroids that don't already have reference elements
	AND relt.AsteroidID IS NULL 
	-- Only return matches where we integrate forward from the quoted epoch to the desired epoch
	-- AND dt0.TimeID <= dt1.TimeID
ORDER BY ast.AsteroidID;

END
$$

DELIMITER ;
