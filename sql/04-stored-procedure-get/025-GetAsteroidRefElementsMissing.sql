DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidRefElementsMissing(
    IN epoch INT,
    IN n0 INT,
    IN n1 INT
)
COMMENT "Get all reference orbital elements for asteroids on the given epoch that have not already been integrated."

BEGIN 

-- Compute TimeID from epoch
SET @TimeID = epoch * 24 * 60;

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
	-- Only take reference elements for asteroids that haven't already been integrated into KS.AsteroidElements
	LEFT JOIN KS.AsteroidElements AS ae ON
		ae.TimeID = dt.TimeID AND
		ae.AsteroidID = ast.AsteroidID
	LEFT JOIN KS.AsteroidVectors AS av ON
		av.TimeID = dt.TimeID AND
		av.AsteroidID = ast.AsteroidID
WHERE
	dt.TimeID = @TimeID AND
	n0 <= ast.AsteroidID AND ast.AsteroidID < n1 AND
	(ae.TimeID IS NULL OR av.TimeID IS NULL)

ORDER BY ast.AsteroidID;

END
$$

DELIMITER ;
