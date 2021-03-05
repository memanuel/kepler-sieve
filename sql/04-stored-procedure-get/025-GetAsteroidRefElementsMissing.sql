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

-- Select only orbital elements at this time that have not been integrated
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
	-- Only take reference elements for asteroids that haven't already been integrated into KS.AsteroidElements
	LEFT JOIN KS.AsteroidElements AS ae ON
		ae.TimeID = @TimeID AND
		ae.AsteroidID = ast.AsteroidID
	LEFT JOIN KS.AsteroidVectors AS av ON
		av.TimeID = @TimeID AND
		av.AsteroidID = ast.AsteroidID
WHERE
    (n0 <= ast.AsteroidID AND ast.AsteroidID < n1) 
    AND
    (elt.TimeID = @TimeID)
    AND
	(ae.TimeID IS NULL OR av.TimeID IS NULL)    
ORDER BY elt.AsteroidID;

END
$$

DELIMITER ;
