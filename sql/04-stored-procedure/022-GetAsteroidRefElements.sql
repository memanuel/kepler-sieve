DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidRefElements(
    IN epoch INT
)
COMMENT "Get all available reference orbital elements for asteroids on the given epoch."

BEGIN 

-- Compute TimeID from epoch
SET @TimeID = epoch * 24 * 60;

-- Select all available orbital elements at this time, sorted by AsteroidID.
SELECT
	elt.AsteroidID,
	elt.TimeID,
	elt.epoch,
	elt.a,
	elt.e,
	elt.inc,
	elt.Omega_node,
	elt.omega_peri,
	elt.f,
	elt.M
FROM
	KS.AsteroidElement_Ref AS elt
WHERE
	elt.TimeID = @TimeID
ORDER BY elt.AsteroidID;

END
$$

DELIMITER ;
