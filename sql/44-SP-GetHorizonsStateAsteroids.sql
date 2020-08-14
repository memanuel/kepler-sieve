DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.GetHorizonsStateAsteroids(
    IN n_ast INT,
    IN epoch INT
)
COMMENT "Get the state vectors (position and velocity) of the first n_ast numbered asteroids at a given epoch (integer MJD only)"

BEGIN 

# Compute MinuteID from epoch
SET @MinuteID = epoch * 24 * 60;

# Get the state vector from JPL.HorizonsVectors and suitable joins
# Use a left join to JPL.MassiveBody, because not all bodies have known masses, just the heavy bodies.
SELECT
	b.BodyID,
	b.BodyName,
	COALESCE(mb.M, 0.0) AS m,
	hv.qx,
	hv.qy,
	hv.qz,
	hv.vx,
	hv.vy,
	hv.vz
FROM
	JPL.SmallBody AS sb
	CROSS JOIN JPL.HorizonsTime AS ht
	INNER JOIN KS.Body AS b ON b.BodyID = sb.BodyID
	INNER JOIN JPL.HorizonsBody AS hb ON hb.BodyID = b.BodyID
	INNER JOIN JPL.HorizonsVectors AS hv ON
		hv.HorizonsBodyID = hb.HorizonsBodyID AND
		hv.HorizonsTimeID = ht.HorizonsTimeID
	LEFT JOIN JPL.MassiveBody AS mb ON mb.HorizonsBodyID = hb.HorizonsBodyID
WHERE
	sb.SmallBodyID <= n_ast AND
	ht.MinuteID = @MinuteID
ORDER BY sb.SmallBodyID;

END
$$

DELIMITER ;
