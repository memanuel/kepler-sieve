DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.GetHorizonsStateCollection(
    IN BodyCollectionName VARCHAR(32),
    IN epoch INT
)
COMMENT "Get the state vectors (position and velocity) of all bodies in the named collection at a given epoch (integer MJD only)"

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
	KS.BodyCollection AS bc
	CROSS JOIN JPL.HorizonsTime AS ht
	INNER JOIN KS.BodyCollectionEntry AS bce ON	bce.BodyCollectionID = bc.BodyCollectionID
	INNER JOIN KS.Body AS b ON b.BodyID = bce.BodyID
	INNER JOIN JPL.HorizonsBody AS hb ON hb.BodyID = b.BodyID
	INNER JOIN JPL.HorizonsVectors AS hv ON
		hv.HorizonsBodyID = hb.HorizonsBodyID AND
		hv.HorizonsTimeID = ht.HorizonsTimeID
	LEFT JOIN JPL.MassiveBody AS mb ON mb.HorizonsBodyID = hb.HorizonsBodyID
WHERE
	bc.BodyCollectionName = BodyCollectionName AND
	ht.MinuteID = @MinuteID
ORDER BY bce.BodyNumber;

END
$$

DELIMITER ;
