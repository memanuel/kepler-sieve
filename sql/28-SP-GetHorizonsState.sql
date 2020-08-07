DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.GetHorizonsState(
    IN HorizonsBodyName VARCHAR(32),
    IN epoch INT
)
COMMENT "Get the state vector (position and velocity) of one object at a given epoch (integer MJD only)"

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
	JPL.HorizonsVectors AS hv
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyID = hv.HorizonsBodyID
	INNER JOIN JPL.HorizonsTime AS ht ON ht.HorizonsTimeID = hv.HorizonsTimeID
	INNER JOIN KS.Body AS b ON b.BodyID = hb.BodyID
	LEFT JOIN JPL.MassiveBody AS mb ON mb.HorizonsBodyID = hv.HorizonsBodyID
WHERE
	hb.HorizonsBodyName = HorizonsBodyName AND ht.MinuteID = @MinuteID;

END
$$

DELIMITER ;
