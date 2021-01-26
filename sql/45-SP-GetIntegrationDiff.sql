DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetIntegrationDiff_Planets(
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get the difference between position and velocities between a Rebound integration and Horizons."

BEGIN 

# Compute TimeID from dates
SET @TimeID_0 = @mjd0 * 24 * 60;
SET @TimeID_1 = @mjd1 * 24 * 60;

# Get the state vector from JPL.HorizonsVectors and suitable joins
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
	hb.HorizonsBodyName = HorizonsBodyName AND ht.TimeID = @TimeID;

END
$$

DELIMITER ;
