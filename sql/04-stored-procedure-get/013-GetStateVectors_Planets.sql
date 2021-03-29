DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetStateVectors_Planets(
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get the state vectors for the massive bodies used in the Planets integration in a date range."

BEGIN 

-- Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

SELECT
	sv.TimeID,
	sv.BodyID,
	sv.mjd,
	sv.qx,
	sv.qy,
	sv.qz,
	sv.vx,
	sv.vy,
	sv.vz
FROM
	KS.StateVectors_Planets AS sv
WHERE 
	sv.TimeID BETWEEN @TimeID_0 AND @TimeID_1
ORDER BY sv.TimeID, sv.BodyID;

END
$$

DELIMITER ;
