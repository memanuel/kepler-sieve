DELIMITER $$

# *****************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetStateVectors_Sun(
    IN mjd0 INT,
    IN mjd1 INT,
    IN stride INT
)
COMMENT "Get the state vectors for Sun in a date range. Stride can be a multiple of 5 minutes."

BEGIN 

-- Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

SELECT
	sv.TimeID,
	sv.mjd,
	sv.qx,
	sv.qy,
	sv.qz,
	sv.vx,
	sv.vy,
	sv.vz
FROM
	KS.StateVectors_Planets AS sv
	INNER JOIN KS.Body AS b ON b.BodyID = sv.BodyID
WHERE
	b.BodyName = 'Sun' 
	AND sv.TimeID BETWEEN @TimeID_0 AND @TimeID_1 
	AND (sv.TimeID % stride) = 0	
ORDER BY sv.TimeID;

END
$$

/*
# *****************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetStateVectors_Sun_Minute(
    IN mjd0 INT,
    IN mjd1 INT,
    IN stride INT
)
COMMENT "Get the state vectors for Sun in a date range. Stride can be a multiple of 1 minute."

BEGIN 

-- Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

SELECT
	sv.TimeID,
	sv.mjd,
	sv.qx,
	sv.qy,
	sv.qz,
	sv.vx,
	sv.vy,
	sv.vz
FROM
	KS.StateVectors_Sun AS sv
WHERE 
	sv.TimeID BETWEEN @TimeID_0 AND @TimeID_1
	AND (sv.TimeID % stride) = 0
ORDER BY sv.TimeID;

END
$$
*/

DELIMITER ;
