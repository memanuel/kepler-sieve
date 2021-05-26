DELIMITER $$

# *****************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetStateVectors_Sun(
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get the state vectors for Sun in a date range. High resolution (one minute)."

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
ORDER BY sv.TimeID;

END
$$

# *****************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetStateVectors_SunAll()
COMMENT "Get all available state vectors for Sun at one minute intervals."

BEGIN 

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
ORDER BY sv.TimeID;

END

$$

DELIMITER ;
