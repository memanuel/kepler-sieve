DELIMITER $$

# *****************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetStateVectors_Earth(
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get the state vectors for Earth in a date range. High resolution (one minute)."

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
	KS.StateVectors_Earth AS sv
WHERE 
	sv.TimeID BETWEEN @TimeID_0 AND @TimeID_1
ORDER BY sv.TimeID;

END
$$

# *****************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetStateVectors_EarthAll()
COMMENT "Get all available state vectors for Earth at one minute intervals."

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
	KS.StateVectors_Earth AS sv
-- WHERE (sv.TimeID % 5) = 0
ORDER BY sv.TimeID;

END

$$

DELIMITER ;
