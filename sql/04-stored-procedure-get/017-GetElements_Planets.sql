DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetElements_Planets(
    IN mjd0 INT,
    IN mjd1 INT,
    IN stride INT
)
COMMENT "Get the orbital elements for the massive bodies used in the Planets integration in a date range. Stride can be a multiple of 5 minutes."

BEGIN 

-- Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

SELECT
	el.TimeID,
	el.BodyID,
	el.mjd,
	el.a,
	el.e,
	el.inc,
	el.Omega_node AS Omega,
	el.omega_peri AS omega,
	el.f,
	el.M
FROM
	KS.OrbitalElements_Planets AS el
WHERE 
	el.TimeID BETWEEN @TimeID_0 AND @TimeID_1
	AND (el.TimeID % stride) = 0
ORDER BY el.TimeID, el.BodyID;

END
$$

DELIMITER ;
