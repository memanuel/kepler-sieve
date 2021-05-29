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

# tau = radians in one circle (2pi)
SET @tau = 2.0 * PI();

SELECT
	el.TimeID,
	el.BodyID,
	el.mjd,
	el.a,
	el.e,
	el.inc,
	el.Omega_node AS Omega,
	el.omega_peri AS omega,
	-- Adjust f and M for the winding number
	el.f + el.WindingNumber  * @tau AS f,
	el.M + el.WindingNumber  * @tau AS M
FROM
	KS.OrbitalElements_Planets AS el
WHERE 
	el.TimeID BETWEEN @TimeID_0 AND @TimeID_1
	AND (el.TimeID % stride) = 0
ORDER BY el.TimeID, el.BodyID;

END
$$

DELIMITER ;
