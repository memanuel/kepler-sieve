DELIMITER $$

-- ****************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetElements_Earth(
    IN mjd0 INT,
    IN mjd1 INT,
    IN stride INT
)
COMMENT "Get the orbital elements for Earth used in the Planets integration in a date range. Stride can be a multiple of 5 minutes."

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
	INNER JOIN KS.Body AS b ON b.BodyID = el.BodyID	
WHERE
	b.BodyName = 'Earth'
	AND el.TimeID BETWEEN @TimeID_0 AND @TimeID_1
	AND (el.TimeID % stride) = 0
ORDER BY el.TimeID;

END
$$

-- ****************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetVectorsElements_Earth(
    IN mjd0 INT,
    IN mjd1 INT,
    IN stride INT
)
COMMENT "Get the state vectors and orbital elements for Earth. Stride can be a multiple of 5 minutes."

BEGIN 

-- Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

SELECT
	-- Key fields
	sv.TimeID,
	sv.BodyID,
	sv.mjd,
	-- State vectors
	sv.qx,
	sv.qy,
	sv.qz,
	sv.vx,
	sv.vy,
	sv.vz,
	-- Orbital elements
	el.a,
	el.e,
	el.inc,
	el.Omega_node AS Omega,
	el.omega_peri AS omega,
	el.f,
	el.M
FROM
	KS.StateVectors_Planets AS sv
	INNER JOIN KS.OrbitalElements_Planets AS el ON
		el.BodyID = sv.BodyID AND
		el.TimeID = sv.TimeID
	INNER JOIN KS.Body AS b ON b.BodyID = sv.BodyID	
WHERE
	b.BodyName = 'Earth'
	AND el.TimeID BETWEEN @TimeID_0 AND @TimeID_1
	AND (el.TimeID % stride) = 0
ORDER BY el.TimeID;

END
$$

DELIMITER ;
