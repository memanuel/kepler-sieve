DELIMITER $$

-- ****************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidElements(
    IN n0 INT,
    IN n1 INT,
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get the orbital elements for a block of asteroids."

BEGIN 

-- Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;
	
# tau = radians in one circle (2pi)
SET @tau = 2.0 * PI();

SELECT
	ae.TimeID,
    ae.AsteroidID,
	ae.mjd,
	ae.a,
	ae.e,
	ae.inc,
	ae.Omega_node AS Omega,
	ae.omega_peri AS omega,
	-- Adjust f and M for the winding number
	ae.f + ae.WindingNumber*@tau AS f,
	ae.M + ae.WindingNumber*@tau AS M,
	ae.WindingNumber
FROM
	KS.AsteroidElements AS ae
WHERE 
	ae.TimeID BETWEEN @TimeID_0 AND @TimeID_1 AND
	ae.AsteroidID BETWEEN n0 AND (n1-1)
ORDER BY ae.AsteroidID, ae.TimeID;

END
$$

-- ****************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidData(
    IN n0 INT,
    IN n1 INT,
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get the state vectors for a block of asteroids."

BEGIN 

-- Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;
	
# tau = radians in one circle (2pi)
SET @tau = 2.0 * PI();

SELECT
	-- Key fields and time
	av.TimeID,
    av.AsteroidID,
	av.mjd,
	-- State vectors
	av.qx,
	av.qy,
	av.qz,
	av.vx,
	av.vy,
	av.vz,
	-- Orbital elements (heliocentric)
	ae.a,
	ae.e,
	ae.inc,
	ae.Omega_node AS Omega,
	ae.omega_peri AS omega,
	-- Adjust f and M for the winding number
	ae.f + ae.WindingNumber*@tau AS f,
	ae.M + ae.WindingNumber*@tau AS M,
	ae.WindingNumber
FROM
	KS.AsteroidVectors AS av 
	INNER JOIN KS.AsteroidElements AS ae ON
		ae.TimeID = av.TimeID AND
		ae.AsteroidID = av.AsteroidID
WHERE
	av.TimeID BETWEEN @TimeID_0 AND @TimeID_1 AND
	av.AsteroidID BETWEEN n0 AND (n1-1)
ORDER BY av.AsteroidID, av.TimeID;

END
$$

DELIMITER ;
