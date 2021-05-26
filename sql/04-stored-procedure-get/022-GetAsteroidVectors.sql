DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidVectors(
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

SELECT
	av.TimeID,
    av.AsteroidID,
	av.mjd,
	av.qx,
	av.qy,
	av.qz,
	av.vx,
	av.vy,
	av.vz
FROM
	KS.AsteroidVectors AS av
WHERE
	av.TimeID BETWEEN @TimeID_0 AND @TimeID_1 AND
	av.AsteroidID BETWEEN n0 AND (n1-1)
ORDER BY av.TimeID, av.AsteroidID;

END
$$

DELIMITER ;
