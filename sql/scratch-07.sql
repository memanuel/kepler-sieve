DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidVectors(
    IN ast0 INT,
    IN ast1 INT
)
COMMENT "Get the state vectors for a block of asteroids."

BEGIN 

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
	av.AsteroidID BETWEEN ast0 AND (ast1-1)
ORDER BY av.TimeID;

END
$$

DELIMITER ;
