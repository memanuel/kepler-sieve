DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidIntegrationDiffByDate()
COMMENT "Get the difference between the MSE asteroid integration and Horizons state vectors."

BEGIN 

-- t1 is pairs of position vectors between MSE and JPL
WITH t1 AS (
SELECT
	av.AsteroidID,
	ast.AsteroidNumber,
	av.MJD,
	(av.qx - hv.qx) AS dx,
	(av.qy - hv.qy) AS dy,
	(av.qz - hv.qz) AS dz
FROM
    -- Start with MSE state vectors
	KS.AsteroidVectors AS av
    -- Join to foreign keys
	INNER JOIN KS.DailyTime AS dt ON dt.TimeID = av.TimeID
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidID = av.AsteroidID
	INNER JOIN KS.Body AS b ON b.BodyID = ast.BodyID
    -- Join to the JPL HorizonsVectors table by way of HorizonsBodyID
	INNER JOIN JPL.HorizonsBody AS hb ON hb.BodyID = b.BodyID
	INNER JOIN JPL.HorizonsVectors AS hv ON
		hv.TimeID = dt.TimeID AND 
		hv.HorizonsBodyID = hb.HorizonsBodyID
), 
-- t2 computes the distance between the MSE and JPL position of each asteroid time snap
-- pick up the semimajor axis a as a proxy for relative error denominator
t2 AS (
SELECT
	t1.AsteroidID,
	t1.MJD,
	sqrt(power(t1.dx,2) + power(t1.dy,2) + power(t1.dz,2)) AS dq,
	elt.a
FROM
	t1 
	INNER JOIN JPL.AsteroidElement AS elt ON elt.AsteroidNumber = t1.AsteroidNumber
)
-- Calculate the average absolute and relative difference by date
SELECT
	t2.MJD,
	avg(t2.dq) AS dq,
	avg(t2.dq / a) AS dq_rel,
	count(t2.AsteroidID) AS match_count
FROM
	t2
GROUP BY t2.MJD
ORDER BY t2.MJD;

END
$$

DELIMITER ;
