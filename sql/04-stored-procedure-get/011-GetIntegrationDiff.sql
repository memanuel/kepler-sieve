DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetIntegrationDiff(
	IN BodyCollectionName VARCHAR(32),
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get the difference between position and velocities between a Rebound integration and Horizons."

BEGIN 

-- Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

WITH id AS(
SELECT
	dt.TimeID,
	dt.mjd,
	b.BodyID,
	SQRT(POW(sv.qx - hv.qx, 2) + POW(sv.qy - hv.qy, 2) + POW(sv.qz - hv.qz, 2)) AS dq,
	SQRT(POW(sv.vx - hv.vx, 2) + POW(sv.vy - hv.vy, 2) + POW(sv.vz - hv.vz, 2)) AS dv
FROM
	-- Start with the body collection we want to check (e.g. Planets)
	KS.BodyCollection AS bc
	-- All the bodies in this collection
	INNER JOIN KS.BodyCollectionEntry AS bce ON bce.BodyCollectionID = bc.BodyCollectionID
	INNER JOIN KS.Body AS b ON b.BodyID = bce.BodyID
	-- Get the HorizonsBodyID associated with this body
	INNER JOIN JPL.HorizonsBody AS hb ON hb.BodyID = b.BodyID
	-- All the daily dates in the input range
	INNER JOIN KS.DailyTime AS dt ON dt.TimeID BETWEEN @TimeID_0 AND @TimeID_1
	-- The position in the KS integration
	INNER JOIN KS.StateVectors_Planets AS sv ON
		sv.TimeID = dt.TimeID AND
		sv.BodyID = b.BodyID
	-- The position according to Horizons
	INNER JOIN JPL.HorizonsVectors AS hv ON
		hv.TimeID = dt.TimeID AND
		hv.HorizonsBodyID = hb.HorizonsBodyID
WHERE
	bc.BodyCollectionName = BodyCollectionName
)
SELECT
	id.TimeID,
	id.BodyID,
	id.mjd,
	id.dq,
	id.dv,
	(id.dq / bv.sd_q) AS dq_rel,
	(id.dv / bv.sd_v) AS dv_rel	
FROM
	id
 	INNER JOIN JPL.BodyVariance AS bv ON bv.BodyID = id.BodyID
ORDER BY id.TimeID, id.BodyID;

END
$$

DELIMITER ;
