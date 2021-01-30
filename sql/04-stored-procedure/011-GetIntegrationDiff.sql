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

# Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

# Query IntegrationDiff table to get difference by TimeID and BodyID for the selected collection
SELECT
	id.TimeID,
	id.BodyID,
	id.MJD,
	id.dq,
	id.dv,
	(id.dq / bv.sd_q) AS dq_rel,
	(id.dv / bv.sd_v) AS dv_rel	
FROM
	KS.BodyCollection AS bc
 	INNER JOIN KS.IntegrationDiff AS id ON
 		id.BodyCollectionID = bc.BodyCollectionID AND
 		id.TimeID BETWEEN @TimeID_0 AND @TimeID_1
 	INNER JOIN JPL.BodyVariance AS bv ON bv.BodyID = id.BodyID
WHERE
	bc.BodyCollectionName = BodyCollectionName
ORDER BY id.TimeID, id.BodyID;

END
$$

DELIMITER ;
