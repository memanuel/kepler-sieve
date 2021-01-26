DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetIntegrationDiffByDate(
	IN BodyCollectionCD VARCHAR(4),
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get the difference between position and velocities between a Rebound integration and Horizons.  Take difference only on planets and group by date."

BEGIN 

# Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

# Query IntegrationDiff table to get difference by TimeID and BodyID for the selected collection
SELECT
	id.TimeID,
	AVG(id.MJD) AS MJD,
	AVG(id.dq / bv.sd_q) AS dq_rel,
	AVG(id.dv / bv.sd_v) AS dv_rel	
FROM
	KS.BodyCollection AS bc
 	INNER JOIN KS.IntegrationDiff AS id ON
 		id.BodyCollectionID = bc.BodyCollectionID AND
 		id.TimeID BETWEEN @TimeID_0 AND @TimeID_1
 	INNER JOIN JPL.BodyVariance AS bv ON bv.BodyID = id.BodyID
 	# Join to planets body collection so DE435 integration only compared on the planets
 	INNER JOIN KS.BodyCollection AS bcp ON bcp.BodyCollectionCD = 'P' 
 	INNER JOIN KS.BodyCollectionEntry AS bce ON
 		bce.BodyCollectionID = bcp.BodyCollectionID AND
 		bce.BodyID = id.BodyID
WHERE
	bc.BodyCollectionCD = BodyCollectionCD
GROUP BY id.TimeID
ORDER BY id.TimeID;

END
$$

DELIMITER ;
