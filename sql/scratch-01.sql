SET @BodyCollectionCD = 'P';
SET @mjd0 = 58000;
SET @mjd1 = 60000;

SET @TimeID_0 = @mjd0 * 24 * 60;
SET @TimeID_1 = @mjd1 * 24 * 60;

SELECT
	id.TimeID,
	id.BodyID,
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
	bc.BodyCollectionCD = @BodyCollectionCD
ORDER BY id.TimeID, id.BodyID;