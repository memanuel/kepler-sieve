explain
SELECT
 	it.TimeID,
 	bw.BodyID,
 	it.MJD,
 	bw.Weight
FROM 
	KS.Body AS b_emb	
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionCD = 'PS3'
	INNER JOIN KS.BarycenterWeight AS bw ON
		bw.BodyCollectionID = bc.BodyCollectionID
	-- INNER JOIN KS.Body AS b ON b.BodyID = bw.BodyID
	INNER JOIN KS.IntegrationTime AS it ON (it.TimeID BETWEEN 59000*24*60 AND 59001*24*60)
	INNER JOIN KS.StateVectors_Planets AS sv ON
		sv.TimeID = it.TimeID AND
		sv.BodyID = bw.BodyID
	
WHERE
	b_emb.BodyName = 'Earth-Moon Barycenter'
GROUP BY it.TimeID
;