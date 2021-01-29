SELECT
 	it.TimeID,
 	b_emb.BodyID,
 	it.MJD,
	SUM(bw.Weight * sv.qx) AS qx,
	SUM(bw.Weight * sv.qy) AS qy,
	SUM(bw.Weight * sv.qz) AS qz,
	SUM(bw.Weight * sv.vx) AS vx,
	SUM(bw.Weight * sv.vy) AS vy,
	SUM(bw.Weight * sv.vz) AS vz
FROM 
	KS.Body AS b_emb	
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionCD = 'PS3'
	INNER JOIN KS.BarycenterWeight AS bw ON
		bw.BodyCollectionID = bc.BodyCollectionID
	INNER JOIN KS.Body AS b ON b.BodyID = bw.BodyID
	INNER JOIN KS.IntegrationTime AS it
	INNER JOIN KS.StateVectors_Planets AS sv ON
		sv.TimeID = it.TimeID AND
		sv.BodyID = b.BodyID
WHERE
	b_emb.BodyName = 'Earth-Moon Barycenter'
GROUP BY it.TimeID
LIMIT 100;