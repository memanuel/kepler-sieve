SELECT
	sv.TimeID,
	sv.BodyID,
	sv.MJD,
	sv.qx,
	sv.qy,
	sv.qz
FROM
	KS.StateVectors_Planets AS sv
WHERE
	sv.TimeID BETWEEN (59000*24*60) AND (59001*24*60)
LIMIT 100;