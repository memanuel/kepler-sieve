SELECT 
	elt.TimeID,
	elt.BodyID,
	elt.MJD,
	elt.a,
	elt.e
FROM KS.OrbitalElements AS elt
WHERE
	elt.TimeID = 59000*24*60 AND
	elt.BodyID IN (301, 399)
;


UPDATE KS.OrbitalElements AS elt
SET elt.a=10
WHERE elt.TimeID = 59000*24*60 AND elt.BodyID=301;
