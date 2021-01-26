CREATE OR REPLACE 
DEFINER = kepler
VIEW KS.IntegrationDiff_DE435
AS
SELECT
	ki.TimeID,
	ki.BodyID,
    ki.MJD,
	sqrt(power(ki.qx - hv.qx, 2) + power(ki.qy - hv.qy, 2) + power(ki.qz - hv.qz, 2)) AS dq,
	sqrt(power(ki.vx - hv.vx, 2) + power(ki.vy - hv.vy, 2) + power(ki.vz - hv.vz, 2)) AS dv
FROM
	# Start with daily snapshots of planets from Horizons
	JPL.HorizonsVectors AS hv
	# Get the BodyID by joining to HorizonsBody
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyID = hv.HorizonsBodyID
	# INNER JOIN KS.Body AS b ON b.BodyID = hb.BodyID
	# Look up corresponding record on integration table
	INNER JOIN KS.Integration_DE435 AS ki ON
		ki.TimeID = hv.TimeID AND
		ki.BodyID = hb.BodyID;
