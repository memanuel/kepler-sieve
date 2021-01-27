DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_IntegrationDiff()
COMMENT "Populate the IntegrationDiff table from views IntegrationDiff_Planets and IntegrationDiff_DE435"

BEGIN 

# Start by emptying out the table
TRUNCATE TABLE KS.IntegrationDiff;
	
# Insert batch of records comparing Integration_Planets to HorizonsVectors
INSERT IGNORE INTO KS.IntegrationDiff
(BodyCollectionID, TimeID, BodyID, MJD, dq, dv)
SELECT
	bc.BodyCollectionID,
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
	# Look up corresponding record on integration table
	INNER JOIN KS.Integration_Planets AS ki ON
		ki.TimeID = hv.TimeID AND
		ki.BodyID = hb.BodyID
	# The body collection that was integrated
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'Planets'
	# Only take bodies in the planets integration!
	INNER JOIN KS.BodyCollection AS bcp ON bcp.BodyCollectionName = 'Planets'
	INNER JOIN KS.BodyCollectionEntry AS bce ON
		bce.BodyCollectionID = bcp.BodyCollectionID AND
		bce.BodyID = hb.BodyID;

# Insert batch of records comparing Integration_DE435 to HorizonsVectors
INSERT IGNORE INTO KS.IntegrationDiff
(BodyCollectionID, TimeID, BodyID, MJD, dq, dv)
SELECT
	bc.BodyCollectionID,
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
	# Look up corresponding record on integration table
	INNER JOIN KS.Integration_DE435 AS ki ON
		ki.TimeID = hv.TimeID AND
		ki.BodyID = hb.BodyID
	# The body collection that was integrated
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'DE435'
	# Only take bodies in the planets integration!
	INNER JOIN KS.BodyCollection AS bcp ON bcp.BodyCollectionName = 'Planets'
	INNER JOIN KS.BodyCollectionEntry AS bce ON
		bce.BodyCollectionID = bcp.BodyCollectionID AND
		bce.BodyID = hb.BodyID;

END
$$

DELIMITER ;
