DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_IntegrationDiff()
COMMENT "Populate the IntegrationDiff table from difference between MSE and JPL solar system integrations."

BEGIN 

-- Start by emptying out the table
TRUNCATE TABLE KS.IntegrationDiff;
	
-- Insert batch of records comparing Integration_Planets to HorizonsVectors
INSERT IGNORE INTO KS.IntegrationDiff
(BodyCollectionID, TimeID, BodyID, MJD, dq, dv)
SELECT
	bc.BodyCollectionID,
	sv.TimeID,
	sv.BodyID,
    sv.MJD,
	sqrt(power(sv.qx - hv.qx, 2) + power(sv.qy - hv.qy, 2) + power(sv.qz - hv.qz, 2)) AS dq,
	sqrt(power(sv.vx - hv.vx, 2) + power(sv.vy - hv.vy, 2) + power(sv.vz - hv.vz, 2)) AS dv
FROM
	-- Start with daily snapshots of planets from Horizons
	JPL.HorizonsVectors AS hv
	-- Get the BodyID by joining to HorizonsBody
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyID = hv.HorizonsBodyID
	-- Look up corresponding record on integration table
	INNER JOIN KS.StateVectors_Planets AS sv ON
		sv.TimeID = hv.TimeID AND
		sv.BodyID = hb.BodyID
	-- The body collection that was integrated
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'Planets'
	-- Only take bodies in the planets integration!
	INNER JOIN KS.BodyCollection AS bcp ON bcp.BodyCollectionName = 'Planets'
	INNER JOIN KS.BodyCollectionEntry AS bce ON
		bce.BodyCollectionID = bcp.BodyCollectionID AND
		bce.BodyID = hb.BodyID;

-- Insert batch of records comparing Integration_DE435 to HorizonsVectors
INSERT IGNORE INTO KS.IntegrationDiff
(BodyCollectionID, TimeID, BodyID, MJD, dq, dv)
SELECT
	bc.BodyCollectionID,
	sv.TimeID,
	sv.BodyID,
    sv.MJD,
	sqrt(power(sv.qx - hv.qx, 2) + power(sv.qy - hv.qy, 2) + power(sv.qz - hv.qz, 2)) AS dq,
	sqrt(power(sv.vx - hv.vx, 2) + power(sv.vy - hv.vy, 2) + power(sv.vz - hv.vz, 2)) AS dv
FROM
	-- Start with daily snapshots of planets from Horizons
	JPL.HorizonsVectors AS hv
	-- Get the BodyID by joining to HorizonsBody
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyID = hv.HorizonsBodyID
	-- Look up corresponding record on integration table
	INNER JOIN KS.StateVectors_DE435 AS sv ON
		sv.TimeID = hv.TimeID AND
		sv.BodyID = hb.BodyID
	-- The body collection that was integrated
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'DE435'
	-- Only take bodies in the planets integration!
	INNER JOIN KS.BodyCollection AS bcp ON bcp.BodyCollectionName = 'Planets'
	INNER JOIN KS.BodyCollectionEntry AS bce ON
		bce.BodyCollectionID = bcp.BodyCollectionID AND
		bce.BodyID = hb.BodyID;

END
$$

DELIMITER ;
