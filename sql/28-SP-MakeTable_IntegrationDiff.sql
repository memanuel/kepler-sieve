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
	bcp.BodyCollectionID,
	idp.TimeID,
	idp.BodyID,
	idp.MJD,
	idp.dq,
	idp.dv
FROM
	KS.IntegrationDiff_Planets AS idp
	INNER JOIN KS.BodyCollection AS bcp ON bcp.BodyCollectionName='Planets';

# Insert batch of records comparing Integration_DE435 to HorizonsVectors
INSERT IGNORE INTO KS.IntegrationDiff
(BodyCollectionID, TimeID, BodyID, MJD, dq, dv)
SELECT
	bcd.BodyCollectionID,
	idd.TimeID,
	idd.BodyID,
	idd.MJD,
	idd.dq,
	idd.dv
FROM
	KS.IntegrationDiff_DE435 AS idd
	INNER JOIN KS.BodyCollection AS bcd ON bcd.BodyCollectionName='DE435';

END
$$

DELIMITER ;
