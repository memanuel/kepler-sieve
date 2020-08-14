DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_MassiveBody()
COMMENT "Populate the MassiveBody table from MassiveBodyImport"
BEGIN 

-- Populate HorizonsBodyName field on MassiveBodyImport table for asteroids
UPDATE 	
	JPL.MassiveBodyImport AS mbi
	INNER JOIN JPL.HorizonsBody AS hb ON hb.SmallBodyID = mbi.AsteroidNumber
SET
	mbi.HorizonsBodyName = hb.HorizonsBodyName;

-- Empty JPL.MassiveBody and populate it from the MassiveBodyImport table
TRUNCATE TABLE JPL.MassiveBody;

INSERT INTO JPL.MassiveBody 
(HorizonsBodyID, M, GM)
SELECT
	hb.HorizonsBodyID,
	mbi.M,
	mbi.GM	
FROM
	JPL.MassiveBodyImport AS mbi
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyName = mbi.HorizonsBodyName;

END
$$

DELIMITER ;
