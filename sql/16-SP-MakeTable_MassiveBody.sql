DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_MassiveBoyd()
COMMENT "Populate the MassiveBody table from MassiveBodyImport"
BEGIN 

-- Populate HorizonsBodyName field on MassiveBodyImport table for asteroids
UPDATE 	
	JPL.MassiveBodyImport AS mbi
	INNER JOIN JPL.HorizonsBody AS hb ON hb.SmallBodyID = mbi.AsteroidNumber
SET
	mbi.HorizonsBodyName = hb.HorizonsBodyName;

INSERT IGNORE INTO JPL.MassiveBody 
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
