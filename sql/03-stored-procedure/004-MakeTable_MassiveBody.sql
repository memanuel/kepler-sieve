DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_MassiveBody()
COMMENT "Populate the MassiveBody table from MassiveBodyImport"
BEGIN 

-- Populate field M from GM
UPDATE
	JPL.MassiveBodyImport AS mbi
	INNER JOIN JPL.MassiveBodyImport AS mbi2 ON mbi2.ParameterName = 'GMS'
SET
	mbi.M = (mbi.GM / mbi2.GM);

-- Populate HorizonsBodyName field on planet barycenters
UPDATE
	JPL.MassiveBodyImport AS mbi
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyID = CAST(RIGHT(mbi.ParameterName, 1) AS INT)
SET
	mbi.HorizonsBodyName = hb.HorizonsBodyName
WHERE (mbi.ParameterName LIKE 'GM_')
AND mbi.ParameterName NOT IN ('GMS', 'GM3', 'GMB', 'GMM');

-- Populate the HorizonsBodyName field on Sun, Earth, EarthMoonBaryCenter and Moon
UPDATE JPL.MassiveBodyImport AS mbi
SET mbi.HorizonsBodyName='LB.Sun'
WHERE mbi.ParameterName = 'GMS';

UPDATE JPL.MassiveBodyImport AS mbi
SET mbi.HorizonsBodyName='LB.Earth-Moon Barycenter'
WHERE mbi.ParameterName = 'GMB';

UPDATE JPL.MassiveBodyImport AS mbi
SET mbi.HorizonsBodyName='LB.Earth'
WHERE mbi.ParameterName = 'GM3';

UPDATE JPL.MassiveBodyImport AS mbi
SET mbi.HorizonsBodyName='LB.Moon'
WHERE mbi.ParameterName = 'GMM';

-- Populate the AsteroidNumber field on MassiveBodyImport from ParameterName field
UPDATE
    JPL.MassiveBodyImport AS mbi
SET 
    mbi.AsteroidNumber = CAST(RIGHT(mbi.ParameterName,4) AS INT) 
WHERE mbi.ParameterName LIKE 'MA____';

-- Populate HorizonsBodyName field on MassiveBodyImport table for asteroids
UPDATE 	
	JPL.MassiveBodyImport AS mbi
	INNER JOIN JPL.HorizonsBody AS hb ON hb.SmallBodyID = mbi.AsteroidNumber
SET
	mbi.HorizonsBodyName = hb.HorizonsBodyName;

-- Empty JPL.MassiveBody and populate it from the MassiveBodyImport table
TRUNCATE TABLE JPL.MassiveBody;

-- Populate JPL.MassiveBody from MassiveBodyImport (of CSVs)
INSERT INTO JPL.MassiveBody 
(HorizonsBodyID, M, GM)
SELECT
	hb.HorizonsBodyID,
	mbi.M,
	mbi.GM	
FROM
	JPL.MassiveBodyImport AS mbi
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyName = mbi.HorizonsBodyName;

-- Populate KS.MassiveBody from JPL.MassiveBody 
INSERT INTO KS.MassiveBody
(BodyID, M, GM)
SELECT
	b.BodyID,
	mb.M,
	mb.GM
FROM
	JPL.MassiveBody AS mb
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyID = mb.HorizonsBodyID 
	INNER JOIN KS.Body AS b ON b.BodyID = hb.BodyID;

END
$$

DELIMITER ;
