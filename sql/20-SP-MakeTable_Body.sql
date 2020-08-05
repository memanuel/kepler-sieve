DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_Body()
COMMENT "Populate the Body table from JPL.HorizonsImport"
BEGIN 

-- Temporarily turn off foreign keys	
SET FOREIGN_KEY_CHECKS=0;
	
-- Truncate body tables to eliminate dependencies and avoid foreign key errors
TRUNCATE TABLE JPL.HorizonsVectors;
TRUNCATE TABLE JPL.HorizonsBody;
TRUNCATE TABLE JPL.LargeBody;
TRUNCATE TABLE JPL.SmallBody;
	
-- Insert the large bodies; these have BodyID = BodyNumber	
INSERT INTO KS.Body
(BodyID, BodyName, BodyTypeID)
SELECT 
	hi.BodyNumber AS BodyID, 
	CONCAT('LB.', hi.BodyName) as BodyName,
	bt.BodyTypeID
FROM 
	JPL.HorizonsImport AS hi 
	INNER JOIN KS.BodyType AS bt ON bt.BodyTypeCD = hi.BodyTypeCD
WHERE bt.IsLargeBody_JPL	
GROUP BY hi.BodyTypeCD, hi.BodyNumber
ORDER BY bt.BodyTypeID, hi.BodyNumber
ON DUPLICATE KEY UPDATE BodyID=BodyID;
	
-- Insert the small bodies; these have BodyID = BodyNumber	+ 1000000
INSERT INTO KS.Body
(BodyID, BodyName, BodyTypeID)
SELECT 
	hi.BodyNumber + 1000000 AS BodyID, 
	CONCAT('SB.', hi.BodyName) as BodyName,
	bt.BodyTypeID
FROM 
	JPL.HorizonsImport AS hi 
	INNER JOIN KS.BodyType AS bt ON bt.BodyTypeCD = hi.BodyTypeCD
WHERE NOT bt.IsLargeBody_JPL	
GROUP BY hi.BodyTypeCD, hi.BodyNumber
ORDER BY bt.BodyTypeID, hi.BodyNumber
ON DUPLICATE KEY UPDATE BodyID=BodyID;


-- Restore foreign keys	
SET FOREIGN_KEY_CHECKS=1;

END
$$

DELIMITER ;
