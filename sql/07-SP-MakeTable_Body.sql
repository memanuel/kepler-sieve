DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KeplerDB.MakeTable_Body()
BEGIN 

-- Truncate the table the inefficient way to avoid permissions problems and restoring foreign keys
DELETE FROM JPL.LargeBody;
DELETE FROM JPL.SmallBody;
DELETE FROM KeplerDB.Body;
	
-- Insert the large bodies; these have BodyID = BodyNumber	
INSERT INTO KeplerDB.Body
(BodyID, BodyName, BodyTypeID)
SELECT 
	hi.BodyNumber AS BodyID, 
	hi.BodyName,
	bt.BodyTypeID
FROM 
	JPL.HorizonsImport AS hi 
	INNER JOIN KeplerDB.BodyType AS bt ON bt.BodyTypeCD = hi.BodyTypeCD
WHERE bt.IsLargeBody_JPL	
GROUP BY hi.BodyTypeCD, hi.BodyNumber
ORDER BY bt.BodyTypeID, hi.BodyNumber;
	
-- Insert the small bodies; these have BodyID = BodyNumber	+ 1000000
INSERT INTO KeplerDB.Body
(BodyID, BodyName, BodyTypeID)
SELECT 
	hi.BodyNumber + 1000000 AS BodyID, 
	hi.BodyName,
	bt.BodyTypeID
FROM 
	JPL.HorizonsImport AS hi 
	INNER JOIN KeplerDB.BodyType AS bt ON bt.BodyTypeCD = hi.BodyTypeCD
WHERE NOT bt.IsLargeBody_JPL	
GROUP BY hi.BodyTypeCD, hi.BodyNumber
ORDER BY bt.BodyTypeID, hi.BodyNumber;

END
$$

DELIMITER ;