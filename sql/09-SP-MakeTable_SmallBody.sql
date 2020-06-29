DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_SmallBody()
BEGIN 

-- Insert into the SmallBody table
INSERT INTO JPL.SmallBody
(SmallBodyID, SmallBodyName, BodyTypeID, BodyID)
SELECT 
	hi.BodyNumber AS LargeBodyID, 
	hi.BodyName AS SmallBodyName,
	bt.BodyTypeID,
	hi.BodyNumber + 1000000 as BodyID
FROM 
	JPL.HorizonsImport AS hi 
	INNER JOIN KeplerDB.BodyType AS bt ON bt.BodyTypeCD = hi.BodyTypeCD
WHERE NOT bt.IsLargeBody_JPL	
GROUP BY hi.BodyTypeCD, hi.BodyNumber
ORDER BY bt.BodyTypeID, hi.BodyNumber;

END
$$

DELIMITER ;