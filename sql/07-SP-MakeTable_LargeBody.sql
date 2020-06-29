DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_LargeBody()
BEGIN 

-- Insert into the LargeBody table
INSERT INTO JPL.LargeBody
(LargeBodyID, LargeBodyName, BodyTypeID, BodyID)
SELECT 
	hi.BodyNumber AS LargeBodyID, 
	hi.BodyName AS LargeBodyName,
	bt.BodyTypeID,
	hi.BodyNumber + 0 as BodyID
FROM 
	JPL.HorizonsImport AS hi 
	INNER JOIN KeplerDB.BodyType AS bt ON bt.BodyTypeCD = hi.BodyTypeCD
WHERE bt.IsLargeBody_JPL	
GROUP BY hi.BodyTypeCD, hi.BodyNumber
ORDER BY bt.BodyTypeID, hi.BodyNumber;

END
$$

DELIMITER ;