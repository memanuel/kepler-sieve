DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_LargeBody()
COMMENT "Populate the LargeBody table from HorizonsImport"
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
	INNER JOIN KS.BodyType AS bt ON bt.BodyTypeCD = hi.BodyTypeCD
WHERE bt.IsLargeBody_JPL	
GROUP BY hi.BodyTypeCD, hi.BodyNumber
ORDER BY bt.BodyTypeID, hi.BodyNumber
ON DUPLICATE KEY UPDATE LargeBodyID = LargeBodyID;

END
$$

DELIMITER ;
