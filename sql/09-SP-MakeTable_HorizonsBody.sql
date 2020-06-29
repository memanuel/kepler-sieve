DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_HorizonsBody()
BEGIN 

-- Insert large bodies into the HorizonsBody table
INSERT INTO JPL.HorizonsBody
(HorizonsBodyID, HorizonsBodyName, BodyTypeID, BodyID, LargeBodyID, SmallBodyID)
SELECT 
	lb.BodyID as HorizonsBodyID,
	CONCAT('LB.', lb.LargeBodyName) as HorizonsBodyName,
	lb.BodyTypeID,
	lb.BodyID,
	lb.LargeBodyID,
	NULL AS SmallBodyID
FROM
	JPL.LargeBody as lb;

-- Insert small bodies into the HorizonsBody table
INSERT INTO JPL.HorizonsBody
(HorizonsBodyID, HorizonsBodyName, BodyTypeID, BodyID, LargeBodyID, SmallBodyID)
SELECT 
	sb.BodyID as HorizonsBodyID,
	CONCAT('SB.', sb.SmallBodyName) as HorizonsBodyName,
	sb.BodyTypeID,
	sb.BodyID,
	NULL AS LargeBodyID,
	sb.SmallBodyID
FROM
	JPL.SmallBody as sb;

END
$$

DELIMITER ;