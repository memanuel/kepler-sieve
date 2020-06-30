DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_HorizonsBody()
COMMENT "Populate the HorizonsBody table from LargeBody and SmallBody"
BEGIN 

-- Insert large bodies into the HorizonsBody table
INSERT INTO JPL.HorizonsBody
(HorizonsBodyID, HorizonsBodyNumber, HorizonsBodyName, BodyTypeID, BodyID, LargeBodyID, SmallBodyID)
SELECT 
	lb.BodyID as HorizonsBodyID,
	lb.LargeBodyID AS BodyNumber,
	CONCAT('LB.', lb.LargeBodyName) as HorizonsBodyName,
	lb.BodyTypeID,
	lb.BodyID,
	lb.LargeBodyID,
	NULL AS SmallBodyID
FROM
	JPL.LargeBody as lb
ON DUPLICATE KEY UPDATE HorizonsBodyID=HorizonsBodyID;

-- Insert small bodies into the HorizonsBody table
INSERT INTO JPL.HorizonsBody
(HorizonsBodyID, HorizonsBodyNumber, HorizonsBodyName, BodyTypeID, BodyID, LargeBodyID, SmallBodyID)
SELECT 
	sb.BodyID as HorizonsBodyID,
	sb.SmallBodyID AS BodyNumber,
	CONCAT('SB.', sb.SmallBodyName) as HorizonsBodyName,
	sb.BodyTypeID,
	sb.BodyID,
	NULL AS LargeBodyID,
	sb.SmallBodyID
FROM
	JPL.SmallBody as sb
ON DUPLICATE KEY UPDATE HorizonsBodyID=HorizonsBodyID;

END
$$

DELIMITER ;
