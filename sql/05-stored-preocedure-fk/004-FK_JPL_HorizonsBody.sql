DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.DropFK_HorizonsBody()
COMMENT "Drop foreign keys on HorizonsBody table"
BEGIN 

ALTER TABLE JPL.HorizonsBody	
	DROP CONSTRAINT IF EXISTS FK_HorizonsBody_BodyTypeID,
	DROP CONSTRAINT IF EXISTS FK_HorizonsBody_BodyID,
	DROP CONSTRAINT IF EXISTS FK_HorizonsBody_LargeBodyID,
	DROP CONSTRAINT IF EXISTS FK_HorizonsBody_SmallBodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.RestoreFK_HorizonsBody()
COMMENT "Restore foreign keys on HorizonsBody table"
BEGIN 

ALTER TABLE JPL.HorizonsBody	
	ADD CONSTRAINT FK_HorizonsBody_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES KS.BodyType(BodyTypeID),
	ADD CONSTRAINT FK_HorizonsBody_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID),
	ADD CONSTRAINT FK_HorizonsBody_LargeBodyID FOREIGN KEY (LargeBodyID) REFERENCES JPL.LargeBody(LargeBodyID),
	ADD CONSTRAINT FK_HorizonsBody_SmallBodyID FOREIGN KEY (SmallBodyID) REFERENCES JPL.SmallBody(SmallBodyID);
END
$$

DELIMITER ;
