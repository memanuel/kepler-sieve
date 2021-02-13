DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.DropFK_LargeBody()
COMMENT "Drop foreign keys on Body table"
BEGIN 

ALTER TABLE JPL.LargeBody	
	DROP CONSTRAINT IF EXISTS FK_LargeBody_BodyTypeID,
	DROP CONSTRAINT IF EXISTS FK_LargeBody_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.RestoreFK_LargeBody()
COMMENT "Restore foreign keys on Body table"
BEGIN 

ALTER TABLE JPL.LargeBody	
	ADD CONSTRAINT FK_LargeBody_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES KS.BodyType(BodyTypeID),
	ADD CONSTRAINT FK_LargeBody_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);

END
$$

DELIMITER ;
