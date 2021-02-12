DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.DropFK_SmallBody()
COMMENT "Drop foreign keys on Body table"
BEGIN 

ALTER TABLE JPL.SmallBody	
	DROP CONSTRAINT IF EXISTS FK_SmallBody_BodyTypeID,
	DROP CONSTRAINT IF EXISTS FK_SmallBody_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.RestoreFK_SmallBody()
COMMENT "Restore foreign keys on Body table"
BEGIN 

ALTER TABLE JPL.SmallBody	
	ADD CONSTRAINT FK_SmallBody_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES KS.BodyType(BodyTypeID),
	ADD CONSTRAINT FK_SmallBody_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);

END
$$

DELIMITER ;
