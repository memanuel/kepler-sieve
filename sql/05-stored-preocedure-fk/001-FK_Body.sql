DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_Body()
COMMENT "Drop foreign keys on Body table"
BEGIN 

ALTER TABLE KS.Body	DROP CONSTRAINT IF EXISTS FK_Body_BodyTypeID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_Body()
COMMENT "Restore foreign keys on Body table"
BEGIN 

ALTER TABLE KS.Body	ADD CONSTRAINT FK_Body_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES BodyType(BodyTypeID);
	
END
$$

DELIMITER ;
