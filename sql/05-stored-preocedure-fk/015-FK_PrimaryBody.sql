DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_PrimaryBody()
COMMENT "Drop foreign keys on PrimaryBody table"
BEGIN 

ALTER TABLE KS.PrimaryBody	
	DROP CONSTRAINT IF EXISTS FK_PrimaryBody_BodyID,
	DROP CONSTRAINT IF EXISTS FK_PrimaryBody_PrimaryBodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_PrimaryBody()
COMMENT "Restore foreign keys on PrimaryBody table"
BEGIN 

ALTER TABLE KS.PrimaryBody
	ADD CONSTRAINT FK_PrimaryBody_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID),
    ADD CONSTRAINT FK_PrimaryBody_PrimaryBodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);

END
$$

DELIMITER ;
