DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_MassiveBody()
COMMENT "Drop foreign keys on MassiveBody table"
BEGIN 

ALTER TABLE KS.MassiveBody	
	DROP CONSTRAINT IF EXISTS FK_MassiveBody_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_MassiveBody()
COMMENT "Restore foreign keys on MassiveBody table"
BEGIN 

ALTER TABLE KS.MassiveBody
	ADD CONSTRAINT FK_MassiveBody_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);
END
$$

DELIMITER ;
