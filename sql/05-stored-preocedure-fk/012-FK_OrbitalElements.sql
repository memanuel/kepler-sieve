DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_OrbitalElements()
COMMENT "Drop foreign keys on OrbitalElements table"
BEGIN 

ALTER TABLE KS.OrbitalElements	
	DROP CONSTRAINT IF EXISTS FK_OrbitalElements_TimeID,
	DROP CONSTRAINT IF EXISTS FK_OrbitalElements_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_OrbitalElements()
COMMENT "Restore foreign keys on OrbitalElements table"
BEGIN 

ALTER TABLE KS.OrbitalElements	
	ADD CONSTRAINT FK_OrbitalElements_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	ADD CONSTRAINT FK_OrbitalElements_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);
END
$$

DELIMITER ;
