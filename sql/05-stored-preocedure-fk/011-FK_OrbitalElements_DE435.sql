DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_OrbitalElements_DE435()
COMMENT "Drop foreign keys on OrbitalElements_DE435 table"
BEGIN 

ALTER TABLE KS.OrbitalElements_DE435	
	DROP CONSTRAINT IF EXISTS FK_OrbitalElements_DE435_TimeID,
	DROP CONSTRAINT IF EXISTS FK_OrbitalElements_DE435_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_OrbitalElements_DE435()
COMMENT "Restore foreign keys on OrbitalElements_DE435 table"
BEGIN 

ALTER TABLE KS.OrbitalElements_DE435	
	ADD CONSTRAINT FK_OrbitalElements_DE435_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	ADD CONSTRAINT FK_OrbitalElements_DE435_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);
END
$$

DELIMITER ;
