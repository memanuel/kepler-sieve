DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_OrbitalElements_Planets()
COMMENT "Drop foreign keys on OrbitalElements_Planets table"
BEGIN 

ALTER TABLE KS.OrbitalElements_Planets	
	DROP CONSTRAINT IF EXISTS FK_OrbitalElements_Planets_TimeID,
	DROP CONSTRAINT IF EXISTS FK_OrbitalElements_Planets_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_OrbitalElements_Planets()
COMMENT "Restore foreign keys on OrbitalElements_Planets table"
BEGIN 

ALTER TABLE KS.OrbitalElements_Planets	
	ADD CONSTRAINT FK_OrbitalElements_Planets_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	ADD CONSTRAINT FK_OrbitalElements_Planets_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);
END
$$

DELIMITER ;
