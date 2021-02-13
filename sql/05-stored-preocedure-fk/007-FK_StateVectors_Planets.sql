DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_StateVectors_Planets()
COMMENT "Drop foreign keys on StateVectors_Planets table"
BEGIN 

ALTER TABLE KS.StateVectors_Planets	
	DROP CONSTRAINT IF EXISTS FK_StateVectors_Planets_TimeID,
	DROP CONSTRAINT IF EXISTS FK_StateVectors_Planets_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_StateVectors_Planets()
COMMENT "Restore foreign keys on StateVectors_Planets table"
BEGIN 

ALTER TABLE KS.StateVectors_Planets	
	ADD CONSTRAINT FK_StateVectors_Planets_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	ADD CONSTRAINT FK_StateVectors_Planets_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);
END
$$

DELIMITER ;
