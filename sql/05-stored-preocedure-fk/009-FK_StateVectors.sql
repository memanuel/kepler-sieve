DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_StateVectors()
COMMENT "Drop foreign keys on StateVectors table"
BEGIN 

ALTER TABLE KS.StateVectors	
	DROP CONSTRAINT IF EXISTS FK_StateVectors_TimeID,
	DROP CONSTRAINT IF EXISTS FK_StateVectors_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_StateVectors()
COMMENT "Restore foreign keys on StateVectors table"
BEGIN 

ALTER TABLE KS.StateVectors	
	ADD CONSTRAINT FK_StateVectors_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	ADD CONSTRAINT FK_StateVectors_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);
END
$$

DELIMITER ;
