DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_StateVectors_DE435()
COMMENT "Drop foreign keys on StateVectors_DE435 table"
BEGIN 

ALTER TABLE KS.StateVectors_DE435	
	DROP CONSTRAINT IF EXISTS FK_StateVectors_DE435_TimeID,
	DROP CONSTRAINT IF EXISTS FK_StateVectors_DE435_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_StateVectors_DE435()
COMMENT "Restore foreign keys on StateVectors_DE435 table"
BEGIN 

ALTER TABLE KS.StateVectors_DE435	
	ADD CONSTRAINT FK_StateVectors_DE435_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	ADD CONSTRAINT FK_StateVectors_DE435_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);
END
$$

DELIMITER ;
