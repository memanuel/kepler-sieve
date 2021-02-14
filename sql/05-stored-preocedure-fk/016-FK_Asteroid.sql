DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_Asteroid()
COMMENT "Drop foreign keys on Asteroid table"
BEGIN 

ALTER TABLE KS.Asteroid	
	DROP CONSTRAINT IF EXISTS FK_Asteroid_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_Asteroid()
COMMENT "Restore foreign keys on Asteroid table"
BEGIN 

ALTER TABLE KS.Asteroid
	ADD CONSTRAINT FK_Asteroid_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);

END
$$

DELIMITER ;
