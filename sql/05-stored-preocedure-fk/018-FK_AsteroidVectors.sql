DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_AsteroidVectors()
COMMENT "Drop foreign keys on AsteroidVectors table"
BEGIN 

ALTER TABLE KS.AsteroidVectors
	DROP CONSTRAINT IF EXISTS FK_AsteroidVectors_AsteroidID,
	DROP CONSTRAINT IF EXISTS FK_AsteroidVectors_TimeID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_AsteroidVectors()
COMMENT "Restore foreign keys on AsteroidVectors table"
BEGIN 

ALTER TABLE KS.AsteroidVectors
    ADD CONSTRAINT FK_AsteroidVectors_AsteroidID FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID),
    ADD CONSTRAINT FK_AsteroidVectors_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID);

END
$$

DELIMITER ;
