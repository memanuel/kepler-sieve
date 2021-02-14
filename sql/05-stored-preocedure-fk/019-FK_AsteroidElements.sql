DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_AsteroidElements()
COMMENT "Drop foreign keys on AsteroidElements table"
BEGIN 

ALTER TABLE KS.AsteroidElements
	DROP CONSTRAINT IF EXISTS FK_AsteroidElements_AsteroidID,
	DROP CONSTRAINT IF EXISTS FK_AsteroidElements_TimeID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_AsteroidElements()
COMMENT "Restore foreign keys on AsteroidElements table"
BEGIN 

ALTER TABLE KS.AsteroidElements
    ADD CONSTRAINT FK_AsteroidElements_AsteroidID FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID),
    ADD CONSTRAINT FK_AsteroidElements_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID);

END
$$

DELIMITER ;
