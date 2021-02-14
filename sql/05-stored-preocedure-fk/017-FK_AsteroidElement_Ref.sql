DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_AsteroidElement_Ref()
COMMENT "Drop foreign keys on AsteroidElement_Ref table"
BEGIN 

ALTER TABLE KS.AsteroidElement_Ref
	DROP CONSTRAINT IF EXISTS FK_AsteroidElement_Ref_AsteroidID,
	DROP CONSTRAINT IF EXISTS FK_AsteroidElement_Ref_TimeID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_AsteroidElement_Ref()
COMMENT "Restore foreign keys on AsteroidElement_Ref table"
BEGIN 

ALTER TABLE KS.AsteroidElement_Ref
    ADD CONSTRAINT FK_AsteroidElement_Ref_AsteroidID FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID),
    ADD CONSTRAINT FK_AsteroidElement_Ref_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID);

END
$$

DELIMITER ;
