DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_BodyCollectionEntry()
COMMENT "Drop foreign keys on BodyCollectionEntry table"
BEGIN 

ALTER TABLE KS.BodyCollectionEntry	
	DROP CONSTRAINT IF EXISTS FK_BodyCollectionEntry_BodyCollectionID,
	DROP CONSTRAINT IF EXISTS FK_BodyCollectionEntry_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_BodyCollectionEntry()
COMMENT "Restore foreign keys on BodyCollectionEntry table"
BEGIN 

ALTER TABLE KS.BodyCollectionEntry	
	ADD CONSTRAINT FK_BodyCollectionEntry_BodyCollectionID
		FOREIGN KEY (BodyCollectionID) REFERENCES KS.BodyCollection(BodyCollectionID),
	ADD CONSTRAINT FK_BodyCollectionEntry_BodyID
		FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);
END
$$

DELIMITER ;
