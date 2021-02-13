DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_BarycenterWeight()
COMMENT "Drop foreign keys on BarycenterWeight table"
BEGIN 

ALTER TABLE KS.BarycenterWeight	
	DROP CONSTRAINT IF EXISTS FK_BarycenterWeight_BodyCollectionID,
	DROP CONSTRAINT IF EXISTS FK_BarycenterWeight_BodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_BarycenterWeight()
COMMENT "Restore foreign keys on BarycenterWeight table"
BEGIN 

ALTER TABLE KS.BarycenterWeight	
	ADD CONSTRAINT FK_BarycenterWeight_BodyCollectionID
		FOREIGN KEY (BodyCollectionID) REFERENCES KS.BodyCollection(BodyCollectionID),
	ADD CONSTRAINT FK_BarycenterWeight_BodyID
		FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);
END
$$

DELIMITER ;
