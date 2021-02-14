DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.DropFK_MassiveBody()
COMMENT "Drop foreign keys on MassiveBody table"
BEGIN 

ALTER TABLE JPL.MassiveBody	
	DROP CONSTRAINT IF EXISTS FK_MassiveBody_HorizonsBodyID;
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.RestoreFK_MassiveBody()
COMMENT "Restore foreign keys on MassiveBody table"
BEGIN 

ALTER TABLE JPL.MassiveBody
	ADD CONSTRAINT FK_MassiveBody_HorizonsBodyID FOREIGN KEY (HorizonsBodyID) REFERENCES JPL.HorizonsBody(HorizonsBodyID);

END
$$

DELIMITER ;
