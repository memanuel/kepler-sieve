DELIMITER $$

-- Drop foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.DropFK_All()
COMMENT "Drop foreign keys all tables"
BEGIN 
    -- JPL Tables
	CALL JPL.DropFK_LargeBody();
	CALL JPL.DropFK_SmallBody();
	CALL JPL.DropFK_HorizonsBody();
	CALL JPL.DropFK_MassiveBody();
    -- KS Tables
	CALL KS.DropFK_Body();
	CALL KS.DropFK_BodyCollectionEntry();
	CALL KS.DropFK_BarycenterWeight();
	CALL KS.DropFK_StateVectors_Planets();
	CALL KS.DropFK_StateVectors_DE435();
	CALL KS.DropFK_StateVectors();
	CALL KS.DropFK_OrbitalElements_Planets();
	CALL KS.DropFK_OrbitalElements_DE435();
	CALL KS.DropFK_OrbitalElements();
	CALL KS.DropFK_MassiveBody();
	CALL KS.DropFK_PrimaryBody();
	CALL KS.DropFK_Asteroid();
	CALL KS.DropFK_AsteroidElement_Ref();
	CALL KS.DropFK_AsteroidVectors();
	CALL KS.DropFK_AsteroidElements();
	
END
$$

-- Restore foreign keys
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.RestoreFK_All()
COMMENT "Restore foreign keys on all tables"
BEGIN 
    -- JPL Tables
	CALL JPL.RestoreFK_LargeBody();
	CALL JPL.RestoreFK_SmallBody();
	CALL JPL.RestoreFK_HorizonsBody();
	CALL JPL.RestoreFK_MassiveBody();
    -- KS Tables
	CALL KS.RestoreFK_Body();
	CALL KS.RestoreFK_BodyCollectionEntry();
	CALL KS.RestoreFK_BarycenterWeight();
	CALL KS.RestoreFK_StateVectors_Planets();
	CALL KS.RestoreFK_StateVectors_DE435();
	CALL KS.RestoreFK_StateVectors();
	CALL KS.RestoreFK_OrbitalElements_Planets();
	CALL KS.RestoreFK_OrbitalElements_DE435();
	CALL KS.RestoreFK_OrbitalElements();
	CALL KS.RestoreFK_MassiveBody();
	CALL KS.RestoreFK_PrimaryBody();
	CALL KS.RestoreFK_Asteroid();
	CALL KS.RestoreFK_AsteroidElement_Ref();
	CALL KS.RestoreFK_AsteroidVectors();
	CALL KS.RestoreFK_AsteroidElements();

END
$$

DELIMITER ;
