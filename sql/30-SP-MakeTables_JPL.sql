DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTables_Horizons()
COMMENT "Regenerate all tables with Horizons data in JPL database."
BEGIN 

CALL KS.MakeTable_Body();
CALL JPL.MakeTable_LargeBody();
CALL JPL.MakeTable_SmallBody();
CALL JPL.MakeTable_HorizonsBody();
CALL JPL.MakeTable_HorizonsTime();
CALL JPL.MakeTable_HorizonsVectors();
CALL JPL.MakeTable_IntegrationDiff();
CALL JPL.MakeTable_BodyVariance()
CALL JPL.MakeTable_MassiveBody();

END $$

DELIMITER ;
