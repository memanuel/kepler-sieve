CREATE OR REPLACE TABLE JPL.MassiveBodyImport(
	ParameterName VARCHAR(6) NOT NULL PRIMARY KEY,
    HorizonsBodyName VARCHAR(32) NULL
        COMMENT "Reference to the HorizonsBody table where applicable",
    AsteroidNumber INT NULL
    	COMMENT "IAU asteroid number where applicable",
	M DOUBLE NOT NULL
		COMMENT "Mass in solar equivalents (sun = 1.0)",
	GM DOUBLE NOT NULL
		COMMENT "Gravitational field strength (G x mass) in AU^3 / Day^2"
)
COMMENT "Mass of heavy objects included in DE 435 integration, sources from technical comments.  Import version of table keyed by HorizonsBodyName.";


CREATE OR REPLACE TABLE JPL.MassiveBody(
    HorizonsBodyID INT NOT NULL PRIMARY KEY
        COMMENT "Reference to the HorizonsBody table",
	M DOUBLE NOT NULL
		COMMENT "Mass in solar equivalents (sun = 1.0)",
	GM DOUBLE NOT NULL
		COMMENT "Gravitational field strength (G x mass) in AU^3 / Day^2",
    CONSTRAINT FK_MassiveBody_HorizonsBody
        FOREIGN KEY (HorizonsBodyID) REFERENCES JPL.HorizonsBody(HorizonsBodyID)
)
COMMENT "Mass of heavy objects included in DE 435 integration, sources from technical comments.";