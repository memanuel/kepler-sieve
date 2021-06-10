-- Create table to import Horizons data files with asteroid directions
CREATE OR REPLACE TABLE JPL.AsteroidDirectionImport(
	RowID INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
	-- Key fields
	AsteroidID INT NOT NULL
		COMMENT "The AsteroidID; foreign key to KS.Asteroid table. Populated by hand.",
	ObservatoryID INT NOT NULL
		COMMENT "Foreign key to ObservatoryID on KS.Observatory; 0 is Geocenter",
	JD DOUBLE NOT NULL
		COMMENT "The Julian Date in the TDB (barycentric dynamical time) frame",
	-- Astrometric RA/DEC
	RA_ast DOUBLE NOT NULL
		COMMENT "Astrometric right ascension in the ICRF; in degrees",
	DEC_ast DOUBLE NOT NULL
		COMMENT "Astrometric declination in the ICRF; in degrees",
	-- Apparent RA/DEC
	RA_app DOUBLE NOT NULL
		COMMENT "Astrometric right ascension in the ICRF; in degrees",
	DEC_app DOUBLE NOT NULL
		COMMENT "Astrometric declination in the ICRF; in degrees",
	-- Magnitude
	Mag DOUBLE NOT NULL
		COMMENT "Apparent magnitude",
	Brightness DOUBLE NOT NULL
		COMMENT "Surface brightness using the H-G model",
	-- Distance between Sun/Asteroid and Earth/Asteroid and derivatives
    r DOUBLE NOT NULL
        COMMENT "Distance from Sun to asteroid; in AU",
    rDot DOUBLE NOT NULL
        COMMENT "Time derivative of r; in AU/day",
    delta DOUBLE NOT NULL
        COMMENT "Distance from Earth observer to asteroid; in AU",
    deltaDot DOUBLE NOT NULL
        COMMENT "Time derivative of delta; in AU/day",
    -- Light Time
    LightTime DOUBLE NOT NULL
        COMMENT "one way down-leg light time from target to observer; in minutes."
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Staging table to import data from files downloaded from Horizons API.";

-- Run this command to load CSV contents into JPL.HorizonsImport table
/*
LOAD DATA INFILE '/ssd1/tmp/mysql/jpl/horizons/asteroid_directions/ast_geocenter.csv'
INTO TABLE JPL.AsteroidDirectionImport
FIELDS TERMINATED BY ","
LINES TERMINATED BY "\r\n"
IGNORE 1 LINES
(AsteroidID, ObservatoryID, JD, RA_ast, DEC_ast, RA_app, DEC_app, Mag, Brightness, r, rDot, delta, deltaDot, LightTime);

LOAD DATA INFILE '/ssd1/tmp/mysql/jpl/horizons/asteroid_directions/ast_palomar.csv'
INTO TABLE JPL.AsteroidDirectionImport
FIELDS TERMINATED BY ","
LINES TERMINATED BY "\r\n"
IGNORE 1 LINES
(AsteroidID, ObservatoryID, JD, RA_ast, DEC_ast, RA_app, DEC_app, Mag, Brightness, r, rDot, delta, deltaDot, LightTime);
*/
