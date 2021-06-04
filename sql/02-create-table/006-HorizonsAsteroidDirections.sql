-- Create table to import Horizons data files with asteroid directions
CREATE OR REPLACE TABLE JPL.AsteroidDirectionImport(
	RowID INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
	AsteroidID INT NOT NULL
		COMMENT "The AsteroidID; foreign key to KS.Asteroid table. Populated by hand.",
	JD DOUBLE NOT NULL
		COMMENT "The Julian Date in the TDB (barycentric dynamical time) frame",
	RA DOUBLE NOT NULL
		COMMENT "Astrometric right ascension in the ICRF; in degrees",
	`DEC` DOUBLE NOT NULL
		COMMENT "Astrometric declination in the ICRF; in degrees",
	dRAxCosD_dt DOUBLE NOT NULL
		COMMENT "Time derivative of RA*cos(D) in arc seconds per hour; multiplied by cos(DEC) to provide linear rate in the plane of sky",
	dDEC_dt DOUBLE NOT NULL
		COMMENT "Time derivative of declination in arc seconds per hour.",
	Mag DOUBLE NOT NULL
		COMMENT "Apparent magnitude",
	Brightness DOUBLE NOT NULL
		COMMENT "Surface brightness using the H-G model",
    r DOUBLE NOT NULL
        COMMENT "Distance from Sun to asteroid; in AU",
    rDot DOUBLE NOT NULL
        COMMENT "Time derivative of r; in AU/day",
    delta DOUBLE NOT NULL
        COMMENT "Distance from Earth observer to asteroid; in AU",
    deltaDot DOUBLE NOT NULL
        COMMENT "Time derivative of delta; in AU/day",
    LightTime DOUBLE NOT NULL
        COMMENT "one way down-leg light time from target to observer; in minutes."
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Staging table to import data from files downloaded from Horizons API.";

-- Run this command to load CSV contents into JPL.HorizonsImport table
/*
LOAD DATA INFILE '/ssd1/tmp/mysql/jpl/horizons/asteroid_directions/ast_0001.csv'
INTO TABLE JPL.AsteroidDirectionImport
FIELDS TERMINATED BY ","
LINES TERMINATED BY "\n"
IGNORE 1 LINES
(AsteroidID, JD, RA, `DEC`, dRAxCosD_dt, dDEC_dt, Mag, Brightness, r, rDot, delta, deltaDot, LightTime);
*/
