-- Create table to import Horizons data files with asteroid directions
CREATE OR REPLACE TABLE JPL.AsteroidDirection(
	-- Key fields
	AsteroidID INT NOT NULL
		COMMENT "The AsteroidID; foreign key to KS.Asteroid table. Populated by hand.",
	ObservatoryID INT NOT NULL
		COMMENT "Foreign key to KS.Observatory table",
	TimeID INT NOT NULL
		COMMENT "Foreign key to IntegrationTime table",
	-- Time as a number (derivable from TimeID)
	mjd DOUBLE NOT NULL,
	-- The quoted astrometric direction
	RA_ast DOUBLE NOT NULL
		COMMENT "The astrometric right ascension quoted by JPL",
	DEC_ast DOUBLE NOT NULL
		COMMENT "The astrometric declination quoted by JPL",
	-- Components of direction derived from astrometric RA and DEC
	ux_ast DOUBLE NOT NULL,
	uy_ast DOUBLE NOT NULL,
	uz_ast DOUBLE NOT NULL,
	-- The quoted apparent direction
	RA_app DOUBLE NOT NULL
		COMMENT "The apparent right ascension quoted by JPL",
	DEC_app DOUBLE NOT NULL
		COMMENT "The apparent declination quoted by JPL",
	-- Components of direction derived from apparent RA and DEC
	ux_app DOUBLE NOT NULL,
	uy_app DOUBLE NOT NULL,
	uz_app DOUBLE NOT NULL,
	-- Distance from Sun
   r DOUBLE NOT NULL
       COMMENT "Distance from Sun to asteroid; in AU",
   rDot DOUBLE NOT NULL
       COMMENT "Time derivative of r; in AU/day",
    -- Distance from Earth
   delta DOUBLE NOT NULL
       COMMENT "Distance from Earth observer to asteroid; in AU",
   deltaDot DOUBLE NOT NULL
       COMMENT "Time derivative of delta; in AU/day",
    -- Evertyghing else
    LightTime DOUBLE NOT NULL
        COMMENT "one way down-leg light time from target to observer; in minutes.",
	Mag DOUBLE NOT NULL
		COMMENT "Apparent magnitude",
	PRIMARY KEY (AsteroidID, ObservatoryID, TimeID),
	UNIQUE KEY UNQ_TimeID_AsteroidID(TimeID, ObservatoryID, AsteroidID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Asteroid directions according to Horizons.";
