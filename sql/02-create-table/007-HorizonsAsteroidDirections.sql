-- Create table to import Horizons data files with asteroid directions
CREATE OR REPLACE TABLE JPL.AsteroidDirection(
	-- Key fields
	AsteroidID INT NOT NULL
		COMMENT "The AsteroidID; foreign key to KS.Asteroid table. Populated by hand.",
	TimeID INT NOT NULL
		COMMENT "Foreign key to IntegrationTime table",
	mjd DOUBLE NOT NULL,
	-- The quoted direction
	RA DOUBLE NOT NULL,
	`DEC` DOUBLE NOT NULL,
	-- Components of direction derived from RA and DEC
	ux DOUBLE NOT NULL,
	uy DOUBLE NOT NULL,
	uz DOUBLE NOT NULL,
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
	PRIMARY KEY (AsteroidID, TimeID),
	UNIQUE KEY UNQ_JPL_AsteroidDirection_TimeID_AsteroidID(TimeID, AsteroidID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Asteroid directions according to Horizons.";
