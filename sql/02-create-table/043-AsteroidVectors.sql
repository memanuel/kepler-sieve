-- Create tables for integrated vectors for asteroids
CREATE OR REPLACE TABLE KS.AsteroidVectors(
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these state vectors; FK to KS.IntegrationTime",
	AsteroidID INT NOT NULL
		COMMENT "The asteroid whose state vectors are described; FK to KS.Asteroid",
	mjd DOUBLE NOT NULL
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame; derivable from TimeID but included for performance.",
	-- Position q = [qx, qy, qz]
	qx DOUBLE NOT NULL
		COMMENT "Position of body (x coordinate) in AU in the barcycentric mean ecliptic frame",
	qy DOUBLE NOT NULL
		COMMENT "Position of body (y coordinate) in AU in the barcycentric mean ecliptic frame",
	qz DOUBLE NOT NULL
		COMMENT "Position of body (z coordinate) in AU in the barcycentric mean ecliptic frame",
	-- Velocity v = [vx, vy, vz]
	vx DOUBLE NOT NULL
		COMMENT "Velocity of body (x coordinate) in AU/day in the barcycentric mean ecliptic frame",
	vy DOUBLE NOT NULL
		COMMENT "Velocity of body (y coordinate) in AU/day in the barcycentric mean ecliptic frame",
	vz DOUBLE NOT NULL
		COMMENT "Velocity of body (z coordinate) in AU/day in the barcycentric mean ecliptic frame",
	-- Keys and constraints
	PRIMARY KEY (TimeID, AsteroidID)
		COMMENT "A state vector is identified by the body and time stamp; use integer time ID for performance.",
	UNIQUE KEY UNQ_AsteroidID_TimeID(AsteroidID, TimeID)
        COMMENT "Allow fast search keyed first by AsteroidID.",
	CONSTRAINT FK_AsteroidVectors_TimeID FOREIGN KEY (TimeID) REFERENCES KS.DailyTime(TimeID),
	CONSTRAINT FK_AsteroidVectors_BodyID	FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "State vectors (position and velocity) for asteroids computed in Rebound using the planets as massive bodies and initial conditions from DE435 at MJD 59000.
Initial conditions for asteroids sourced from AsteroidElement_Ref"
