-- Create tables for integrated directions from Earth to asteroids
CREATE OR REPLACE TABLE KS.AsteroidDirections(
	AsteroidID INT NOT NULL
		COMMENT "The asteroid whose state vectors are described; FK to KS.Asteroid",
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp when light is LEAVING the asteroid (t_ast); FK to KS.IntegrationTime",
	mjd DOUBLE NOT NULL
		COMMENT "The MJD in the TDB frame when when light leaves the asteroid (t_ast); denormalized from TimeID.",
	-- Direction u = [ux, uy, uz]
	ux DOUBLE NOT NULL
		COMMENT "Position of body (x coordinate) in AU in the barcycentric mean ecliptic frame",
	uy DOUBLE NOT NULL
		COMMENT "Position of body (y coordinate) in AU in the barcycentric mean ecliptic frame",
	uz DOUBLE NOT NULL
		COMMENT "Position of body (z coordinate) in AU in the barcycentric mean ecliptic frame",
	-- Light time
	LightTime DOUBLE NOT NULL
		COMMENT "Time for light leaving asteroid to reach Earth in MINUTES not days.",
	-- Keys and constraints
	PRIMARY KEY (AsteroidID, TimeID)
		COMMENT "A state vector is identified by the body and time stamp; use integer time ID for performance.",
	UNIQUE KEY UNQ_TimeID_AsteroidID(TimeID, AsteroidID)
        COMMENT "Allow fast search keyed first by TimeID.",
	CONSTRAINT FK_AsteroidDirections_TimeID FOREIGN KEY (TimeID) REFERENCES KS.DailyTime(TimeID),
	CONSTRAINT FK_AsteroidDirections_AsteroidID	FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Directions from Earth to asteroid for all known asteroids computed in rebound using the planets as massive bodies and initial conditions from DE435 at MJD 59000."