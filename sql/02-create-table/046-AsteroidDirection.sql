-- Create tables for integrated directions from Earth to asteroids
CREATE OR REPLACE TABLE KS.AsteroidDirections(
	AsteroidID INT NOT NULL
		COMMENT "The asteroid whose state vectors are described; FK to KS.Asteroid",
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp when light is ARRIVING on Earth (tObs); FK to KS.IntegrationTime",
	tObs DOUBLE NOT NULL
		COMMENT "The MJD in the TDB frame when when light arrives at Earth geocenter; denormalized from TimeID.",
	-- Distance from asteroid to earth
	-- delta DOUBLE NOT NULL,
	-- delta_dot DOUBLE NOT NULL,
	-- Direction u = [ux, uy, uz]
	ux DOUBLE NOT NULL
		COMMENT "Direction of body (x coordinate) from Earth geocenter in the BME frame",
	uy DOUBLE NOT NULL
		COMMENT "Direction of body (y coordinate) from Earth geocenter in the BME frame",
	uz DOUBLE NOT NULL
		COMMENT "Direction of body (z coordinate) from Earth geocenter in the BME frame",
	-- Light time
	LightTime DOUBLE NOT NULL
		COMMENT "Time for light leaving asteroid to reach Earth in MINUTES not days; tAst = tObs - LightTime / 1440.0.",
	-- Keys and constraints
	PRIMARY KEY (AsteroidID, TimeID)
		COMMENT "A state vector is identified by the body and time stamp; use integer time ID for performance.",
	UNIQUE KEY UNQ_TimeID_AsteroidID(TimeID, AsteroidID)
        COMMENT "Allow fast search keyed first by TimeID.",
	CONSTRAINT FK_AsteroidDirections_TimeID FOREIGN KEY (TimeID) REFERENCES KS.DailyTime(TimeID),
	CONSTRAINT FK_AsteroidDirections_AsteroidID	FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Directions from Earth geocenter to asteroid; keyed by observation time. Calculated from MSE solar system integration done in rebound.."
