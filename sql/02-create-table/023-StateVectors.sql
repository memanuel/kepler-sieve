-- Create tables for integrated vectors
-- May include asteroids and derived quantities e.g. Earth-Moon Barycenter
CREATE OR REPLACE TABLE KS.StateVectors(
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these state vectors; FK to KS.IntegrationTime",
	BodyID INT NOT NULL
		COMMENT "The Body whose state vectors are described; FK to JS.Body",
	MJD DOUBLE NOT NULL
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame; derivable from IntegrationTimeID but included for performance.",
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
	PRIMARY KEY (TimeID, BodyID)
		COMMENT "A state vector is identified by the body and time stamp; use integer time ID for performance.",
	UNIQUE KEY UNQ_BodyID_TimeID(BodyID, TimeID)
		COMMENT "Allow fast search keyed first by BodyID.",
	CONSTRAINT FK_StateVectors_TimeID
		FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	CONSTRAINT FK_StateVectors_BodyID
		FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
ENGINE='aria' TRANSACTIONAL=0
COMMENT "State vectors (position and velocity) for Solar Systems bodies computed in Rebound using the planets as massive bodies and initial conditions from DE435 at MJD 59000. 
Not limited to bodies used in the integration.  May include asteroids and derived quantities (e.g. Earth-Moon Barycenter)."
